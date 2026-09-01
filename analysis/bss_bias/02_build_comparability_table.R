# ==============================================================================
# 02_build_comparability_table.R
#
# Purpose:
#   Turn the raw per-fishery comparability rows written by 01_fit_bss_bias.R
#   (bss_b_comparability_raw.csv -- written BEFORE Stan fitting, for every
#   fishery-year that made it through data prep) into the standalone
#   comparability deliverable: a basin x year x named-fishery table of season
#   duration and timing, spatial coverage, and survey-design flags, with an
#   explicit comparability tier per row.
#
#   SEASON DURATION AND TIMING are shown alongside the other comparability
#   factors, not just as a month span. They vary far more than the month span
#   suggests -- Skagit fall salmon runs 91-140 days across years, Snohomish
#   fall salmon 33-91, and Stillaguamish 2023-24 was 12 days (season cut
#   short on Chinook impacts). Some of the interannual variation in `b` is
#   therefore a window artifact rather than a change in vehicle-index bias,
#   which matters directly when interpreting how variable `b` is across
#   years. start_shift_days / pct_of_ref_days quantify it per fishery-year;
#   both are REPORTED ONLY and do not affect the tier (see the case_when()).
#
#   NO MCMC REQUIRED. This script only needs bss_b_comparability_raw.csv to
#   exist with at least one row -- it does not read bss_b_summary.csv or any
#   fit output. Per README.md's go/no-go section, THIS TABLE ALONE is a
#   legitimate meeting deliverable even if every Stan fit is still running,
#   partially failed, or was never attempted: it answers "how many
#   comparable fishery-years do we even have," which the group needs to
#   agree on before any b-vs-year number matters.
#
# Usage:
#   Rscript analysis/bss_bias/02_build_comparability_table.R
#
# Outputs (analysis/bss_bias/outputs/):
#   bss_b_comparability.csv   -- the flagged/tiered table
#   bss_b_comparability.html  -- gt-rendered version (grouped by basin, tiers colored)
# ==============================================================================

library(tidyverse)
library(gt)
library(cli)
library(here)

OUT_DIR <- here::here("analysis", "bss_bias", "outputs")
raw_path <- file.path(OUT_DIR, "bss_b_comparability_raw.csv")

if (!file.exists(raw_path)) {
  cli::cli_abort(
    "{.file {raw_path}} not found -- run 01_fit_bss_bias.R first (it writes \\
     this file per fishery-year during data prep, before any Stan fitting)."
  )
}

raw <- read_csv(raw_path, show_col_types = FALSE) |> distinct(fishery_name, .keep_all = TRUE)

# ------------------------------------------------------------------------------
# Basin / display-label derivation. Prefers a hand-maintained lookup
# (analysis/bss_bias/lookup/fishery_labels.csv: fishery_name, basin,
# fishery_label, season_label) when a row exists there; otherwise derives
# basin from the name by the same regex Phase 0 uses, and falls back to the
# bare fishery_name as the label. Edit the lookup file to fix a mis-derived
# label rather than editing this script.
# ------------------------------------------------------------------------------

lookup_path <- here::here("analysis", "bss_bias", "lookup", "fishery_labels.csv")
lookup <- if (file.exists(lookup_path)) {
  read_csv(lookup_path, show_col_types = FALSE)
} else {
  tibble(fishery_name = character(), basin = character(), fishery_label = character(), season_label = character())
}

derive_basin <- function(x) {
  case_when(
    str_detect(x, regex("skagit", ignore_case = TRUE)) ~ "Skagit",
    str_detect(x, regex("snohomish", ignore_case = TRUE)) ~ "Snohomish",
    str_detect(x, regex("stillaguamish", ignore_case = TRUE)) ~ "Stillaguamish",
    TRUE ~ "Other"
  )
}

comp <- raw |>
  left_join(lookup, by = "fishery_name") |>
  mutate(
    basin = coalesce(basin, derive_basin(fishery_name)),
    fishery_label = coalesce(fishery_label, fishery_name),
    season_label = coalesce(season_label, str_extract(fishery_name, "\\d{4}(-\\d{2,4})?")),
    year_start = suppressWarnings(as.integer(str_extract(season_label, "^\\d{4}")))
  ) |>
  arrange(basin, fishery_label, year_start)

# ------------------------------------------------------------------------------
# Per-series (basin x fishery_label) comparison against the most recent year
# ------------------------------------------------------------------------------

# Tolerances for the season duration/timing comparison. Explicit and
# renegotiable, like the tier rules below -- state them in the room.
# Chosen to be loose enough that ordinary year-to-year scheduling jitter does
# not flag, but tight enough to catch a season that is materially shorter or
# starts in a different part of the run:
#   Skagit fall salmon runs 91-140 days across years; Snohomish 33-91.
#   Stillaguamish 2023-24 was 12 days (season cut short on Chinook impacts) --
#   a real management action, and exactly the kind of thing that should be
#   visible in this table rather than buried in the raw CSV.
DURATION_TOL_PCT <- 0.25   # window length within +/-25% of the reference year
TIMING_TOL_DAYS  <- 14     # season start within +/-14 days of the reference year's start

comp <- comp |>
  group_by(basin, fishery_label) |>
  arrange(desc(year_start), .by_group = TRUE) |>
  mutate(
    is_reference = row_number() == 1,
    ref_crc_areas          = crc_areas[is_reference][1],
    ref_section_nums       = section_nums[is_reference][1],
    ref_season_month_span  = season_month_span[is_reference][1],
    ref_count_types        = count_types_present[is_reference][1],
    ref_p_TI_bank           = p_TI_bank[is_reference][1],
    ref_p_TI_boat            = p_TI_boat[is_reference][1],
    ref_study_design          = study_design[is_reference][1],
    ref_n_days                 = n_days_in_window[is_reference][1],
    ref_start_doy               = as.integer(format(date_start[is_reference][1], "%j")),

    # DURATION and TIMING, relative to the reference year. season_month_span
    # alone misses both: it treats Aug 14 and Aug 30 as the same start, and a
    # 91-day and a 140-day "Sep-Nov" window as identical.
    start_doy          = as.integer(format(date_start, "%j")),
    start_shift_days   = start_doy - ref_start_doy,
    pct_of_ref_days    = n_days_in_window / ref_n_days,

    flag_crc_changed      = !is_reference & (crc_areas != ref_crc_areas),
    flag_sections_changed  = !is_reference & (section_nums != ref_section_nums),
    flag_window_shifted     = !is_reference & (season_month_span != ref_season_month_span),
    flag_counttypes_changed  = !is_reference & (count_types_present != ref_count_types),
    flag_pTI_changed          = !is_reference & (
      abs(coalesce(p_TI_bank, -999) - coalesce(ref_p_TI_bank, -999)) > 1e-6 |
      abs(coalesce(p_TI_boat, -999) - coalesce(ref_p_TI_boat, -999)) > 1e-6
    ),
    flag_design_changed        = !is_reference & (study_design != ref_study_design),
    flag_no_vehicle_counts       = !has_vehicle_counts,
    flag_thin_intA                 = FALSE,  # filled in below once bss_b_stan_dims.csv is joined, if available

    # REPORTED AND FLAGGED, but deliberately NOT wired into comparability_tier
    # below -- surfacing season duration/timing is a different decision from
    # letting it downgrade a fishery-year, and that second one belongs to the
    # group. To make them count, add these two to the "comparable-with-caveat"
    # line of the case_when().
    flag_duration_changed = !is_reference & !is.na(pct_of_ref_days) &
      abs(pct_of_ref_days - 1) > DURATION_TOL_PCT,
    flag_timing_shifted   = !is_reference & !is.na(start_shift_days) &
      abs(start_shift_days) > TIMING_TOL_DAYS
  ) |>
  ungroup() |>
  arrange(basin, fishery_label, desc(year_start))

# Join stan-dims (IntA) if 01 has produced it yet -- optional, not required.
dims_path <- file.path(OUT_DIR, "bss_b_stan_dims.csv")
if (file.exists(dims_path)) {
  dims <- read_csv(dims_path, show_col_types = FALSE) |> select(fishery_name, IntA, b_weakly_informed)
  comp <- comp |>
    left_join(dims, by = "fishery_name") |>
    mutate(flag_thin_intA = coalesce(b_weakly_informed, FALSE))
}

# ------------------------------------------------------------------------------
# Comparability tier -- explicit, renegotiable rules (state them in the room)
# ------------------------------------------------------------------------------

comp <- comp |>
  mutate(
    comparability_tier = case_when(
      flag_no_vehicle_counts ~ "not-estimable",
      is_reference ~ "reference",
      flag_design_changed | flag_counttypes_changed | flag_pTI_changed ~ "not-comparable",
      flag_crc_changed | flag_sections_changed | flag_window_shifted | flag_thin_intA ~ "comparable-with-caveat",
      TRUE ~ "comparable"
    )
  )

write_csv(comp, file.path(OUT_DIR, "bss_b_comparability.csv"))
cli::cli_alert_success("Wrote {nrow(comp)} rows to bss_b_comparability.csv")
comp |> count(basin, comparability_tier) |> print(n = 50)

# ------------------------------------------------------------------------------
# gt render, grouped by basin, tier-colored
# ------------------------------------------------------------------------------

tier_colors <- c(
  "reference"              = "#e8e8e8",
  "comparable"              = "#d7f0d1",
  "comparable-with-caveat"   = "#fdeecb",
  "not-comparable"            = "#f7d6d6",
  "not-estimable"               = "#e2e2e2"
)

display <- comp |>
  select(basin, fishery_label, season_label, comparability_tier,
         # Season duration and timing, up front rather than derived-and-hidden:
         # these were already computed into bss_b_comparability_raw.csv but
         # never surfaced, so a 12-day season and a 140-day one looked alike.
         date_start, date_end, n_days_in_window, n_days_open,
         start_shift_days, pct_of_ref_days, season_month_span,
         chosen_est_cg, crc_areas, section_nums, study_design, count_types_present,
         p_TI_bank, p_TI_boat,
         flag_duration_changed, flag_timing_shifted,
         flag_crc_changed, flag_sections_changed, flag_window_shifted,
         flag_pTI_changed, flag_design_changed, flag_counttypes_changed,
         flag_no_vehicle_counts, flag_thin_intA)

# Render booleans as plain glyphs rather than relying on a specific gt version's
# fmt_icon() (icon-font support varies by gt version) -- keeps this robust for
# an unattended overnight/next-morning run. Formatted BEFORE the boolean
# conversion so the numeric columns keep their types.
display <- display |>
  mutate(
    pct_of_ref_days  = if_else(is.na(pct_of_ref_days), NA_character_,
                                paste0(round(100 * pct_of_ref_days), "%")),
    start_shift_days = case_when(
      is.na(start_shift_days)  ~ NA_character_,
      start_shift_days == 0     ~ "—",
      start_shift_days > 0      ~ paste0("+", start_shift_days, "d"),
      TRUE                       ~ paste0(start_shift_days, "d")
    )
  ) |>
  mutate(across(starts_with("flag_"), ~ if_else(.x, "changed", "—")))

gt_tbl <- display |>
  gt(groupname_col = "basin", rowname_col = "season_label") |>
  tab_header(title = "BSS effort-bias (b) comparability table",
             subtitle = "Season duration and timing shown alongside the other comparability factors -- see README.md") |>
  data_color(columns = comparability_tier,
             fn = scales::col_factor(palette = unname(tier_colors), domain = names(tier_colors))) |>
  tab_spanner(label = "Season duration & timing",
              columns = c(date_start, date_end, n_days_in_window, n_days_open,
                          start_shift_days, pct_of_ref_days, season_month_span)) |>
  tab_spanner(label = "Spatial & design",
              columns = c(crc_areas, section_nums, study_design, count_types_present,
                          p_TI_bank, p_TI_boat)) |>
  cols_label(fishery_label = "Fishery", season_month_span = "Month span",
             date_start = "Start", date_end = "End",
             n_days_in_window = "Days", n_days_open = "Open days",
             start_shift_days = "Start vs ref", pct_of_ref_days = "Length vs ref",
             crc_areas = "CRC area(s)", section_nums = "Section(s)",
             comparability_tier = "Tier") |>
  tab_footnote(
    footnote = paste0(
      "Start vs ref / Length vs ref compare each year to the most recent year in its series. ",
      "Flagged beyond +/-", TIMING_TOL_DAYS, " days or +/-", round(100 * DURATION_TOL_PCT),
      "% respectively. These are REPORTED ONLY -- they do not currently affect the Tier."
    ),
    locations = cells_column_labels(columns = c(start_shift_days, pct_of_ref_days))
  ) |>
  opt_row_striping()

gt::gtsave(gt_tbl, file.path(OUT_DIR, "bss_b_comparability.html"))
cli::cli_alert_success("Wrote bss_b_comparability.html")
