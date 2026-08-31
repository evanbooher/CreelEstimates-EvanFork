# ==============================================================================
# 02_build_comparability_table.R
#
# Purpose:
#   Turn the raw per-fishery comparability rows written by 01_fit_bss_bias.R
#   (bss_b_comparability_raw.csv -- written BEFORE Stan fitting, for every
#   fishery-year that made it through data prep) into the standalone
#   comparability deliverable: a basin x year x named-fishery table of date
#   windows, spatial coverage, and survey-design flags, with an explicit
#   comparability tier per row.
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
    flag_thin_intA                 = FALSE  # filled in below once bss_b_stan_dims.csv is joined, if available
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
  select(basin, fishery_label, season_label, chosen_est_cg, season_month_span,
         crc_areas, section_nums, study_design, count_types_present,
         p_TI_bank, p_TI_boat, comparability_tier,
         flag_crc_changed, flag_sections_changed, flag_window_shifted,
         flag_pTI_changed, flag_design_changed, flag_counttypes_changed,
         flag_no_vehicle_counts, flag_thin_intA)

# Render booleans as plain glyphs rather than relying on a specific gt version's
# fmt_icon() (icon-font support varies by gt version) -- keeps this robust for
# an unattended overnight/next-morning run.
display <- display |>
  mutate(across(starts_with("flag_"), ~ if_else(.x, "changed", "—")))

gt_tbl <- display |>
  gt(groupname_col = "basin", rowname_col = "season_label") |>
  tab_header(title = "BSS effort-bias (b) comparability table",
             subtitle = "Reviewed BEFORE the b-vs-year plots -- see README.md") |>
  data_color(columns = comparability_tier,
             fn = scales::col_factor(palette = unname(tier_colors), domain = names(tier_colors))) |>
  cols_label(fishery_label = "Fishery", season_month_span = "Season window",
             crc_areas = "CRC area(s)", section_nums = "Section(s)",
             comparability_tier = "Tier") |>
  opt_row_striping()

gt::gtsave(gt_tbl, file.path(OUT_DIR, "bss_b_comparability.html"))
cli::cli_alert_success("Wrote bss_b_comparability.html")
