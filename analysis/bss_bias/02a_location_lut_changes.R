# ==============================================================================
# 02a_location_lut_changes.R
#
# Purpose:
#   Answer, precisely, what changed about a fishery's location definition
#   across years -- and what did not. Reads only the committed capture from
#   00c_probe_location_lut.R, so it needs no VPN, no database and no network.
#
# ------------------------------------------------------------------------------
# WHICH SIGNAL TO TRUST, AND IN WHAT ORDER
#
#   The lookup carries three descriptions of where a fishery happened, and they
#   are NOT equally reliable:
#
#   1. SURVEY SITES (location_id) -- the primary signal. A site is a physical,
#      named place: a boat ramp, an access point, a hole. If the same sites are
#      surveyed in two years, the fishery covered the same water, however the
#      paperwork around them changed. This is the test that answers "is this
#      the same fishery in aggregate".
#
#   2. SECTIONS -- administrative groupings of sites, and they get redrawn. The
#      Skagit split sections at Hwy 9 to give closures a boundary: the spatial
#      units changed while the footprint did not. So a section change is a
#      PROMPT TO CHECK THE CLOSURES, not a finding on its own -- and the
#      closures themselves would have changed fishing dynamics, which is the
#      part that could actually move `b`.
#
#   3. RIVER MILES -- reported here, but treat as unvalidated. They need
#      checking in their own right, so nothing in this script's conclusions
#      rests on them. bss_b_lut_rm_consistency.csv is the check: the same
#      site should carry the same river mile every year, and where it does not,
#      the river miles are not yet trustworthy for that fishery.
#
#   Water body is nonetheless part of the grain throughout, because a section
#   number only identifies a place within one, and every tributary restarts at
#   river mile 0 at its confluence.
#
# Usage:
#   Rscript analysis/bss_bias/02a_location_lut_changes.R
#
# Inputs:
#   analysis/bss_bias/lookup/fishery_location_lut.csv  (from 00c; tracked)
#
# Outputs (analysis/bss_bias/outputs/):
#   bss_b_lut_site_stability.csv    -- PRIMARY: per year, sites retained/added/dropped vs prior and first year
#   bss_b_lut_site_moves.csv        -- sites that persisted but changed section (the re-partitioning signature)
#   bss_b_lut_location_changes.csv  -- per site: years present, persistent/added/dropped/intermittent
#   bss_b_lut_section_year.csv      -- fishery_type x year x water body x section
#   bss_b_lut_water_body_year.csv   -- fishery_type x year x water body
#   bss_b_lut_changes.csv           -- per fishery_type: which attributes varied, with per-year values
#   bss_b_lut_stability.csv         -- per fishery_type: which attributes held constant
#   bss_b_lut_rm_consistency.csv    -- sites whose river mile moved between years (RM validation)
#   bss_b_lut_anomalies.csv         -- case-variant labels, split sections, zero-length spans
#   bss_b_lut_changes.html          -- gt rendering of the change table
#   figures/fig9_lut_site_persistence.{png,pdf}
#   figures/fig10_lut_site_turnover.{png,pdf}
#   figures/fig11_lut_section_year.{png,pdf}
# ==============================================================================

library(tidyverse)
library(gt)
library(cli)
library(here)

source(here::here("analysis", "bss_bias", "common.R"))

OUT_DIR <- here::here("analysis", "bss_bias", "outputs")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

LUT_PATH <- here::here("analysis", "bss_bias", "lookup", "fishery_location_lut.csv")
if (!file.exists(LUT_PATH)) {
  cli::cli_abort(c(
    "{.file {LUT_PATH}} not found.",
    "i" = "Run {.file analysis/bss_bias/00c_probe_location_lut.R} (VPN required) and commit the result."
  ))
}

lut_raw <- read_csv(LUT_PATH, show_col_types = FALSE)
cli::cli_alert_info("Read {nrow(lut_raw)} row(s) covering {n_distinct(lut_raw$fishery_name)} fishery-year(s).")

# ------------------------------------------------------------------------------
# Normalisation
# ------------------------------------------------------------------------------
# location_type and survey_type carry case-variant duplicates in the source
# ("Site"/"site", "Census"/"census"). Left alone they split a group in two, and
# the pipeline's exact-match tests -- prep_dwg_census_expan()'s
# `location_type == "Site"`, for one -- silently drop the odd-cased rows.

title_case <- function(x) str_replace(str_to_lower(x), "^(.)", toupper)

case_anomalies <- lut_raw |>
  filter(location_type != title_case(location_type) | survey_type != title_case(survey_type)) |>
  transmute(anomaly = "case-variant type label", fishery_name,
            section_num = as.integer(section_num), water_body_code, location_code,
            detail = paste0("location_type=", location_type, ", survey_type=", survey_type))

lut <- lut_raw |>
  mutate(
    location_type = title_case(location_type),
    survey_type   = title_case(survey_type),
    fishery_type  = fishery_type_from_name(fishery_name),
    # fishery_start_year, not a regex over the name: the database's own answer,
    # and it resolves the multi-season names ("2022-23" -> 2022) the same way
    # the b-series parity anchor does.
    year          = as.integer(fishery_start_year),
    section_num   = as.integer(section_num),
    across(c(upper_rm, lower_rm), as.numeric)
  )

basin_of <- function(x) {
  case_when(
    str_detect(x, regex("skagit",        ignore_case = TRUE)) ~ "Skagit",
    str_detect(x, regex("snohomish",     ignore_case = TRUE)) ~ "Snohomish",
    str_detect(x, regex("stillaguamish", ignore_case = TRUE)) ~ "Stillaguamish",
    TRUE ~ "Other"
  )
}

# Water bodies take a fixed colour each, assigned once and used everywhere --
# colour follows the entity, so the Skykomish is the same hue in every figure.
# The eight CAT hues are used in their documented order; that palette passes the
# categorical checks at all eight slots.
WATER_BODIES <- sort(unique(lut$water_body_code))
WB_COLORS <- setNames(unname(CAT)[seq_along(WATER_BODIES)], WATER_BODIES)

finite_or_na <- function(x) if_else(is.finite(x), x, NA_real_)

# location_id is the key, not location_code: 141 ids carry only 138 distinct
# codes, so codes are not unique. Codes are for display.
sites <- lut |>
  distinct(basin = basin_of(fishery_type), fishery_type, year,
           location_id, location_code, section_num, water_body_code, survey_type)

# ==============================================================================
# PRIMARY: is it the same set of sites?
# ==============================================================================

jaccard <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

site_sets <- sites |>
  group_by(basin, fishery_type, year) |>
  summarise(ids = list(sort(unique(location_id))), .groups = "drop") |>
  arrange(basin, fishery_type, year)

site_stability <- site_sets |>
  group_by(fishery_type) |>
  mutate(
    # Explicit prior-row shift rather than lag(): lag() on a list column fills
    # the first element with NA, and intersect(x, NA) is character(0) -- a
    # first year would report 0 retained instead of "no prior year".
    prev_ids  = c(list(NULL), ids[-n()]),
    first_ids = rep(ids[1], n()),
    n_sites   = map_int(ids, length),
    # vs the previous year: the year-on-year turnover
    n_retained_prev = map2_int(ids, prev_ids, ~ if (is.null(.y)) NA_integer_ else length(intersect(.x, .y))),
    n_added_prev    = map2_int(ids, prev_ids, ~ if (is.null(.y)) NA_integer_ else length(setdiff(.x, .y))),
    n_dropped_prev  = map2_int(ids, prev_ids, ~ if (is.null(.y)) NA_integer_ else length(setdiff(.y, .x))),
    jaccard_prev    = map2_dbl(ids, prev_ids, ~ if (is.null(.y)) NA_real_ else jaccard(.x, .y)),
    # vs the fishery's FIRST year: cumulative drift, which year-on-year hides
    jaccard_first   = map2_dbl(ids, first_ids, jaccard),
    n_shared_first  = map2_int(ids, first_ids, ~ length(intersect(.x, .y)))
  ) |>
  ungroup() |>
  mutate(
    # The judgement this table exists to support. "Same sites" is the aggregate
    # footprint test; it deliberately says nothing about how they were grouped.
    site_verdict = case_when(
      is.na(jaccard_prev)   ~ "first year",
      jaccard_prev == 1     ~ "same sites as prior year",
      jaccard_prev >= 0.9   ~ "near-identical (>=90% shared)",
      jaccard_prev >= 0.7   ~ "mostly shared (70-90%)",
      TRUE                  ~ "substantially different (<70%)"
    )
  ) |>
  select(basin, fishery_type, year, n_sites, n_retained_prev, n_added_prev, n_dropped_prev,
         jaccard_prev, n_shared_first, jaccard_first, site_verdict)

write_csv(site_stability, file.path(OUT_DIR, "bss_b_lut_site_stability.csv"))
cli::cli_alert_success("Wrote bss_b_lut_site_stability.csv ({nrow(site_stability)} fishery-years).")

cli::cli_h2("Site-set stability (the aggregate-footprint test)")
site_stability |>
  select(fishery_type, year, n_sites, n_added_prev, n_dropped_prev, jaccard_prev, site_verdict) |>
  print(n = 100)

# ------------------------------------------------------------------------------
# Sites that stayed but were re-grouped
# ------------------------------------------------------------------------------
# The Hwy 9 signature: the site is still surveyed, the section it belongs to
# changed. Separating this from "the site went away" is the whole point -- one
# is a paperwork change, the other is a footprint change. A re-grouping is
# still worth following up, because sections are what closures are written
# against, and a closure change WOULD move fishing dynamics.

site_moves <- sites |>
  group_by(basin, fishery_type, location_id, location_code) |>
  filter(n_distinct(paste(water_body_code, section_num)) > 1) |>
  summarise(
    n_years     = n_distinct(year),
    assignments = paste0(year, ": ", water_body_code, " s", section_num, collapse = "; "),
    n_sections  = n_distinct(section_num),
    changed_water_body = n_distinct(water_body_code) > 1,
    .groups = "drop"
  ) |>
  arrange(basin, fishery_type, location_code)

write_csv(site_moves, file.path(OUT_DIR, "bss_b_lut_site_moves.csv"))
cli::cli_alert_success("Wrote bss_b_lut_site_moves.csv ({nrow(site_moves)} sites re-grouped between years).")
if (nrow(site_moves) > 0) {
  cli::cli_h3("Sites kept but re-assigned -- check the closures for these years")
  site_moves |> count(fishery_type, name = "n_sites_regrouped") |> print()
}

# ------------------------------------------------------------------------------
# Per-site presence
# ------------------------------------------------------------------------------

fishery_year_span <- sites |>
  distinct(fishery_type, year) |>
  group_by(fishery_type) |>
  summarise(all_years = list(sort(unique(year))), .groups = "drop")

location_changes <- sites |>
  group_by(basin, fishery_type, location_id, location_code) |>
  summarise(
    years_present = list(sort(unique(year))),
    sections      = paste(sort(unique(section_num)), collapse = "|"),
    water_bodies  = paste(sort(unique(water_body_code)), collapse = "|"),
    survey_types  = paste(sort(unique(survey_type)), collapse = "|"),
    .groups = "drop"
  ) |>
  left_join(fishery_year_span, by = "fishery_type") |>
  mutate(
    n_years_present = map_int(years_present, length),
    n_years_fishery = map_int(all_years, length),
    first_year      = map_int(years_present, ~ .x[1]),
    last_year       = map_int(years_present, ~ tail(.x, 1)),
    years           = map_chr(years_present, ~ paste(.x, collapse = ",")),
    status = case_when(
      n_years_present == n_years_fishery                                       ~ "persistent",
      map2_lgl(years_present, all_years, ~ identical(.x, tail(.y, length(.x)))) ~ "added",
      map2_lgl(years_present, all_years, ~ identical(.x, head(.y, length(.x)))) ~ "dropped",
      TRUE                                                                     ~ "intermittent"
    )
  ) |>
  select(basin, fishery_type, location_code, sections, water_bodies, survey_types,
         status, years, n_years_present, n_years_fishery, first_year, last_year) |>
  arrange(basin, fishery_type, status, location_code)

write_csv(location_changes, file.path(OUT_DIR, "bss_b_lut_location_changes.csv"))
cli::cli_alert_success("Wrote bss_b_lut_location_changes.csv ({nrow(location_changes)} sites).")

cli::cli_h3("Site turnover by fishery")
location_changes |>
  count(fishery_type, status) |>
  pivot_wider(names_from = status, values_from = n, values_fill = 0) |>
  print()

# ==============================================================================
# SUPPORTING: how the sites were grouped
# ==============================================================================
# Water body is part of the grain, not a label hung off it: Stillaguamish 2023
# carries section 5 on both the North and South Forks, so (year, section_num)
# alone does not identify a row.

section_year <- lut |>
  group_by(basin = basin_of(fishery_type), fishery_type, year, water_body_code, section_num) |>
  summarise(
    n_locations = n_distinct(location_id),
    n_index     = n_distinct(location_id[survey_type == "Index"]),
    n_census    = n_distinct(location_id[survey_type == "Census"]),
    rm_lower    = suppressWarnings(min(lower_rm, na.rm = TRUE)),
    rm_upper    = suppressWarnings(max(upper_rm, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  # min/max over an all-NA river-mile group returns +/-Inf, not NA.
  mutate(across(c(rm_lower, rm_upper), finite_or_na), rm_span = rm_upper - rm_lower) |>
  arrange(basin, fishery_type, year, water_body_code, section_num)

write_csv(section_year, file.path(OUT_DIR, "bss_b_lut_section_year.csv"))

water_body_year <- lut |>
  group_by(basin = basin_of(fishery_type), fishery_type, year, water_body_code, water_body_desc) |>
  summarise(
    n_sections  = n_distinct(section_num),
    n_locations = n_distinct(location_id),
    n_index     = n_distinct(location_id[survey_type == "Index"]),
    n_census    = n_distinct(location_id[survey_type == "Census"]),
    rm_lower    = suppressWarnings(min(lower_rm, na.rm = TRUE)),
    rm_upper    = suppressWarnings(max(upper_rm, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  mutate(across(c(rm_lower, rm_upper), finite_or_na)) |>
  arrange(basin, fishery_type, year, water_body_code)

write_csv(water_body_year, file.path(OUT_DIR, "bss_b_lut_water_body_year.csv"))
cli::cli_alert_success("Wrote bss_b_lut_section_year.csv and bss_b_lut_water_body_year.csv.")

# ------------------------------------------------------------------------------
# Attribute-level diff
# ------------------------------------------------------------------------------

fishery_year_attrs <- lut |>
  group_by(basin = basin_of(fishery_type), fishery_type, year) |>
  summarise(
    n_sites          = n_distinct(location_id),
    n_index_sites    = n_distinct(location_id[survey_type == "Index"]),
    n_census_units   = n_distinct(location_id[survey_type == "Census"]),
    n_sections       = n_distinct(section_num),
    sections         = paste(sort(unique(section_num)), collapse = ","),
    n_water_bodies   = n_distinct(water_body_code),
    water_bodies     = paste(sort(unique(water_body_code)), collapse = ","),
    n_surveyors      = n_distinct(na.omit(surveyor_num)),
    fishery_start    = min(fishery_start_date),
    fishery_end      = max(fishery_end_date),
    .groups = "drop"
  )

attr_long <- fishery_year_attrs |>
  mutate(across(everything(), as.character)) |>
  pivot_longer(-c(basin, fishery_type, year), names_to = "attribute", values_to = "value")

attr_summary <- attr_long |>
  # `values` pastes year=value in row order, so the per-year list has to be
  # chronological rather than whatever order the pivot left.
  arrange(basin, fishery_type, attribute, year) |>
  group_by(basin, fishery_type, attribute) |>
  summarise(
    n_years  = n_distinct(year),
    n_values = n_distinct(value),
    values   = paste0(year, "=", value, collapse = "; "),
    .groups = "drop"
  ) |>
  mutate(changed = n_values > 1)

lut_changes <- attr_summary |> filter(changed)  |> arrange(basin, fishery_type, desc(n_values), attribute)
lut_stable  <- attr_summary |> filter(!changed) |> arrange(basin, fishery_type, attribute)

write_csv(lut_changes, file.path(OUT_DIR, "bss_b_lut_changes.csv"))
write_csv(lut_stable,  file.path(OUT_DIR, "bss_b_lut_stability.csv"))
cli::cli_alert_success(
  "Wrote bss_b_lut_changes.csv ({nrow(lut_changes)} varying) and \\
   bss_b_lut_stability.csv ({nrow(lut_stable)} constant)."
)

# ==============================================================================
# River-mile validation
# ==============================================================================
# Not a result -- a check on whether the river miles can be used at all. The
# same site should carry the same river mile every year. Where it does not,
# either the site moved or the river miles were re-measured, and until that is
# resolved no distance-based comparison for that fishery means anything.

rm_consistency <- lut |>
  filter(!is.na(lower_rm) | !is.na(upper_rm)) |>
  group_by(basin = basin_of(fishery_type), fishery_type, location_id, location_code) |>
  summarise(
    n_years      = n_distinct(year),
    n_rm_values  = n_distinct(paste(lower_rm, upper_rm)),
    rm_by_year   = paste0(year, ": ", lower_rm, "-", upper_rm, collapse = "; "),
    .groups = "drop"
  ) |>
  filter(n_years > 1, n_rm_values > 1) |>
  arrange(basin, fishery_type, location_code)

write_csv(rm_consistency, file.path(OUT_DIR, "bss_b_lut_rm_consistency.csv"))
if (nrow(rm_consistency) > 0) {
  cli::cli_alert_warning(
    "{nrow(rm_consistency)} site(s) carry more than one river mile across years -- \\
     river miles are NOT yet reliable for {.val {sort(unique(rm_consistency$fishery_type))}}."
  )
} else {
  cli::cli_alert_success("Every site carries a consistent river mile across years.")
}

# ------------------------------------------------------------------------------
# Anomalies
# ------------------------------------------------------------------------------

split_sections <- section_year |>
  group_by(fishery_type, year, section_num) |>
  filter(n_distinct(water_body_code) > 1) |>
  summarise(wbs = paste(sort(unique(water_body_code)), collapse = " + "), .groups = "drop") |>
  transmute(anomaly = "section spans >1 water body",
            fishery_name = paste(fishery_type, year), section_num,
            water_body_code = wbs, location_code = NA_character_,
            detail = "(year, section_num) is not a unique key for this fishery-year")

zero_spans <- section_year |>
  filter(!is.na(rm_span), rm_span == 0) |>
  transmute(anomaly = "zero-length river-mile span",
            fishery_name = paste(fishery_type, year), section_num, water_body_code,
            location_code = NA_character_,
            detail = paste0("rm ", rm_lower, "-", rm_upper))

anomalies <- bind_rows(case_anomalies, split_sections, zero_spans)
if (nrow(anomalies) > 0) {
  write_csv(anomalies, file.path(OUT_DIR, "bss_b_lut_anomalies.csv"))
  cli::cli_alert_warning("{nrow(anomalies)} anomal{?y/ies} in the lookup:")
  print(anomalies, width = Inf)
} else {
  cli::cli_alert_success("No lookup anomalies found.")
}

# ==============================================================================
# Figures
# ==============================================================================

# ------------------------------------------------------------------------------
# Fig 9 -- site persistence, coloured by section
# ------------------------------------------------------------------------------
# One chart, two questions. A row that is filled in every year is a site the
# fishery kept; a gap is a site added or dropped. A row that stays filled but
# CHANGES COLOUR is a site kept and re-grouped -- the Hwy 9 case, where the
# footprint held and only the paperwork moved.
#
# Section number takes the sequential ramp, not categorical hues: sections are
# ordinal (numbered downstream to upstream), and there are up to nine of them,
# past where a categorical palette should be asked to go.

site_grid <- sites |>
  group_by(fishery_type, location_id) |>
  mutate(order_key = min(section_num)) |>
  ungroup() |>
  arrange(fishery_type, order_key, location_code) |>
  mutate(site_label = fct_rev(fct_inorder(paste0(location_code, "  [", water_body_code, "]"))))

fig9 <- site_grid |>
  ggplot(aes(x = factor(year), y = site_label, fill = section_num)) +
  geom_tile(colour = SURFACE, linewidth = 0.8) +
  facet_wrap(~fishery_type, scales = "free_y", ncol = 2) +
  scale_fill_gradientn(colours = SEQ_RAMP, name = "Section", breaks = scales::breaks_width(2)) +
  labs(
    title = "Which survey sites each fishery used, by year",
    subtitle = "A gap is a site added or dropped. A row that stays filled but changes colour is a site kept and re-assigned to a different section -- the footprint held, the grouping moved.",
    x = NULL, y = NULL
  ) +
  theme_bss() +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 6),
        legend.position = "right")

save_fig(fig9, "fig9_lut_site_persistence", width = 13, height = 14)
cli::cli_alert_success("Wrote fig9_lut_site_persistence.")

# ------------------------------------------------------------------------------
# Fig 10 -- year-on-year site turnover
# ------------------------------------------------------------------------------
# The summary of fig 9: how many sites carried over, how many are new, how many
# went away. Retained/added/dropped is a status-like distinction, so it takes
# the reserved status hues with a legend rather than arbitrary series colours.

turnover <- site_stability |>
  filter(!is.na(n_retained_prev)) |>
  select(fishery_type, year, Retained = n_retained_prev, Added = n_added_prev, Dropped = n_dropped_prev) |>
  pivot_longer(c(Retained, Added, Dropped), names_to = "change", values_to = "n") |>
  mutate(change = factor(change, levels = c("Retained", "Added", "Dropped")))

fig10 <- turnover |>
  ggplot(aes(x = factor(year), y = n, fill = change)) +
  geom_col(width = 0.7, colour = SURFACE, linewidth = 1) +
  facet_wrap(~fishery_type, scales = "free_y", ncol = 3) +
  scale_fill_manual(
    values = c(Retained = BASELINE_COL, Added = STATUS[["good"]], Dropped = STATUS[["serious"]]),
    name = NULL
  ) +
  labs(
    title = "Site turnover against the prior year",
    subtitle = "Retained sites are the fishery's continuity. A year that is mostly grey covered the same water as the one before it, whatever happened to the section numbering.",
    x = NULL, y = "Sites"
  ) +
  theme_bss() +
  theme(legend.position = "top")

save_fig(fig10, "fig10_lut_site_turnover", width = 12, height = 7)
cli::cli_alert_success("Wrote fig10_lut_site_turnover.")

# ------------------------------------------------------------------------------
# Fig 11 -- section x year, by water body
# ------------------------------------------------------------------------------
# The supporting view: how the sites were partitioned. Read this AFTER fig 9 --
# a change here with no change there is a re-grouping, not a footprint change.

fig11 <- section_year |>
  mutate(row_label = fct_rev(fct_inorder(paste0(water_body_code, "  s", section_num)))) |>
  ggplot(aes(x = factor(year), y = row_label, fill = n_locations)) +
  geom_tile(colour = SURFACE, linewidth = 1.2) +
  geom_text(aes(label = n_locations), colour = INK, size = 3) +
  facet_wrap(~fishery_type, scales = "free_y", ncol = 2) +
  scale_fill_gradientn(colours = SEQ_RAMP, name = "Sites") +
  labs(
    title = "How sites were grouped into sections, by year",
    subtitle = "Supporting view. A change here without a matching change in site persistence is a re-drawn boundary, not a change in the water fished.",
    x = NULL, y = NULL
  ) +
  theme_bss() +
  theme(panel.grid = element_blank(), legend.position = "right")

save_fig(fig11, "fig11_lut_section_year", width = 12, height = 10)
cli::cli_alert_success("Wrote fig11_lut_section_year.")

# ------------------------------------------------------------------------------
# gt rendering
# ------------------------------------------------------------------------------

gt_tbl <- site_stability |>
  select(basin, fishery_type, year, n_sites, n_added_prev, n_dropped_prev,
         jaccard_prev, jaccard_first, site_verdict) |>
  mutate(across(c(jaccard_prev, jaccard_first), ~ round(.x, 3))) |>
  gt(groupname_col = "fishery_type") |>
  tab_header(
    title = "Is it the same fishery, year to year?",
    subtitle = "Measured on the set of survey sites, which are physical places -- not on sections, which are administrative and get redrawn."
  ) |>
  cols_label(n_sites = "Sites", n_added_prev = "Added", n_dropped_prev = "Dropped",
             jaccard_prev = "Overlap vs prior", jaccard_first = "Overlap vs first yr",
             site_verdict = "Verdict") |>
  tab_footnote(
    footnote = "Overlap is |shared| / |union| of site sets. 1.00 means exactly the same sites.",
    locations = cells_column_labels(columns = jaccard_prev)
  ) |>
  opt_row_striping()

gtsave(gt_tbl, file.path(OUT_DIR, "bss_b_lut_changes.html"))
cli::cli_alert_success("Wrote bss_b_lut_changes.html")
