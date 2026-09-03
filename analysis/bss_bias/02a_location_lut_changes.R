# ==============================================================================
# 02a_location_lut_changes.R
#
# Purpose:
#   Answer, precisely, what changed about a fishery's LOCATION DEFINITION
#   across years -- and what did not. Reads only the committed capture from
#   00c_probe_location_lut.R, so it needs no VPN, no database and no network.
#
#   This is not bookkeeping. The capture shows the spatial extent of these
#   fisheries moving substantially between years: Snohomish fall salmon goes
#   from one section on the Snohomish in 2021 to four sections spanning the
#   Snohomish, Skykomish and Snoqualmie in 2024. Comparing `b` across those
#   years is comparing different stretches of river, not the same fishery in
#   two years. Any statement about interannual variability in `b` has to be
#   read against this table.
#
# ------------------------------------------------------------------------------
# TWO THINGS THAT MAKE THE NAIVE COMPARISON WRONG
#
#   1. RIVER MILE IS ONLY MEANINGFUL WITHIN A WATER BODY. Every tributary
#      restarts at RM 0 at its confluence, so pooling river miles across water
#      bodies invents a span nobody surveyed. Skagit spring Chinook upper does
#      not reach down to RM 0 in 2024 -- it ADDS the Cascade River (section 2,
#      Cascade RM 0-0.9, at the confluence near Marblemount) while the Skagit
#      part stays exactly where it was, 67.7-78.2. Read without water body,
#      that looks like the fishery growing 68 river miles downstream.
#
#   2. SECTION NUMBER IS A LABEL, NOT A PLACE. Skagit fall salmon section 2 is
#      RM 11.4-28.8 in 2021 and RM 11.1-22.5 from 2022 on, because section 2
#      was split and everything above it renumbered. Tracking "section 2"
#      across years therefore tracks different water.
#
#   So the comparable unit is (water_body, river-mile interval), and every
#   summary here is grouped by water body. `river_miles` is the length of the
#   UNION of a year's section intervals within a water body, which is
#   section-numbering-independent and does not double-count where sections
#   abut or overlap.
#
# Usage:
#   Rscript analysis/bss_bias/02a_location_lut_changes.R
#
# Inputs:
#   analysis/bss_bias/lookup/fishery_location_lut.csv  (from 00c; tracked)
#
# Outputs (analysis/bss_bias/outputs/):
#   bss_b_lut_section_year.csv      -- fishery_type x year x water body x section, counts + RM span
#   bss_b_lut_water_body_year.csv   -- fishery_type x year x water body: sections, locations, river miles
#   bss_b_lut_changes.csv           -- per fishery_type: which attributes varied, with per-year values
#   bss_b_lut_stability.csv         -- per fishery_type: which attributes held constant
#   bss_b_lut_section_renumbering.csv -- section numbers whose water changed between years
#   bss_b_lut_location_changes.csv  -- per location: present in which years, and added/dropped/persistent
#   bss_b_lut_anomalies.csv         -- case-variant labels, split sections, zero-length spans
#   bss_b_lut_changes.html          -- gt rendering of the change table
#   figures/fig9_lut_section_year.{png,pdf}
#   figures/fig10_lut_river_extent.{png,pdf}
#   figures/fig11_lut_river_miles.{png,pdf}
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
# ("Site" and "site", "Census" and "census"). Left alone they split a group in
# two, and any downstream exact-match test -- the pipeline is full of them, e.g.
# prep_dwg_census_expan()'s `location_type == "Site"` -- silently drops the
# odd-cased rows. Fold the case here, and WRITE OUT which rows needed it: a
# label entered two ways is a data-entry finding in its own right, not just an
# inconvenience for this script.

title_case <- function(x) str_replace(str_to_lower(x), "^(.)", toupper)

case_anomalies <- lut_raw |>
  filter(location_type != title_case(location_type) | survey_type != title_case(survey_type)) |>
  transmute(anomaly = "case-variant type label", fishery_name, section_num,
            water_body_code, location_code,
            detail = paste0("location_type=", location_type, ", survey_type=", survey_type))

lut <- lut_raw |>
  mutate(
    location_type = title_case(location_type),
    survey_type   = title_case(survey_type),
    fishery_type  = fishery_type_from_name(fishery_name),
    # fishery_start_year, not a regex over the name: it is the database's own
    # answer, and it resolves the multi-season names ("2022-23" -> 2022)
    # the same way the b-series parity anchor does.
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

# Water bodies take a FIXED colour each, assigned once in sorted order and used
# in every figure -- colour follows the entity, so the Skykomish is the same
# hue wherever it appears. The eight CAT hues are used in their documented
# order; that palette passes the categorical checks at all eight slots.
WATER_BODIES <- sort(unique(lut$water_body_code))
WB_COLORS <- setNames(unname(CAT)[seq_along(WATER_BODIES)], WATER_BODIES)

# Length of the UNION of a set of river-mile intervals. Sections abut (one ends
# where the next begins) and occasionally overlap, so summing (upper - lower)
# double-counts and max(upper) - min(lower) invents water between disjoint
# pieces. Neither is the river miles actually covered.
interval_union_length <- function(lo, hi) {
  ok <- !is.na(lo) & !is.na(hi) & hi >= lo
  if (!any(ok)) return(NA_real_)
  lo <- lo[ok]; hi <- hi[ok]
  o  <- order(lo); lo <- lo[o]; hi <- hi[o]
  total <- 0; cur_lo <- lo[1]; cur_hi <- hi[1]
  for (i in seq_along(lo)[-1]) {
    if (lo[i] <= cur_hi) {
      cur_hi <- max(cur_hi, hi[i])
    } else {
      total <- total + (cur_hi - cur_lo); cur_lo <- lo[i]; cur_hi <- hi[i]
    }
  }
  total + (cur_hi - cur_lo)
}

finite_or_na <- function(x) if_else(is.finite(x), x, NA_real_)

# ------------------------------------------------------------------------------
# T1: fishery_type x year x water body x section
# ------------------------------------------------------------------------------
# Water body is part of the GRAIN, not a label hung off it: Stillaguamish 2023
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
  mutate(across(c(rm_lower, rm_upper), finite_or_na),
         rm_span = rm_upper - rm_lower) |>
  arrange(basin, fishery_type, year, water_body_code, section_num)

write_csv(section_year, file.path(OUT_DIR, "bss_b_lut_section_year.csv"))
cli::cli_alert_success("Wrote bss_b_lut_section_year.csv ({nrow(section_year)} rows).")

# ------------------------------------------------------------------------------
# T1b: fishery_type x year x water body
# ------------------------------------------------------------------------------
# The section-numbering-independent view. Two years cover the same water if
# their river_miles agree per water body, however the sections were labelled.

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
  left_join(
    section_year |>
      group_by(fishery_type, year, water_body_code) |>
      summarise(river_miles = interval_union_length(rm_lower, rm_upper), .groups = "drop"),
    by = c("fishery_type", "year", "water_body_code")
  ) |>
  arrange(basin, fishery_type, year, water_body_code)

write_csv(water_body_year, file.path(OUT_DIR, "bss_b_lut_water_body_year.csv"))
cli::cli_alert_success("Wrote bss_b_lut_water_body_year.csv ({nrow(water_body_year)} rows).")

cli::cli_h3("River miles by water body and year")
water_body_year |>
  select(fishery_type, year, water_body_code, n_sections, river_miles) |>
  print(n = 100)

# ------------------------------------------------------------------------------
# T1c: section renumbering
# ------------------------------------------------------------------------------
# A section number that denotes different water in different years is the thing
# that makes a naive year-on-year "section 2 vs section 2" comparison wrong.

section_renumbering <- section_year |>
  group_by(basin, fishery_type, section_num) |>
  summarise(
    n_years      = n_distinct(year),
    water_bodies = paste(sort(unique(water_body_code)), collapse = "|"),
    spans        = paste0(year, ": ", water_body_code, " ",
                          round(rm_lower, 1), "-", round(rm_upper, 1), collapse = "; "),
    n_distinct_spans = n_distinct(paste(water_body_code, rm_lower, rm_upper)),
    .groups = "drop"
  ) |>
  mutate(section_moved = n_distinct_spans > 1) |>
  arrange(basin, fishery_type, section_num)

write_csv(section_renumbering, file.path(OUT_DIR, "bss_b_lut_section_renumbering.csv"))
cli::cli_alert_success(
  "Wrote bss_b_lut_section_renumbering.csv ({sum(section_renumbering$section_moved)} of \
   {nrow(section_renumbering)} section numbers denote different water in different years)."
)

# ------------------------------------------------------------------------------
# T2/T3: what varied, what held
# ------------------------------------------------------------------------------
# One row per fishery_type x year x attribute, then split on whether the
# attribute took more than one value across that fishery's years. Reporting the
# constant attributes as well as the varying ones is deliberate: "these six
# things did not change" is the harder half to see and the more reassuring one.

fishery_year_attrs <- lut |>
  group_by(basin = basin_of(fishery_type), fishery_type, year) |>
  summarise(
    n_sections       = n_distinct(section_num),
    sections         = paste(sort(unique(section_num)), collapse = ","),
    n_locations      = n_distinct(location_id),
    n_index_sites    = n_distinct(location_id[survey_type == "Index"]),
    n_census_units   = n_distinct(location_id[survey_type == "Census"]),
    n_water_bodies   = n_distinct(water_body_code),
    water_bodies     = paste(sort(unique(water_body_code)), collapse = ","),
    n_surveyors      = n_distinct(na.omit(surveyor_num)),
    fishery_start    = min(fishery_start_date),
    fishery_end      = max(fishery_end_date),
    .groups = "drop"
  ) |>
  # River miles enter the diff PER WATER BODY, never pooled -- a single
  # rm_lower/rm_upper across water bodies is a span nobody surveyed.
  left_join(
    water_body_year |>
      mutate(wb_miles = paste0(water_body_code, "=", round(river_miles, 1))) |>
      group_by(fishery_type, year) |>
      summarise(river_miles_by_water_body = paste(sort(wb_miles), collapse = "; "),
                river_miles_total = round(sum(river_miles, na.rm = TRUE), 1),
                .groups = "drop"),
    by = c("fishery_type", "year")
  )

attr_long <- fishery_year_attrs |>
  mutate(across(everything(), as.character)) |>
  pivot_longer(-c(basin, fishery_type, year), names_to = "attribute", values_to = "value")

attr_summary <- attr_long |>
  # Explicit: `values` pastes year=value in row order, so the per-year list has
  # to be chronological rather than whatever order the pivot happened to leave.
  arrange(basin, fishery_type, attribute, year) |>
  group_by(basin, fishery_type, attribute) |>
  summarise(
    n_years   = n_distinct(year),
    n_values  = n_distinct(value),
    values    = paste0(year, "=", value, collapse = "; "),
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

# ------------------------------------------------------------------------------
# Anomalies
# ------------------------------------------------------------------------------
# Collected rather than quietly handled. Each is a data-entry finding for
# whoever maintains the lookup, and each has a way of biting silently:
#   * a case-variant label is dropped by the pipeline's exact-match tests, e.g.
#     prep_dwg_census_expan()'s `location_type == "Site"`;
#   * one section number spanning two water bodies breaks (year, section) as a
#     key, here and in any downstream join;
#   * a zero-length river-mile span is a section covering no river.

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

# ------------------------------------------------------------------------------
# T4: location-level presence
# ------------------------------------------------------------------------------
# The finest-grain answer: which individual survey locations came and went.
# `status` classifies the shape of the change, since "dropped after 2022" and
# "used in 2021 and 2024 only" are different stories about a program.

loc_years <- lut |> distinct(basin = basin_of(fishery_type), fishery_type, location_code, section_num, survey_type, year)

fishery_year_span <- loc_years |>
  distinct(fishery_type, year) |>
  group_by(fishery_type) |>
  summarise(all_years = list(sort(unique(year))), .groups = "drop")

location_changes <- loc_years |>
  group_by(basin, fishery_type, location_code, section_num, survey_type) |>
  summarise(years_present = list(sort(unique(year))), .groups = "drop") |>
  left_join(fishery_year_span, by = "fishery_type") |>
  mutate(
    n_years_present = map_int(years_present, length),
    n_years_fishery = map_int(all_years, length),
    first_year      = map_int(years_present, ~ .x[1]),
    last_year       = map_int(years_present, ~ tail(.x, 1)),
    years           = map_chr(years_present, ~ paste(.x, collapse = ",")),
    status = case_when(
      n_years_present == n_years_fishery                       ~ "persistent",
      map2_lgl(years_present, all_years, ~ identical(.x, tail(.y, length(.x)))) ~ "added",
      map2_lgl(years_present, all_years, ~ identical(.x, head(.y, length(.x)))) ~ "dropped",
      TRUE                                                     ~ "intermittent"
    )
  ) |>
  select(basin, fishery_type, location_code, section_num, survey_type,
         status, years, n_years_present, n_years_fishery, first_year, last_year) |>
  arrange(basin, fishery_type, status, section_num, location_code)

write_csv(location_changes, file.path(OUT_DIR, "bss_b_lut_location_changes.csv"))
cli::cli_alert_success("Wrote bss_b_lut_location_changes.csv ({nrow(location_changes)} locations).")

cli::cli_h3("Location turnover by fishery")
location_changes |> count(fishery_type, status) |> pivot_wider(names_from = status, values_from = n, values_fill = 0) |> print()

# ------------------------------------------------------------------------------
# Fig 9 -- section x year coverage, split by water body
# ------------------------------------------------------------------------------
# Form: a grid, because the question is "which cells exist" and an ABSENT cell
# is the finding. Fill is a magnitude (locations in that section), so it takes
# the single-hue sequential ramp; an absent section stays surface-coloured
# rather than becoming another category. Water body is a row band, not a
# colour, because a section number only means something inside one.

fig9 <- section_year |>
  mutate(
    year = factor(year),
    # section_year is already sorted by water body then section number, so
    # fct_inorder() below preserves upstream order -- factor()'s alphabetical
    # default would put section 10 between 1 and 2.
    row_label = paste0(water_body_code, "  s", section_num)
  ) |>
  ggplot(aes(x = year, y = fct_rev(fct_inorder(row_label)), fill = n_locations)) +
  geom_tile(colour = SURFACE, linewidth = 1.2) +
  geom_text(aes(label = n_locations), colour = INK, size = 3) +
  facet_wrap(~fishery_type, scales = "free_y", ncol = 2) +
  scale_fill_gradientn(colours = SEQ_RAMP, name = "Locations") +
  labs(
    title = "Which water each fishery covered, by year",
    subtitle = "One row per water body and section; number = survey locations in it. A blank cell is water the fishery did not cover that year.",
    x = NULL, y = NULL
  ) +
  theme_bss() +
  theme(panel.grid = element_blank(), legend.position = "right")

save_fig(fig9, "fig9_lut_section_year", width = 12, height = 10)
cli::cli_alert_success("Wrote fig9_lut_section_year.")

# ------------------------------------------------------------------------------
# Fig 10 -- river-mile extent, one panel per water body
# ------------------------------------------------------------------------------
# Form: an interval per year, because the question is about a SPAN moving.
# Faceted by fishery x water body and given a free x scale, since river mile
# restarts at 0 in every tributary -- one shared axis would put the Cascade's
# RM 0 next to the Skagit's and imply they are the same place. Colour is the
# water body, fixed per entity across every figure; the facet strip carries the
# same identity in text, so nothing depends on colour alone.

extent <- section_year |>
  filter(!is.na(rm_lower), !is.na(rm_upper)) |>
  mutate(panel = paste0(fishery_type, "\n", water_body_code))

fig10 <- extent |>
  ggplot(aes(y = fct_rev(factor(year)), colour = water_body_code)) +
  geom_segment(aes(x = rm_lower, xend = rm_upper, yend = fct_rev(factor(year))),
               linewidth = 2.5, lineend = "butt", alpha = 0.85) +
  geom_point(aes(x = rm_lower), size = 1.6) +
  geom_point(aes(x = rm_upper), size = 1.6) +
  facet_wrap(~panel, scales = "free_x", ncol = 3) +
  scale_colour_manual(values = WB_COLORS, name = NULL) +
  labs(
    title = "River miles surveyed, by year and water body",
    subtitle = "One segment per section. Each panel has its own river-mile axis: every tributary restarts at RM 0 at its confluence, so miles are comparable only within a water body.",
    x = "River mile", y = NULL
  ) +
  theme_bss() +
  theme(legend.position = "top")

save_fig(fig10, "fig10_lut_river_extent", width = 12, height = 11)
cli::cli_alert_success("Wrote fig10_lut_river_extent.")

# ------------------------------------------------------------------------------
# Fig 11 -- how much water, of which kind
# ------------------------------------------------------------------------------
# The section-numbering-independent summary: river miles covered per year,
# stacked by water body. This is the one to read when asking whether two years
# of `b` describe comparable fisheries -- it does not care how the sections
# were labelled or split.

fig11 <- water_body_year |>
  filter(!is.na(river_miles)) |>
  ggplot(aes(x = factor(year), y = river_miles, fill = water_body_code)) +
  geom_col(width = 0.7, colour = SURFACE, linewidth = 1) +
  facet_wrap(~fishery_type, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = WB_COLORS, name = NULL) +
  labs(
    title = "River miles covered per year",
    subtitle = "Length of the union of a year's section intervals within each water body -- unaffected by how sections were numbered or split.",
    x = NULL, y = "River miles"
  ) +
  theme_bss() +
  theme(legend.position = "top")

save_fig(fig11, "fig11_lut_river_miles", width = 12, height = 7)
cli::cli_alert_success("Wrote fig11_lut_river_miles.")

# ------------------------------------------------------------------------------
# gt rendering of the change table
# ------------------------------------------------------------------------------

gt_tbl <- lut_changes |>
  select(basin, fishery_type, attribute, n_values, values) |>
  gt(groupname_col = "fishery_type") |>
  tab_header(
    title = "What changed in the location lookup, by fishery",
    subtitle = "Only attributes that took more than one value across a fishery's years. Constants are in bss_b_lut_stability.csv."
  ) |>
  cols_label(n_values = "Distinct values", values = "By year") |>
  opt_row_striping()

gtsave(gt_tbl, file.path(OUT_DIR, "bss_b_lut_changes.html"))
cli::cli_alert_success("Wrote bss_b_lut_changes.html")

cli::cli_h2("Most-changed attributes")
lut_changes |> count(attribute, sort = TRUE) |> print(n = 20)
