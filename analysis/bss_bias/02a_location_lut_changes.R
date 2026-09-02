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
#   from 1 section / 12 locations / RM 0-20.5 in 2021 to 4 sections / 45
#   locations / RM 0-49.6 in 2024. Comparing `b` across those years is
#   comparing different stretches of river, not the same fishery in two
#   years. Any statement about interannual variability in `b` has to be read
#   against this table.
#
# Usage:
#   Rscript analysis/bss_bias/02a_location_lut_changes.R
#
# Inputs:
#   analysis/bss_bias/lookup/fishery_location_lut.csv  (from 00c; tracked)
#
# Outputs (analysis/bss_bias/outputs/):
#   bss_b_lut_section_year.csv      -- fishery_type x year x section, counts + river-mile span
#   bss_b_lut_changes.csv           -- per fishery_type: which attributes varied, with per-year values
#   bss_b_lut_stability.csv         -- per fishery_type: which attributes held constant
#   bss_b_lut_location_changes.csv  -- per location: present in which years, and added/dropped/persistent
#   bss_b_lut_case_anomalies.csv    -- rows whose type labels differ only in case (written only if any)
#   bss_b_lut_changes.html          -- gt rendering of the change table
#   figures/fig9_lut_section_year.{png,pdf}
#   figures/fig10_lut_river_extent.{png,pdf}
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
  select(fishery_name, section_num, location_code, location_type, survey_type)

if (nrow(case_anomalies) > 0) {
  cli::cli_alert_warning("{nrow(case_anomalies)} row(s) have case-variant type labels:")
  print(case_anomalies)
  write_csv(case_anomalies, file.path(OUT_DIR, "bss_b_lut_case_anomalies.csv"))
}

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

# ------------------------------------------------------------------------------
# T1: fishery_type x year x section
# ------------------------------------------------------------------------------

section_year <- lut |>
  group_by(basin = basin_of(fishery_type), fishery_type, year, section_num) |>
  summarise(
    n_locations = n_distinct(location_id),
    n_index     = n_distinct(location_id[survey_type == "Index"]),
    n_census    = n_distinct(location_id[survey_type == "Census"]),
    rm_lower    = suppressWarnings(min(lower_rm, na.rm = TRUE)),
    rm_upper    = suppressWarnings(max(upper_rm, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  # min/max over an all-NA river-mile group returns +/-Inf, not NA.
  mutate(across(c(rm_lower, rm_upper), ~ if_else(is.finite(.x), .x, NA_real_))) |>
  arrange(basin, fishery_type, year, section_num)

write_csv(section_year, file.path(OUT_DIR, "bss_b_lut_section_year.csv"))
cli::cli_alert_success("Wrote bss_b_lut_section_year.csv ({nrow(section_year)} rows).")

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
    water_bodies     = paste(sort(unique(water_body_code)), collapse = ","),
    rm_lower         = suppressWarnings(min(lower_rm, na.rm = TRUE)),
    rm_upper         = suppressWarnings(max(upper_rm, na.rm = TRUE)),
    n_surveyors      = n_distinct(na.omit(surveyor_num)),
    fishery_start    = min(fishery_start_date),
    fishery_end      = max(fishery_end_date),
    .groups = "drop"
  ) |>
  mutate(across(c(rm_lower, rm_upper), ~ if_else(is.finite(.x), .x, NA_real_)))

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
# Fig 9 -- section x year coverage
# ------------------------------------------------------------------------------
# Form: a grid, because the question is "which cells exist", and an ABSENT cell
# is the finding. Fill is a magnitude (locations in that section), so it takes
# the single-hue sequential ramp; an absent section is left as surface, not
# given a colour, so gaps read as gaps rather than as another category.

fig9 <- section_year |>
  # Numeric levels, not factor()'s default alphabetical sort -- section 10 must
  # not land between 1 and 2.
  mutate(year = factor(year),
         section_num = factor(section_num, levels = sort(unique(section_num)))) |>
  ggplot(aes(x = year, y = fct_rev(section_num), fill = n_locations)) +
  geom_tile(colour = SURFACE, linewidth = 1.2) +
  geom_text(aes(label = n_locations), colour = INK, size = 3) +
  facet_wrap(~fishery_type, scales = "free_y", ncol = 2) +
  scale_fill_gradientn(colours = SEQ_RAMP, name = "Locations") +
  labs(
    title = "Which sections a fishery covered, by year",
    subtitle = "Cell = section present in the location lookup; number = survey locations in it. A blank cell is a section the fishery did not cover that year.",
    x = NULL, y = "Section"
  ) +
  theme_bss() +
  theme(panel.grid = element_blank(), legend.position = "right")

save_fig(fig9, "fig9_lut_section_year", width = 11, height = 8)
cli::cli_alert_success("Wrote fig9_lut_section_year.")

# ------------------------------------------------------------------------------
# Fig 10 -- river-mile extent
# ------------------------------------------------------------------------------
# Form: an interval per year, because the question is about a SPAN moving, and a
# span is what a segment shows. Two series (index vs census), so the first two
# categorical hues in their documented order, plus a legend -- identity is never
# carried by colour alone.

extent <- lut |>
  filter(!is.na(lower_rm) | !is.na(upper_rm)) |>
  group_by(basin = basin_of(fishery_type), fishery_type, year, survey_type) |>
  summarise(
    rm_lower = suppressWarnings(min(lower_rm, na.rm = TRUE)),
    rm_upper = suppressWarnings(max(upper_rm, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  filter(is.finite(rm_lower), is.finite(rm_upper))

fig10 <- extent |>
  ggplot(aes(y = fct_rev(factor(year)), colour = survey_type)) +
  geom_segment(aes(x = rm_lower, xend = rm_upper, yend = fct_rev(factor(year))),
               linewidth = 2, lineend = "round",
               position = position_dodge(width = 0.5)) +
  facet_wrap(~fishery_type, scales = "free_x", ncol = 2) +
  scale_colour_manual(values = c(Index = CAT[["blue"]], Census = CAT[["orange"]]), name = NULL) +
  labs(
    title = "River-mile extent surveyed, by year",
    subtitle = "Span from the lowest to the highest river mile in the location lookup. A span that moves between years is a fishery that is not spatially comparable across them.",
    x = "River mile", y = NULL
  ) +
  theme_bss() +
  theme(legend.position = "top")

save_fig(fig10, "fig10_lut_river_extent", width = 11, height = 7)
cli::cli_alert_success("Wrote fig10_lut_river_extent.")

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
