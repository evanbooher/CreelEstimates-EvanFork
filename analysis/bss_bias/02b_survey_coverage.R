# ==============================================================================
# 02b_survey_coverage.R
#
# Purpose:
#   What was actually SURVEYED, by date and section, for every fishery-year --
#   and how that evolved through the season across years.
#
#   The location lookup (02a) says what was DEFINED. It cannot answer this,
#   because sections do not all open at once: lower sections generally open
#   first and upper sections later, so a section defined for a fishery-year may
#   carry no survey at all until weeks into the season. Two years with an
#   identical defined footprint can have quite different realised coverage.
#
#   The question this answers, per fishery: on a given date, was the same
#   section surveyed as in other years -- and where did the season's coverage
#   ramp up, hold, and fall away?
#
# ------------------------------------------------------------------------------
# WHAT COUNTS AS "SURVEYED", AND WHY THE DISTINCTION MATTERS FOR `b`
#
#   INDEX counts   -- dwg$effort, tie_in_indicator == 0. Counted at sites.
#   CENSUS counts  -- dwg$effort, tie_in_indicator == 1. Totals the whole
#                     section-block, scheduled against an index count.
#   INTERVIEWS     -- dwg$interview. Can occur anywhere in the block.
#
#   `b` is a single pooled scalar (vector[G] b in the Stan model, no section
#   index) carried by the index-count likelihood and NOT by the census one, so
#   the census counts are what set the scale of the latent effort the index
#   counts are compared against. Since a census count is dropped unless an
#   index count exists for the same section on the same day, the tie_in_day
#   flag is effectively the count of census ANCHORS -- and it is the headline
#   here rather than a derived afterthought: a fishery-year with plenty of
#   survey days but few paired ones has little pinning `b` down, however busy
#   its calendar looks.
#
# ------------------------------------------------------------------------------
# TWO THINGS THIS TABLE IS CAREFUL ABOUT
#
#   OPEN IS ASSUME-OPEN. prep_days() builds expand_grid(..., open = TRUE) and
#   flips to FALSE only where a closure record exists, so a missing closure and
#   a genuinely open day are indistinguishable. `open` here is reported the
#   same way and should be read as "not recorded closed".
#
#   SECTION NUMBERS MOVE BETWEEN YEARS. They are labels on blocks, and blocks
#   get redrawn (see 02a). Every output here therefore carries water body
#   alongside the number, and the across-year comparison is offered on both
#   calendar date and day-of-season so a shifted window does not masquerade as
#   a coverage change.
#
# Usage:
#   Rscript analysis/bss_bias/02b_survey_coverage.R
#   No VPN. Reads the captured windows in lookup/fishery_params.csv and caches
#   each fetch under outputs/cache/dwg/, so re-runs are fast.
#
# Outputs (analysis/bss_bias/outputs/):
#   bss_b_survey_days.csv           -- ATOMIC: fishery x date x section, all counts and flags
#   bss_b_survey_by_fishery.csv     -- per fishery-year: days surveyed, tie-in days, coverage
#   bss_b_survey_by_section.csv     -- per fishery-year x section: first/last survey, spread
#   bss_b_survey_date_match.csv     -- per fishery x month-day x section: which years surveyed it
#   bss_b_fishery_composition.{csv,html} -- one line per fishery-year: defined vs surveyed
#   bss_b_common_window_cost.csv    -- the window all years share, and the tie-in days a refit there costs
#   bss_b_parity_confound.csv       -- whether pink parity is separable from window start, per fishery
#   figures/fig13_survey_calendar.{png,pdf}
#   figures/fig14_season_ramp.{png,pdf}
#   figures/fig15_tie_in_days.{png,pdf}
#   figures/fig16_fishery_composition.{png,pdf}
# ==============================================================================

library(tidyverse)
library(cli)
library(here)
library(creelutils)
library(gt)
library(zoo)   # rollmean, for the season-ramp smoothing

source(here::here("analysis", "bss_bias", "common.R"))
source(here::here("analysis", "bss_bias", "fishery_data.R"))

OUT_DIR <- here::here("analysis", "bss_bias", "outputs")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

on.exit(disconnect_fishery_data(), add = TRUE)

# Same source of truth for scope as 01: whatever is queued for fitting is what
# gets a coverage row, so the two describe the same set of fishery-years.
discovery_path <- file.path(OUT_DIR, "fishery_discovery_target.csv")
if (!file.exists(discovery_path)) {
  cli::cli_abort("Run 00_discover_fisheries.R first -- {.file {discovery_path}} not found.")
}
target_fisheries <- read_csv(discovery_path, show_col_types = FALSE) |>
  filter(basin_match == "target", include_in_run) |>
  pull(fishery_name_raw) |>
  unique()

cli::cli_alert_info("Building survey coverage for {length(target_fisheries)} fishery-year(s).")

# Water body per (fishery_name, section_num), from the committed lookup. A bare
# section number does not identify a place, so it is carried everywhere.
LUT_PATH <- here::here("analysis", "bss_bias", "lookup", "fishery_location_lut.csv")
section_water_body <- if (file.exists(LUT_PATH)) {
  read_csv(LUT_PATH, show_col_types = FALSE) |>
    mutate(section_num = as.double(section_num)) |>
    distinct(fishery_name, section_num, water_body_code) |>
    # A section straddling two water bodies (Stillaguamish 2023-24 section 5)
    # gets both, joined as one label rather than duplicating the survey row.
    group_by(fishery_name, section_num) |>
    summarise(water_body = paste(sort(unique(water_body_code)), collapse = " + "), .groups = "drop")
} else {
  cli::cli_alert_warning("{.file {LUT_PATH}} not found -- sections will carry no water body.")
  tibble(fishery_name = character(), section_num = double(), water_body = character())
}

# ------------------------------------------------------------------------------
# Per-fishery coverage grid
# ------------------------------------------------------------------------------
# Built on the FULL date x section grid, not just the rows that carry data, so
# a day with no survey is a row saying so rather than an absence to be inferred.

survey_grid_for <- function(fishery_name) {
  est_dates <- resolve_window(fishery_name)
  if (is.null(est_dates)) {
    cli::cli_alert_warning("  {.val {fishery_name}}: {window_skip_reason()}")
    return(NULL)
  }
  d_start <- as.Date(est_dates$est_date_start)
  d_end   <- as.Date(est_dates$est_date_end)

  dwg <- tryCatch(fetch_fishery_dwg(fishery_name, est_dates), error = function(e) {
    cli::cli_alert_danger("  {.val {fishery_name}}: fetch failed -- {conditionMessage(e)}")
    NULL
  })
  if (is.null(dwg)) return(NULL)

  # Every zero-row path has to return the SAME COLUMNS as the populated one.
  # A tibble that is merely empty still joins and filters correctly; one that is
  # missing a column fails later, far from here, on whatever first refers to it.
  EFF_COLS <- tibble(event_date = as.Date(character()), section_num = double(),
                     tie_in_indicator = double())
  INT_COLS <- tibble(event_date = as.Date(character()), section_num = double())

  in_window <- function(d, template) {
    if (is.null(d) || nrow(d) == 0 || !all(names(template) %in% names(d))) return(template)
    d |> filter(between(event_date, d_start, d_end))
  }
  eff <- in_window(dwg$effort, EFF_COLS)
  int <- in_window(dwg$interview, INT_COLS)

  # Sections from the UNION of every table that names one, so a section that
  # was surveyed but never defined (or vice versa) still appears as a row.
  sections <- sort(unique(na.omit(as.double(c(
    dwg$fishery_manager$section_num, eff$section_num, int$section_num
  )))))
  if (length(sections) == 0) {
    cli::cli_alert_warning("  {.val {fishery_name}}: no section numbers anywhere.")
    return(NULL)
  }

  # recorded_closed must be present in BOTH branches: without it the empty case
  # joins in no such column and the `open` mutate below has nothing to coalesce.
  closures <- dwg$closures
  closed <- if (is.null(closures) || nrow(closures) == 0) {
    tibble(section_num = double(), event_date = as.Date(character()),
           recorded_closed = logical())
  } else {
    closures |>
      mutate(event_date = as.Date(substr(as.character(event_date), 1, 10)),
             section_num = as.double(section_num)) |>
      filter(between(event_date, d_start, d_end)) |>
      distinct(section_num, event_date) |>
      mutate(recorded_closed = TRUE)
  }

  tally <- function(d, name) {
    if (nrow(d) == 0) return(tibble(event_date = as.Date(character()), section_num = double(),
                                    !!name := integer()))
    d |>
      mutate(section_num = as.double(section_num)) |>
      count(event_date, section_num, name = name)
  }

  idx <- tally(eff |> filter(tie_in_indicator == 0), "n_index_counts")
  cen <- tally(eff |> filter(tie_in_indicator == 1), "n_census_counts")
  ints <- tally(int, "n_interviews")

  expand_grid(event_date = seq(d_start, d_end, by = "day"), section_num = sections) |>
    left_join(idx,  by = c("event_date", "section_num")) |>
    left_join(cen,  by = c("event_date", "section_num")) |>
    left_join(ints, by = c("event_date", "section_num")) |>
    left_join(closed, by = c("event_date", "section_num")) |>
    mutate(
      fishery_name = fishery_name,
      across(c(n_index_counts, n_census_counts, n_interviews), ~ coalesce(.x, 0L)),
      # "Not recorded closed" -- absence of a closure record is not evidence of
      # an open day. See the header.
      open       = !coalesce(recorded_closed, FALSE),
      surveyed   = n_index_counts > 0 | n_census_counts > 0 | n_interviews > 0,
      # The quantity that actually informs `b`: an index count and the census
      # count it pairs with, in the same block on the same day.
      tie_in_day = n_index_counts > 0 & n_census_counts > 0,
      status = case_when(
        !open                                    ~ "closed",
        tie_in_day                               ~ "index + census",
        n_index_counts  > 0                      ~ "index only",
        n_census_counts > 0                      ~ "census only",
        n_interviews    > 0                      ~ "interviews only",
        TRUE                                     ~ "open, no survey"
      )
    ) |>
    select(-recorded_closed)
}

grids <- map(target_fisheries, function(fn) {
  cli::cli_alert_info("{.val {fn}}")
  survey_grid_for(fn)
})
survey_days <- compact(grids) |> bind_rows()

if (nrow(survey_days) == 0) cli::cli_abort("No fishery produced a coverage grid.")

survey_days <- survey_days |>
  left_join(section_water_body, by = c("fishery_name", "section_num")) |>
  mutate(
    water_body   = coalesce(water_body, "(not in lookup)"),
    fishery_type = fishery_type_from_name(fishery_name),
    year         = as.integer(str_extract(fishery_name, "\\d{4}")),
    month_day    = format(event_date, "%m-%d"),
    # Day of year, so years stack on one axis. Guarded against a window that
    # crosses 31 December -- none currently do, but a wrapped doy would fold
    # January back before September and silently scramble the panel.
    doy          = as.integer(format(event_date, "%j")),
    doy          = if (max(doy) - min(doy) > 300) if_else(doy < 100L, doy + 365L, doy) else doy,
    # Day of season, so a fishery whose window shifts between years can still
    # be compared on "how far into the season are we".
    day_of_season = as.integer(event_date - min(event_date)) + 1L,
    .by = fishery_name
  ) |>
  relocate(fishery_name, fishery_type, year, event_date, day_of_season, month_day,
           water_body, section_num)

write_csv(survey_days, file.path(OUT_DIR, "bss_b_survey_days.csv"))
cli::cli_alert_success("Wrote bss_b_survey_days.csv ({nrow(survey_days)} date x section rows).")

# ------------------------------------------------------------------------------
# Per fishery-year
# ------------------------------------------------------------------------------

by_fishery <- survey_days |>
  group_by(fishery_type, year, fishery_name) |>
  summarise(
    n_days_window     = n_distinct(event_date),
    n_sections        = n_distinct(section_num),
    n_days_surveyed   = n_distinct(event_date[surveyed]),
    n_days_tie_in     = n_distinct(event_date[tie_in_day]),
    n_section_days    = n(),
    n_sd_surveyed     = sum(surveyed),
    n_sd_tie_in       = sum(tie_in_day),
    n_sd_closed       = sum(!open),
    pct_days_surveyed = round(100 * n_distinct(event_date[surveyed]) / n_distinct(event_date), 1),
    pct_sd_tie_in     = round(100 * sum(tie_in_day) / n(), 1),
    .groups = "drop"
  ) |>
  arrange(fishery_type, year)

write_csv(by_fishery, file.path(OUT_DIR, "bss_b_survey_by_fishery.csv"))
cli::cli_h2("Coverage per fishery-year")
by_fishery |>
  select(fishery_type, year, n_days_window, n_sections, n_days_surveyed,
         n_days_tie_in, pct_days_surveyed, pct_sd_tie_in) |>
  print(n = 100)

# ------------------------------------------------------------------------------
# Per fishery-year x section -- when each section came into the season
# ------------------------------------------------------------------------------
# This is the sequential-opening pattern, measured: first and last surveyed
# date per section, as a day of season so years are comparable.

by_section <- survey_days |>
  filter(surveyed) |>
  group_by(fishery_type, year, water_body, section_num) |>
  summarise(
    first_survey      = min(event_date),
    last_survey       = max(event_date),
    first_doy_season  = min(day_of_season),
    last_doy_season   = max(day_of_season),
    n_days_surveyed   = n_distinct(event_date),
    n_days_tie_in     = n_distinct(event_date[tie_in_day]),
    .groups = "drop"
  ) |>
  arrange(fishery_type, year, water_body, section_num)

write_csv(by_section, file.path(OUT_DIR, "bss_b_survey_by_section.csv"))
cli::cli_alert_success("Wrote bss_b_survey_by_section.csv.")

cli::cli_h2("When each section entered the season (day of season of first survey)")
by_section |>
  select(fishery_type, year, water_body, section_num, first_doy_season, last_doy_season, n_days_surveyed) |>
  print(n = 200)

# ------------------------------------------------------------------------------
# On a given date, was the same section surveyed across years?
# ------------------------------------------------------------------------------
# The direct answer. Keyed on calendar month-day, since that is the question as
# asked -- a section on 15 September, across every year the fishery ran.

date_match <- survey_days |>
  group_by(fishery_type, month_day, water_body, section_num) |>
  summarise(
    n_years_in_window  = n_distinct(year),
    n_years_surveyed   = n_distinct(year[surveyed]),
    n_years_tie_in     = n_distinct(year[tie_in_day]),
    years_in_window    = paste(sort(unique(year)), collapse = ","),
    years_surveyed     = paste(sort(unique(year[surveyed])), collapse = ","),
    .groups = "drop"
  ) |>
  mutate(
    agreement = case_when(
      n_years_surveyed == 0                  ~ "never surveyed on this date",
      n_years_surveyed == n_years_in_window  ~ "surveyed in every year in window",
      TRUE                                   ~ "surveyed in some years only"
    )
  ) |>
  arrange(fishery_type, month_day, water_body, section_num)

write_csv(date_match, file.path(OUT_DIR, "bss_b_survey_date_match.csv"))
cli::cli_alert_success("Wrote bss_b_survey_date_match.csv.")

cli::cli_h2("Date-by-date agreement across years")
date_match |>
  filter(n_years_in_window > 1) |>
  count(fishery_type, agreement) |>
  pivot_wider(names_from = agreement, values_from = n, values_fill = 0) |>
  print()

# ------------------------------------------------------------------------------
# What each year actually covered
# ------------------------------------------------------------------------------
# The one-line-per-fishery-year answer to "how has this fishery changed", with
# what was DEFINED and what was SURVEYED side by side. They are not the same:
# Stillaguamish 2022-23 defined nine sections and surveyed eight, and its
# section 7 was surveyed on exactly one day. A defined footprint is a plan; the
# survey record is what the estimate rests on.
#
# The last column is the one that governs `b`. Everything else describes the
# fishery; paired index + census days describe the evidence.

lut_defined <- if (file.exists(LUT_PATH)) {
  read_csv(LUT_PATH, show_col_types = FALSE) |>
    mutate(fishery_type = fishery_type_from_name(fishery_name),
           year = as.integer(fishery_start_year)) |>
    group_by(fishery_type, year) |>
    summarise(
      water_bodies       = paste(sort(unique(str_replace(water_body_code, "^Stillaguamish - ", ""))),
                                 collapse = "+"),
      sections_defined   = n_distinct(section_num),
      index_sites_defined = n_distinct(location_id[str_to_lower(location_type) == "site"]),
      .groups = "drop"
    )
} else NULL

composition <- by_fishery |>
  select(fishery_type, year, n_days_window, n_sections, n_days_surveyed, n_days_tie_in) |>
  rename(sections_surveyed = n_sections)

if (!is.null(lut_defined)) {
  composition <- composition |>
    left_join(lut_defined, by = c("fishery_type", "year")) |>
    relocate(water_bodies, sections_defined, index_sites_defined, .after = year)
}

composition <- composition |>
  mutate(parity = year_parity_label(year)) |>
  arrange(fishery_type, year)

write_csv(composition, file.path(OUT_DIR, "bss_b_fishery_composition.csv"))
cli::cli_alert_success("Wrote bss_b_fishery_composition.csv.")

cli::cli_h2("What each year covered")
print(composition, n = 40, width = Inf)

if (requireNamespace("gt", quietly = TRUE) && !is.null(lut_defined)) {
  gt_comp <- composition |>
    gt::gt(groupname_col = "fishery_type") |>
    gt::tab_header(
      title = "What each fishery-year covered",
      subtitle = "Defined comes from the location lookup, surveyed from the effort and interview record. The last column is what informs b."
    ) |>
    gt::tab_spanner("Defined", columns = c(water_bodies, sections_defined, index_sites_defined)) |>
    gt::tab_spanner("Surveyed", columns = c(n_days_window, sections_surveyed, n_days_surveyed, n_days_tie_in)) |>
    gt::cols_label(water_bodies = "Water", sections_defined = "Sections",
                   index_sites_defined = "Index sites", n_days_window = "Window days",
                   sections_surveyed = "Sections", n_days_surveyed = "Days surveyed",
                   n_days_tie_in = "Paired index+census days", parity = "Year") |>
    gt::opt_row_striping()
  gt::gtsave(gt_comp, file.path(OUT_DIR, "bss_b_fishery_composition.html"))
  cli::cli_alert_success("Wrote bss_b_fishery_composition.html")
}

# ------------------------------------------------------------------------------
# The common window, and what restricting to it would cost
# ------------------------------------------------------------------------------
# `b` is a SINGLE POOLED SCALAR -- vector[G] b, no day or section index -- so it
# is an effort-weighted average of the day-and-section ratios inside whatever
# window the fishery-year ran. Two years with different windows therefore
# average over different periods, and a difference between their `b` values is
# not attributable to anything else about those years.
#
# That bites hardest on the pink-year question. In Skagit fall salmon the odd
# (pink) years start mid-August and the even years start 1 September, so parity
# and window start are perfectly aligned across the five years available. The
# same holds for Snohomish. Odd-year `b` averages in ~2.5 weeks that even-year
# `b` never sees -- weeks with different species composition, a different
# spatial spread of effort, and fewer sections open. A pink-vs-even difference
# in `b` cannot be separated from that.
#
# Note this is NOT the "pink years have more effort" objection. `b` multiplies
# the expected count given latent effort, so a uniform change in effort leaves
# it untouched. `b` moves when effort redistributes relative to the index
# sites -- which is a real and plausible pink mechanism, and exactly the one
# the window confound would masquerade as.
#
# The common window is the intersection of every year's window for a fishery,
# on day-of-year. Refitting there is the test that separates the two. This
# table is what that test would cost in the currency that matters: paired
# index + census days.

# survey_days$doy, not a fresh format() call -- that column carries the
# year-boundary guard and these must agree with the calendar figure.
common_window <- survey_days |>
  group_by(fishery_type, year) |>
  summarise(doy_start = min(doy), doy_end = max(doy), .groups = "drop") |>
  group_by(fishery_type) |>
  summarise(common_doy_start = max(doy_start),
            common_doy_end   = min(doy_end), .groups = "drop") |>
  mutate(common_n_days = common_doy_end - common_doy_start + 1,
         common_start = format(as.Date(common_doy_start - 1, origin = "2001-01-01"), "%b %d"),
         common_end   = format(as.Date(common_doy_end   - 1, origin = "2001-01-01"), "%b %d"))

common_window_cost <- survey_days |>
  left_join(common_window, by = "fishery_type") |>
  mutate(in_common = doy >= common_doy_start & doy <= common_doy_end) |>
  group_by(fishery_type, year, common_start, common_end, common_n_days) |>
  summarise(
    n_days_window     = n_distinct(event_date),
    tie_in_days_all   = n_distinct(event_date[tie_in_day]),
    tie_in_days_common = n_distinct(event_date[tie_in_day & in_common]),
    .groups = "drop"
  ) |>
  mutate(tie_in_days_lost = tie_in_days_all - tie_in_days_common,
         parity = year_parity_label(year)) |>
  arrange(fishery_type, year)

write_csv(common_window_cost, file.path(OUT_DIR, "bss_b_common_window_cost.csv"))
cli::cli_alert_success("Wrote bss_b_common_window_cost.csv.")

cli::cli_h2("Common window across years, and the tie-in days a refit there would cost")
common_window_cost |>
  select(fishery_type, year, parity, common_start, common_end,
         n_days_window, tie_in_days_all, tie_in_days_common, tie_in_days_lost) |>
  print(n = 40)

# Is parity confounded with window start? Where every odd year starts earlier
# than every even year, the two cannot be told apart with these data.
parity_confound <- survey_days |>
  group_by(fishery_type, year) |>
  summarise(doy_start = min(doy), .groups = "drop") |>
  mutate(parity = year_parity_label(year)) |>
  group_by(fishery_type) |>
  summarise(
    odd_starts  = paste(sort(doy_start[parity == "odd (pink)"]), collapse = ","),
    even_starts = paste(sort(doy_start[parity == "even"]), collapse = ","),
    separable = !(all(doy_start[parity == "odd (pink)"] < min(doy_start[parity == "even"])) |
                  all(doy_start[parity == "odd (pink)"] > max(doy_start[parity == "even"]))),
    .groups = "drop"
  ) |>
  mutate(verdict = if_else(separable,
                           "start dates overlap -- parity is separable",
                           "every odd year starts on one side of every even year -- CONFOUNDED"))

write_csv(parity_confound, file.path(OUT_DIR, "bss_b_parity_confound.csv"))
cli::cli_h2("Is a pink-year effect separable from window start?")
print(parity_confound, n = 20, width = Inf)

# ==============================================================================
# Figures
# ==============================================================================
# Survey status is a small ordered set of STATES, not a magnitude and not an
# identity, so it takes the reserved status hues plus two neutrals -- and every
# figure that uses it carries a legend, since a state must never be read from
# colour alone. The order runs worst-to-best for `b`: a closed day, an open day
# nobody worked, a day worked but not paired, and a paired day.
STATUS_LEVELS <- c("closed", "open, no survey", "interviews only",
                   "census only", "index only", "index + census")
# Two neutrals for the states with no survey, then four hues in an order the
# palette validator passes at every adjacent pair -- the previous version put
# CAT yellow next to STATUS warning, two yellows a reader cannot separate.
STATUS_COLORS <- c(
  "closed"          = GRID_COLOR,
  "open, no survey" = BASELINE_COL,
  "interviews only" = CAT[["violet"]],
  "census only"     = CAT[["orange"]],
  "index only"      = CAT[["aqua"]],
  "index + census"  = CAT[["blue"]]
)

survey_days <- survey_days |> mutate(status = factor(status, levels = STATUS_LEVELS))

# Water bodies keep a fixed hue each, assigned in sorted order exactly as 02a
# does, so a water body is the same colour in both scripts' figures.
WB_LEVELS <- sort(unique(survey_days$water_body))
WB_COLORS_SURVEY <- setNames(rep_len(unname(CAT), length(WB_LEVELS)), WB_LEVELS)

# ------------------------------------------------------------------------------
# Fig 13 -- the survey calendar
# ------------------------------------------------------------------------------
# One panel per fishery-year, rows are sections, columns are calendar days. This
# is where sequential opening is visible directly: a section whose row starts
# blank and fills in later opened later. Reading down the panels of one fishery
# compares years on the same calendar axis.

# Day of year on the x axis, not calendar date: on a date axis the four
# Stillaguamish seasons occupy four thin slivers of a 2022-2025 span, and the
# whole point is to compare them. Stacked on day of year they sit directly
# above one another.
doy_label <- function(x) format(as.Date(x - 1, origin = "2001-01-01"), "%b %d")

plot_calendar <- function(df, ftype) {
  rng <- range(df$doy)
  df |>
    mutate(
      # The fishery name is in the title; repeating the basin in every row
      # label costs a third of the plot width.
      #
      # Ordered on SECTION NUMBER, not water body. Sections run downstream to
      # upstream, so the number is the spatial order; sorting by water body
      # first scatters them -- 2023 Snohomish reads Skykomish s3, Snohomish s1,
      # Snohomish s2 down the panel. fct_rev puts s1 at the top.
      section_label = fct_rev(fct_inorder(paste0(
        str_remove(water_body, "^[A-Za-z]+ - "), "  s", section_num
      )))
    ) |>
    ggplot(aes(x = doy, y = section_label, fill = status)) +
    geom_tile() +
    # space = "free_y" so a year with 3 sections gets a third the height of one
    # with 9, instead of the same band padded out with whitespace.
    facet_grid(year ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_manual(values = STATUS_COLORS, drop = FALSE, name = NULL) +
    scale_x_continuous(
      limits = rng, expand = c(0, 0),
      breaks = scales::breaks_width(14), labels = doy_label
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      title = paste0(ftype, ": what was surveyed, by date and section"),
      subtitle = "Fishery years stacked on a shared day-of-year axis. 'index + census' is pairing anchoring b term.",
      x = NULL, y = NULL
    ) +
    theme_bss() +
    theme(
      panel.grid = element_blank(),
      panel.spacing.y = unit(3, "pt"),
      legend.position = "top",
      legend.key.size = unit(10, "pt"),
      strip.placement = "outside",
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size = 8)
    )
}

for (ft in sort(unique(survey_days$fishery_type))) {
  # section_num first: fct_inorder() inside plot_calendar takes its level order
  # from this arrange, and the section number is the spatial order.
  df <- survey_days |> filter(fishery_type == ft) |> arrange(section_num, water_body)
  # Height from the TOTAL rows across years, since space = "free_y" makes each
  # panel proportional rather than equal.
  n_rows <- df |> distinct(year, water_body, section_num) |> nrow()
  save_fig(plot_calendar(df, ft),
           paste0("fig13_survey_calendar_", safe_name(ft)),
           width = 11, height = max(4, 0.22 * n_rows + 1.8))
}
cli::cli_alert_success("Wrote fig13_survey_calendar_* ({n_distinct(survey_days$fishery_type)} fisheries).")

# ------------------------------------------------------------------------------
# Fig 14 -- season ramp, years overlaid
# ------------------------------------------------------------------------------
# The across-year comparison in one panel per fishery: how many sections were
# actually being surveyed as the season progressed. Lower sections opening
# first and upper ones later shows up as a rising line; a year that ramps on a
# different schedule separates from the others visibly.
#
# Plotted against DAY OF SEASON rather than calendar date, so a fishery whose
# window start moves between years (Snohomish 2022 began 29 October, 2023 on 1
# September) is compared on how far into its own season it was.

smooth_days <- function(x, k = 7) {
  k <- min(k, length(x))
  if (k < 3) return(as.numeric(x))          # too short to smooth; show it raw
  if (k %% 2 == 0) k <- k - 1               # rollmean(align = "center") wants odd
  zoo::rollmean(x, k = k, fill = NA, align = "center")
}

ramp <- survey_days |>
  group_by(fishery_type, year, day_of_season) |>
  summarise(n_sections_surveyed = n_distinct(section_num[surveyed]),
            n_sections_tie_in   = n_distinct(section_num[tie_in_day]),
            .groups = "drop") |>
  # A daily line over a survey schedule is mostly zeros; a rolling mean shows
  # the schedule rather than the sampling calendar's sawtooth. The window
  # shrinks for a short season -- Stillaguamish 2023-24 ran 12 days -- because
  # rollmean() errors outright when k exceeds the series length.
  group_by(fishery_type, year) |>
  arrange(day_of_season, .by_group = TRUE) |>
  mutate(across(c(n_sections_surveyed, n_sections_tie_in),
                ~ smooth_days(.x), .names = "{.col}_7d")) |>
  ungroup()

fig14 <- ramp |>
  ggplot(aes(x = day_of_season, y = n_sections_surveyed_7d, colour = factor(year))) +
  # na.rm: the rolling mean is undefined for the half-window at each end of a
  # season. Those NAs are expected, so drop them quietly rather than emitting a
  # "removed N rows" warning that reads like a data problem.
  geom_line(linewidth = 1, na.rm = TRUE) +
  facet_wrap(~fishery_type, scales = "free", ncol = 2) +
  scale_colour_manual(values = colorRampPalette(SEQ_RAMP)(n_distinct(ramp$year)),
                      name = NULL) +
  labs(
    title = "How coverage built through the season, by year",
    subtitle = "Sections under survey, 7-day rolling mean, against day of season rather than calendar date so a shifted window is not read as a coverage change. A rising line is sections opening in sequence.",
    x = "Day of season", y = "Sections surveyed"
  ) +
  theme_bss() +
  theme(legend.position = "top")

save_fig(fig14, "fig14_season_ramp", width = 12, height = 9)
cli::cli_alert_success("Wrote fig14_season_ramp.")

# ------------------------------------------------------------------------------
# Fig 15 -- paired index+census days
# ------------------------------------------------------------------------------
# The `b`-relevant summary. Census counts are what set the scale of the latent
# effort that index counts are compared against, and prep_dwg_effort_census()
# drops a census count with no same-day index count in the same section -- so
# this is effectively the number of census ANCHORS available to a year's
# estimate. Index counts elsewhere in the section still contribute, through the
# season-level effort mean they share, but with no census in a section there is
# nothing to compare them against and `b` falls back on its prior.

fig15 <- survey_days |>
  filter(tie_in_day) |>
  count(fishery_type, year, water_body, section_num, name = "n_tie_in_days") |>
  mutate(section_label = paste0(water_body, "  s", section_num)) |>
  ggplot(aes(x = factor(year), y = n_tie_in_days, fill = water_body)) +
  geom_col(width = 0.7, colour = SURFACE, linewidth = 0.8) +
  facet_wrap(~fishery_type, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = WB_COLORS_SURVEY, name = NULL) +
  labs(
    title = "Paired index + census days, the observations behind b",
    subtitle = "Census counts set the scale of the effort that index counts are measured against, and are recorded only alongside a same-day index count. A year with many survey days but few paired ones has few anchors for b.",
    x = NULL, y = "Section-days with both counts"
  ) +
  theme_bss() +
  theme(legend.position = "top")

save_fig(fig15, "fig15_tie_in_days", width = 12, height = 7)
cli::cli_alert_success("Wrote fig15_tie_in_days.")

# ------------------------------------------------------------------------------
# Fig 16 -- how each fishery changed, year to year
# ------------------------------------------------------------------------------
# The composition table read visually. Four measures on four rows sharing one
# year axis, because they are on wildly different scales -- days in the window
# run to 140, paired index+census days to 18 -- and forcing them onto one axis
# would flatten the measure that matters most. Each row keeps its own scale and
# says what it is.
#
# The rows run from what was planned to what the estimate rests on:
#   sections defined -> sections surveyed -> days surveyed -> paired days.
# A gap opening between consecutive rows is where a plan stopped becoming data.

MEASURE_LEVELS <- c(
  "Sections defined"          = "sections_defined",
  "Sections surveyed"         = "sections_surveyed",
  "Days surveyed"             = "n_days_surveyed",
  "Paired index+census days"  = "n_days_tie_in"
)
have <- MEASURE_LEVELS[MEASURE_LEVELS %in% names(composition)]

comp_long <- composition |>
  select(fishery_type, year, all_of(unname(have))) |>
  pivot_longer(-c(fishery_type, year), names_to = "measure", values_to = "value") |>
  mutate(measure = factor(measure, levels = unname(have), labels = names(have)))

# One series per panel, so no colour encoding and no legend: the facet strip
# already names the fishery, and a hue here would carry no information.
fig16 <- comp_long |>
  ggplot(aes(x = year, y = value, group = fishery_type)) +
  geom_line(linewidth = 0.7, colour = CAT[["blue"]], alpha = 0.85) +
  geom_point(size = 2.2, colour = CAT[["blue"]]) +
  geom_text(aes(label = value), vjust = -0.9, size = 2.7, colour = INK_SECOND) +
  facet_grid(measure ~ fishery_type, scales = "free", switch = "y") +
  scale_x_continuous(breaks = sort(unique(comp_long$year))) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.22))) +
  labs(
    title = "How each fishery changed, year to year",
    subtitle = "Top row is what the location lookup defined; the rest is what the survey record holds. The bottom row -- section-days carrying both an index and a census count -- is the only one that informs b.",
    x = NULL, y = NULL
  ) +
  theme_bss() +
  theme(strip.placement = "outside", panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

save_fig(fig16, "fig16_fishery_composition", width = 14, height = 8)
cli::cli_alert_success("Wrote fig16_fishery_composition.")

cli::cli_h2("Done")
cli::cli_alert_info(
  "Read {.file bss_b_survey_date_match.csv} for the date-by-date answer, and \\
   fig13 per fishery for the season shape."
)
