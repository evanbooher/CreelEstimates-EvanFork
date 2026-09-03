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
#   `b` is the ratio linking an index count to the census count it is paired
#   with, and it is a single pooled scalar (vector[G] b in the Stan model, no
#   section index). So the quantity that actually informs `b` is not "was this
#   section surveyed" but "did this section-day carry BOTH an index and a
#   census count". That is the tie_in_day flag, and it is the headline here
#   rather than a derived afterthought: a fishery-year with plenty of survey
#   days but few paired ones has a weakly-informed `b` however busy its
#   calendar looks.
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
#   figures/fig13_survey_calendar.{png,pdf}
#   figures/fig14_season_ramp.{png,pdf}
#   figures/fig15_tie_in_days.{png,pdf}
# ==============================================================================

library(tidyverse)
library(cli)
library(here)
library(creelutils)
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

  in_window <- function(d) {
    if (is.null(d) || nrow(d) == 0 || !"event_date" %in% names(d)) return(d[0, , drop = FALSE])
    d |> filter(between(event_date, d_start, d_end))
  }
  eff <- in_window(dwg$effort)
  int <- in_window(dwg$interview)

  # Sections from the UNION of every table that names one, so a section that
  # was surveyed but never defined (or vice versa) still appears as a row.
  sections <- sort(unique(na.omit(as.double(c(
    dwg$fishery_manager$section_num, eff$section_num, int$section_num
  )))))
  if (length(sections) == 0) {
    cli::cli_alert_warning("  {.val {fishery_name}}: no section numbers anywhere.")
    return(NULL)
  }

  closures <- dwg$closures
  closed <- if (is.null(closures) || nrow(closures) == 0) {
    tibble(section_num = double(), event_date = as.Date(character()))
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
STATUS_COLORS <- c(
  "closed"          = GRID_COLOR,
  "open, no survey" = BASELINE_COL,
  "interviews only" = STATUS[["warning"]],
  "census only"     = CAT[["orange"]],
  "index only"      = CAT[["yellow"]],
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

plot_calendar <- function(df, ftype) {
  df |>
    mutate(section_label = fct_rev(fct_inorder(paste0(water_body, "  s", section_num)))) |>
    ggplot(aes(x = event_date, y = section_label, fill = status)) +
    geom_tile(colour = NA) +
    facet_wrap(~year, ncol = 1, scales = "free_y", strip.position = "right") +
    scale_fill_manual(values = STATUS_COLORS, drop = FALSE, name = NULL) +
    scale_x_date(date_labels = "%b %d", date_breaks = "2 weeks") +
    labs(
      title = paste0(ftype, ": what was surveyed, by date and section"),
      subtitle = "Each row is a section, each column a day of the estimation window. A row that starts blank and fills in later is a section that opened later. 'Index + census' is the pairing that informs b.",
      x = NULL, y = NULL
    ) +
    theme_bss() +
    theme(panel.grid = element_blank(), legend.position = "top",
          axis.text.y = element_text(size = 7))
}

for (ft in sort(unique(survey_days$fishery_type))) {
  df <- survey_days |> filter(fishery_type == ft) |> arrange(water_body, section_num)
  n_rows <- df |> distinct(year, water_body, section_num) |> nrow()
  save_fig(plot_calendar(df, ft),
           paste0("fig13_survey_calendar_", safe_name(ft)),
           width = 13, height = max(5, 0.30 * n_rows + 2))
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
  geom_line(linewidth = 1) +
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
# The `b`-relevant summary. `b` is informed only by section-days carrying BOTH
# an index and a census count, so this is the count of observations actually
# behind each year's estimate -- which a survey-day total can overstate badly.

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
    subtitle = "b is informed only by section-days carrying both an index count and the census count it pairs with. A year with many survey days but few paired ones has a weakly-informed b.",
    x = NULL, y = "Section-days with both counts"
  ) +
  theme_bss() +
  theme(legend.position = "top")

save_fig(fig15, "fig15_tie_in_days", width = 12, height = 7)
cli::cli_alert_success("Wrote fig15_tie_in_days.")

cli::cli_h2("Done")
cli::cli_alert_info(
  "Read {.file bss_b_survey_date_match.csv} for the date-by-date answer, and \\
   fig13 per fishery for the season shape."
)
