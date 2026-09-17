# ==============================================================================
# 13_census_index_ratio.R -- daily paired census-to-index ratios
#
# The empirical quantity behind `b`, one point per section per paired day,
# rather than the single seasonal scalar the model reports. This is what shows
# WHY a `b` is what it is: the mainstem 2025 vehicle ratios climb across the
# three anchor dates while 2024's sit flat, and one scalar has to reconcile
# anchors that disagree systematically.
#
# Reads the cached DWG through the same scope rules the fits use, so the days
# plotted are exactly the days that informed `b`. No fits required, no VPN.
#
#   Rscript analysis/bss_bias/13_census_index_ratio.R
#
# Outputs (analysis/bss_bias/outputs/):
#   bss_b_census_index_ratio_daily.csv
#   figures/fig21_census_index_ratio_vehicle.png
#   figures/fig21_census_index_ratio_trailer.png
#
# NOT plotted: the fitted `b` as a reference line. The raw ratio is census
# anglers per index OBJECT; `b` is the bias on the index after expansion by
# anglers-per-vehicle. They differ by that expansion factor, so drawing them on
# one axis would invite reading a gap that is not there. Within a year the
# expansion is a constant, so the SHAPE of these series is the shape of `b`'s
# evidence -- which is the point of the figure.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(cli); library(here)
})

OUT_DIR <- here::here("analysis", "bss_bias", "outputs")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

source(here::here("analysis", "bss_bias", "common.R"))
source(here::here("analysis", "bss_bias", "fishery_data.R"))
source(here::here("analysis", "bss_bias", "scope_rules.R"))

FITS <- tibble::tribble(
  ~scope, ~fishery_name,
  "MS",   "Stillaguamish salmon and gamefish 2024-25",
  "MS",   "Stillaguamish salmon and gamefish 2025-26",
  "NF",   "Stillaguamish salmon and gamefish 2024-25",
  "NF",   "Stillaguamish salmon and gamefish 2025-26"
) |>
  mutate(year = as.integer(str_extract(fishery_name, "\\d{4}")),
         fork = if_else(scope == "MS", "Mainstem", "North Fork"))

bank_or_boat <- function(d) {
  if ("angler_final" %in% names(d)) return(as.character(d$angler_final))
  lbl <- as.character(d$count_type %||% rep(NA_character_, nrow(d)))
  dplyr::if_else(str_starts(lbl, "Boat"), "boat", "bank")
}

# ------------------------------------------------------------------------------
# One scope-resolved pull per fit
# ------------------------------------------------------------------------------
ratios_one <- function(scope_tag, fishery_name, year, fork) {
  RUN_SCOPE <<- SCOPE_PRESETS[[scope_tag]]

  est_dates <- resolve_window(fishery_name)
  win <- fishery_window_limit(fishery_name)
  if (!is.null(win)) {
    est_dates$est_date_start <- win$est_date_start
    est_dates$est_date_end   <- win$est_date_end
  }
  keep <- fishery_section_limit(fishery_name)

  eff <- fetch_fishery_dwg(fishery_name, est_dates)$effort |>
    mutate(event_date  = as.Date(event_date),
           section_num = suppressWarnings(as.double(section_num))) |>
    filter(between(event_date,
                   as.Date(est_dates$est_date_start),
                   as.Date(est_dates$est_date_end)))
  if (!is.null(keep)) eff <- eff |> filter(section_num %in% keep)
  eff$gear <- bank_or_boat(eff)

  is_census <- eff$tie_in_indicator %in% c(1, TRUE, "TRUE", "true")

  cen <- eff[is_census, , drop = FALSE] |>
    group_by(section_num, event_date) |>
    summarise(census_anglers_all  = sum(count_quantity, na.rm = TRUE),
              census_anglers_boat = sum(count_quantity[gear == "boat"], na.rm = TRUE),
              .groups = "drop")

  idx <- eff[!is_census, , drop = FALSE] |>
    group_by(section_num, event_date) |>
    summarise(index_vehicles = sum(count_quantity[count_type == "Vehicle Only"], na.rm = TRUE),
              index_trailers = sum(count_quantity[count_type == "Trailers Only"], na.rm = TRUE),
              .groups = "drop")

  # A paired day: the section carries BOTH counts. Census alone is dropped by
  # the model, index alone has no anchor to be measured against.
  inner_join(cen, idx, by = c("section_num", "event_date")) |>
    mutate(scope = scope_tag, fishery_name = fishery_name, year = year, fork = fork,
           .before = 1)
}

cli_h1("13 -- daily census-to-index ratios")
.orig_scope <- RUN_SCOPE
daily <- pmap_dfr(list(FITS$scope, FITS$fishery_name, FITS$year, FITS$fork), ratios_one)
RUN_SCOPE <- .orig_scope

# ------------------------------------------------------------------------------
# Stable reach identity, and the two ratio channels
# ------------------------------------------------------------------------------
# SECTION NUMBERS ARE NOT STABLE BETWEEN YEARS. The same mainstem water is
# sections 2 and 3 in 2024 and 1 and 2 in 2025, so colouring or faceting by
# section_num would compare different reaches across years. Numbering runs
# upstream, so rank within a fork-year recovers the correspondence.
daily <- daily |>
  group_by(fork, year) |>
  mutate(reach = {
    r <- dense_rank(section_num)
    if (max(r) == 1) "Whole fork"
    else if (max(r) == 2) c("Lower", "Upper")[r]
    else paste("Reach", r)          # never indexes past a two-element vector
  }) |>
  ungroup() |>
  mutate(
    # Denominator of zero is "not measured", not a ratio of infinity.
    ratio_vehicle = census_anglers_all  / na_if(index_vehicles, 0),
    ratio_trailer = census_anglers_boat / na_if(index_trailers, 0),
    # One axis for two years: the anchor dates are the same calendar dates in
    # both years, so a common-year date aligns them without faking a trend.
    season_date = as.Date(format(event_date, "2000-%m-%d"))
  )

write_csv(daily, file.path(OUT_DIR, "bss_b_census_index_ratio_daily.csv"))
cli_alert_success("Wrote bss_b_census_index_ratio_daily.csv ({nrow(daily)} paired section-days).")

# ------------------------------------------------------------------------------
# Figures -- one per index channel, so either can be shared on its own
# ------------------------------------------------------------------------------
# Two hues, validated for colour-vision deficiency (worst adjacent pair dE 24.7
# protan, 33.6 normal vision). Year carries colour; reach carries shape and
# line type, so identity is never colour alone.
YEAR_COLS <- c("2024" = "#2a78d6", "2025" = "#eb6834")

ratio_plot <- function(d, ycol, ylab) {
  d <- d |> filter(!is.na(.data[[ycol]]))
  if (nrow(d) == 0) return(NULL)
  ggplot(d, aes(x = season_date, y = .data[[ycol]],
                colour = factor(year), shape = reach, linetype = reach,
                group = interaction(year, reach))) +
    geom_hline(yintercept = 1, linetype = "dotted", colour = "grey55") +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.6) +
    facet_wrap(~fork, nrow = 1, scales = "free_y") +
    scale_colour_manual(values = YEAR_COLS, name = NULL) +
    scale_shape_manual(values = c(16, 17, 15), name = NULL) +
    scale_linetype_manual(values = c("solid", "longdash", "dotted"), name = NULL) +
    scale_x_date(date_labels = "%b %d") +
    labs(x = NULL, y = ylab) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom",
          strip.background = element_rect(fill = "white", colour = "black"))
}

walk2(
  c("ratio_vehicle", "ratio_trailer"),
  c("Census anglers per vehicle counted", "Census boat anglers per trailer counted"),
  function(ycol, ylab) {
    p <- ratio_plot(daily, ycol, ylab)
    channel <- str_remove(ycol, "ratio_")
    if (is.null(p)) {
      cli_alert_warning("No {channel} ratios to plot -- every denominator was zero.")
      return(invisible(NULL))
    }
    f <- file.path(FIG_DIR, paste0("fig21_census_index_ratio_", channel, ".png"))
    ggsave(f, p, width = 8, height = 4, dpi = 300, bg = "white")
    cli_alert_success("Wrote {.file {basename(f)}}.")
  }
)

cli_alert_info("Ratios are census anglers per index OBJECT, not `b` -- they differ \\
                by the anglers-per-vehicle expansion, which is constant within a \\
                year. Compare shapes across years, not levels against `b`.")
