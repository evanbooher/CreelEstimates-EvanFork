# ==============================================================================
# 03_plot_b_series.R
#
# Purpose:
#   The top-priority deliverable for tomorrow's meeting: visuals that make
#   the b[1] (vehicle-count bias) temporal-stability story legible on their
#   own. Six figures, built per the plan in analysis/bss_bias/README.md.
#
# On the `dataviz` skill: that skill's procedure (form -> color-by-job ->
# validate -> marks -> hover -> accessibility) targets interactive HTML/CSS
# charts. This deliverable is static R/ggplot2 output (PNG/PDF) for a
# printed/projected meeting deck, so the HOVER/INTERACTION and CSS-token
# steps don't apply -- but the TRANSFERABLE principles are followed
# throughout and called out inline where used:
#   - form chosen by the data's job (magnitude-over-time -> point + interval,
#     not a bar chart; distribution -> density, not a boxplot)
#   - color assigned BY JOB, not decoration: categorical (fishery identity),
#     sequential (year, one hue light->dark), status (comparability tier) --
#     never a rainbow, never color standing in for more than one job at once
#   - the skill's validated categorical palette and status palette (hex
#     values below, from references/palette.md) used AS-IS rather than
#     picked by eye, since they're already verified colorblind-safe in
#     their documented order
#   - identity is never color-alone: shape/linetype double-encode
#     comparability tier alongside its status color
#   - one axis; a legend for every series (color/shape symbology, not direct
#     text labels -- tried end-labels for Skagit's <=4 concurrent fisheries
#     first, but real data has series close enough in year/value that
#     ggrepel-placed labels became illegible; a legend scales better);
#     recessive gridlines/axes; thin lines
#
# Usage:
#   Rscript analysis/bss_bias/03_plot_b_series.R
#   Requires bss_b_summary.csv (01) and bss_b_comparability.csv (02) to exist.
#
# Outputs (analysis/bss_bias/outputs/figures/):
#   fig1_b_vehicle_by_year.png/.pdf   -- the hero figure
#   fig1b_b_trailer_by_year.png/.pdf  -- same chart for b[2] (trailer bias); secondary, not in the composite
#   fig2_comparability_strip.png/.pdf
#   fig1_fig2_composite.png/.pdf      -- Figure 1 + 2 stacked, aligned x-axis (the slide)
#   fig3_posterior_densities.png/.pdf         -- b[1] (vehicle)
#   fig3b_posterior_densities_trailer.png/.pdf -- b[2] (trailer)
#   fig4_vehicle_vs_trailer_bias.png/.pdf
#   fig5_information_diagnostic.png/.pdf
# ==============================================================================

library(tidyverse)
library(patchwork)
library(here)
library(cli)

OUT_DIR <- here::here("analysis", "bss_bias", "outputs")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

summary_path <- file.path(OUT_DIR, "bss_b_summary.csv")
comp_path    <- file.path(OUT_DIR, "bss_b_comparability.csv")
if (!file.exists(summary_path)) cli::cli_abort("{.file {summary_path}} not found -- run 01_fit_bss_bias.R first.")
if (!file.exists(comp_path))    cli::cli_abort("{.file {comp_path}} not found -- run 02_build_comparability_table.R first.")

b_summary <- read_csv(summary_path, show_col_types = FALSE)
comp      <- read_csv(comp_path, show_col_types = FALSE)

# ------------------------------------------------------------------------------
# Palette (from the dataviz skill's references/palette.md, used verbatim --
# this is the pre-validated categorical order, not a freehand choice)
# ------------------------------------------------------------------------------

CAT <- c(blue = "#2a78d6", orange = "#eb6834", aqua = "#1baf7a", yellow = "#eda100",
         magenta = "#e87ba4", green = "#008300", violet = "#4a3aa7", red = "#e34948")

STATUS <- c(good = "#0ca30c", warning = "#fab219", serious = "#ec835a", critical = "#d03b3b")

INK          <- "#0b0b0b"
INK_SECOND   <- "#52514e"
INK_MUTED    <- "#898781"
GRID_COLOR   <- "#e1e0d9"
BASELINE_COL <- "#c3c2b7"
SURFACE      <- "#fcfcfb"

TIER_COLORS <- c(
  "reference"              = INK_MUTED,
  "comparable"              = STATUS[["good"]],
  "comparable-with-caveat"   = STATUS[["warning"]],
  "not-comparable"            = STATUS[["critical"]],
  "not-estimable"               = INK_MUTED
)
TIER_SHAPES <- c(
  "reference"              = 16,  # filled circle
  "comparable"              = 16,  # filled circle
  "comparable-with-caveat"   = 21,  # filled circle w/ ring (fill + stroke)
  "not-comparable"            = 1,   # hollow circle
  "not-estimable"               = 4    # x
)

theme_bss <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.major = element_line(color = GRID_COLOR, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = BASELINE_COL, linewidth = 0.4),
      axis.text = element_text(color = INK_SECOND),
      axis.title = element_text(color = INK_SECOND),
      plot.title = element_text(color = INK, face = "bold"),
      plot.subtitle = element_text(color = INK_SECOND),
      plot.caption = element_text(color = INK_MUTED, size = 8, hjust = 0),
      strip.text = element_text(color = INK, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(color = INK_SECOND),
      legend.text = element_text(color = INK_SECOND),
      plot.background = element_rect(fill = SURFACE, color = NA),
      panel.background = element_rect(fill = SURFACE, color = NA)
    )
}

save_fig <- function(plot, name, width = 9, height = 6) {
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), plot, width = width, height = height, dpi = 300, bg = SURFACE)
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), plot, width = width, height = height, device = cairo_pdf, bg = SURFACE)
}

# ------------------------------------------------------------------------------
# Assemble a plotting frame for ONE bias_type ("vehicle" or "trailer"),
# joined to comparability -- shared by Figure 1 (vehicle) and Figure 1b
# (trailer) below.
# ------------------------------------------------------------------------------

build_bias_plot_df <- function(bias_type_val) {
  b_summary |>
    filter(bias_type == bias_type_val) |>
    left_join(
      comp |> select(fishery_name, basin, fishery_label, season_label, year_start, comparability_tier),
      by = "fishery_name"
    ) |>
    mutate(
      comparability_tier = replace_na(comparability_tier, "not-estimable"),
      basin = factor(basin, levels = c("Skagit", "Snohomish", "Stillaguamish"))
    )
}

plot_df <- build_bias_plot_df("vehicle")
if (nrow(plot_df) == 0) {
  cli::cli_abort("No vehicle-bias rows to plot -- bss_b_summary.csv may be empty (no fits completed yet).")
}

# Colors assigned from EVERY fishery in b_summary (not just vehicle rows), so
# a given fishery is the same color in both the vehicle and trailer figures
# -- a reader matching a series across the two plots shouldn't have to
# re-learn the color key. facet_wrap(~basin) already separates basins, so
# Snohomish/Stillaguamish's single series each still reads unambiguously
# despite not sharing one "neutral" color the way Skagit's four don't either.
all_series    <- b_summary |> left_join(comp |> select(fishery_name, fishery_label), by = "fishery_name") |>
  distinct(fishery_label) |> pull(fishery_label) |> sort()
series_colors <- setNames(unname(CAT)[((seq_along(all_series) - 1) %% length(CAT)) + 1], all_series)

# ------------------------------------------------------------------------------
# Figure 1 / 1b -- b[1] (vehicle) and b[2] (trailer) vs year, by basin,
# comparability-tier-encoded. Same chart for both parameters so they read as
# a matched pair -- only the bias_type filter, axis label, and title differ.
# ------------------------------------------------------------------------------

make_b_series_fig <- function(plot_df, param_label, y_lab) {
  ggplot(plot_df, aes(x = year_start, y = median, group = fishery_label, color = fishery_label)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = INK_MUTED, linewidth = 0.4) +
    geom_segment(aes(xend = year_start, y = q2.5, yend = q97.5), linewidth = 0.5, alpha = 0.55, lineend = "round") +
    geom_segment(aes(xend = year_start, y = q10, yend = q90), linewidth = 1.4, lineend = "round") +
    geom_line(linewidth = 0.6, alpha = 0.7) +
    geom_point(aes(shape = comparability_tier, fill = fishery_label), size = 3, stroke = 1) +
    scale_color_manual(values = series_colors, name = "Fishery") +
    scale_fill_manual(values = series_colors, guide = "none") +
    scale_shape_manual(values = TIER_SHAPES, name = "Comparability tier") +
    facet_wrap(~ basin, ncol = 1, scales = "free_y") +
    labs(
      title = paste0("BSS ", param_label, " across years"),
      subtitle = "Thick segment = 80% credible interval (q10–q90); thin segment = 95% (q2.5–q97.5). Dashed line = b=1 (\"no bias\").",
      x = NULL, y = y_lab,
      caption = paste0(
        "Model: ", unique(b_summary$model_file)[1], " | fit config: ", unique(b_summary$fit_config)[1],
        " | prior: lognormal(0, ", unique(b_summary$prior_sigma)[1], ") | est_cg: fixed per-fishery-name target ",
        "(Chinook/Coho/Sockeye harvest -- see README.md's Catch-group selection).\n",
        "Marker shape encodes comparability tier (see Figure 2 / bss_b_comparability.csv) -- open circle = not-comparable, x = not-estimable."
      )
    ) +
    theme_bss() +
    theme(legend.position = "top")
}

fig1 <- make_b_series_fig(plot_df, "vehicle-count bias (b₁)", expression(b[1]~"(vehicle bias)"))
save_fig(fig1, "fig1_b_vehicle_by_year", width = 9, height = 9)

plot_df_trailer <- build_bias_plot_df("trailer")
if (nrow(plot_df_trailer) == 0) {
  cli::cli_alert_warning("No trailer-bias rows to plot -- skipping Figure 1b.")
} else {
  fig1_trailer <- make_b_series_fig(plot_df_trailer, "trailer-count bias (b₂)", expression(b[2]~"(trailer bias)"))
  save_fig(fig1_trailer, "fig1b_b_trailer_by_year", width = 9, height = 9)
}

# ------------------------------------------------------------------------------
# Figure 2 -- comparability strip, x-aligned under Figure 1
# ------------------------------------------------------------------------------

strip_df <- comp |>
  mutate(basin = factor(basin, levels = c("Skagit", "Snohomish", "Stillaguamish")))

fig2 <- ggplot(strip_df, aes(x = year_start, y = fishery_label, fill = comparability_tier)) +
  geom_tile(color = SURFACE, linewidth = 1) +
  scale_fill_manual(values = TIER_COLORS, name = "Comparability tier") +
  facet_grid(basin ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Year", y = NULL, title = "Comparability tier by fishery-year") +
  theme_bss() +
  theme(panel.grid = element_blank(), strip.text.y = element_text(angle = 0))

save_fig(fig2, "fig2_comparability_strip", width = 9, height = 4)

# Composite: Figure 1 over Figure 2, aligned x-axis -- the actual slide.
composite <- fig1 / fig2 + plot_layout(heights = c(3, 1))
save_fig(composite, "fig1_fig2_composite", width = 9, height = 12)

# ------------------------------------------------------------------------------
# Figure 3 / 3b -- posterior densities OVERLAID per fishery TYPE (name with
# the year stripped out, so all years of e.g. "Skagit fall salmon" share one
# panel), alpha transparency, sequential year ramp (one hue, light->dark;
# most recent year emphasized). Year needs a real legend -- unlike a
# ridge-per-row layout, color/fill is the ONLY thing distinguishing years
# once they're overlaid in shared (x, density) space. b_draws/<safe_name>.rds
# holds BOTH b[1] and b[2] columns (saved together in 01_fit_bss_bias.R), so
# one parameterized builder produces both figures rather than only ever
# reading b[1].
# ------------------------------------------------------------------------------

make_posterior_density_fig <- function(draws_files, param_col, param_symbol, x_lab) {
  draws_long <- map_dfr(draws_files, function(f) {
    fn <- tools::file_path_sans_ext(basename(f))
    d <- readRDS(f)
    if (!(param_col %in% names(d))) return(NULL)
    tibble(fishery_name_safe = fn, val = d[[param_col]])
  })

  # match safe_name() back to fishery_name via bss_b_summary's own safe-name mapping
  name_map <- b_summary |> distinct(fishery_name) |>
    mutate(fishery_name_safe = stringr::str_replace_all(fishery_name, "[^[:alnum:]]", "_"))

  draws_long <- draws_long |>
    left_join(name_map, by = "fishery_name_safe") |>
    left_join(comp |> select(fishery_name, basin, fishery_label, year_start), by = "fishery_name")

  if (nrow(draws_long) == 0) return(NULL)

  # The year token sits in different positions across naming conventions
  # ("...salmon 2022" vs. "...Chinook 2024 upper"), so replace (not strip)
  # it with a space and squish, rather than assuming it's a trailing token.
  draws_long <- draws_long |>
    mutate(fishery_type = stringr::str_squish(
      stringr::str_replace(fishery_name, "\\d{4}(-\\d{2,4})?", " ")
    ))

  seq_ramp <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
  year_levels <- sort(unique(draws_long$year_start))
  ramp_fun <- colorRampPalette(seq_ramp)
  year_colors <- setNames(ramp_fun(length(year_levels)), year_levels)

  ggplot(draws_long, aes(x = val, fill = factor(year_start), color = factor(year_start))) +
    geom_density(alpha = 0.45, linewidth = 0.4) +
    geom_vline(xintercept = 1, linetype = "dashed", color = INK_MUTED, linewidth = 0.4) +
    scale_fill_manual(values = year_colors, name = "Year") +
    scale_color_manual(values = year_colors, guide = "none") +
    facet_wrap(~ fishery_type, scales = "free") +
    labs(title = paste0(param_symbol, " posterior densities, overlaid by year within each fishery"),
         subtitle = "Darker = more recent year. Each panel is one fishery across all its years.",
         x = x_lab, y = "Density") +
    theme_bss() +
    theme(legend.position = "top")
}

draws_files <- list.files(file.path(OUT_DIR, "b_draws"), pattern = "\\.rds$", full.names = TRUE)
if (length(draws_files) > 0) {
  fig3 <- make_posterior_density_fig(draws_files, "b[1]", "b₁", expression(b[1]))
  if (!is.null(fig3)) {
    save_fig(fig3, "fig3_posterior_densities", width = 10, height = 8)
  } else {
    cli::cli_alert_warning("No draws matched fishery names for Figure 3 (vehicle) -- skipping.")
  }

  fig3_trailer <- make_posterior_density_fig(draws_files, "b[2]", "b₂", expression(b[2]))
  if (!is.null(fig3_trailer)) {
    save_fig(fig3_trailer, "fig3b_posterior_densities_trailer", width = 10, height = 8)
  } else {
    cli::cli_alert_warning("No draws matched fishery names for Figure 3b (trailer) -- skipping.")
  }
} else {
  cli::cli_alert_warning("No files in outputs/b_draws/ yet -- skipping Figures 3/3b (run after some fits complete).")
}

# ------------------------------------------------------------------------------
# Figure 4 -- b[1] (vehicle) vs b[2] (trailer)
# ------------------------------------------------------------------------------

wide_bias <- b_summary |>
  select(fishery_name, bias_type, median) |>
  pivot_wider(names_from = bias_type, values_from = median) |>
  left_join(comp |> select(fishery_name, basin, comparability_tier), by = "fishery_name") |>
  filter(!is.na(vehicle), !is.na(trailer)) |>
  mutate(comparability_tier = replace_na(comparability_tier, "not-estimable"),
         basin = factor(basin, levels = c("Skagit", "Snohomish", "Stillaguamish")))

if (nrow(wide_bias) > 0) {
  lims <- range(c(wide_bias$vehicle, wide_bias$trailer, 1), na.rm = TRUE)
  fig4 <- ggplot(wide_bias, aes(x = vehicle, y = trailer, color = basin, shape = comparability_tier)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = INK_MUTED) +
    geom_hline(yintercept = 1, linetype = "dashed", color = INK_MUTED, linewidth = 0.3) +
    geom_vline(xintercept = 1, linetype = "dashed", color = INK_MUTED, linewidth = 0.3) +
    geom_point(size = 3) +
    scale_color_manual(values = c(Skagit = CAT[["blue"]], Snohomish = CAT[["orange"]], Stillaguamish = CAT[["aqua"]])) +
    scale_shape_manual(values = TIER_SHAPES) +
    coord_equal(xlim = lims, ylim = lims) +
    labs(title = "Vehicle bias (b₁) vs. trailer bias (b₂)",
         subtitle = "Dotted line = 1:1; dashed lines = b=1. Points off the diagonal: the two count types disagree.",
         x = expression(b[1]~"(vehicle)"), y = expression(b[2]~"(trailer)"), color = "Basin", shape = "Comparability") +
    theme_bss()
  save_fig(fig4, "fig4_vehicle_vs_trailer_bias", width = 7, height = 7)
}

# ------------------------------------------------------------------------------
# Figure 5 -- information diagnostic: prior_contraction vs b[1] median
# (the honest-broker figure: low-information years cluster near the prior)
# ------------------------------------------------------------------------------

info_df <- b_summary |> filter(bias_type == "vehicle") |>
  left_join(comp |> select(fishery_name, basin), by = "fishery_name") |>
  mutate(basin = factor(basin, levels = c("Skagit", "Snohomish", "Stillaguamish")))

if (nrow(info_df) > 0) {
  fig5 <- ggplot(info_df, aes(x = prior_contraction, y = median, color = basin)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = INK_MUTED, linewidth = 0.3) +
    geom_vline(xintercept = 0.10, linetype = "dotted", color = STATUS[["warning"]], linewidth = 0.4) +
    geom_point(aes(size = ess_bulk), alpha = 0.85) +
    scale_color_manual(values = c(Skagit = CAT[["blue"]], Snohomish = CAT[["orange"]], Stillaguamish = CAT[["aqua"]])) +
    labs(
      title = "Information diagnostic: how much did the data move b₁ from its prior?",
      subtitle = "prior_contraction near 0 = posterior ≈ prior (b is not informed by this fishery-year's data). Dotted line = the 0.10 flag threshold used upstream.",
      x = "Prior contraction (1 − posterior var / prior var)", y = expression(b[1]~"median"),
      color = "Basin", size = "ESS (bulk)"
    ) +
    theme_bss()
  save_fig(fig5, "fig5_information_diagnostic", width = 8, height = 6)
}

# ------------------------------------------------------------------------------
# Figure 6 -- Figure 1's series with the four candidate options drawn as
# reference bands, so the choice is visible AGAINST the data it summarizes.
# Optional: only runs if 04_candidate_options.R has already produced its
# output (re-run this script after 04 to get this figure).
# ------------------------------------------------------------------------------

options_path <- file.path(OUT_DIR, "bss_b_candidate_options.csv")
if (file.exists(options_path)) {
  # NOTE: opts already carries its own `basin` column (from 04's
  # group_by(basin, fishery_label, bias_type)) -- do NOT re-join `basin` from
  # `comp` here, that would produce basin.x/basin.y and silently break
  # facet_wrap(~ basin)'s per-panel routing of these hlines (inherited from fig1).
  opts <- read_csv(options_path, show_col_types = FALSE) |>
    filter(bias_type == "vehicle") |>
    mutate(basin = factor(basin, levels = c("Skagit", "Snohomish", "Stillaguamish"))) |>
    filter(!is.na(estimate))

  option_labels <- c(time_series_mean = "Mean", most_recent_year = "Most recent",
                      precision_weighted = "Precision-weighted", hierarchical_shrinkage = "Hierarchical")
  option_shapes  <- c(time_series_mean = "solid", most_recent_year = "dotted",
                       precision_weighted = "dashed", hierarchical_shrinkage = "dotdash")

  fig6 <- fig1 +
    geom_hline(
      data = opts, aes(yintercept = estimate, linetype = option),
      color = INK_SECOND, linewidth = 0.5, inherit.aes = FALSE
    ) +
    scale_linetype_manual(values = option_shapes, labels = option_labels, name = "Candidate option") +
    labs(title = "b₁ across years, with candidate bias-correction options overlaid")

  save_fig(fig6, "fig6_candidate_options_overlay", width = 9, height = 9)
  cli::cli_alert_success("Figure 6 written (candidate options overlaid).")
} else {
  cli::cli_alert_info("bss_b_candidate_options.csv not found yet -- run 04_candidate_options.R then re-run this script for Figure 6.")
}

cli::cli_alert_success("Figures written to analysis/bss_bias/outputs/figures/")
cli::cli_alert_info("Open fig1_fig2_composite.png first -- that's the slide.")
