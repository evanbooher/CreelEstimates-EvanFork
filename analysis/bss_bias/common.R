# ==============================================================================
# common.R -- shared palette, theme, figure saving, and label derivations
#
# Sourced by 03_plot_b_series.R, 05_predictive_toy.R, 06_variability_analysis.R
# and the Quarto report, so a colour or theme change lands everywhere at once.
# Previously copy-pasted per script, which is exactly how a "consistent" set of
# figures quietly stops matching.
#
# Palette values are the dataviz skill's validated categorical and status
# palettes, used as-is rather than picked by eye -- they are already verified
# colourblind-safe in their documented order.
# ==============================================================================

CAT <- c(blue = "#2a78d6", orange = "#eb6834", aqua = "#1baf7a", yellow = "#eda100",
         magenta = "#e87ba4", green = "#008300", violet = "#4a3aa7", red = "#e34948")

STATUS <- c(good = "#0ca30c", warning = "#fab219", serious = "#ec835a", critical = "#d03b3b")

INK          <- "#0b0b0b"
INK_SECOND   <- "#52514e"
INK_MUTED    <- "#898781"
GRID_COLOR   <- "#e1e0d9"
BASELINE_COL <- "#c3c2b7"
SURFACE      <- "#fcfcfb"

# Sequential ramp for year (light -> dark, one hue): "more recent = darker".
SEQ_RAMP <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")

# Pink (odd) vs non-pink (even) year. Deliberately NOT from the sequential
# ramp -- parity is a categorical distinction, not a magnitude.
PARITY_COLORS <- c("even" = CAT[["blue"]], "odd (pink)" = CAT[["magenta"]])

theme_bss <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = GRID_COLOR, linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = BASELINE_COL, linewidth = 0.4),
      axis.text = ggplot2::element_text(color = INK_SECOND),
      axis.title = ggplot2::element_text(color = INK_SECOND),
      plot.title = ggplot2::element_text(color = INK, face = "bold"),
      plot.subtitle = ggplot2::element_text(color = INK_SECOND),
      plot.caption = ggplot2::element_text(color = INK_MUTED, size = 8, hjust = 0),
      strip.text = ggplot2::element_text(color = INK, face = "bold"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(color = INK_SECOND),
      legend.text = ggplot2::element_text(color = INK_SECOND),
      plot.background = ggplot2::element_rect(fill = SURFACE, color = NA),
      panel.background = ggplot2::element_rect(fill = SURFACE, color = NA)
    )
}

# Writes PNG (for slides/inline) and PDF (vector, for print). FIG_DIR is
# expected to be defined by the calling script.
save_fig <- function(plot, name, width = 9, height = 6) {
  ggplot2::ggsave(file.path(FIG_DIR, paste0(name, ".png")), plot,
                  width = width, height = height, dpi = 300, bg = SURFACE)
  ggplot2::ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), plot,
                  width = width, height = height, device = cairo_pdf, bg = SURFACE)
}

# Fishery TYPE = fishery name with the year token removed, e.g.
#   "Skagit fall salmon 2022"        -> "Skagit fall salmon"
#   "Skagit spring Chinook 2024 upper" -> "Skagit spring Chinook upper"
#   "Stillaguamish salmon and gamefish 2022-23" -> "Stillaguamish salmon and gamefish"
# This is the unit a b-series belongs to: all years of one named fishery.
#
# REPLACES the year token with a space rather than stripping it, because the
# year sits mid-string for the upper/lower Skagit fisheries -- stripping it
# would run the surrounding words together ("...Chinookupper").
fishery_type_from_name <- function(x) {
  stringr::str_squish(stringr::str_replace(x, "\\d{4}(-\\d{2,4})?", " "))
}

# Puget Sound pink salmon return in ODD years. year_start is the correct
# anchor even for multi-season names: "Stillaguamish salmon and gamefish
# 2022-23" is a fall-2022 fishery, so it takes 2022's parity (confirmed with
# the fishery lead).
year_parity_label <- function(year_start) {
  dplyr::if_else(year_start %% 2 == 1, "odd (pink)", "even")
}

# ------------------------------------------------------------------------------
# Series colours, keyed on TARGET SPECIES rather than one arbitrary hue per
# series.
#
# Six unrelated hues for six series reads as a rainbow and makes a reader learn
# six arbitrary pairings. Colouring by target species instead means a reader
# learns three -- and the grouping is the analytically meaningful one, because
# the target species IS the catch group each fit ran against
# (fishery_target_catch_group() in 01: Chinook / Coho / Sockeye). So "orange =
# a Coho-target fishery" holds across every basin panel.
#
# Where one basin has several series on the same target (Skagit spring Chinook
# lower and upper), they take light/dark shades of that target's hue -- related
# fisheries look related, which is the whole point.
# ------------------------------------------------------------------------------

TARGET_COLORS <- c(
  Chinook = CAT[["blue"]],
  Coho    = CAT[["orange"]],
  Sockeye = CAT[["violet"]],
  Other   = INK_MUTED
)

# Blend a hex colour toward white (amount > 0) or black (amount < 0).
# Base R only -- not worth a colorspace dependency for this.
shade_color <- function(hex, amount) {
  v <- grDevices::col2rgb(hex)[, 1]
  target <- if (amount >= 0) c(255, 255, 255) else c(0, 0, 0)
  out <- v + (target - v) * abs(amount)
  grDevices::rgb(out[1], out[2], out[3], maxColorValue = 255)
}

# est_cg is built as species_lifestage_finmark_fate, so the target species is
# everything before the first underscore.
target_species_from_est_cg <- function(est_cg) {
  sp <- stringr::str_extract(as.character(est_cg), "^[^_]+")
  dplyr::if_else(sp %in% names(TARGET_COLORS), sp, "Other")
}

# Named colour vector over fishery_type, shaded within each target species.
# `df` needs fishery_type and target_species columns.
series_palette_by_target <- function(df) {
  key <- df |>
    dplyr::distinct(fishery_type, target_species) |>
    dplyr::filter(!is.na(fishery_type)) |>
    dplyr::arrange(target_species, fishery_type)

  shaded <- purrr::map_dfr(split(key, key$target_species), function(g) {
    base <- TARGET_COLORS[[g$target_species[1]]]
    n <- nrow(g)
    # One series on this target keeps the base hue; several spread light->dark
    # so they stay recognisably the same species.
    amounts <- if (n == 1) 0 else seq(0.35, -0.25, length.out = n)
    g$color <- vapply(amounts, function(a) shade_color(base, a), character(1))
    g
  })

  stats::setNames(shaded$color, shaded$fishery_type)
}
