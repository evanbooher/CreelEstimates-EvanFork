# ==============================================================================
# 06_variability_analysis.R
#
# Purpose:
#   The FIRST-ORDER QUESTION for the DFW/tribal technical-staff meeting:
#   how variable are the BSS effort-bias terms across years, rivers and
#   fisheries -- and what does that variability imply for predicting `b` in a
#   fishery-year where we cannot estimate it (no census / tie-in counts)?
#
#   Per direction: extract the bias terms BY LIKELIHOOD TYPE (b[1] vehicle,
#   b[2] trailer) BY RIVER AND YEAR, assess interannual variability, and look
#   for a pink-year (odd/even) effect in the point estimates in hand.
#
#   The conceptual model and a known-truth demonstration of the predictive
#   approach live in 05_predictive_toy.R -- read that first if the ladder
#   (random effects -> pink meta-regression -> multilevel) is unfamiliar.
#   This script applies the same ladder to the real estimates.
#
# Usage:
#   Rscript analysis/bss_bias/06_variability_analysis.R
#   Requires bss_b_summary.csv (01) and bss_b_comparability.csv (02).
#   No DB, no VPN, no Stan -- reads the pipeline's CSV outputs only.
#
# ------------------------------------------------------------------------------
# SCOPE DECISIONS (agreed, and worth restating in the room)
#
# [V1] Variance components are computed over ALL fishery-years, not only those
#      the comparability table calls `comparable`/`reference`. That is the
#      agreed default. It means the variability reported here INCLUDES
#      variation driven by changing season windows and spatial coverage, not
#      only changes in true vehicle-index bias. The comparability tier and the
#      season duration/timing columns from 02 are how you see which years
#      carry that caveat. A tier-restricted sensitivity run is one filter away
#      (see INCLUDE_TIERS below).
#
# [V2] Everything is modelled on the LOG scale. `b` is strictly positive with
#      a lognormal(0, sigma) prior, so log(b) is the natural symmetric scale
#      and effects read as proportional changes. SEs come via the delta
#      method, se_log ~= sd / median.
#
# [V3] LOW VARIANCE IN b[2] MAY MEAN "NO DATA", NOT "STABLE". Trailer
#      estimates are frequently prior-dominated -- a fishery-year with no
#      trailer index counts still yields a b[2] posterior, it just reproduces
#      the prior. A tight tau^2 across such years is the prior talking, not
#      evidence of a stable bias term. Every table below therefore reports
#      n_informed alongside n_years, and T2 carries an informed-only
#      sensitivity column. Read those before quoting a b[2] number.
#
# [V4] The PINK-YEAR analysis (T3, Figure 7) uses only the FALL-TIMED
#      fisheries. Pinks return in the fall of odd years, so a fishery has to
#      be fishing while they are in the river for a pink year to move its
#      effort. Skagit spring Chinook (Apr-Jul) and summer sockeye (Jun-Jul)
#      have no in-season overlap and are excluded -- including them would not
#      test the hypothesis, only dilute it. T1 and T2 still cover every
#      fishery. See PINK_RELEVANT_TARGETS below.
#
# Outputs (analysis/bss_bias/outputs/ and outputs/figures/):
#   bss_b_T1_inventory.csv / .html        -- every bias term: river x fishery x year x likelihood type
#   bss_b_T2_variability.csv / .html      -- per-series tau^2, I^2, Q, and the PREDICTION interval
#   bss_b_T3_pink_effect.csv / .html      -- odd/even contrast + multilevel variance decomposition
#   fig7_pink_year_effect.png/.pdf        -- b by year with pink years marked
#   fig8_predictive_distribution.png/.pdf -- predictive distribution for an unobserved year vs. observed years
# ==============================================================================

library(tidyverse)
library(gt)
library(cli)
library(here)

if (!requireNamespace("metafor", quietly = TRUE)) {
  cli::cli_abort(c(
    "Package {.pkg metafor} is required.",
    "i" = "Install with: {.code install.packages(\"metafor\")}"
  ))
}

OUT_DIR <- here::here("analysis", "bss_bias", "outputs")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

source(here::here("analysis", "bss_bias", "common.R"))

summary_path <- file.path(OUT_DIR, "bss_b_summary.csv")
comp_path    <- file.path(OUT_DIR, "bss_b_comparability.csv")
if (!file.exists(summary_path)) cli::cli_abort("{.file {summary_path}} not found -- run 01_fit_bss_bias.R first.")
if (!file.exists(comp_path))    cli::cli_abort("{.file {comp_path}} not found -- run 02_build_comparability_table.R first.")

# NULL = all tiers (the [V1] default). Set to e.g.
# c("reference", "comparable", "comparable-with-caveat") for a sensitivity run.
INCLUDE_TIERS <- NULL

# Below this prior contraction, the posterior has barely moved off its prior
# and the estimate is carrying little data -- see [V3]. Matches the threshold
# 01_fit_bss_bias.R uses for informed_flag == "prior-dominated".
INFORMED_MIN_CONTRACTION <- 0.10

b_summary <- read_csv(summary_path, show_col_types = FALSE)
comp      <- read_csv(comp_path, show_col_types = FALSE)

# A bss_b_comparability.csv written by an earlier version of 02 lacks
# fishery_type (and the duration/timing columns). Say so plainly rather than
# failing several steps later with "object 'fishery_type' not found".
comp_needed <- c("fishery_type", "n_days_in_window", "start_shift_days", "pct_of_ref_days")
comp_absent <- setdiff(comp_needed, names(comp))
if (length(comp_absent) > 0) {
  cli::cli_abort(c(
    "{.file {comp_path}} is missing column{?s} {.val {comp_absent}}.",
    "x" = "It was written by an earlier version of 02_build_comparability_table.R.",
    "i" = "Re-run {.file analysis/bss_bias/02_build_comparability_table.R}, then re-run this script."
  ))
}

# ------------------------------------------------------------------------------
# Analysis frame
# ------------------------------------------------------------------------------

dat <- b_summary |>
  left_join(
    comp |> select(fishery_name, basin, fishery_label, fishery_type, year_start,
                   comparability_tier, n_days_in_window, start_shift_days, pct_of_ref_days),
    by = "fishery_name"
  ) |>
  mutate(
    parity   = year_parity_label(year_start),
    is_pink  = as.integer(year_start %% 2 == 1),
    log_b    = log(median),
    se_log   = sd / median,           # delta method, see [V2]
    vi       = se_log^2,
    informed = !is.na(prior_contraction) & prior_contraction >= INFORMED_MIN_CONTRACTION &
                 (is.na(informed_flag) | informed_flag != "unconverged")
  ) |>
  filter(!is.na(log_b), is.finite(log_b), !is.na(vi), is.finite(vi), vi > 0)

if (!is.null(INCLUDE_TIERS)) {
  n_before <- nrow(dat)
  dat <- dat |> filter(comparability_tier %in% INCLUDE_TIERS)
  cli::cli_alert_info("INCLUDE_TIERS set -- {nrow(dat)} of {n_before} estimate(s) retained.")
}

if (nrow(dat) == 0) cli::cli_abort("No usable bias estimates after filtering -- is bss_b_summary.csv populated?")

# Unconverged fits are excluded from the MODELS (they are not estimates of
# anything) but stay visible in the T1 inventory so nothing silently vanishes.
model_dat <- dat |> filter(is.na(informed_flag) | informed_flag != "unconverged")

# ------------------------------------------------------------------------------
# [V4] PINK-YEAR ANALYSIS IS RESTRICTED TO FALL-TIMED FISHERIES.
#
# Puget Sound pinks return in the fall of odd years. A fishery has to be
# fishing while they are in the river for a pink year to change its effort,
# and therefore its vehicle:angler relationship. That holds for the fall
# fisheries and not for the others:
#
#   Skagit fall salmon         Aug-Nov   overlaps the pink return    IN
#   Snohomish fall salmon      Aug-Nov   overlaps                    IN
#   Stillaguamish salmon+gf    Sep-Nov   overlaps                    IN
#   Skagit spring Chinook      Apr-Jul   no overlap                  OUT
#   Skagit summer sockeye      Jun-Jul   no overlap                  OUT
#
# Including the spring and summer fisheries would not test the hypothesis --
# there is no in-season angler basis for pinks to affect them -- it would only
# add years whose odd/even split is noise, diluting the contrast and making a
# real effect harder to see.
#
# Selected via target species, which is the same key the catch groups use:
# the Coho-target fisheries ARE the fall-timed ones.
# ------------------------------------------------------------------------------

PINK_RELEVANT_TARGETS <- "Coho"

pink_dat <- model_dat |>
  mutate(target_species = target_species_from_est_cg(est_cg)) |>
  filter(target_species %in% PINK_RELEVANT_TARGETS)

cli::cli_alert_info(
  "Pink-year analysis restricted to fall-timed ({.val {PINK_RELEVANT_TARGETS}}-target) fisheries: \\
   {n_distinct(pink_dat$fishery_type)} series, {nrow(pink_dat)} estimate(s) \\
   (of {nrow(model_dat)} overall) -- see [V4]."
)

cli::cli_alert_info(
  "{nrow(dat)} bias estimate(s) across {n_distinct(dat$fishery_type)} fishery series, \\
   {n_distinct(dat$basin)} basin(s), {n_distinct(dat$year_start)} year(s); \\
   {sum(!dat$informed)} weakly-informed (see [V3])."
)

# ------------------------------------------------------------------------------
# T1 -- the inventory. The literal ask: bias terms by likelihood type, by
# river and year.
# ------------------------------------------------------------------------------

T1 <- dat |>
  transmute(
    basin, fishery_type, fishery_name, year_start, parity,
    likelihood = bias_type,                       # "vehicle" = b[1], "trailer" = b[2]
    b_median = median, b_sd = sd,
    q10, q90, q2.5, q97.5,
    prior_contraction, informed_flag, informed,
    rhat, ess_bulk,
    comparability_tier, n_days_in_window
  ) |>
  arrange(basin, fishery_type, likelihood, year_start)

write_csv(T1, file.path(OUT_DIR, "bss_b_T1_inventory.csv"))

T1 |>
  mutate(across(c(b_median, b_sd, q10, q90, q2.5, q97.5, prior_contraction), ~ round(.x, 3))) |>
  select(basin, fishery_type, year_start, parity, likelihood, b_median,
         q10, q90, prior_contraction, informed_flag, comparability_tier) |>
  gt(groupname_col = "basin") |>
  tab_header(title = "T1. Effort-bias terms by likelihood type, river and year",
             subtitle = "b[1] = vehicle-count bias, b[2] = trailer-count bias. q10-q90 is the 80% credible interval.") |>
  cols_label(fishery_type = "Fishery", year_start = "Year", parity = "Year type",
             likelihood = "Likelihood", b_median = "b (median)",
             prior_contraction = "Prior contraction", informed_flag = "Flag",
             comparability_tier = "Tier") |>
  tab_footnote(
    footnote = "Prior contraction near 0 means the posterior barely moved off its prior -- the estimate carries little data. See [V3] in the script header.",
    locations = cells_column_labels(columns = prior_contraction)
  ) |>
  opt_row_striping() |>
  gt::gtsave(file.path(OUT_DIR, "bss_b_T1_inventory.html"))

cli::cli_alert_success("T1 written ({nrow(T1)} rows).")

# ------------------------------------------------------------------------------
# T2 -- per-series interannual variability, and the predictive distribution
# for an unobserved year.
#
# tau^2 is the interannual variance AFTER accounting for each estimate's own
# uncertainty -- i.e. real year-to-year movement in b, not noise in measuring
# it. That distinction is the whole point of using a random-effects fit rather
# than the SD of the point estimates.
#
# pi_lb / pi_ub are the PREDICTION interval: where a NEW year's b is expected
# to fall. It is wider than the confidence interval on the series mean, and it
# is the interval that answers "what should we use for a year we cannot
# measure". Reporting the CI in its place would understate uncertainty.
# ------------------------------------------------------------------------------

fit_series <- function(df, label) {
  n_years <- nrow(df)
  base <- tibble(
    n_years      = n_years,
    n_informed   = sum(df$informed),
    mean_b       = exp(mean(df$log_b)),
    sd_log_b_raw = if (n_years > 1) sd(df$log_b) else NA_real_
  )
  if (n_years < 2) {
    return(bind_cols(base, tibble(
      tau2 = NA_real_, tau = NA_real_, I2 = NA_real_, Q_p = NA_real_,
      pooled_b = exp(df$log_b[1]), ci_lb = NA_real_, ci_ub = NA_real_,
      pi_lb = NA_real_, pi_ub = NA_real_, pi_width_ratio = NA_real_,
      tau2_informed_only = NA_real_,
      note = "single year -- no interannual variance estimable"
    )))
  }
  f <- try(metafor::rma(yi = log_b, vi = vi, data = df, method = "REML"), silent = TRUE)
  if (inherits(f, "try-error")) {
    return(bind_cols(base, tibble(
      tau2 = NA_real_, tau = NA_real_, I2 = NA_real_, Q_p = NA_real_,
      pooled_b = NA_real_, ci_lb = NA_real_, ci_ub = NA_real_,
      pi_lb = NA_real_, pi_ub = NA_real_, pi_width_ratio = NA_real_,
      tau2_informed_only = NA_real_,
      note = "model did not converge"
    )))
  }
  p <- predict(f)

  # [V3] sensitivity: the same tau^2 over informed estimates only. If this is
  # much larger than the all-estimates tau^2, the apparent stability was
  # prior-dominated years pulling everything toward b = 1.
  inf_df <- df |> filter(informed)
  tau2_inf <- if (nrow(inf_df) >= 2) {
    fi <- try(metafor::rma(yi = log_b, vi = vi, data = inf_df, method = "REML"), silent = TRUE)
    if (inherits(fi, "try-error")) NA_real_ else fi$tau2
  } else NA_real_

  bind_cols(base, tibble(
    tau2 = f$tau2, tau = sqrt(f$tau2), I2 = f$I2, Q_p = f$QEp,
    pooled_b = exp(as.numeric(f$b)),
    ci_lb = exp(p$ci.lb), ci_ub = exp(p$ci.ub),
    pi_lb = exp(p$pi.lb), pi_ub = exp(p$pi.ub),
    # How much wider the predictive interval is than the CI on the mean --
    # the cost, in uncertainty, of predicting a new year rather than
    # describing past ones.
    pi_width_ratio = (p$pi.ub - p$pi.lb) / (p$ci.ub - p$ci.lb),
    tau2_informed_only = tau2_inf,
    note = NA_character_
  ))
}

T2 <- model_dat |>
  group_by(basin, fishery_type, bias_type) |>
  group_modify(~ fit_series(.x, .y)) |>
  ungroup() |>
  arrange(bias_type, basin, fishery_type)

write_csv(T2, file.path(OUT_DIR, "bss_b_T2_variability.csv"))

T2 |>
  mutate(across(c(mean_b, tau, pooled_b, ci_lb, ci_ub, pi_lb, pi_ub, pi_width_ratio,
                  tau2, tau2_informed_only), ~ round(.x, 3)),
         I2 = round(I2, 1), Q_p = signif(Q_p, 3)) |>
  select(basin, fishery_type, bias_type, n_years, n_informed, pooled_b,
         ci_lb, ci_ub, pi_lb, pi_ub, tau, I2, Q_p, tau2_informed_only) |>
  gt(groupname_col = "bias_type") |>
  tab_header(title = "T2. Interannual variability in b, and the predictive distribution for a new year",
             subtitle = "tau = SD of true year-to-year variation, after removing each estimate's own uncertainty.") |>
  tab_spanner(label = "Series mean (95% CI)", columns = c(pooled_b, ci_lb, ci_ub)) |>
  tab_spanner(label = "Prediction interval for a NEW year", columns = c(pi_lb, pi_ub)) |>
  cols_label(fishery_type = "Fishery", bias_type = "Likelihood", n_years = "N yrs",
             n_informed = "N informed", pooled_b = "b", ci_lb = "lo", ci_ub = "hi",
             pi_lb = "lo", pi_ub = "hi", Q_p = "Q p",
             tau2_informed_only = "tau² (informed only)") |>
  tab_footnote(
    footnote = "The PREDICTION interval, not the CI -- it is where a new year's b is expected to fall, and is deliberately wider than the interval on the mean of past years.",
    locations = cells_column_spanners(spanners = "Prediction interval for a NEW year")
  ) |>
  tab_footnote(
    footnote = "If tau² (informed only) is much larger than tau², the apparent stability came from prior-dominated years sitting near b = 1, not from a genuinely stable bias term.",
    locations = cells_column_labels(columns = tau2_informed_only)
  ) |>
  opt_row_striping() |>
  gt::gtsave(file.path(OUT_DIR, "bss_b_T2_variability.html"))

cli::cli_alert_success("T2 written ({nrow(T2)} series).")

# ------------------------------------------------------------------------------
# T3 -- pink-year effect and the variance decomposition across rivers/years.
#
# One multilevel fit per likelihood type:
#   log b ~ pink + (1 | basin/fishery_type) + (1 | observation)
# The moderator gives the pink effect; the random-effect variances give
# Thomas's "across years and rivers" decomposition. Fitting with and without
# the moderator shows how much of the interannual variance pink explains.
# ------------------------------------------------------------------------------

fit_decomposition <- function(df) {
  df <- df |> mutate(obs_id = row_number())
  rand <- list(~ 1 | basin, ~ 1 | fishery_type, ~ 1 | obs_id)

  m0 <- try(metafor::rma.mv(yi = log_b, V = vi, random = rand, data = df, method = "REML"), silent = TRUE)
  m1 <- try(metafor::rma.mv(yi = log_b, V = vi, mods = ~ is_pink, random = rand, data = df, method = "REML"), silent = TRUE)
  if (inherits(m0, "try-error") || inherits(m1, "try-error")) return(NULL)

  # sigma2 follows the order of `rand`
  tibble(
    n_estimates      = nrow(df),
    n_series         = n_distinct(df$fishery_type),
    n_basins         = n_distinct(df$basin),
    var_basin        = m0$sigma2[1],
    var_fishery      = m0$sigma2[2],
    var_year         = m0$sigma2[3],
    var_measurement  = mean(df$vi),
    pink_log_effect  = as.numeric(m1$b[2]),
    pink_ratio       = exp(as.numeric(m1$b[2])),   # multiplicative effect on b
    pink_ci_lo       = exp(m1$ci.lb[2]),
    pink_ci_hi       = exp(m1$ci.ub[2]),
    pink_p           = m1$pval[2],
    var_year_no_pink   = m0$sigma2[3],
    var_year_with_pink = m1$sigma2[3],
    pct_year_var_explained = 100 * (1 - m1$sigma2[3] / m0$sigma2[3])
  )
}

T3 <- pink_dat |>
  group_by(bias_type) |>
  group_modify(~ fit_decomposition(.x) %||% tibble()) |>
  ungroup()

if (nrow(T3) == 0) {
  cli::cli_alert_warning("Multilevel decomposition did not converge for any likelihood type -- skipping T3.")
} else {
  write_csv(T3, file.path(OUT_DIR, "bss_b_T3_pink_effect.csv"))

  T3 |>
    mutate(across(c(var_basin, var_fishery, var_year, var_measurement,
                    pink_ratio, pink_ci_lo, pink_ci_hi), ~ round(.x, 4)),
           pink_p = signif(pink_p, 3),
           pct_year_var_explained = round(pct_year_var_explained, 1)) |>
    select(bias_type, n_estimates, n_series, n_basins,
           var_basin, var_fishery, var_year, var_measurement,
           pink_ratio, pink_ci_lo, pink_ci_hi, pink_p, pct_year_var_explained) |>
    gt() |>
    tab_header(title = "T3. Where the variability lives, and the pink-year effect",
               subtitle = "Variance components on the log scale; pink effect as a multiplicative factor on b.") |>
    tab_spanner(label = "Variance components", columns = starts_with("var_")) |>
    tab_spanner(label = "Pink (odd) year effect", columns = c(pink_ratio, pink_ci_lo, pink_ci_hi, pink_p)) |>
    cols_label(bias_type = "Likelihood", n_estimates = "N est", n_series = "N series", n_basins = "N basins",
               var_basin = "Between rivers", var_fishery = "Between fisheries",
               var_year = "Between years", var_measurement = "Estimation error",
               pink_ratio = "x b", pink_ci_lo = "lo", pink_ci_hi = "hi", pink_p = "p",
               pct_year_var_explained = "% of year variance explained by pink") |>
    tab_footnote(
      footnote = "A pink ratio of 1.0 means no effect. 1.25 would mean pink years run 25% higher.",
      locations = cells_column_labels(columns = pink_ratio)
    ) |>
    opt_row_striping() |>
    gt::gtsave(file.path(OUT_DIR, "bss_b_T3_pink_effect.html"))

  cli::cli_alert_success("T3 written.")
  cli::cli_h2("Variance decomposition")
  print(T3 |> select(bias_type, var_basin, var_fishery, var_year, var_measurement,
                     pink_ratio, pink_ci_lo, pink_ci_hi, pct_year_var_explained))
}

# ------------------------------------------------------------------------------
# Figure 7 -- the pink-year picture
# ------------------------------------------------------------------------------

make_pink_fig <- function(bt, param_symbol) {
  d <- pink_dat |> filter(bias_type == bt)
  if (nrow(d) == 0) return(NULL)

  parity_means <- d |> group_by(fishery_type, parity) |>
    summarise(mean_b = exp(mean(log_b)), .groups = "drop")

  ggplot(d, aes(x = year_start, y = median)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = INK_MUTED, linewidth = 0.4) +
    geom_hline(data = parity_means, aes(yintercept = mean_b, color = parity),
               linetype = "dotted", linewidth = 0.5) +
    geom_linerange(aes(ymin = q10, ymax = q90, color = parity), linewidth = 1.2, alpha = 0.8) +
    geom_point(aes(color = parity, shape = informed), size = 3) +
    scale_color_manual(values = PARITY_COLORS, name = "Year type") +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), name = "Informed by data",
                       labels = c("TRUE" = "yes", "FALSE" = "weakly / prior-dominated")) +
    facet_wrap(~ fishery_type, scales = "free_y") +
    labs(
      title = paste0(param_symbol, ": is there a pink-year (odd) effect?"),
      subtitle = "Bars = 80% credible interval. Dotted lines = mean of odd vs. even years within each fishery.\nPuget Sound pink salmon return in odd years.",
      x = NULL, y = param_symbol,
      caption = "Hollow points are weakly informed -- their position reflects the prior more than the data, so they should not drive the read."
    ) +
    theme_bss()
}

fig7 <- make_pink_fig("vehicle", "b₁ (vehicle bias)")
if (!is.null(fig7)) save_fig(fig7, "fig7_pink_year_effect", width = 10, height = 7)

fig7b <- make_pink_fig("trailer", "b₂ (trailer bias)")
if (!is.null(fig7b)) save_fig(fig7b, "fig7b_pink_year_effect_trailer", width = 10, height = 7)

# ------------------------------------------------------------------------------
# Figure 8 -- the money figure: what would we use for a year we cannot measure?
#
# Observed years, with the predictive distribution for an unobserved year
# drawn behind them as a band. If the band is tight, historical b transfers
# well; if it is wide, it does not, and that is the answer to the meeting's
# first-order question.
# ------------------------------------------------------------------------------

make_prediction_fig <- function(bt, param_symbol) {
  obs <- model_dat |> filter(bias_type == bt)
  pred <- T2 |> filter(bias_type == bt, !is.na(pi_lb))
  if (nrow(obs) == 0 || nrow(pred) == 0) return(NULL)

  ggplot() +
    geom_rect(data = pred,
              aes(xmin = -Inf, xmax = Inf, ymin = pi_lb, ymax = pi_ub),
              fill = CAT[["blue"]], alpha = 0.13) +
    geom_hline(data = pred, aes(yintercept = pooled_b), color = CAT[["blue"]], linewidth = 0.5) +
    geom_hline(yintercept = 1, linetype = "dashed", color = INK_MUTED, linewidth = 0.4) +
    geom_linerange(data = obs, aes(x = year_start, ymin = q10, ymax = q90),
                   color = INK_SECOND, linewidth = 1.1, alpha = 0.75) +
    geom_point(data = obs, aes(x = year_start, y = median, fill = parity),
               shape = 21, size = 3, color = SURFACE, stroke = 0.8) +
    scale_fill_manual(values = PARITY_COLORS, name = "Year type") +
    facet_wrap(~ fishery_type, scales = "free_y") +
    labs(
      title = paste0(param_symbol, ": predictive distribution for a year we cannot measure"),
      subtitle = "Shaded band = 95% prediction interval for a NEW year, from that fishery's own history.\nLine = series mean. Points = observed years with their 80% intervals.",
      x = NULL, y = param_symbol,
      caption = paste0(
        "The band is the answer to \"what do we use when there is no census data this year, and how uncertain is it?\" ",
        "It is wider than a confidence interval on the mean by construction. Single-year series have no band -- nothing to predict from."
      )
    ) +
    theme_bss()
}

fig8 <- make_prediction_fig("vehicle", "b₁ (vehicle bias)")
if (!is.null(fig8)) save_fig(fig8, "fig8_predictive_distribution", width = 10, height = 7)

fig8b <- make_prediction_fig("trailer", "b₂ (trailer bias)")
if (!is.null(fig8b)) save_fig(fig8b, "fig8b_predictive_distribution_trailer", width = 10, height = 7)

cli::cli_alert_success("Figures 7/7b/8/8b written to {.file {FIG_DIR}}.")
cli::cli_alert_info("Tables: bss_b_T1_inventory, bss_b_T2_variability, bss_b_T3_pink_effect (.csv and .html).")
