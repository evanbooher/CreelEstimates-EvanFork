# ==============================================================================
# 05_predictive_toy.R
#
# Purpose:
#   A CONCEPTUAL MODEL and TOY EXAMPLE for predicting the BSS effort-bias term
#   `b` in a fishery-year where we CANNOT estimate it from that year's own
#   data (no census / tie-in counts). Built for the DFW/tribal technical-staff
#   meeting, per direction to bring a conceptual model and worked toy example
#   rather than only empirical results.
#
#   Runs standalone in seconds. NO database, NO VPN, NO Stan, NO dependency on
#   any pipeline output -- it simulates a world where the truth is KNOWN, so
#   the method can be demonstrated and interrogated before anyone argues about
#   the real numbers. That property is deliberate: it lets tomorrow's
#   conversation be about the approach.
#
# ------------------------------------------------------------------------------
# THE PROBLEM
#
#   For a fishery-year with census data, BSS estimates `b` directly. For a
#   fishery-year WITHOUT it, `b` is unidentified and the posterior just
#   reproduces its lognormal(0, sigma) prior -- i.e. "no bias, wide
#   uncertainty", which is an assumption, not an estimate. The question is
#   what to use instead, and how much uncertainty it should carry.
#
#   Everything below works on the LOG scale. `b` is strictly positive and
#   lognormal-priored, so log(b) is the natural symmetric scale and
#   differences on it read as proportional changes.
#
# ------------------------------------------------------------------------------
# THE LADDER (each rung nests inside the next)
#
#   Rung 0 -- "use last year's b". No pooling, no uncertainty propagation.
#             The status quo, included as the comparison baseline.
#
#   Rung 1 -- RANDOM EFFECTS / PARTIAL POOLING, per fishery series:
#               log b_hat[i,y] ~ Normal(log theta[i,y], se[i,y]^2)   <- what we measured
#               log theta[i,y] ~ Normal(mu[i], tau^2)                <- true year-to-year variation
#             tau^2 IS the "how variable are the bias terms" answer. And the
#             PREDICTIVE DISTRIBUTION for an unmeasured year is exactly the
#             random-effects PREDICTION interval:
#               log theta[i,new] ~ Normal(mu_hat[i], tau_hat^2 + SE(mu_hat[i])^2)
#             Note this is wider than the confidence interval on mu -- it has
#             to be, because it describes a NEW year, not the mean of past
#             ones. Conflating the two understates uncertainty, which is the
#             most common way this kind of analysis goes wrong.
#
#   Rung 2 -- META-REGRESSION with a pink-year covariate:
#               log theta[i,y] ~ Normal(mu[i] + beta * pink[y], tau^2)
#             Puget Sound pink salmon return in ODD years, bringing large
#             effort influxes that plausibly shift the vehicle:angler
#             relationship. beta estimates that effect; the DROP in tau^2 from
#             rung 1 to rung 2 measures how much of the interannual
#             variability pink year explains. Prediction for a new year then
#             conditions on whether that year is a pink year.
#
#   Rung 3 -- POOLING ACROSS RIVERS, so a data-poor fishery borrows strength:
#               mu[i] ~ Normal(M, sigma_river^2)
#             Gives the across-river vs. across-year variance decomposition,
#             and yields a usable predictive distribution even for a fishery
#             with little or no history of its own.
#
#   Rung 4 -- NOT BUILT HERE. Re-specify `b` hierarchically INSIDE the Stan
#             model so posterior uncertainty propagates exactly instead of
#             being approximated in two stages. This is the README's deferred
#             "(d2)". Rungs 1-3 are a two-stage plug-in approximation to it:
#             they treat each se[i,y] as fixed and known, which is very
#             standard practice but is an approximation worth stating out loud.
#
# Usage:
#   Rscript analysis/bss_bias/05_predictive_toy.R
#
# Outputs (analysis/bss_bias/outputs/figures/ and outputs/):
#   figA_toy_simulated_world.png/.pdf   -- the known-truth world we simulated
#   figB_toy_prediction_check.png/.pdf  -- held-out year vs. its predictive distribution
#   figC_toy_variance_components.png/.pdf -- where the variability lives; tau^2 before/after pink
#   bss_b_toy_recovery.csv              -- true vs. recovered parameter values
#
# Objects left in the global environment for the report to reuse after
# source()-ing this file: toy_sim, toy_fits, toy_recovery, figA/figB/figC.
# ==============================================================================

library(tidyverse)
library(cli)
library(here)

if (!requireNamespace("metafor", quietly = TRUE)) {
  cli::cli_abort(c(
    "Package {.pkg metafor} is required for the predictive-model toy example.",
    "i" = "Install with: {.code install.packages(\"metafor\")}",
    "i" = "It provides rma()/rma.mv(), which give multilevel variance components \\
           and prediction intervals directly -- hand-rolling those is error-prone."
  ))
}

OUT_DIR <- here::here("analysis", "bss_bias", "outputs")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Palette / theme shared with 03_plot_b_series.R. Kept in sync by hand for now
# (both are small); if a third script needs them, factor them out properly.
CAT <- c(blue = "#2a78d6", orange = "#eb6834", aqua = "#1baf7a", yellow = "#eda100",
         magenta = "#e87ba4", green = "#008300", violet = "#4a3aa7", red = "#e34948")
INK          <- "#0b0b0b"
INK_SECOND   <- "#52514e"
INK_MUTED    <- "#898781"
GRID_COLOR   <- "#e1e0d9"
BASELINE_COL <- "#c3c2b7"
SURFACE      <- "#fcfcfb"

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
# 1. Simulate a world where the truth is known
#
# Truth values are chosen to resemble what the real fits have been producing
# (b roughly 0.6-2.0, i.e. log b roughly -0.5 to +0.7, with per-estimate
# posterior SEs on log scale in the 0.05-0.20 range). They are NOT calibrated
# to the real data -- this is a demonstration of method, not a result.
# ------------------------------------------------------------------------------

TOY_TRUTH <- list(
  M            = log(1.10),  # grand mean log-bias across rivers
  sigma_river  = 0.20,       # river-to-river SD in mean log-bias
  tau          = 0.15,       # year-to-year SD WITHIN a river, after pink is accounted for
  beta_pink    = 0.25,       # pink (odd) year effect on log-bias; exp(0.25) ~ +28%
  se_range     = c(0.05, 0.20)
)

simulate_toy_world <- function(rivers = c("River A", "River B", "River C"),
                               years = 2019:2025,
                               truth = TOY_TRUTH,
                               seed = 42) {
  set.seed(seed)

  river_means <- tibble(
    river  = rivers,
    mu_i   = rnorm(length(rivers), truth$M, truth$sigma_river)
  )

  tidyr::expand_grid(river = rivers, year = years) |>
    left_join(river_means, by = "river") |>
    mutate(
      pink = as.integer(year %% 2 == 1),  # Puget Sound pinks return in ODD years
      # true year-specific log-bias
      log_theta = rnorm(n(), mu_i + truth$beta_pink * pink, truth$tau),
      # what a BSS fit for that year would report: a noisy read on log_theta
      se        = runif(n(), truth$se_range[1], truth$se_range[2]),
      log_b_hat = rnorm(n(), log_theta, se),
      vi        = se^2,
      b_hat     = exp(log_b_hat),
      theta     = exp(log_theta)
    )
}

toy_sim <- simulate_toy_world()

# ------------------------------------------------------------------------------
# 2. Climb the ladder
# ------------------------------------------------------------------------------

# Rung 1: per-river random effects. tau^2 = interannual variability.
fit_rung1 <- function(dat) {
  dat |>
    group_by(river) |>
    group_modify(~ {
      f <- metafor::rma(yi = log_b_hat, vi = vi, data = .x, method = "REML")
      p <- predict(f)
      tibble(
        mu_hat = as.numeric(f$b), se_mu = f$se, tau2 = f$tau2,
        I2 = f$I2,
        # CONFIDENCE interval on the mean vs. PREDICTION interval for a new
        # year -- the second is the one that answers "what should we use next
        # year", and it is meaningfully wider.
        ci_lb = p$ci.lb, ci_ub = p$ci.ub,
        pi_lb = p$pi.lb, pi_ub = p$pi.ub
      )
    }) |>
    ungroup()
}

# Rung 2: add the pink-year covariate. Compare tau^2 against rung 1.
fit_rung2 <- function(dat) {
  dat |>
    group_by(river) |>
    group_modify(~ {
      f <- metafor::rma(yi = log_b_hat, vi = vi, mods = ~ pink, data = .x, method = "REML")
      tibble(
        intercept   = as.numeric(f$b[1]),
        beta_pink   = as.numeric(f$b[2]),
        beta_se     = f$se[2],
        beta_ci_lb  = f$ci.lb[2],
        beta_ci_ub  = f$ci.ub[2],
        tau2_resid  = f$tau2
      )
    }) |>
    ungroup()
}

# Rung 3: multilevel -- pool across rivers, with a pink moderator. Gives the
# across-river vs. across-year decomposition in one fit.
fit_rung3 <- function(dat) {
  dat <- dat |> mutate(obs_id = row_number())
  metafor::rma.mv(
    yi = log_b_hat, V = vi,
    mods = ~ pink,
    random = list(~ 1 | river, ~ 1 | obs_id),
    data = dat, method = "REML"
  )
}

toy_fits <- list(
  rung1 = fit_rung1(toy_sim),
  rung2 = fit_rung2(toy_sim),
  rung3 = fit_rung3(toy_sim)
)

# ------------------------------------------------------------------------------
# 3. Did it recover the truth? (the point of a toy example)
# ------------------------------------------------------------------------------

r3 <- toy_fits$rung3
toy_recovery <- tibble::tribble(
  ~parameter,                          ~truth,                    ~recovered,
  "Grand mean log-bias (M)",           TOY_TRUTH$M,               as.numeric(r3$b[1]),
  "Pink-year effect (beta)",           TOY_TRUTH$beta_pink,       as.numeric(r3$b[2]),
  "Between-river SD (sigma_river)",    TOY_TRUTH$sigma_river,     sqrt(r3$sigma2[1]),
  "Within-river year SD (tau)",        TOY_TRUTH$tau,             sqrt(r3$sigma2[2])
) |>
  mutate(
    difference = recovered - truth,
    across(c(truth, recovered, difference), ~ round(.x, 3))
  )

write_csv(toy_recovery, file.path(OUT_DIR, "bss_b_toy_recovery.csv"))

cli::cli_h2("Toy example: parameter recovery")
print(toy_recovery)

# Mean tau^2 with and without the pink covariate -- how much of the
# interannual variability does pink year explain?
tau2_no_pink   <- mean(toy_fits$rung1$tau2)
tau2_with_pink <- mean(toy_fits$rung2$tau2_resid)
cli::cli_alert_info(
  "Mean tau^2 without pink covariate: {round(tau2_no_pink, 4)}; \\
   with pink: {round(tau2_with_pink, 4)} \\
   ({round(100 * (1 - tau2_with_pink / tau2_no_pink))}% of interannual variance explained)."
)

# ------------------------------------------------------------------------------
# 4. Held-out check: predict a year we pretended not to observe
#
# This is the whole use case in miniature -- a fishery-year with no census
# data. Drop the most recent year, fit on the rest, and ask whether the
# predictive distribution actually covers the value we hid.
# ------------------------------------------------------------------------------

holdout_year <- max(toy_sim$year)
train <- toy_sim |> filter(year != holdout_year)
test  <- toy_sim |> filter(year == holdout_year)

holdout_check <- train |>
  group_by(river) |>
  group_modify(~ {
    f <- metafor::rma(yi = .x$log_b_hat, vi = .x$vi, mods = ~ pink, data = .x, method = "REML")
    # predict for the held-out year's pink status
    pink_new <- test$pink[test$river == .y$river][1]
    p <- predict(f, newmods = pink_new)
    tibble(pred = p$pred, pi_lb = p$pi.lb, pi_ub = p$pi.ub)
  }) |>
  ungroup() |>
  left_join(test |> select(river, log_theta, log_b_hat), by = "river") |>
  mutate(
    covered = log_theta >= pi_lb & log_theta <= pi_ub,
    across(c(pred, pi_lb, pi_ub, log_theta, log_b_hat), ~ round(.x, 3))
  )

cli::cli_h2("Toy example: held-out year {holdout_year} vs. its predictive distribution")
print(holdout_check)

# ------------------------------------------------------------------------------
# 5. Figures
# ------------------------------------------------------------------------------

figA <- ggplot(toy_sim, aes(x = year, y = b_hat)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = INK_MUTED, linewidth = 0.4) +
  geom_segment(aes(xend = year, y = exp(log_b_hat - 1.96 * se), yend = exp(log_b_hat + 1.96 * se)),
               color = INK_MUTED, linewidth = 0.5, alpha = 0.6) +
  geom_line(aes(y = theta), color = CAT[["aqua"]], linewidth = 0.7, alpha = 0.9) +
  geom_point(aes(fill = factor(pink, levels = c(0, 1), labels = c("even", "odd (pink)"))),
             shape = 21, size = 3, color = SURFACE, stroke = 0.8) +
  scale_fill_manual(values = c("even" = CAT[["blue"]], "odd (pink)" = CAT[["magenta"]]), name = "Year type") +
  facet_wrap(~ river, ncol = 1) +
  labs(
    title = "Toy world: what we simulated",
    subtitle = "Green line = TRUE year-specific bias. Points = what a BSS fit would report, with its 95% interval.\nPink (odd) years sit systematically higher by construction.",
    x = NULL, y = "b (bias term)",
    caption = "Simulated data with known truth -- a demonstration of method, not a result. See TOY_TRUTH in 05_predictive_toy.R."
  ) +
  theme_bss()

save_fig(figA, "figA_toy_simulated_world", width = 9, height = 8)

figB <- holdout_check |>
  ggplot(aes(x = river)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = INK_MUTED, linewidth = 0.4) +
  geom_linerange(aes(ymin = pi_lb, ymax = pi_ub), color = CAT[["blue"]], linewidth = 6, alpha = 0.25) +
  geom_point(aes(y = pred), color = CAT[["blue"]], size = 3.5) +
  geom_point(aes(y = log_theta, shape = covered), color = CAT[["red"]], size = 3.5, stroke = 1.2) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 4),
                     name = paste0("Truth inside\nprediction interval")) +
  coord_flip() +
  labs(
    title = paste0("Held-out check: can we predict ", holdout_year, " without observing it?"),
    subtitle = "Blue band = 95% PREDICTIVE distribution from earlier years only. Red = the true value we hid.\nThis is the real use case: a fishery-year with no census data.",
    x = NULL, y = "log(b)",
    caption = "Prediction interval, not confidence interval -- it describes a NEW year, so it is deliberately wider than the CI on the mean."
  ) +
  theme_bss()

save_fig(figB, "figB_toy_prediction_check", width = 9, height = 5)

var_comp <- tibble(
  component = c("Between rivers", "Between years (within river)", "Within-year estimation error"),
  variance  = c(r3$sigma2[1], r3$sigma2[2], mean(toy_sim$vi))
) |>
  mutate(
    component = factor(component, levels = rev(component)),
    share = variance / sum(variance)
  )

figC <- ggplot(var_comp, aes(x = share, y = component, fill = component)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = scales::percent(share, accuracy = 1)),
            hjust = -0.15, color = INK_SECOND, size = 3.5) +
  scale_fill_manual(values = unname(CAT)[1:3], guide = "none") +
  scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Where does the variability live?",
    subtitle = paste0(
      "Adding the pink-year covariate cuts mean interannual variance from ",
      round(tau2_no_pink, 4), " to ", round(tau2_with_pink, 4),
      " (", round(100 * (1 - tau2_with_pink / tau2_no_pink)), "% explained)."
    ),
    x = "Share of total variance", y = NULL,
    caption = "This decomposition is the first-order question: if 'between years' dominates, last year's b is a poor stand-in for this year's."
  ) +
  theme_bss()

save_fig(figC, "figC_toy_variance_components", width = 9, height = 4.5)

cli::cli_alert_success("Toy example complete. Figures A-C written to {.file {FIG_DIR}}.")
cli::cli_alert_info("Recovery table: {.file {file.path(OUT_DIR, 'bss_b_toy_recovery.csv')}}")
