# ==============================================================================
# 07_catch_sensitivity.R
#
# Purpose:
#   How much does an error in the effort-index bias term `b` move the quantity
#   anyone actually cares about -- ESTIMATED CATCH?
#
#   The answer is fixed by the model's structure, not by simulation. In the BSS
#   Stan model:
#
#     line 205  c[a] ~ neg_binomial_2(lambda_C_S[...] * h[a], r_C)
#               CPUE is fit from interviews alone. No `b`, no effort term.
#
#     line 187  V_I[i] ~ poisson((lambda_E_S_I[...] * p_TI * R_V[1] + ...) * b[1])
#               `b` multiplies the Poisson mean for the index counts. With the
#               counts V_I fixed as data, lambda_E is proportional to 1/b.
#               R_V cannot absorb it -- interviews pin R_V separately (line 210).
#
#     line 240  lambda_Ctot_S[s][d,g] = lambda_E_S[s][d,g] * L[d] * lambda_C_S[s][d,g]
#               Catch = effort x trip length x CPUE.
#
#   So `b` reaches catch through exactly ONE channel (effort), and:
#
#         catch is proportional to 1 / b        (elasticity = -1)
#
#   b > 1 (index over-counts relative to anglers) revises effort and catch DOWN.
#   b < 1 revises them UP.
#
#   WHAT BREAKS THE -1: census. Line 201, `E_s[e] ~ poisson(lambda_E_S_I * p_TI)`,
#   carries no `b` and pins effort directly. Census is the ONLY thing in the
#   model that separates `b` from effort. In a census-free year the confound is
#   total and -1 is exact. See T4 below, which reports per fishery-year which
#   index channels were live and whether census was present, rather than
#   asserting -1 everywhere.
#
#   A COROLLARY WORTH STATING IN THE ROOM: the PERCENTAGE effect is identical
#   for every catch group. Chinook encounters and Coho harvest move by exactly
#   the same percent for a given error in `b`, because both are effort x CPUE
#   and only effort moves. Only the absolute number of fish differs. One ratio
#   analysis, N translations.
#
# Usage:
#   Rscript analysis/bss_bias/07_catch_sensitivity.R
#   No DB, no VPN, no Stan. Reads the pipeline's CSVs + saved b draws only.
#
# ------------------------------------------------------------------------------
# SCOPE DECISIONS
#
# [S1] THE BACKTEST IS THE HEADLINE, NOT THE LADDER. A ladder of round b values
#      answers "what if b were 2?" -- a question nobody is asking. The decision
#      on the table is "use a predicted b instead of measuring one", so T6 asks
#      exactly that, retrospectively: for every fishery-year that HAS a measured
#      b, what would catch have been had we imported the series prediction?
#
# [S2] LEAVE-ONE-OUT, ALWAYS. The pooled_b in T2 includes the year being tested.
#      Scoring a year against a mean that contains it flatters the method. Every
#      backtest row therefore refits the meta-analysis WITHOUT that year and
#      predicts it from the others. This is why T6's intervals are wider than
#      T2's, and the wider ones are the honest ones.
#
# [S3] THE BACKTEST IS AN UPPER BOUND ON ERROR. Elasticity -1 assumes the year
#      has no census. Historical years mostly DID have census, which would have
#      partially corrected an imported b. So T6 answers "what if we imported b
#      AND collected no census" -- which is precisely this year's proposal, but
#      it must be said aloud or someone will correctly object that the
#      historical years were not run that way.
#
# [S4] A MEASURED b IS NOT TRUTH, IT IS AN ESTIMATE. The ratio distribution
#      draws from the measured posterior (all saved draws) rather than using a
#      point, so the yardstick carries its own uncertainty. Prior-dominated and
#      unconverged years are excluded outright -- a b that reproduces its prior
#      is not something to score a prediction against.
#
# [S5] EXCHANGEABILITY IS THE WEAK LINK, AND T6 TESTS IT RATHER THAN ASSUMING
#      IT. The meta-analysis assumes a new year is drawn from the same
#      distribution as past years -- exactly what the comparability work
#      (shifting sections, moving season starts, pink parity) puts in doubt. The
#      calibration line printed at the end is the empirical check: if ~95% of
#      held-out years land inside their own 95% prediction interval, the
#      assumption is earning its keep. Materially less, and it is failing.
#
# Outputs (analysis/bss_bias/outputs/ and outputs/figures/):
#   bss_b_T4_pass_through.csv   -- per fishery-year: live index channels, census presence, pass-through class
#   bss_b_T5_ladder.csv         -- deterministic b -> catch multiplier, empirical tiers + round ladder
#   bss_b_T6_loo_backtest.csv   -- LOO backtest: ratio distribution per fishery-year x bias type
#   bss_b_T6_calibration.csv    -- per-series PI coverage, the [S5] check
#   bss_b_T7_direct_sensitivity.csv -- vary b on a real fitted dataset, read off
#                                      the catch: multiplier, % change, and fish
#                                      where a C_sum baseline exists
#   fig17_catch_ladder.png/.pdf
#   fig18_loo_backtest.png/.pdf
#   fig19_exposure.png/.pdf        -- how far catch could move per series, ranked
# ==============================================================================

library(tidyverse)
library(gt)
library(cli)
library(here)

if (!requireNamespace("metafor", quietly = TRUE)) {
  cli::cli_abort(c("Package {.pkg metafor} is required.",
                   "i" = "Install with: {.code install.packages(\"metafor\")}"))
}

OUT_DIR    <- here::here("analysis", "bss_bias", "outputs")
FIG_DIR    <- file.path(OUT_DIR, "figures")
DRAWS_DIR  <- file.path(OUT_DIR, "b_draws")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

source(here::here("analysis", "bss_bias", "common.R"))

# Same definition as fishery_data.R's. Copied rather than sourced: that file
# pulls in creelutils and the fishery-params lookup, and this script is meant to
# run with no DB, no VPN and no pipeline dependencies beyond the CSVs.
safe_name <- function(x) stringr::str_replace_all(x, "[^[:alnum:]]", "_")

set.seed(20260914)     # the ratio distributions are sampled; keep runs reproducible
N_SIM <- 8000          # draws per fishery-year x bias type for the ratio distribution

# Matches 01_fit_bss_bias.R's informed_flag threshold and 06's [V3].
INFORMED_MIN_CONTRACTION <- 0.10

# A leave-one-out fit needs enough years LEFT OVER to estimate tau^2 at all.
# Holding one out of 3 leaves 2, which is the bare minimum metafor will take.
LOO_MIN_SERIES_YEARS <- 3

# The round ladder, for intuition only -- the empirical tiers are the result.
# 0.1 is deliberately absent: it implies a 10x catch revision and is far
# outside anything observed (the lowest PI bound across all 12 series is 0.17),
# so it would dominate any figure it appeared in and invite an argument about a
# value we have never seen. Added back by editing this line if asked.
ROUND_LADDER <- c(0.5, 0.8, 1.0, 1.2, 1.5, 2.0)

# ------------------------------------------------------------------------------
# Inputs
# ------------------------------------------------------------------------------

req <- function(path, who) {
  if (!file.exists(path)) cli::cli_abort("{.file {path}} not found -- run {.file {who}} first.")
  read_csv(path, show_col_types = FALSE)
}

b_summary <- req(file.path(OUT_DIR, "bss_b_summary.csv"),      "01_fit_bss_bias.R")
comp      <- req(file.path(OUT_DIR, "bss_b_comparability.csv"), "02_build_comparability_table.R")
dims      <- req(file.path(OUT_DIR, "bss_b_stan_dims.csv"),     "01_fit_bss_bias.R")

# Written by 01b_fit_catch_groups.R. Optional on purpose: without it every
# result below is still complete, just expressed as a percentage rather than in
# fish. Do not make this a hard dependency -- the MCMC that produces it takes
# hours and the analysis must not wait on it.
# Globbed, not a single path: 01b runs as several concurrent processes and each
# writes its own tagged file (see OUTPUT_TAG in 01), because append_csv_row() is
# a read-modify-write that concurrent runs corrupt. distinct() on the key
# absorbs the one overlap that can occur -- the same fishery-year fitted under
# the same catch group by two jobs whose fishery filters overlap.
catch_files <- list.files(OUT_DIR, pattern = "^bss_catch_baseline.*\\.csv$", full.names = TRUE)
catch_base <- if (length(catch_files) > 0) {
  map_dfr(catch_files, ~read_csv(.x, show_col_types = FALSE, col_types = cols(.default = col_guess()))) |>
    distinct(fishery_name, est_cg, .keep_all = TRUE)
} else NULL
if (length(catch_files) > 1) {
  cli::cli_alert_info("Merged {length(catch_files)} catch-baseline files from parallel 01b runs.")
}
if (is.null(catch_base)) {
  cli::cli_alert_info(
    "No {.file bss_catch_baseline.csv} -- results will be in percent only. \\
     Run {.file 09_read_production_estimates.R} to pull season totals from the \\
     production pipeline's fishery_analyses/ output."
  )
}

# A fishery-year with ZERO records for a catch group is a RESULT, not a gap.
# 01 skips those at its "no matching records" guard, so they never reach the
# baseline file -- and in a Chinook-impact table an absent row reads as
# "not examined" when it means "examined, and the answer was none". That
# distinction matters most for exactly the impact-limited groups where a zero
# is the finding. Carried over from 00d's inventory with catch = 0.
inv_path <- file.path(OUT_DIR, "bss_catch_inventory.csv")
if (file.exists(inv_path) && !is.null(catch_base)) {
  zero_rows <- read_csv(inv_path, show_col_types = FALSE) |>
    filter(!is.na(catch_group), !is.na(n_records), n_records == 0) |>
    transmute(fishery_name, est_cg, C_sum_median = 0,
              baseline_source = "zero encounters (no records)")
  if (nrow(zero_rows) > 0) {
    catch_base <- catch_base |>
      mutate(baseline_source = "fitted") |>
      bind_rows(anti_join(zero_rows, catch_base, by = c("fishery_name", "est_cg")))
    cli::cli_alert_info(
      "Carried {nrow(zero_rows)} zero-encounter fishery-year x group combination{?s} \\
       from the inventory as catch = 0, so they read as examined rather than missing."
    )
  }
} else if (!is.null(catch_base)) {
  catch_base <- mutate(catch_base, baseline_source = "fitted")
}

dat <- b_summary |>
  left_join(
    comp |> select(fishery_name, basin, fishery_label, fishery_type, year_start),
    by = "fishery_name"
  ) |>
  mutate(
    log_b    = log(median),
    se_log   = sd / median,          # delta method, as in 06 [V2]
    vi       = se_log^2,
    informed = !is.na(prior_contraction) & prior_contraction >= INFORMED_MIN_CONTRACTION &
                 (is.na(informed_flag) | informed_flag != "unconverged")
  )

usable <- dat |>
  filter(!is.na(median), median > 0, !is.na(vi), vi > 0, is.finite(vi), informed)

cli::cli_alert_info(
  "{nrow(usable)} usable fishery-year x bias-type estimates out of {nrow(dat)} \\
   (informed, converged, positive posterior median)."
)

# ------------------------------------------------------------------------------
# T4 -- pass-through diagnostic
#
# Elasticity -1 is exact only under stated conditions. Rather than assert it,
# report per fishery-year what the data actually supports:
#
#   E_n > 0        census present -> census pins effort directly, so an error in
#                  an imported b is partly absorbed. Pass-through is DAMPED,
#                  somewhere strictly between -1 and 0.
#   A_n > 0        angler index counts present. These carry no b (line ~196) but
#                  run through p_I, a free beta(0.5,0.5) parameter confounded
#                  with effort the same way b is -- so they resist an imported b
#                  only weakly, via that prior. Flagged, not counted as an anchor.
#   V_n, T_n       vehicle / trailer index channels. Full -1 pass-through
#                  requires importing the b for EVERY live channel. Import b[1]
#                  alone into a fishery that also runs trailer counts and the
#                  trailer channel pulls back -- pass-through lands between
#                  -1 and 0 even with no census.
# ------------------------------------------------------------------------------

T4 <- dims |>
  select(fishery_name, V_n, T_n, A_n, E_n, IntA, IntC) |>
  left_join(comp |> select(fishery_name, basin, fishery_type, year_start), by = "fishery_name") |>
  mutate(
    has_vehicle_index = V_n > 0,
    has_trailer_index = T_n > 0,
    has_angler_index  = A_n > 0,
    has_census        = E_n > 0,
    n_b_channels      = as.integer(has_vehicle_index) + as.integer(has_trailer_index),
    pass_through_class = case_when(
      has_census                      ~ "damped (census anchors effort)",
      has_angler_index                ~ "near-full (angler counts resist weakly via p_I)",
      n_b_channels >= 1               ~ "full (-1): b and effort perfectly confounded",
      TRUE                            ~ "no index channel -- check inputs"
    ),
    # The scenario actually proposed for this year: no census collected. Under
    # that counterfactual the census column above is irrelevant and what governs
    # is whether every live b-channel gets an imported value.
    census_free_pass_through = case_when(
      has_angler_index  ~ "near-full",
      n_b_channels >= 1 ~ "full (-1)",
      TRUE              ~ "undefined"
    )
  ) |>
  arrange(basin, fishery_type, year_start)

write_csv(T4, file.path(OUT_DIR, "bss_b_T4_pass_through.csv"))
cli::cli_alert_success("T4 pass-through written ({nrow(T4)} fishery-years).")

# ------------------------------------------------------------------------------
# T5 -- the deterministic ladder
#
# catch(b_alt) / catch(b_ref) = b_ref / b_alt. Two reference framings:
#   "empirical"  b_ref = the series pooled_b; tiers = that series' own T2
#                prediction-interval bounds. What the predicted b actually
#                implies for this year, per fishery.
#   "round"      b_ref = 1 (no bias correction); tiers = ROUND_LADDER. Basin-
#                agnostic intuition for the shape of 1/b.
# ------------------------------------------------------------------------------

t2_path <- file.path(OUT_DIR, "bss_b_T2_variability.csv")
T2 <- if (file.exists(t2_path)) read_csv(t2_path, show_col_types = FALSE) else NULL
if (is.null(T2)) cli::cli_abort("{.file {t2_path}} not found -- run 06_variability_analysis.R first.")

ladder_round <- tidyr::expand_grid(
  basin = NA_character_, fishery_type = "(any)", bias_type = "(any)",
  framing = "round", b_ref = 1, b_alt = ROUND_LADDER
) |>
  mutate(tier = case_when(b_alt < 1 ~ "below 1", b_alt == 1 ~ "reference", TRUE ~ "above 1"))

ladder_emp <- T2 |>
  filter(!is.na(pooled_b), !is.na(pi_lb), !is.na(pi_ub)) |>
  select(basin, fishery_type, bias_type, pooled_b, pi_lb, pi_ub) |>
  pivot_longer(c(pi_lb, pooled_b, pi_ub), names_to = "tier", values_to = "b_alt") |>
  mutate(
    framing = "empirical",
    tier = recode(tier, pi_lb = "PI low", pooled_b = "predicted (pooled)", pi_ub = "PI high")
  ) |>
  left_join(T2 |> select(basin, fishery_type, bias_type, b_ref = pooled_b),
            by = c("basin", "fishery_type", "bias_type"))

T5 <- bind_rows(ladder_emp, ladder_round) |>
  mutate(
    catch_multiplier = b_ref / b_alt,
    pct_change_catch = 100 * (catch_multiplier - 1),
    direction = case_when(
      abs(catch_multiplier - 1) < 1e-9 ~ "no change",
      catch_multiplier > 1             ~ "catch revised UP",
      TRUE                             ~ "catch revised DOWN"
    )
  ) |>
  select(framing, basin, fishery_type, bias_type, tier, b_ref, b_alt,
         catch_multiplier, pct_change_catch, direction) |>
  arrange(framing, basin, fishery_type, bias_type, b_alt)

write_csv(T5, file.path(OUT_DIR, "bss_b_T5_ladder.csv"))
cli::cli_alert_success("T5 ladder written ({nrow(T5)} rows).")

# ------------------------------------------------------------------------------
# T7 -- THE DIRECT SENSITIVITY: vary b on a real dataset, read off the catch
#
# The plainest form of the question. Take a fishery-year that has actually been
# fitted, hold everything else at what that fit produced, and swap `b` for a
# range of values. Because catch goes exactly as 1/b, no re-fitting is needed:
#
#     catch(b_alt) = C_sum(fitted) * b_fitted / b_alt
#
# The anchor is that fishery-year's OWN fitted b, not 1 -- C_sum came out of a
# fit that used b_fitted, so that is the point the curve passes through.
#
# Tiers are both the round ladder (readable, basin-agnostic) and the series'
# own T2 prediction bounds (what b could actually plausibly be here).
#
# CAVEAT, and it is the same one T4 records: this scales the b for ONE
# likelihood type. The clean 1/b holds when every live index channel's b moves
# together, or when only one channel is in use. Scale b[1] alone in a fishery
# that also runs trailer counts and the trailer channel resists, so the true
# effect is somewhat smaller than shown. Read T4 alongside this.
# ------------------------------------------------------------------------------

b_fitted <- dat |>
  filter(!is.na(median), median > 0) |>
  select(fishery_name, basin, fishery_type, year_start, bias_type,
         b_fitted = median, informed, informed_flag)

series_bounds <- T2 |>
  select(basin, fishery_type, bias_type, pi_lb, pooled_b, pi_ub)

# One row per fishery-year x bias type x tier.
tier_grid <- b_fitted |>
  left_join(series_bounds, by = c("basin", "fishery_type", "bias_type")) |>
  mutate(round_tiers = list(ROUND_LADDER)) |>
  rowwise() |>
  mutate(b_values = list(c(
    stats::setNames(ROUND_LADDER, paste0("b = ", ROUND_LADDER)),
    stats::setNames(c(pi_lb, pooled_b, pi_ub),
                    c("series PI low", "series predicted", "series PI high")),
    stats::setNames(b_fitted, "as fitted")
  ))) |>
  ungroup() |>
  select(-round_tiers) |>
  mutate(tier = map(b_values, names), b_alt = map(b_values, unname)) |>
  select(-b_values) |>
  unnest(c(tier, b_alt)) |>
  filter(!is.na(b_alt), b_alt > 0)

T7 <- tier_grid |>
  mutate(
    tier_kind = case_when(
      tier == "as fitted"          ~ "anchor",
      str_starts(tier, "series")   ~ "empirical",
      TRUE                          ~ "round"
    ),
    catch_multiplier = b_fitted / b_alt,
    pct_change_catch = 100 * (catch_multiplier - 1),
    direction = case_when(
      abs(catch_multiplier - 1) < 1e-9 ~ "no change",
      catch_multiplier > 1             ~ "catch UP",
      TRUE                             ~ "catch DOWN"
    )
  )

# In fish, wherever a C_sum baseline exists. One ratio row becomes one row per
# catch group, because the multiplier is identical across groups.
if (!is.null(catch_base) && all(c("fishery_name", "est_cg", "C_sum_median") %in% names(catch_base))) {
  T7 <- T7 |>
    left_join(catch_base |> select(fishery_name, est_cg, C_sum_median, baseline_source),
              by = "fishery_name", relationship = "many-to-many") |>
    # A zero baseline stays zero at every b: 0 * anything is 0. That is correct
    # and worth seeing -- no value of b turns an unobserved encounter into one.
    mutate(catch_estimate = C_sum_median * catch_multiplier)
} else {
  T7 <- T7 |> mutate(est_cg = NA_character_, C_sum_median = NA_real_,
                     baseline_source = NA_character_, catch_estimate = NA_real_)
}

T7 <- T7 |>
  select(basin, fishery_type, fishery_name, year_start, bias_type,
         tier, tier_kind, b_fitted, b_alt, catch_multiplier, pct_change_catch,
         direction, est_cg, catch_baseline = C_sum_median, baseline_source,
         catch_estimate, informed_flag) |>
  arrange(basin, fishery_type, year_start, bias_type, b_alt)

write_csv(T7, file.path(OUT_DIR, "bss_b_T7_direct_sensitivity.csv"))
cli::cli_alert_success("T7 direct sensitivity written ({nrow(T7)} rows).")

# ------------------------------------------------------------------------------
# T6 -- the leave-one-out backtest  [S1] [S2]
#
# For each series (basin x fishery_type x bias_type) and each year in it:
#   1. refit the random-effects model on THE OTHER YEARS ONLY
#   2. draw from that model's predictive for the held-out year
#          log b_pred ~ Normal(mu_loo, tau^2_loo + SE(mu_loo)^2)
#      (the same arithmetic T2's prediction interval uses -- verified: metafor
#      takes z = 1.96, and PI half-width = 1.96 * sqrt(tau^2 + SE^2))
#   3. draw from the MEASURED posterior for that year (the saved b draws)
#   4. ratio = b_measured / b_pred, draw by draw
#
# ratio = catch_predicted / catch_measured, because catch goes as 1/b.
#   ratio > 1  importing b OVERSTATES catch  (predicted b came in too low)
#   ratio < 1  importing b UNDERSTATES catch (predicted b came in too high)
#
# The two uncertainties are independent sources and are sampled independently:
# the predictive is P(b_new | other years); the measured posterior is this
# year's own estimate of that same quantity.
# ------------------------------------------------------------------------------

# Saved draws are per fishery-year, columns "b[1]" (vehicle) / "b[2]" (trailer).
read_b_draws <- function(fishery_name, bias_type) {
  path <- file.path(DRAWS_DIR, paste0(safe_name(fishery_name), ".rds"))
  if (!file.exists(path)) return(NULL)
  d <- try(readRDS(path), silent = TRUE)
  if (inherits(d, "try-error")) return(NULL)
  col <- switch(bias_type, vehicle = "b[1]", trailer = "b[2]", NA_character_)
  if (is.na(col) || !col %in% names(d)) return(NULL)
  v <- as.numeric(d[[col]])
  v <- v[is.finite(v) & v > 0]
  if (length(v) < 50) return(NULL)
  v
}

loo_one <- function(series_df, i) {
  held <- series_df[i, ]
  rest <- series_df[-i, ]
  if (nrow(rest) < 2) return(NULL)

  f <- try(metafor::rma(yi = log_b, vi = vi, data = rest, method = "REML"), silent = TRUE)
  if (inherits(f, "try-error")) return(NULL)

  mu_loo  <- as.numeric(f$b)
  se_loo  <- as.numeric(f$se)
  sd_pred <- sqrt(f$tau2 + se_loo^2)       # predictive SD for a NEW year, log scale

  # Measured side: the full posterior, not a point [S4]. Falls back to the
  # lognormal implied by the summary if the draws file is missing, so a series
  # is never silently dropped for a file-management reason.
  meas <- read_b_draws(held$fishery_name, held$bias_type)
  meas_source <- if (is.null(meas)) "summary (lognormal approx)" else "posterior draws"
  log_meas <- if (is.null(meas)) {
    rnorm(N_SIM, held$log_b, held$se_log)
  } else {
    log(sample(meas, N_SIM, replace = TRUE))
  }

  log_pred  <- rnorm(N_SIM, mu_loo, sd_pred)
  log_ratio <- log_meas - log_pred          # = log(catch_pred / catch_meas)
  ratio     <- exp(log_ratio)

  # Two different intervals, and conflating them biases the calibration check.
  #
  # sd_pred covers a new year's TRUE b: sqrt(tau^2 + SE^2). That is the right
  # interval for the value you would IMPORT, and it is what pred_pi_* reports.
  #
  # The calibration check compares that prediction against a MEASURED b, which
  # carries its own estimation error v_i on top of the true value. So the
  # interval a held-out observation should fall inside is
  # sqrt(tau^2 + SE^2 + v_i) -- wider. Checking against sd_pred alone asks the
  # prediction to hit a noisy target exactly and reports too many misses; an
  # earlier version of this did that and put coverage at 72% when the honest
  # figure is higher.
  sd_obs <- sqrt(f$tau2 + se_loo^2 + held$vi)
  pi_lo  <- mu_loo - 1.96 * sd_pred     # prediction for this year's true b
  pi_hi  <- mu_loo + 1.96 * sd_pred
  obs_lo <- mu_loo - 1.96 * sd_obs      # where a measured b should land
  obs_hi <- mu_loo + 1.96 * sd_obs

  tibble(
    basin = held$basin, fishery_type = held$fishery_type, bias_type = held$bias_type,
    fishery_name = held$fishery_name, year_start = held$year_start,
    n_years_in_loo_fit = nrow(rest),
    measured_b   = held$median,
    predicted_b  = exp(mu_loo),
    pred_pi_lb   = exp(pi_lo),
    pred_pi_ub   = exp(pi_hi),
    tau_loo      = sqrt(f$tau2),
    obs_pi_lb    = exp(obs_lo),
    obs_pi_ub    = exp(obs_hi),
    inside_95_pi = held$log_b >= obs_lo && held$log_b <= obs_hi,
    # The stricter question, kept alongside: would the measured value have
    # fallen inside the interval we would actually QUOTE for this year's b?
    inside_95_pred = held$log_b >= pi_lo && held$log_b <= pi_hi,
    ratio_point  = held$median / exp(mu_loo),
    ratio_median = median(ratio),
    ratio_q10    = unname(quantile(ratio, 0.10)),
    ratio_q90    = unname(quantile(ratio, 0.90)),
    ratio_q2.5   = unname(quantile(ratio, 0.025)),
    ratio_q97.5  = unname(quantile(ratio, 0.975)),
    pct_err_point   = 100 * (held$median / exp(mu_loo) - 1),
    abs_pct_err_med = 100 * median(abs(ratio - 1)),
    p_overstates    = mean(ratio > 1),
    measured_source = meas_source
  )
}

series_list <- usable |>
  group_by(basin, fishery_type, bias_type) |>
  filter(n() >= LOO_MIN_SERIES_YEARS) |>
  group_split()

if (length(series_list) == 0) {
  cli::cli_abort(
    "No series has {LOO_MIN_SERIES_YEARS}+ usable years -- a leave-one-out \\
     backtest is not estimable. Check informed_flag / prior_contraction in \\
     {.file bss_b_summary.csv}."
  )
}

T6 <- map_dfr(series_list, function(sdf) {
  sdf <- arrange(sdf, year_start)
  map_dfr(seq_len(nrow(sdf)), ~loo_one(sdf, .x))
})

if (nrow(T6) == 0) cli::cli_abort("Leave-one-out produced no rows -- every series fit failed.")

# In-fish translation, when 01b has run. The ratio is unitless and identical
# across catch groups, so this is a pure join-and-multiply: one ratio row
# becomes one row per catch group.
if (!is.null(catch_base)) {
  needed <- c("fishery_name", "est_cg", "C_sum_median")
  missing_cols <- setdiff(needed, names(catch_base))
  if (length(missing_cols) > 0) {
    cli::cli_warn("{.file bss_catch_baseline.csv} lacks column{?s} {.val {missing_cols}} -- skipping the in-fish join.")
  } else {
    T6 <- T6 |>
      left_join(catch_base |> select(all_of(needed)), by = "fishery_name",
                relationship = "many-to-many") |>
      mutate(
        catch_measured  = C_sum_median,
        catch_predicted = C_sum_median * ratio_point,
        fish_error      = catch_predicted - catch_measured,
        fish_err_lo     = C_sum_median * (ratio_q2.5  - 1),
        fish_err_hi     = C_sum_median * (ratio_q97.5 - 1)
      )
  }
}

write_csv(T6, file.path(OUT_DIR, "bss_b_T6_loo_backtest.csv"))
cli::cli_alert_success("T6 backtest written ({nrow(T6)} held-out fishery-years).")

# --- Calibration: the [S5] empirical check on exchangeability ------------------

T6_cal <- T6 |>
  distinct(basin, fishery_type, bias_type, fishery_name, year_start, inside_95_pi,
           inside_95_pred, abs_pct_err_med, pct_err_point) |>
  group_by(basin, fishery_type, bias_type) |>
  summarise(
    n_years_tested     = n(),
    n_inside_95_pi     = sum(inside_95_pi),
    pct_inside_95_pi   = 100 * mean(inside_95_pi),
    pct_inside_95_pred = 100 * mean(inside_95_pred),
    median_abs_pct_err = median(abs(pct_err_point)),
    worst_pct_err      = pct_err_point[which.max(abs(pct_err_point))],
    .groups = "drop"
  ) |>
  arrange(bias_type, basin, fishery_type)

write_csv(T6_cal, file.path(OUT_DIR, "bss_b_T6_calibration.csv"))

cov_df <- distinct(T6, fishery_name, bias_type, inside_95_pi, inside_95_pred)
overall_cov  <- 100 * mean(cov_df$inside_95_pi)
overall_pred <- 100 * mean(cov_df$inside_95_pred)
cli::cli_alert_info(
  "Calibration [S5]: {round(overall_cov)}% of held-out years fell inside their \\
   leave-one-out 95% interval for a MEASURED b (nominal 95%, n = {nrow(cov_df)})."
)
cli::cli_alert_info(
  "Stricter: {round(overall_pred)}% fell inside the narrower interval we would \\
   QUOTE for this year's true b. The gap between the two is the measurement \\
   error in the yardstick, not a failure of the prediction."
)

# ------------------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------------------

# fig17 -- the 1/b curve. ONE deterministic relationship, so one line and no
# legend box (the title names it); the series tiers ride on it as points. Log-log,
# where 1/b is straight -- on a linear axis it reads as if low b were punished
# harder than high b, which is an artifact of the ruler, not the model.
curve_df <- tibble(b_alt = exp(seq(log(0.3), log(3), length.out = 400)),
                   catch_multiplier = 1 / b_alt)

emp_pts <- T5 |>
  filter(framing == "empirical") |>
  mutate(series = paste0(fishery_type, " (", bias_type, ")"),
         catch_multiplier_v1 = 1 / b_alt)

fig17 <- ggplot(curve_df, aes(b_alt, catch_multiplier)) +
  geom_hline(yintercept = 1, color = BASELINE_COL, linewidth = 0.4) +
  geom_vline(xintercept = 1, color = BASELINE_COL, linewidth = 0.4) +
  geom_line(color = INK_SECOND, linewidth = 0.9) +
  geom_point(
    data = emp_pts, aes(y = catch_multiplier_v1, color = bias_type, shape = tier),
    size = 2.6, stroke = 0.9
  ) +
  scale_x_log10(breaks = c(0.3, 0.5, 0.8, 1, 1.25, 1.5, 2, 3)) +
  scale_y_log10(breaks = c(0.33, 0.5, 0.8, 1, 1.25, 2, 3)) +
  scale_color_manual(values = c(vehicle = CAT[["blue"]], trailer = CAT[["orange"]]), name = NULL) +
  scale_shape_manual(values = c("PI low" = 1, "predicted (pooled)" = 19, "PI high" = 2), name = NULL) +
  labs(
    title = "Estimated catch moves as 1 / b",
    x = "Effort-index bias term b (log scale)",
    y = "Catch multiplier relative to b = 1 (log scale)"
  ) +
  theme_bss()

save_fig(fig17, "fig17_catch_ladder", width = 9, height = 6)

# fig18 -- the backtest. Interval per held-out year on a log ratio axis, so
# over- and under-prediction are visually symmetric. Reference line at 1.
fig18_df <- T6 |>
  distinct(basin, fishery_type, bias_type, fishery_name, year_start,
           ratio_median, ratio_q10, ratio_q90, ratio_q2.5, ratio_q97.5, inside_95_pi) |>
  mutate(series = paste0(fishery_type, "\n(", bias_type, ")"))

fig18 <- ggplot(fig18_df, aes(x = factor(year_start), y = ratio_median, color = bias_type)) +
  geom_hline(yintercept = 1, color = INK_SECOND, linewidth = 0.5, linetype = "22") +
  geom_linerange(aes(ymin = ratio_q2.5, ymax = ratio_q97.5), linewidth = 0.5, alpha = 0.55) +
  geom_linerange(aes(ymin = ratio_q10,  ymax = ratio_q90),  linewidth = 1.6, alpha = 0.9) +
  geom_point(size = 2.4, color = SURFACE) +
  geom_point(size = 1.7) +
  facet_wrap(~series, scales = "free_x", ncol = 4) +
  scale_y_log10() +
  scale_color_manual(values = c(vehicle = CAT[["blue"]], trailer = CAT[["orange"]]), name = NULL) +
  labs(
    title = "What importing a predicted b would have cost, year by year",
    x = NULL, y = "catch (imported b) / catch (measured b), log scale"
  ) +
  theme_bss()

save_fig(fig18, "fig18_loo_backtest", width = 11, height = 7)

# fig19 -- HOW EXPOSED IS EACH FISHERY to getting b wrong?
#
# An earlier version drew catch against b, one line per fishery-year. That was
# the wrong form: catch = 1/b is deterministic and identical everywhere, so
# every panel was the same fixed-slope line shifted sideways, and fig17 already
# shows that curve once. Nothing about a fishery was visible in it.
#
# What actually differs between fisheries is how much of the curve is IN PLAY
# -- how wide that series' plausible b range is, and therefore how far catch
# could move. That is a range comparison across categories, so: one row per
# series, a segment spanning the catch outcomes implied by its own prediction
# interval, ordered by exposure. Reading down the axis ranks the fisheries by
# how much a missed b would cost.
exposure <- T2 |>
  filter(!is.na(pooled_b), !is.na(pi_lb), !is.na(pi_ub), pi_lb > 0) |>
  transmute(
    basin, fishery_type, bias_type,
    series = paste0(fishery_type, "  (", bias_type, ")"),
    # b LOW -> catch UP, and vice versa: the multiplier flips the bounds.
    mult_hi = pooled_b / pi_lb,
    mult_lo = pooled_b / pi_ub,
    span    = mult_hi / mult_lo
  ) |>
  arrange(span) |>
  mutate(series = factor(series, levels = series))

mult_breaks <- c(0.25, 0.5, 0.67, 1, 1.5, 2, 3, 4)
mult_labels <- ifelse(
  abs(mult_breaks - 1) < 1e-9, "no change",
  sprintf("%+.0f%%", 100 * (mult_breaks - 1))
)

fig19 <- ggplot(exposure, aes(y = series, color = bias_type)) +
  geom_vline(xintercept = 1, color = INK_SECOND, linewidth = 0.5, linetype = "22") +
  geom_linerange(aes(xmin = mult_lo, xmax = mult_hi), linewidth = 2.4, alpha = 0.9) +
  geom_point(aes(x = mult_lo), size = 2.2, shape = 18) +
  geom_point(aes(x = mult_hi), size = 2.2, shape = 18) +
  geom_text(
    aes(x = mult_hi, label = sprintf("%.1fx span", span)),
    hjust = -0.25, size = 3, color = INK_SECOND, show.legend = FALSE
  ) +
  scale_x_log10(breaks = mult_breaks, labels = mult_labels,
                expand = expansion(mult = c(0.05, 0.22))) +
  scale_color_manual(values = c(vehicle = CAT[["blue"]], trailer = CAT[["orange"]]),
                     name = NULL) +
  labs(
    title = "How much could estimated catch move, if b is wrong?",
    x = "Change in estimated catch (log scale)", y = NULL
  ) +
  theme_bss() +
  theme(panel.grid.major.y = element_blank())

save_fig(fig19, "fig19_exposure", width = 10, height = 5.5)

cli::cli_alert_success("Figures written to {.path {FIG_DIR}}.")

cli::cli_h2("T7 -- what changing b does to catch, per fishery (most recent year)")
T7 |>
  filter(tier_kind != "anchor") |>
  group_by(fishery_type, bias_type) |>
  filter(year_start == max(year_start)) |>
  ungroup() |>
  filter(tier_kind == "empirical") |>
  select(fishery_type, bias_type, tier, b_alt, pct_change_catch, catch_estimate) |>
  arrange(fishery_type, bias_type, b_alt) |>
  print(n = Inf)

# ------------------------------------------------------------------------------
# Console summary -- the numbers to walk into the meeting with
# ------------------------------------------------------------------------------

cli::cli_h2("Backtest summary by series")
print(T6_cal, n = Inf)

cli::cli_h2("Worst single-year errors")
T6 |>
  distinct(fishery_name, bias_type, year_start, measured_b, predicted_b, pct_err_point) |>
  arrange(desc(abs(pct_err_point))) |>
  head(10) |>
  print(n = Inf)

cli::cli_h2("Pass-through classes present")
T4 |> count(pass_through_class, census_free_pass_through) |> print(n = Inf)
