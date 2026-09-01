# ==============================================================================
# 04_candidate_options.R
#
# Purpose:
#   Compute, per fishery-series (basin x fishery_type -- the fishery name
#   WITHOUT its year, so all years of one fishery form a series) and separately for
#   b[1] (vehicle) and b[2] (trailer), the four candidate bias-correction
#   options the meeting needs to choose between:
#     (a) time-series mean
#     (b) most recent year only
#     (c) precision-weighted average (inverse posterior variance), with
#         Cochran's Q / I^2 / tau^2 heterogeneity diagnostics
#     (d1) hierarchical / partial-pooling shrinkage, POST-HOC over the
#          per-year posterior summaries (log scale, DerSimonian-Laird
#          random-effects meta-analysis) -- NOT a re-fit of the Stan model.
#          (d2), a proper in-model hierarchical re-specification of `b`
#          itself, is deliberately NOT built here -- see README.md's
#          "explicitly deferred" note. (d1) is a same-day stand-in that
#          argues for the same conclusion without a multi-day model rewrite.
#
# Inclusion rule (toggle these two constants live in the meeting if the group
# wants to see it change the numbers):
#   comparability_tier %in% INCLUDED_TIERS
#   informed_flag != "unconverged"
#
# Usage:
#   Rscript analysis/bss_bias/04_candidate_options.R
#   Requires bss_b_summary.csv and bss_b_comparability.csv.
#
# Outputs (analysis/bss_bias/outputs/):
#   bss_b_candidate_options.csv
#   bss_b_candidate_options.html  (gt render)
# ==============================================================================

library(tidyverse)
library(gt)
library(cli)
library(here)

OUT_DIR <- here::here("analysis", "bss_bias", "outputs")
summary_path <- file.path(OUT_DIR, "bss_b_summary.csv")
comp_path    <- file.path(OUT_DIR, "bss_b_comparability.csv")
if (!file.exists(summary_path)) cli::cli_abort("{.file {summary_path}} not found -- run 01_fit_bss_bias.R first.")
if (!file.exists(comp_path))    cli::cli_abort("{.file {comp_path}} not found -- run 02_build_comparability_table.R first.")

INCLUDED_TIERS <- c("reference", "comparable", "comparable-with-caveat")

b_summary <- read_csv(summary_path, show_col_types = FALSE)
comp      <- read_csv(comp_path, show_col_types = FALSE)

dat <- b_summary |>
  left_join(comp |> select(fishery_name, basin, fishery_type, year_start, comparability_tier),
            by = "fishery_name") |>
  filter(comparability_tier %in% INCLUDED_TIERS, informed_flag != "unconverged") |>
  mutate(
    log_b  = log(median),
    se_log = sd / median   # delta method: Var(log X) ~= Var(X) / X^2, so SE(log X) ~= SE(X)/X
  )

if (nrow(dat) == 0) {
  cli::cli_abort("No rows survive the inclusion filter (INCLUDED_TIERS / informed_flag != 'unconverged'). Nothing to compute options over yet.")
}

# ------------------------------------------------------------------------------
# DerSimonian-Laird random-effects meta-analysis, on the log(b) scale.
# Implemented from the standard formulas (no metafor dependency, to keep this
# runnable with a minimal package set on a laptop tonight).
# ------------------------------------------------------------------------------

dl_meta <- function(y, v) {
  k <- length(y)
  w_fe <- 1 / v
  theta_fe <- sum(w_fe * y) / sum(w_fe)
  var_fe <- 1 / sum(w_fe)

  Q  <- sum(w_fe * (y - theta_fe)^2)
  df <- k - 1
  p_Q <- if (df > 0) pchisq(Q, df, lower.tail = FALSE) else NA_real_
  I2 <- if (df > 0) max(0, (Q - df) / Q) * 100 else NA_real_

  tau2 <- if (df > 0) {
    denom <- sum(w_fe) - sum(w_fe^2) / sum(w_fe)
    max(0, (Q - df) / denom)
  } else 0

  w_re <- 1 / (v + tau2)
  theta_re <- sum(w_re * y) / sum(w_re)
  var_re <- 1 / sum(w_re)

  # per-year shrinkage estimates (empirical Bayes): pull each y_i toward theta_re
  # in proportion to its own precision vs. the between-year precision (1/tau2)
  shrunk <- if (tau2 > 0) {
    w_within <- 1 / v
    w_tau <- 1 / tau2
    (w_within * y + w_tau * theta_re) / (w_within + w_tau)
  } else {
    rep(theta_re, k)
  }

  list(theta_fe = theta_fe, var_fe = var_fe, Q = Q, df = df, p_Q = p_Q, I2 = I2, tau2 = tau2,
       theta_re = theta_re, var_re = var_re, shrunk = shrunk)
}

# ------------------------------------------------------------------------------
# Per-series, per-bias-type option computation
# ------------------------------------------------------------------------------

compute_options_for_series <- function(df_series) {
  df_series <- df_series |> arrange(desc(year_start))
  n_years <- nrow(df_series)

  # (a) time-series mean -- unweighted mean of posterior medians, two uncertainty views
  mean_est <- mean(df_series$median)
  se_between <- if (n_years > 1) sd(df_series$median) / sqrt(n_years) else NA_real_
  se_within  <- sqrt(mean(df_series$post_var) / n_years)

  # (b) most recent year only
  recent <- df_series |> slice(1)

  out <- tibble::tribble(
    ~option, ~estimate, ~lower, ~upper, ~notes,
    "time_series_mean", mean_est, mean_est - 1.96*se_between, mean_est + 1.96*se_between,
      paste0("SE(between-year) shown; SE(within-year, pooled)=", round(se_within, 4), " also computed -- see notes"),
    "most_recent_year", recent$median, recent$q2.5, recent$q97.5,
      paste0("Year ", recent$year_start, "; ", n_years - 1, " earlier comparable year(s) discarded")
  )

  if (n_years >= 2) {
    dl <- dl_meta(df_series$log_b, df_series$se_log^2)
    pw_est <- exp(dl$theta_fe); pw_lo <- exp(dl$theta_fe - 1.96*sqrt(dl$var_fe)); pw_hi <- exp(dl$theta_fe + 1.96*sqrt(dl$var_fe))
    re_est <- exp(dl$theta_re); re_lo <- exp(dl$theta_re - 1.96*sqrt(dl$var_re)); re_hi <- exp(dl$theta_re + 1.96*sqrt(dl$var_re))
    out <- bind_rows(out, tibble::tribble(
      ~option, ~estimate, ~lower, ~upper, ~notes,
      "precision_weighted", pw_est, pw_lo, pw_hi,
        paste0("Fixed-effect (inverse-variance) on log scale. Q=", round(dl$Q,2), ", df=", dl$df,
               ", p=", signif(dl$p_Q,3), ", I2=", round(dl$I2,1), "%, tau2=", round(dl$tau2,4),
               ". High I2 => precision-weighting understates real between-year spread -- see hierarchical row."),
      "hierarchical_shrinkage", re_est, re_lo, re_hi,
        paste0("Post-hoc DerSimonian-Laird random-effects pooling (log scale), NOT an in-model re-fit. tau2=",
               round(dl$tau2, 4), " (", if (dl$tau2 < 1e-6) "effectively full pooling" else "meaningful between-year variance", ").")
    ))
  } else {
    out <- bind_rows(out, tibble::tribble(
      ~option, ~estimate, ~lower, ~upper, ~notes,
      "precision_weighted", NA_real_, NA_real_, NA_real_, "Needs >=2 comparable years; not computed.",
      "hierarchical_shrinkage", NA_real_, NA_real_, NA_real_, "Needs >=2 comparable years; not computed."
    ))
  }

  out |> mutate(n_years_used = n_years)
}

candidate_options <- dat |>
  group_by(basin, fishery_type, bias_type) |>
  group_modify(~ compute_options_for_series(.x)) |>
  ungroup()

write_csv(candidate_options, file.path(OUT_DIR, "bss_b_candidate_options.csv"))
cli::cli_alert_success("Wrote {nrow(candidate_options)} rows to bss_b_candidate_options.csv")
candidate_options |> filter(bias_type == "vehicle") |> select(basin, fishery_type, option, estimate, lower, upper, n_years_used) |> print(n = 100)

# ------------------------------------------------------------------------------
# gt render (vehicle bias only, the headline parameter)
# ------------------------------------------------------------------------------

gt_tbl <- candidate_options |>
  filter(bias_type == "vehicle") |>
  mutate(across(c(estimate, lower, upper), ~ round(.x, 3))) |>
  select(basin, fishery_type, option, estimate, lower, upper, n_years_used, notes) |>
  gt(groupname_col = "basin", rowname_col = "fishery_type") |>
  tab_header(title = "Candidate bias-correction options (b₁, vehicle)",
             subtitle = paste0("Included tiers: ", paste(INCLUDED_TIERS, collapse = ", "))) |>
  cols_label(option = "Option", estimate = "Estimate", lower = "Lower 95%", upper = "Upper 95%",
             n_years_used = "N years") |>
  opt_row_striping()

gt::gtsave(gt_tbl, file.path(OUT_DIR, "bss_b_candidate_options.html"))
cli::cli_alert_success("Wrote bss_b_candidate_options.html")

cli::cli_alert_info("Figure 6 (candidate options overlaid on the b-vs-year series) is generated by re-running 03_plot_b_series.R after this script -- see README.md.")
