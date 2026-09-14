# ==============================================================================
# 08_diagnose_catch_fit.R
#
# Purpose:
#   Fit ONE fishery-year and ONE catch group, keep the stanfit, and put the
#   model's implied CPUE next to the CPUE actually observed in interviews.
#
#   Written because bss_catch_baseline.csv returned a season Chinook total of
#   40,426 for Skagit fall salmon 2021 against 25 interviewed fish. The
#   suspicion is that C_sum collapses toward E_sum -- catch is effort x CPUE
#   (Stan line 240), so a CPUE stuck near exp(0) = 1 fish/hour makes the catch
#   total equal the effort total. 1 fish/hour is an absurd CPUE for these
#   fisheries; it is seen only in exceptional cases like a pink year. If the
#   model reports it, the catch half of the model is not being informed.
#
#   THE DECISIVE NUMBER is cpue_ratio at the bottom: model-implied CPUE divided
#   by observed CPUE. Near 1 means the catch model is working. Orders of
#   magnitude above 1 means it is not, and no season total from that fit means
#   anything.
#
# Usage (R Console):
#   source("analysis/bss_bias/08_diagnose_catch_fit.R")
#
#   Defaults to Snohomish fall salmon 2021, coho harvest, "quick". Override
#   before sourcing:
#     DIAG_FISHERY <- "Skagit fall salmon 2021"
#     DIAG_GROUP   <- "chinook_all"        # chinook_all | coho_harvest
#     DIAG_CONFIG  <- "prod"               # smoke | lite | quick | prod
#
#   Keeps the full stanfit AND the Stan input list in outputs/fits/ (gitignored),
#   available afterwards as `diag_fit` and `diag_inp`. Observed CPUE is taken
#   from diag_inp$c and diag_inp$h -- the exact vectors the model was handed --
#   so the comparison cannot drift from what was actually fitted.
# ==============================================================================

library(tidyverse)
library(cli)
library(here)

if (!exists("DIAG_FISHERY", inherits = FALSE)) DIAG_FISHERY <- "Snohomish fall salmon 2021"
if (!exists("DIAG_GROUP",   inherits = FALSE)) DIAG_GROUP   <- "coho_harvest"
if (!exists("DIAG_CONFIG",  inherits = FALSE)) DIAG_CONFIG  <- "quick"

source(here::here("analysis", "bss_bias", "catch_groups.R"))
source(here::here("analysis", "bss_bias", "fishery_data.R"))

if (!DIAG_GROUP %in% names(CATCH_GROUPS)) {
  cli::cli_abort("Unknown DIAG_GROUP {.val {DIAG_GROUP}}. Available: {.val {names(CATCH_GROUPS)}}")
}
grp <- CATCH_GROUPS[[DIAG_GROUP]]
ecg <- catch_group_label(grp)

cli::cli_h1("08 -- catch-fit diagnostic")
cli::cli_alert_info("Fishery: {.val {DIAG_FISHERY}}")
cli::cli_alert_info("Group:   {.val {ecg}}")
cli::cli_alert_info("Config:  {.val {DIAG_CONFIG}}")

# ------------------------------------------------------------------------------
# 1. Fit, reusing 01's machinery exactly -- same preps, same preflight, same
#    priors. Nothing here is a parallel implementation.
# ------------------------------------------------------------------------------

CATCH_BASELINE_ONLY <- TRUE
RUN_CATCH_GROUP     <- grp
ONLY_FISHERIES      <- DIAG_FISHERY
FIT_CONFIG_NAME     <- DIAG_CONFIG
SAVE_FITS           <- TRUE        # keeps the stanfit AND the Stan inputs

source(here::here("analysis", "bss_bias", "01_fit_bss_bias.R"), local = FALSE)

FITS_DIR_D <- here::here("analysis", "bss_bias", "outputs", "fits")
fit_path   <- file.path(FITS_DIR_D, paste0(safe_name(DIAG_FISHERY), ".rds"))
inp_path   <- file.path(FITS_DIR_D, paste0(safe_name(DIAG_FISHERY), "__inputs.rds"))
if (!file.exists(fit_path)) {
  cli::cli_abort("No stanfit at {.file {fit_path}} -- the fit did not complete. Read the ledger.")
}
diag_fit  <- readRDS(fit_path)
diag_inp  <- if (file.exists(inp_path)) readRDS(inp_path) else NULL

# ------------------------------------------------------------------------------
# 2. Observed CPUE -- from the EXACT vectors the model was handed
#
# `c` is fish per interview and `h` is person-hours (fishing_time *
# person_count_final; prep_inputs_bss line 155), so sum(c)/sum(h) is the
# observed catch rate in the same units as lambda_C. Taking these from the
# saved inputs rather than recomputing from dwg$interview matters: the raw
# interview table carries fishing_start_time / fishing_end_time, not a fishing
# time, and reproducing the prep chain here would risk a number that differs
# from what the model actually fitted -- which is precisely the comparison
# being made.
# ------------------------------------------------------------------------------

if (is.null(diag_inp) || is.null(diag_inp$c) || is.null(diag_inp$h)) {
  cli::cli_abort(c(
    "No saved Stan inputs at {.file {inp_path}}.",
    "i" = "01_fit_bss_bias.R writes these when SAVE_FITS is TRUE -- re-pull and re-run."
  ))
}
obs_fish  <- sum(diag_inp$c, na.rm = TRUE)
obs_hours <- sum(diag_inp$h, na.rm = TRUE)
obs_cpue  <- obs_fish / obs_hours

cli::cli_h2("Observed, from the interviews the model was given")
cli::cli_alert_info("Interviews with CPUE data (IntC): {diag_inp$IntC}")
cli::cli_alert_info("Fish of this group:  {round(obs_fish)}")
cli::cli_alert_info("Person-hours:        {round(obs_hours)}")
cli::cli_alert_info("Observed CPUE:       {signif(obs_cpue, 4)} fish/person-hour")

# ------------------------------------------------------------------------------
# 3. What the model says
# ------------------------------------------------------------------------------

d <- posterior::as_draws_df(diag_fit)

grab <- function(pat) {
  cols <- grep(pat, names(d), value = TRUE)
  if (length(cols) == 0) return(NULL)
  as.matrix(d[cols])
}

C <- grab("^C_sum$"); E <- grab("^E_sum$"); muC <- grab("^mu_C\\[")

C_med <- median(C, na.rm = TRUE)
E_med <- median(E, na.rm = TRUE)
model_cpue <- C_med / E_med

cli::cli_h2("What the model produced")
cli::cli_alert_info("C_sum median: {round(C_med)}")
cli::cli_alert_info("E_sum median: {round(E_med)}  (angler hours)")
cli::cli_alert_info("C_sum / E_sum = {signif(model_cpue, 3)} fish/hour  <- the model's CPUE")

if (!is.null(muC)) {
  cli::cli_alert_info(
    "exp(mu_C) by gear x section, median: {paste(signif(exp(apply(muC, 2, median)), 3), collapse = ', ')}"
  )
  cli::cli_alert_info("mu_C is the season-long catch-rate intercept; exp() is fish/hour.")
  cli::cli_alert_info("Prior on mu_mu_C is log(0.02), i.e. 0.02 fish/hour. init = \"0\" starts it at 1.0.")
}

# Convergence, so a bad number cannot be blamed on the wrong thing.
conv <- posterior::summarise_draws(
  posterior::subset_draws(d, variable = c("C_sum", "E_sum", "mu_C")),
  median = ~median(.x, na.rm = TRUE), rhat = posterior::rhat, ess_bulk = posterior::ess_bulk
)
cli::cli_h2("Convergence")
print(as.data.frame(conv), row.names = FALSE)

nan_frac <- mean(!is.finite(C))
if (nan_frac > 0) {
  cli::cli_alert_warning(
    "{round(100 * nan_frac)}% of C_sum draws are non-finite -- poisson_rng overflows past a rate of 2^30, \\
     which is itself evidence the catch rate is being driven far too high."
  )
}

# ------------------------------------------------------------------------------
# 4. The verdict
# ------------------------------------------------------------------------------

cpue_ratio <- model_cpue / obs_cpue

cli::cli_h2("Verdict")
cat(sprintf(
  "  observed CPUE   %10.4f fish/person-hour   (%d fish / %d person-hours)\n",
  obs_cpue, round(obs_fish), round(obs_hours)
))
cat(sprintf("  model CPUE      %10.4f fish/person-hour   (C_sum / E_sum)\n", model_cpue))
cat(sprintf("  ratio           %10.1fx\n\n", cpue_ratio))

if (!is.finite(cpue_ratio)) {
  cli::cli_alert_danger("Ratio is not finite -- see the non-finite C_sum warning above.")
} else if (cpue_ratio > 5) {
  cli::cli_alert_danger(
    "Model CPUE is {round(cpue_ratio)}x the observed rate. The catch half of the model is NOT \\
     being informed by the interviews, and NO season total from this fit is usable."
  )
  if (abs(model_cpue - 1) < 0.35) {
    cli::cli_alert_danger(
      "It is also sitting near 1.0 fish/hour, which is exp(0) -- the init value. mu_C has not \\
       moved off its starting point."
    )
  }
} else if (cpue_ratio > 2 || cpue_ratio < 0.5) {
  cli::cli_alert_warning(
    "Model CPUE is {signif(cpue_ratio, 2)}x the observed rate. Some divergence is expected -- the \\
     model expands to unsampled days and sections, and observed CPUE is a raw ratio with no \\
     weighting -- but this is more than that should explain."
  )
} else {
  cli::cli_alert_success(
    "Model CPUE is {signif(cpue_ratio, 2)}x the observed rate. The catch model is tracking the \\
     interviews, and the season total is behaving as it should."
  )
}

cli::cli_alert_info("The stanfit is in {.code diag_fit} for anything this did not print.")
