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
#   Keeps the full stanfit in outputs/fits/ (gitignored) so the object is
#   available afterwards as `diag_fit` for anything this script does not print.
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
# 1. Observed CPUE, straight from the interviews -- no model involved
# ------------------------------------------------------------------------------

est_dates <- resolve_window(DIAG_FISHERY)
if (is.null(est_dates)) cli::cli_abort("No estimation window for {.val {DIAG_FISHERY}}.")
dwg <- fetch_fishery_dwg(DIAG_FISHERY, est_dates)

hit <- match_catch_group(dwg$catch, grp)
obs_fish <- sum(suppressWarnings(as.numeric(as.character(hit$fish_count))), na.rm = TRUE)

# Interview fishing time. Column naming has varied, so take the first that is
# present rather than assuming one.
time_col <- intersect(c("fishing_time_total", "fishing_time", "trip_time_total"),
                      names(dwg$interview))
if (length(time_col) == 0) {
  cli::cli_abort("No fishing-time column found in dwg$interview: {.val {names(dwg$interview)}}")
}
obs_hours <- sum(suppressWarnings(as.numeric(dwg$interview[[time_col[1]]])), na.rm = TRUE)
obs_cpue  <- obs_fish / obs_hours

cli::cli_h2("Observed, from interviews only")
cli::cli_alert_info("Fish of this group:   {round(obs_fish)}")
cli::cli_alert_info("Angler hours ({time_col[1]}): {round(obs_hours)}")
cli::cli_alert_info("Observed CPUE:        {signif(obs_cpue, 3)} fish/hour")

# ------------------------------------------------------------------------------
# 2. Fit, reusing 01's machinery exactly -- same preps, same preflight, same
#    priors. Nothing here is a parallel implementation.
# ------------------------------------------------------------------------------

CATCH_BASELINE_ONLY <- TRUE
RUN_CATCH_GROUP     <- grp
ONLY_FISHERIES      <- DIAG_FISHERY
FIT_CONFIG_NAME     <- DIAG_CONFIG
SAVE_FITS           <- TRUE        # keep the stanfit so it can be interrogated

source(here::here("analysis", "bss_bias", "01_fit_bss_bias.R"), local = FALSE)

fit_path <- file.path(here::here("analysis", "bss_bias", "outputs", "fits"),
                      paste0(safe_name(DIAG_FISHERY), ".rds"))
if (!file.exists(fit_path)) {
  cli::cli_abort("No stanfit at {.file {fit_path}} -- the fit did not complete. Read the ledger.")
}
diag_fit <- readRDS(fit_path)

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
  "  observed CPUE   %10.4f fish/hour   (%d fish / %d hours, interviews only)\n",
  obs_cpue, round(obs_fish), round(obs_hours)
))
cat(sprintf("  model CPUE      %10.4f fish/hour   (C_sum / E_sum)\n", model_cpue))
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
