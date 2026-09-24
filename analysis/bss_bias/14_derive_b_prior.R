# ==============================================================================
# 14_derive_b_prior.R -- the b priors used in census-free fishery-years, and
# where each number came from
#
# One entry per prior in PRIORS below. The script reads the named fits from
# bss_b_summary.csv, applies the method described next, and writes every
# derivation -- inputs, intermediate quantities, result, and the value
# actually applied -- to bss_b_prior_ledger.csv. The fw_creel_bprior_*.Rmd
# files carry the applied values in their YAML and cite a prior_id from that
# ledger, so a number in a render can always be traced back to fits.
#
#   Rscript analysis/bss_bias/14_derive_b_prior.R
#
# METHOD -- random-effects predictive prior on log(b):
#
#   y_i   = log(median b) for each historical fit i
#   SE_i  = (log q97.5 - log q2.5) / (2 * 1.96)    within-year SE, log scale
#   tau^2 = max(0, var(y) - mean(SE^2))            between-year variance,
#                                                  sampling error removed
#   mu    = mean(y)
#   SE_mu = sqrt((tau^2 + mean(SE^2)) / n)         uncertainty in the centre
#   sigma = sqrt(tau^2 + SE_mu^2) * sigma_widen    predictive SD for a new year
#
#   prior: b ~ lognormal(mu, sigma)
#
# SE_i comes from the 95% interval because bss_b_summary.csv is all that is
# committed -- draws are not saved by default. For a roughly lognormal
# posterior this matches the SD of log(b) draws closely.
#
# tau_fits can differ from center_fits. That is how a data-poor reach borrows
# a between-year variance from a better-measured one: its centre and its own
# within-year SEs, but the neighbour's tau. The Stillaguamish North Fork does
# this -- its own tau^2 truncates to zero because its within-year SEs swamp
# the between-year difference, which says its data cannot measure its
# variability, not that it has none.
#
# n = 2 FOR EVERY PRIOR BELOW. With two years, tau^2 rests on one degree of
# freedom. The formula carries that uncertainty through SE_mu, but it cannot
# know whether the two years happened to fall close together. The script says
# so on every run rather than letting a two-point sigma pass as a measured
# one, and `contains_center_intervals` checks the result against the fits it
# came from: a prior whose 95% range does not contain a year's own posterior
# interval is narrower than that year's evidence alone.
#
# sigma_widen stays 1 unless a widening is agreed, and then it is recorded
# here, per prior, rather than typed into a YAML.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(cli); library(here)
})

OUT_DIR      <- here::here("analysis", "bss_bias", "outputs")
SUMMARY_PATH <- file.path(OUT_DIR, "bss_b_summary.csv")
LEDGER_PATH  <- file.path(OUT_DIR, "bss_b_prior_ledger.csv")

Z95 <- qnorm(0.975)

# ------------------------------------------------------------------------------
# The fits each prior is built from. Exact fishery_name keys from
# bss_b_summary.csv, scope tag included -- not a pattern, so the inputs are
# readable here without running anything.
# ------------------------------------------------------------------------------
STILLY_NF <- c("Stillaguamish salmon and gamefish 2024-25 [NF]",
               "Stillaguamish salmon and gamefish 2025-26 [NF]")
STILLY_MS <- c("Stillaguamish salmon and gamefish 2024-25 [MS]",
               "Stillaguamish salmon and gamefish 2025-26 [MS]")
SNO_MAIN  <- c("Snohomish fall salmon 2023 [SN_MAIN]",
               "Snohomish fall salmon 2024 [SN_MAIN]")

# applied_mu / applied_sigma: the value a render actually uses. NULL means
# "the derived value, rounded to 2 dp". Set explicitly when the applied value
# is a decision rather than the formula's output, and say why in `decision`.
PRIORS <- list(
  list(prior_id    = "stilly_nf_2026_vehicle",
       target      = "Stillaguamish salmon and gamefish 2026, North Fork",
       bias_type   = "vehicle",
       center_fits = STILLY_NF,
       tau_fits    = STILLY_MS,
       sigma_widen = 1,
       applied_mu  = 0.20, applied_sigma = 0.43,
       decision    = paste(
         "Centre and within-year SEs from the North Fork's own fits; tau borrowed",
         "from the mainstem, since the NF's own tau^2 truncates to 0. Applied",
         "values were derived before this script existed, from the fits as they",
         "stood then; later refits moved the result by about 0.01, so the applied",
         "pair is kept rather than silently replaced.")),

  list(prior_id    = "stilly_nf_2026_trailer",
       target      = "Stillaguamish salmon and gamefish 2026, North Fork",
       bias_type   = "trailer",
       center_fits = STILLY_NF,
       tau_fits    = STILLY_MS,
       sigma_widen = 1,
       applied_mu  = 0, applied_sigma = 1,
       decision    = paste(
         "Model default kept. The derived prior (mainstem trailer tau; MS went",
         "3.03 -> 0.86 between years) is wider in the upper tail than",
         "lognormal(0, 1), so history is no better than uninformative for this",
         "channel. Census boat anglers were under 20 in every fit behind it.")),

  list(prior_id    = "sno_main_2026_vehicle",
       target      = "Snohomish fall salmon 2026, mainstem",
       bias_type   = "vehicle",
       center_fits = SNO_MAIN,
       tau_fits    = SNO_MAIN,
       sigma_widen = 1,
       applied_mu  = NULL, applied_sigma = NULL,
       decision    = paste(
         "Derived as-is, no widening, pending co-manager review. The mainstem's",
         "own tau is used: 657 and 175 census boat anglers, prior_contraction",
         "0.997-0.999 in both fits.")),

  list(prior_id    = "sno_main_2026_trailer",
       target      = "Snohomish fall salmon 2026, mainstem",
       bias_type   = "trailer",
       center_fits = SNO_MAIN,
       tau_fits    = SNO_MAIN,
       sigma_widen = 1,
       applied_mu  = NULL, applied_sigma = NULL,
       decision    = paste(
         "Derived as-is, no widening, pending co-manager review. Unlike the",
         "Stillaguamish trailer term, the derived prior is far narrower than",
         "lognormal(0, 1), so history is informative here."))
)

# ------------------------------------------------------------------------------
# Derivation
# ------------------------------------------------------------------------------

if (!file.exists(SUMMARY_PATH)) {
  cli_abort("{.file {SUMMARY_PATH}} not found -- run 01_fit_bss_bias.R first.")
}
b_all <- read_csv(SUMMARY_PATH, show_col_types = FALSE)
summary_md5 <- unname(tools::md5sum(SUMMARY_PATH))

pick_fits <- function(fits, bias_type, prior_id) {
  d <- b_all |> filter(.data$fishery_name %in% fits, .data$bias_type == !!bias_type)
  missing <- setdiff(fits, d$fishery_name)
  if (length(missing) > 0) {
    cli_abort(c("{.val {prior_id}}: no {bias_type} row in bss_b_summary.csv for:",
                set_names(missing, rep("x", length(missing)))))
  }
  # One row per fit. A second would mean two runs of the same scope were
  # appended rather than replaced, and averaging them is not the method.
  if (anyDuplicated(d$fishery_name)) {
    cli_abort("{.val {prior_id}}: more than one {bias_type} row for a fit -- dedupe bss_b_summary.csv.")
  }
  d |>
    arrange(match(fishery_name, fits)) |>
    mutate(y  = log(median),
           se = (log(q97.5) - log(q2.5)) / (2 * Z95))
}

between_year_tau2 <- function(d) max(0, var(d$y) - mean(d$se^2))

derive_one <- function(p) {
  cen   <- pick_fits(p$center_fits, p$bias_type, p$prior_id)
  tau_d <- pick_fits(p$tau_fits,    p$bias_type, p$prior_id)
  n <- nrow(cen)
  if (n < 2 || nrow(tau_d) < 2) {
    cli_abort("{.val {p$prior_id}}: needs at least two fits for both the centre and tau.")
  }
  if (n < 3) {
    cli_alert_warning(
      "{.val {p$prior_id}}: centre from n = {n} years -- tau^2 rests on {n - 1} \\
       degree{?s} of freedom. Read sigma as provisional."
    )
  }

  tau2  <- between_year_tau2(tau_d)
  mu    <- mean(cen$y)
  se_mu <- sqrt((tau2 + mean(cen$se^2)) / n)
  sigma <- sqrt(tau2 + se_mu^2) * p$sigma_widen

  lo <- exp(mu - Z95 * sigma)
  hi <- exp(mu + Z95 * sigma)

  applied_source <- if (is.null(p$applied_mu)) "derived" else "set"
  applied_mu     <- p$applied_mu    %||% round(mu, 2)
  applied_sigma  <- p$applied_sigma %||% round(sigma, 2)

  tibble(
    prior_id       = p$prior_id,
    target         = p$target,
    bias_type      = p$bias_type,
    method         = "random-effects predictive, log scale",
    center_fits    = paste(p$center_fits, collapse = "; "),
    center_medians = paste(sprintf("%.4f", cen$median), collapse = "; "),
    center_se_log  = paste(sprintf("%.4f", cen$se), collapse = "; "),
    n_center       = n,
    tau_fits       = paste(p$tau_fits, collapse = "; "),
    tau_borrowed   = !setequal(p$center_fits, p$tau_fits),
    tau2_own       = between_year_tau2(cen),
    tau2           = tau2,
    tau            = sqrt(tau2),
    se_mu          = se_mu,
    sigma_widen    = p$sigma_widen,
    mu             = mu,
    sigma          = sigma,
    prior_median_b = exp(mu),
    prior_lo95_b   = lo,
    prior_hi95_b   = hi,
    # TRUE only if every centre fit's own 95% posterior interval sits inside
    # the prior's 95% range.
    contains_center_intervals = all(lo <= cen$q2.5 & cen$q97.5 <= hi),
    applied_mu     = applied_mu,
    applied_sigma  = applied_sigma,
    applied_source = applied_source,
    decision       = p$decision,
    summary_md5    = summary_md5,
    derived_at     = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )
}

cli_h1("14 -- b priors for census-free fishery-years")
ledger <- map_dfr(PRIORS, derive_one)

write_csv(ledger, LEDGER_PATH)

cli_h2("Derived vs applied")
ledger |>
  transmute(prior_id,
            derived = sprintf("lognormal(%.2f, %.2f)", mu, sigma),
            range95 = sprintf("%.2f-%.2f", prior_lo95_b, prior_hi95_b),
            tau     = sprintf("%.3f%s", tau, if_else(tau_borrowed, " (borrowed)", "")),
            applied = sprintf("lognormal(%.2f, %.2f) [%s]", applied_mu, applied_sigma, applied_source),
            contains = contains_center_intervals) |>
  print(n = Inf, width = Inf)

narrow <- ledger |> filter(!contains_center_intervals)
if (nrow(narrow) > 0) {
  cli_alert_warning(
    "Prior 95% range does not contain every source fit's own 95% interval: \\
     {.val {narrow$prior_id}}."
  )
}

# Applied values that are "set" rather than derived, and that the formula no
# longer lands near. Not an error -- NF vehicle is set on purpose -- but it
# should be seen on every run, not discovered later.
drift <- ledger |>
  filter(applied_source == "set",
         abs(applied_mu - mu) > 0.05 | abs(applied_sigma - sigma) > 0.05)
if (nrow(drift) > 0) {
  cli_alert_info(
    "Set values differ from today's derivation by more than 0.05 for \\
     {.val {drift$prior_id}} -- see `decision` in the ledger for why."
  )
}

cli_alert_success("Ledger written: {.file {LEDGER_PATH}} ({nrow(ledger)} prior{?s}).")
