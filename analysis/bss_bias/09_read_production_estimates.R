# ==============================================================================
# 09_read_production_estimates.R
#
# Purpose:
#   Read BSS estimates produced by the PRODUCTION pipeline
#   (template_scripts/fw_creel.Rmd, writing into fishery_analyses/) and hand
#   them to this analysis in the shape 06/07 already expect.
#
#   This replaces 01b_fit_catch_groups.R as the source of season totals. The
#   fits are run by the real pipeline, so the numbers are the agency's own --
#   not a reimplementation that has to be argued for.
#
# What it reads, per analysis folder:
#   outputs/estimates_bss.rds   named list, one element per est_cg, each with
#                                 $draws          extract(stanfit) -- C_sum,
#                                                 E_sum, b, and everything else
#                                 $season_results catch/effort quantiles
#                                 $summary        summary(fit)$summary, carries Rhat
#   inputs/inputs_bss.rds       named list per est_cg; `c` (fish per interview)
#                                 and `h` (person-hours) if present
#
# Folder layout (see setup_analysis_structure.R):
#   fishery_analyses/<project>/<fishery_name>/<identifier>_<id4>_.../
#       inputs/  outputs/  figures/
#
#   The leaf folder name carries a per-run uuid fragment, so a fishery-year
#   re-run lands in a NEW folder rather than overwriting. This script keeps the
#   most recently modified folder per (fishery_name, est_cg) and reports the
#   ones it passed over, rather than silently picking one.
#
# Usage (R Console):
#   source("analysis/bss_bias/09_read_production_estimates.R")
#
#   Override before sourcing:
#     PROD_ROOT    <- "some/other/fishery_analyses"
#     PROD_PROJECT <- "bss_bias"      # NULL (default) = every project folder
#     PROD_ALL_RUNS <- TRUE           # keep every run, not just the newest
#
# Outputs:
#   bss_catch_baseline.csv        <- 07_catch_sensitivity.R reads this
#   bss_production_index.csv      what was found, where, and when it was written
#   b_draws_production/<name>.rds b[1]/b[2] draws from the production fits,
#                                 for checking against outputs/b_draws/
# ==============================================================================

library(tidyverse)
library(cli)
library(here)

OUT_DIR    <- here::here("analysis", "bss_bias", "outputs")
PROD_DRAWS <- file.path(OUT_DIR, "b_draws_production")
dir.create(OUT_DIR,    recursive = TRUE, showWarnings = FALSE)
dir.create(PROD_DRAWS, recursive = TRUE, showWarnings = FALSE)

if (!exists("PROD_ROOT",     inherits = FALSE)) PROD_ROOT     <- here::here("fishery_analyses")
if (!exists("PROD_PROJECT",  inherits = FALSE)) PROD_PROJECT  <- NULL
if (!exists("PROD_ALL_RUNS", inherits = FALSE)) PROD_ALL_RUNS <- FALSE

safe_name <- function(x) stringr::str_replace_all(x, "[^[:alnum:]]", "_")

cli::cli_h1("09 -- read production BSS estimates")
if (!dir.exists(PROD_ROOT)) {
  cli::cli_abort(c(
    "{.path {PROD_ROOT}} does not exist.",
    "i" = "Run template_scripts/fw_creel.Rmd first, or set PROD_ROOT."
  ))
}
cli::cli_alert_info("Root: {.path {PROD_ROOT}}")

# ------------------------------------------------------------------------------
# Locate analysis folders
# ------------------------------------------------------------------------------

est_files <- list.files(PROD_ROOT, pattern = "^estimates_bss\\.rds$",
                        recursive = TRUE, full.names = TRUE)
if (length(est_files) == 0) {
  cli::cli_abort("No {.file estimates_bss.rds} found anywhere under {.path {PROD_ROOT}}.")
}

# .../<project>/<fishery_name>/<analysis_folder>/outputs/estimates_bss.rds
meta <- tibble(est_path = est_files) |>
  mutate(
    analysis_folder = basename(dirname(dirname(est_path))),
    fishery_name    = basename(dirname(dirname(dirname(est_path)))),
    project_name    = basename(dirname(dirname(dirname(dirname(est_path))))),
    inputs_path     = file.path(dirname(dirname(est_path)), "inputs", "inputs_bss.rds"),
    modified        = file.info(est_path)$mtime
  )

if (!is.null(PROD_PROJECT)) meta <- filter(meta, project_name %in% PROD_PROJECT)
if (nrow(meta) == 0) cli::cli_abort("No analyses matched PROD_PROJECT {.val {PROD_PROJECT}}.")

cli::cli_alert_info("Found {nrow(meta)} analysis folder{?s} across {n_distinct(meta$fishery_name)} fishery-year{?s}.")

# ------------------------------------------------------------------------------
# Pull one analysis folder
# ------------------------------------------------------------------------------

summarise_vec <- function(x) {
  x <- as.numeric(x)
  fin <- is.finite(x)
  if (!any(fin)) {
    return(list(mean = NA_real_, sd = NA_real_, median = NA_real_,
                q2.5 = NA_real_, q97.5 = NA_real_, n_finite_frac = 0, n_draws = length(x)))
  }
  list(
    mean   = mean(x[fin]), sd = stats::sd(x[fin]), median = stats::median(x[fin]),
    q2.5   = unname(stats::quantile(x[fin], 0.025)),
    q97.5  = unname(stats::quantile(x[fin], 0.975)),
    n_finite_frac = mean(fin), n_draws = length(x)
  )
}

# Rhat lives in summary(fit)$summary, rownames = parameter. Returned as NA
# rather than dropped when absent, so a missing diagnostic is visible.
rhat_of <- function(summ, par) {
  if (is.null(summ) || is.null(rownames(summ)) || !"Rhat" %in% colnames(summ)) return(NA_real_)
  if (!par %in% rownames(summ)) return(NA_real_)
  unname(summ[par, "Rhat"])
}

read_one <- function(row) {
  est <- try(readRDS(row$est_path), silent = TRUE)
  if (inherits(est, "try-error") || !is.list(est) || length(est) == 0) {
    cli::cli_alert_warning("Unreadable or empty: {.path {row$est_path}}")
    return(NULL)
  }
  inp <- if (file.exists(row$inputs_path)) {
    try(readRDS(row$inputs_path), silent = TRUE)
  } else NULL
  if (inherits(inp, "try-error")) inp <- NULL

  imap_dfr(est, function(e, ecg) {
    draws <- e$draws
    if (is.null(draws) || is.null(draws$C_sum)) {
      cli::cli_alert_warning("{row$fishery_name} / {ecg}: no C_sum draws -- skipped.")
      return(NULL)
    }
    cs <- summarise_vec(draws$C_sum)
    es <- if (!is.null(draws$E_sum)) summarise_vec(draws$E_sum) else summarise_vec(NA_real_)

    # Observed CPUE from the vectors the model was handed, when the inputs were
    # saved. sum(c)/sum(h) is fish per person-hour, the same units as lambda_C,
    # which makes it directly comparable to C_sum/E_sum.
    i <- if (!is.null(inp) && !is.null(inp[[ecg]])) inp[[ecg]] else NULL
    obs_fish  <- if (!is.null(i$c)) sum(i$c, na.rm = TRUE) else NA_real_
    obs_hours <- if (!is.null(i$h)) sum(i$h, na.rm = TRUE) else NA_real_

    # b from the production fit, kept so these can be checked against the
    # b series in outputs/b_draws/ rather than assumed to agree.
    b <- draws$b
    if (!is.null(b)) {
      bdf <- tibble(`b[1]` = as.numeric(b[, 1]),
                    `b[2]` = if (ncol(b) >= 2) as.numeric(b[, 2]) else NA_real_)
      saveRDS(bdf, file.path(PROD_DRAWS, paste0(safe_name(row$fishery_name), "__", safe_name(ecg), ".rds")))
    }

    tibble(
      fishery_name = row$fishery_name, est_cg = ecg, project_name = row$project_name,
      analysis_folder = row$analysis_folder, modified = row$modified,
      C_sum_mean = cs$mean, C_sum_sd = cs$sd, C_sum_median = cs$median,
      C_sum_q2.5 = cs$q2.5, C_sum_q97.5 = cs$q97.5,
      C_sum_n_finite_frac = cs$n_finite_frac,
      E_sum_mean = es$mean, E_sum_sd = es$sd, E_sum_median = es$median,
      E_sum_q2.5 = es$q2.5, E_sum_q97.5 = es$q97.5,
      E_sum_n_finite_frac = es$n_finite_frac,
      n_draws = cs$n_draws,
      C_sum_rhat = rhat_of(e$summary, "C_sum"),
      E_sum_rhat = rhat_of(e$summary, "E_sum"),
      b1_median = if (!is.null(b)) stats::median(as.numeric(b[, 1]), na.rm = TRUE) else NA_real_,
      b2_median = if (!is.null(b) && ncol(b) >= 2) stats::median(as.numeric(b[, 2]), na.rm = TRUE) else NA_real_,
      obs_fish = obs_fish, obs_person_hours = obs_hours,
      obs_cpue = obs_fish / obs_hours,
      model_cpue = cs$median / es$median
    )
  })
}

all_rows <- map_dfr(seq_len(nrow(meta)), ~read_one(meta[.x, ]))
if (nrow(all_rows) == 0) cli::cli_abort("No usable estimates found.")

# ------------------------------------------------------------------------------
# One row per (fishery_name, est_cg): newest run wins, older ones reported
# ------------------------------------------------------------------------------

if (!PROD_ALL_RUNS) {
  dupes <- all_rows |> count(fishery_name, est_cg) |> filter(n > 1)
  baseline <- all_rows |>
    arrange(fishery_name, est_cg, desc(modified)) |>
    distinct(fishery_name, est_cg, .keep_all = TRUE)
  if (nrow(dupes) > 0) {
    cli::cli_alert_warning(
      "{nrow(dupes)} fishery-year x catch group combination{?s} had more than one run; kept the newest."
    )
    all_rows |>
      semi_join(dupes, by = c("fishery_name", "est_cg")) |>
      arrange(fishery_name, est_cg, desc(modified)) |>
      select(fishery_name, est_cg, analysis_folder, modified, C_sum_median) |>
      print(n = Inf)
  }
} else {
  baseline <- all_rows
}

write_csv(baseline, file.path(OUT_DIR, "bss_catch_baseline.csv"))
write_csv(all_rows, file.path(OUT_DIR, "bss_production_index.csv"))
cli::cli_alert_success("Wrote {nrow(baseline)} baseline row{?s} to {.file bss_catch_baseline.csv}.")

# ------------------------------------------------------------------------------
# Sanity view -- the numbers worth looking at before anything downstream runs
# ------------------------------------------------------------------------------

cli::cli_h2("Season totals")
baseline |>
  mutate(group = str_extract(est_cg, "^[^_]+")) |>
  select(fishery_name, group, C_sum_median, E_sum_median,
         obs_fish, obs_cpue, model_cpue, C_sum_rhat, n_draws) |>
  arrange(fishery_name, group) |>
  print(n = Inf)

bad_rhat <- baseline |> filter(!is.na(C_sum_rhat), C_sum_rhat > 1.05)
if (nrow(bad_rhat) > 0) {
  cli::cli_h2("Rhat above 1.05 -- do not quote these")
  bad_rhat |> select(fishery_name, est_cg, C_sum_rhat, E_sum_rhat) |> print(n = Inf)
}

nonfinite <- baseline |> filter(C_sum_n_finite_frac < 1)
if (nrow(nonfinite) > 0) {
  cli::cli_h2("Non-finite C_sum draws present")
  cli::cli_alert_info("poisson_rng overflows past a rate of 2^30; the median is over the finite draws only.")
  nonfinite |> select(fishery_name, est_cg, C_sum_n_finite_frac) |> print(n = Inf)
}

cpue_off <- baseline |> filter(!is.na(obs_cpue), obs_cpue > 0, model_cpue / obs_cpue > 2)
if (nrow(cpue_off) > 0) {
  cli::cli_h2("Model CPUE more than 2x the observed rate")
  cli::cli_alert_info("Some divergence is expected -- the model expands to unsampled days and sections.")
  cpue_off |>
    mutate(ratio = model_cpue / obs_cpue) |>
    select(fishery_name, est_cg, obs_cpue, model_cpue, ratio) |>
    print(n = Inf)
}

cli::cli_alert_info("Next: {.code Rscript analysis/bss_bias/07_catch_sensitivity.R}")
