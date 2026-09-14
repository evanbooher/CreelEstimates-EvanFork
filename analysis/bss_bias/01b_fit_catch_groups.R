# ==============================================================================
# 01b_fit_catch_groups.R
#
# Purpose:
#   Re-fit the BSS for a NAMED CATCH GROUP so that 07_catch_sensitivity.R can
#   express its results in fish instead of percent.
#
#   WHY THIS IS A SEPARATE SCRIPT FROM THE ANALYSIS:
#   `b` is INVARIANT to the catch group. Two independent reasons, both checked
#   against the code rather than assumed:
#
#     1. The interview set does not change. prep_dwg_interview_catch()
#        replicates every interview for each catch group and left-joins the
#        counts with replace_na(fish_count, 0) -- interviews with no fish of
#        that group are kept at zero, not dropped. So IntA, IntC and every
#        effort input are identical across groups.
#
#     2. The effort and catch sub-models share no parameters. In the Stan model
#        block, `b`, R_V, R_T, p_I, mu_E and B1 sit on one side; mu_C, omega_C
#        and r_C on the other. B1 appears in lambda_E_S only, not lambda_C_S.
#        They meet solely in generated quantities (line 240). The posterior
#        factorises.
#
#   So a re-fit reproduces the same `b` up to Monte Carlo error, and buys
#   exactly one thing: C_sum in fish. That is why this run sets
#   CATCH_BASELINE_ONLY = TRUE -- writing its `b` into bss_b_summary.csv would
#   duplicate every series and double the year counts in 06's T2.
#
#   The re-fit does give a free check on the claim: bss_b_invariance_check.csv
#   holds each run's `b`, to be compared against the stored value. A material
#   difference falsifies the reasoning above and should stop the analysis.
#
# Usage -- SERIAL, from the RStudio Console (the simple path, live output):
#
#   source("analysis/bss_bias/01b_fit_catch_groups.R")
#
#   That runs both catch groups against Snohomish + Stillaguamish, most recent
#   year each, one fit at a time with everything printing to the Console. To
#   narrow it, set any of these first (they persist between sources -- reset
#   them or restart R to change a run):
#
#     GROUP_KEY  <- "chinook_all"          # chinook_all | coho_harvest | all
#     FISHERY_RE <- "Snohomish"            # regex over fishery_name
#     YEARS_MODE <- "latest"               # latest | all
#
# Usage -- from a shell, optionally several at once:
#   Rscript analysis/bss_bias/01b_fit_catch_groups.R <group> [fishery-regex] [years]
#
#   Each invocation is one OS process. FIT_CONFIGS$quick uses 2 chains on 2
#   cores, so on an 8-core machine run FOUR invocations concurrently and no
#   more -- see the run block at the bottom of this file.
#
#   CONCURRENCY: 01's append_csv_row() is a read-modify-write over the whole
#   CSV, which is not safe for several processes at once. Each run therefore
#   sets OUTPUT_TAG so it writes its own files (..._<group>_<fisheries>.csv);
#   07_catch_sensitivity.R globs and merges them. Do not remove the tag to
#   "tidy up" the output directory -- concurrent runs will corrupt each
#   other's reads silently before they fail loudly.
#
# Outputs (appended, one row per fishery-year x catch group):
#   bss_catch_baseline.csv        -- C_sum / E_sum posterior summaries  <- 07 reads this
#   bss_b_invariance_check.csv    -- the `b` this run produced, for the check above
#   bss_b_fit_ledger.csv          -- as usual, via 01's own ledger
# ==============================================================================

library(cli)
library(here)

# Settings come from, in order of precedence: variables already in the global
# environment, then command-line arguments, then the defaults. The first branch
# is what makes this sourceable straight from the RStudio Console --
#
#   source("analysis/bss_bias/01b_fit_catch_groups.R")           # everything, serial
#   GROUP_KEY <- "coho_harvest"; FISHERY_RE <- "Snohomish"
#   source("analysis/bss_bias/01b_fit_catch_groups.R")           # one combination
#
# NOTE: those variables PERSIST between sources. Reset them (or restart R)
# before a run you want to use different settings.
args <- commandArgs(trailingOnly = TRUE)
if (!exists("GROUP_KEY",  inherits = FALSE)) GROUP_KEY  <- if (length(args) >= 1) args[[1]] else "all"
if (!exists("FISHERY_RE", inherits = FALSE)) FISHERY_RE <- if (length(args) >= 2) args[[2]] else "Snohomish|Stillaguamish"
if (!exists("YEARS_MODE", inherits = FALSE)) YEARS_MODE <- if (length(args) >= 3) args[[3]] else "latest"

group_key  <- GROUP_KEY
fishery_re <- FISHERY_RE
years_mode <- YEARS_MODE

# ------------------------------------------------------------------------------
# The catch groups. Each field is a str_detect PATTERN, so alternation works and
# "NA" matches the literal string prep_dwg_interview_catch() coerces NA to.
# ------------------------------------------------------------------------------

CATCH_GROUPS <- list(
  # Any encounter of Chinook -- adults and jacks, every mark status including
  # unknown and unrecorded, kept and released.
  chinook_all = list(
    species    = "Chinook",
    life_stage = "Adult|Jack",
    fin_mark   = "UM|AD|UNK|NA",
    fate       = "Released|Kept"
  ),
  # Coho harvest. Note this is NOT the pipeline default
  # (Coho_Adult_AD|UM_Kept) -- it adds jacks, so it is a genuinely new fit.
  coho_harvest = list(
    species    = "Coho",
    life_stage = "Adult|Jack",
    fin_mark   = "UM|AD",
    fate       = "Kept"
  )
)

keys <- if (identical(group_key, "all")) names(CATCH_GROUPS) else group_key
bad  <- setdiff(keys, names(CATCH_GROUPS))
if (length(bad) > 0) {
  cli::cli_abort(c("Unknown catch group {.val {bad}}.",
                   "i" = "Available: {.val {names(CATCH_GROUPS)}}"))
}

# ------------------------------------------------------------------------------
# Which fishery-years. Read from the same discovery CSV 01 uses, so this script
# never invents a fishery name that the pipeline does not recognise.
# ------------------------------------------------------------------------------

OUT_DIR <- here::here("analysis", "bss_bias", "outputs")
disc_path <- file.path(OUT_DIR, "fishery_discovery_target.csv")
if (!file.exists(disc_path)) {
  cli::cli_abort("{.file {disc_path}} not found -- run 00_discover_fisheries.R first.")
}

disc <- readr::read_csv(disc_path, show_col_types = FALSE) |>
  dplyr::filter(basin_match == "target", include_in_run,
                stringr::str_detect(fishery_name, fishery_re))

if (nrow(disc) == 0) {
  cli::cli_abort("No fishery in the discovery CSV matches {.val {fishery_re}}.")
}

# `latest` keeps only the most recent year of each fishery TYPE -- enough to
# anchor this year's scale in fish, and four fits instead of eighteen. Extend
# with `all` once the clock allows.
if (identical(years_mode, "latest")) {
  disc <- disc |>
    dplyr::mutate(ftype = stringr::str_squish(
      stringr::str_replace(fishery_name, "\\d{4}(-\\d{2,4})?", " "))) |>
    dplyr::group_by(ftype) |>
    dplyr::slice_max(order_by = fishery_name, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
}

targets <- sort(unique(disc$fishery_name))

# Announce this process's PID so 01b_launch_jobs.R's job_kill() can stop these
# runs specifically, rather than taskkill-ing every Rscript on the machine.
# Written here rather than captured by the launcher because system2(wait =
# FALSE) returns no PID. Removed on clean completion; a stale file left by a
# crash is harmless, job_kill() skips PIDs that are no longer running.
JOB_PID_FILE <- file.path(
  here::here("analysis", "bss_bias", "outputs", "logs"),
  sprintf("01b_%s.pid", gsub("[^[:alnum:]]+", "_", paste(group_key, fishery_re, sep = "_")))
)
dir.create(dirname(JOB_PID_FILE), recursive = TRUE, showWarnings = FALSE)
writeLines(as.character(Sys.getpid()), JOB_PID_FILE)

cli::cli_h1("01b -- catch-group baseline fits")
cli::cli_alert_info("Catch group{?s}: {.val {keys}}")
# Two calls on purpose: a cli string may carry only ONE quantity when it also
# carries a {?s} plural marker, and "{length(targets)} ... {targets}" gives it
# two -- which aborts with "Multiple quantities for pluralization" rather than
# degrading. Same trap as the section-restriction message in 01.
cli::cli_alert_info("Queued {length(targets)} fishery-year{?s}.")
cli::cli_alert_info("{.val {targets}}")

# ------------------------------------------------------------------------------
# Run. Setting the knobs BEFORE sourcing 01 is what makes this work: 01 defines
# each of them with the "only if the caller has not" idiom, so these survive.
# ------------------------------------------------------------------------------

CATCH_BASELINE_ONLY <- TRUE
ONLY_FISHERIES      <- targets
FIT_CONFIG_NAME     <- "quick"

for (k in keys) {
  cli::cli_h2("Catch group: {k}")
  RUN_CATCH_GROUP <- CATCH_GROUPS[[k]]
  # Every CSV this process writes gets its own name. append_csv_row() in 01 is
  # a read-modify-write and is not safe for concurrent processes -- without
  # this, two jobs racing on bss_b_comparability_raw.csv leave one of them
  # reading a truncated file and failing with "object 'fishery_name' not
  # found". Merged back together by 07_catch_sensitivity.R, which globs.
  # Per-process output files exist only to stop CONCURRENT runs corrupting each
  # other's read-modify-write in append_csv_row(). A serial run from the Console
  # has no such problem, so it writes the normal filenames and leaves the output
  # directory tidy. 07 reads either shape.
  OUTPUT_TAG <- if (interactive()) "" else gsub("[^[:alnum:]]+", "_", paste(k, fishery_re, sep = "_"))
  source(here::here("analysis", "bss_bias", "01_fit_bss_bias.R"), local = FALSE)
}

unlink(JOB_PID_FILE)

cli::cli_alert_success(
  "Done. 07_catch_sensitivity.R will pick up {.file bss_catch_baseline.csv} automatically."
)

# ------------------------------------------------------------------------------
# RUN BLOCK -- four concurrent processes on 8 cores (2 chains x 2 cores each).
#
#   PowerShell:
#     Start-Process Rscript "analysis/bss_bias/01b_fit_catch_groups.R chinook_all  Snohomish      latest"
#     Start-Process Rscript "analysis/bss_bias/01b_fit_catch_groups.R coho_harvest Snohomish      latest"
#     Start-Process Rscript "analysis/bss_bias/01b_fit_catch_groups.R chinook_all  Stillaguamish  latest"
#     Start-Process Rscript "analysis/bss_bias/01b_fit_catch_groups.R coho_harvest Stillaguamish  latest"
#
#   bash / macOS / Linux:
#     Rscript analysis/bss_bias/01b_fit_catch_groups.R chinook_all  Snohomish     latest &
#     Rscript analysis/bss_bias/01b_fit_catch_groups.R coho_harvest Snohomish     latest &
#     Rscript analysis/bss_bias/01b_fit_catch_groups.R chinook_all  Stillaguamish latest &
#     Rscript analysis/bss_bias/01b_fit_catch_groups.R coho_harvest Stillaguamish latest &
#     wait
#
# Concurrent appends are handled by OUTPUT_TAG (see the header): each process
# writes its own CSVs and 07 merges them. The per-fishery DWG cache is read-only
# by the time these run -- if a fishery-year has never been fetched, run it once
# on its own first so the two catch-group jobs for it do not race on writing the
# same cache file.
# ------------------------------------------------------------------------------
