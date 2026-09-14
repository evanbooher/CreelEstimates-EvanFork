# ==============================================================================
# 01b_launch_jobs.R
#
# Purpose:
#   Launch the 01b_fit_catch_groups.R fits as parallel background processes
#   FROM THE R CONSOLE, so no shell is involved.
#
#   RStudio's Terminal tab is Git Bash on most Windows installs, not
#   PowerShell, and `Rscript` is frequently not on PATH there. Both problems
#   disappear when R launches the processes itself: R always knows where its
#   own Rscript lives (R.home("bin")).
#
#   Source this file, or select-all and run it, from the R Console.
#
# Concurrency:
#   FIT_CONFIGS$quick in 01_fit_bss_bias.R uses 2 chains on 2 cores, so each
#   job occupies 2 cores. MAX_CONCURRENT is set for an 8-core machine. Raise it
#   only if you have the cores -- oversubscribing makes every fit slower, not
#   just the extra ones.
#
# Output:
#   Each job writes its own log under outputs/logs/, so four runs do not
#   interleave into one unreadable stream. Watch any of them with
#   job_tail() below.
# ==============================================================================

JOBS <- list(
  c(group = "chinook_all",  fisheries = "Snohomish"),
  c(group = "coho_harvest", fisheries = "Snohomish"),
  c(group = "chinook_all",  fisheries = "Stillaguamish"),
  c(group = "coho_harvest", fisheries = "Stillaguamish")
)

# "latest" = most recent year per fishery (4 fits, ~30 min wall clock).
# "all"    = every year in the discovery CSV (~18 fits, hours). Switch to "all"
#            once the first pass has proven itself.
YEARS_MODE     <- "latest"
MAX_CONCURRENT <- 4

# ------------------------------------------------------------------------------

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required. install.packages(\"here\")")
}

ROOT    <- here::here()
SCRIPT  <- here::here("analysis", "bss_bias", "01b_fit_catch_groups.R")
LOG_DIR <- here::here("analysis", "bss_bias", "outputs", "logs")
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(SCRIPT)) stop("Not found: ", SCRIPT)

# R always knows its own Rscript, whatever PATH says. file.path() +
# .Platform$r_arch keeps this correct on Windows multi-arch installs.
RSCRIPT <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
if (!file.exists(RSCRIPT)) stop("Could not locate Rscript at: ", RSCRIPT)

job_log <- function(job) {
  file.path(LOG_DIR, sprintf("01b_%s_%s.log", job[["group"]], job[["fisheries"]]))
}

launch_one <- function(job) {
  log_path <- job_log(job)
  # Truncate any previous log so a re-run is not read as still-running output.
  cat("", file = log_path)
  system2(
    RSCRIPT,
    args   = c(shQuote(SCRIPT), job[["group"]], shQuote(job[["fisheries"]]), YEARS_MODE),
    wait   = FALSE,          # background: the Console returns immediately
    stdout = log_path,
    stderr = log_path        # one file per job, both streams
  )
  message(sprintf("launched: %-13s %-14s -> %s",
                  job[["group"]], job[["fisheries"]], basename(log_path)))
  invisible(log_path)
}

if (length(JOBS) > MAX_CONCURRENT) {
  warning(sprintf(
    "%d jobs but MAX_CONCURRENT is %d. All are launched at once -- this script does not queue. Split the list by hand if that matters.",
    length(JOBS), MAX_CONCURRENT
  ))
}

message("Rscript:  ", RSCRIPT)
message("Project:  ", ROOT)
message("Years:    ", YEARS_MODE)
message("Logs:     ", LOG_DIR, "\n")

invisible(lapply(JOBS, launch_one))

message("\nAll jobs launched in the background. The Console is free.")
message("Watch one with:   job_tail(1)")
message("Check them all:   job_status()")
message("\n07_catch_sensitivity.R does NOT wait on these -- run it now:")
message("  Rscript analysis/bss_bias/07_catch_sensitivity.R")

# ------------------------------------------------------------------------------
# Helpers, left in the global environment for use after this script returns
# ------------------------------------------------------------------------------

job_tail <- function(i = 1, n = 25) {
  p <- job_log(JOBS[[i]])
  if (!file.exists(p)) return(message("no log yet: ", p))
  cat(tail(readLines(p, warn = FALSE), n), sep = "\n")
}

# A job is "done" once 01b prints its closing line; anything else with recent
# output is still running. Deliberately crude -- it reads the log, it does not
# track the OS process, so a crashed job shows as "running" until you look.
job_status <- function() {
  data.frame(
    group     = vapply(JOBS, function(j) j[["group"]], ""),
    fisheries = vapply(JOBS, function(j) j[["fisheries"]], ""),
    lines     = vapply(JOBS, function(j) {
      p <- job_log(j); if (file.exists(p)) length(readLines(p, warn = FALSE)) else 0L
    }, integer(1)),
    done      = vapply(JOBS, function(j) {
      p <- job_log(j)
      if (!file.exists(p)) return(FALSE)
      any(grepl("will pick up", readLines(p, warn = FALSE), fixed = TRUE))
    }, logical(1)),
    modified  = vapply(JOBS, function(j) {
      p <- job_log(j)
      if (file.exists(p)) format(file.info(p)$mtime, "%H:%M:%S") else NA_character_
    }, ""),
    stringsAsFactors = FALSE
  )
}
