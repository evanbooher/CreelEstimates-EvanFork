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

# One sanitised key per job, shared by the log and PID filenames -- and by
# 01b_fit_catch_groups.R, which derives the same name for its own PID file.
# Sanitised because a fisheries filter is a REGEX: "Snohomish|Stillaguamish"
# is a legitimate value and a pipe is not a legal Windows filename character.
job_key <- function(job) {
  gsub("[^[:alnum:]]+", "_", paste(job[["group"]], job[["fisheries"]], sep = "_"))
}

job_log <- function(job) file.path(LOG_DIR, paste0("01b_", job_key(job), ".log"))
job_pid_file <- function(job) file.path(LOG_DIR, paste0("01b_", job_key(job), ".pid"))

launch_one <- function(job) {
  log_path <- job_log(job)
  # Truncate any previous log so a re-run is not read as still-running output,
  # and clear any stale PID file from a crashed run.
  cat("", file = log_path)
  unlink(job_pid_file(job))
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
message("Stop them with:   job_kill()")
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

# A job is "done" once 01b prints its closing line. `running` comes from the
# PID file and an actual liveness check, so a job that crashed shows as neither
# done nor running -- which is the case worth noticing, and the one an
# earlier log-only version of this could not distinguish.
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
    running   = vapply(seq_along(JOBS), function(k) pid_alive(job_pid(k)), logical(1)),
    modified  = vapply(JOBS, function(j) {
      p <- job_log(j)
      if (file.exists(p)) format(file.info(p)$mtime, "%H:%M:%S") else NA_character_
    }, ""),
    stringsAsFactors = FALSE
  )
}

# Is this PID actually alive? Used by both job_status() and job_kill() so a
# stale PID file from a crashed run is never reported as running, and never
# passed to a kill command that would then error.
pid_alive <- function(pid) {
  if (is.na(pid)) return(FALSE)
  if (.Platform$OS.type == "windows") {
    out <- suppressWarnings(system2(
      "tasklist", c("/FI", shQuote(sprintf("PID eq %s", pid)), "/NH"),
      stdout = TRUE, stderr = NULL
    ))
    any(grepl(as.character(pid), out, fixed = TRUE))
  } else {
    !inherits(try(tools::pskill(pid, 0), silent = TRUE), "try-error")
  }
}

job_pid <- function(i) {
  p <- job_pid_file(JOBS[[i]])
  if (!file.exists(p)) return(NA_integer_)
  suppressWarnings(as.integer(readLines(p, warn = FALSE)[1]))
}

# Stop the jobs THIS launcher started -- not every Rscript on the machine.
# `taskkill /IM Rscript.exe` would also take out any unrelated R job, and on a
# machine running a long backfill that is an expensive mistake.
#
#   job_kill()     all jobs
#   job_kill(2)    just job 2
job_kill <- function(i = seq_along(JOBS)) {
  for (k in i) {
    job <- JOBS[[k]]
    pid <- job_pid(k)
    label <- sprintf("%s / %s", job[["group"]], job[["fisheries"]])
    if (is.na(pid)) {
      message(sprintf("no PID file  %-30s (finished, or never started)", label))
      next
    }
    if (!pid_alive(pid)) {
      message(sprintf("not running  %-30s (pid %s -- stale file, removing)", label, pid))
      unlink(job_pid_file(job))
      next
    }
    ok <- if (.Platform$OS.type == "windows") {
      system2("taskkill", c("/F", "/PID", pid), stdout = FALSE, stderr = FALSE)
    } else {
      system2("kill", c("-9", pid), stdout = FALSE, stderr = FALSE)
    }
    if (identical(as.integer(ok), 0L)) {
      message(sprintf("killed       %-30s (pid %s)", label, pid))
      unlink(job_pid_file(job))
    } else {
      message(sprintf("KILL FAILED  %-30s (pid %s) -- end it in Task Manager", label, pid))
    }
  }
  invisible(NULL)
}
