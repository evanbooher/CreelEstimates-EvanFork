# ==============================================================================
# fishery_data.R
#
# The shared DATA LAYER for analysis/bss_bias/: everything between "a fishery
# name" and "the raw dwg tables for its estimation window", with nothing about
# Stan, priors or model inputs.
#
# Sourced by 01_fit_bss_bias.R (which fits `b`) and by 02b_survey_coverage.R
# (which only needs the raw effort/interview record). It exists so those two
# cannot drift apart on which window a fishery-year uses or how its dates are
# typed -- a difference there would silently make the coverage tables describe
# a different span than the fits.
#
# Requires, from the sourcing script: tidyverse, cli, here, creelutils, and
# creelutils' resolve_dates() when DATA_SOURCE is "internal".
#
# Provides:
#   DATA_SOURCE, DB_CONN           -- where data comes from, and the one connection
#   FISHERY_PARAMS, resolve_window(), window_skip_reason()
#   safe_name(), normalize_dwg_dates()
#   fetch_fishery_dwg()            -- fetch + normalise, cached to disk
#   disconnect_fishery_data()      -- closes DB_CONN if one was opened
# ==============================================================================

# DEFAULTS TO "external" (data.wa.gov) ON PURPOSE, not "internal": it makes
# this entire analysis reproducible by anyone with the repo and plain
# internet access -- no VPN, no WDFW DB credentials. Meeting participants
# (DFW + tribal technical staff) can re-run start to finish and get the same
# numbers, which is the difference between "here are my results" and "here
# is the analysis." Verified equivalent for BSS input purposes in check 0.A
# above. Switch to "internal" only if you specifically need DB-only data.
DATA_SOURCE <- "external"

# creelutils::fetch_data()'s documented conn = NULL auto-connect-when-internal
# path does not work against the currently installed creelutils version --
# confirmed locally: fetch_data(fishery_name=, data_source="internal") errors
# "`conn` is required" even after a successful DB connection message. Fix:
# open ONE connection here (matches the pattern in chore/multi-fishery-trip-
# summary's multi_fishery_creel_summary.R) and pass it explicitly into every
# fetch_data() call below, instead of reconnecting per fishery-year.
DB_CONN <- if (identical(DATA_SOURCE, "internal")) {
  cli::cli_alert_info("Connecting to internal DB...")
  creelutils::connect_creel_db()
} else {
  NULL  # ignored by fetch_data() when data_source == "external"
}

# ------------------------------------------------------------------------------
# Estimation window: the one remaining VPN-only dependency, and how an
# external run gets around it.
#
# resolve_dates() reads creelutils::fishery_lut(), which is internal-DB only --
# there is no data.wa.gov equivalent, and the public Socrata fishery registry
# (vkjc-s5u8) carries fishery_name but no start/end dates. So under
# DATA_SOURCE = "external" we read the windows captured by
# 00b_capture_fishery_params.R instead of querying.
#
# That is not only a VPN workaround. resolve_dates() computes
#   resolved_end = min(fishery_end_date, Sys.Date() - 1)
# so for an IN-SEASON fishery the window depends on what day the script runs.
# Reading a committed capture pins the window, making two runs on different
# days comparable -- which the live lookup cannot guarantee.
# ------------------------------------------------------------------------------

FISHERY_PARAMS_PATH <- here::here("analysis", "bss_bias", "lookup", "fishery_params.csv")

FISHERY_PARAMS <- if (file.exists(FISHERY_PARAMS_PATH)) {
  read_csv(FISHERY_PARAMS_PATH, show_col_types = FALSE)
} else {
  NULL
}

if (!identical(DATA_SOURCE, "internal")) {
  if (is.null(FISHERY_PARAMS)) {
    cli::cli_abort(c(
      "DATA_SOURCE is {.val {DATA_SOURCE}} but {.file {FISHERY_PARAMS_PATH}} does not exist.",
      "x" = "Estimation windows cannot be resolved without it (fishery_lut is VPN-only).",
      "i" = "Run {.file analysis/bss_bias/00b_capture_fishery_params.R} once with VPN and commit the result,",
      "i" = "or set {.code DATA_SOURCE <- \"internal\"} if you have DB access."
    ))
  }
  cli::cli_alert_info(
    "Using estimation windows captured {.val {FISHERY_PARAMS$captured_at[1]}} \\
     ({nrow(FISHERY_PARAMS)} fishery-year(s)) -- no VPN required."
  )
}

# Single entry point for "what window does this fishery-year use?", so the
# internal and external paths cannot drift apart. Returns NULL (not an error)
# when unresolvable; callers decide whether that is a skip.
resolve_window <- function(fishery_name) {
  if (identical(DATA_SOURCE, "internal")) {
    return(tryCatch(resolve_dates(fishery_name, "", "", conn = DB_CONN), error = function(e) NULL))
  }
  row <- FISHERY_PARAMS |>
    filter(fishery_name == .env$fishery_name, resolve_status == "ok")
  if (nrow(row) != 1) return(NULL)
  list(est_date_start = as.character(row$est_date_start[1]),
       est_date_end   = as.character(row$est_date_end[1]))
}

# Message for the skip when resolve_window() comes back empty -- names the
# actual remedy rather than just reporting the symptom.
window_skip_reason <- function() {
  if (identical(DATA_SOURCE, "internal")) {
    "Could not resolve estimation window from fishery_lut."
  } else {
    paste0("No captured estimation window for this fishery in lookup/fishery_params.csv. ",
           "Re-run 00b_capture_fishery_params.R with VPN to capture it, or set ",
           "DATA_SOURCE <- \"internal\" if you have DB access.")
  }
}

safe_name <- function(x) stringr::str_replace_all(x, "[^[:alnum:]]", "_")

# Coerce event_date to Date on every dwg table that has one, whatever the
# fetch handed back (Date, POSIXct, character, or an untyped column from a
# zero-row result). substr() to 10 characters first so an ISO timestamp
# ("2021-05-19T00:00:00.000") parses as cleanly as a bare date; the whole
# thing is a no-op on a column that is already Date.
normalize_dwg_dates <- function(dwg) {
  for (tbl in names(dwg)) {
    if (is.data.frame(dwg[[tbl]]) && "event_date" %in% names(dwg[[tbl]])) {
      dwg[[tbl]]$event_date <- as.Date(substr(as.character(dwg[[tbl]]$event_date), 1, 10))
    }
  }
  dwg
}


# ------------------------------------------------------------------------------
# Cached fetch
# ------------------------------------------------------------------------------
# One fishery-year's raw tables are a few MB and the external fetch is the
# slowest step that is not Stan. Caching them to disk means the coverage pass
# can iterate over all 28 fishery-years repeatedly -- while a figure is being
# tuned, say -- without re-hitting data.wa.gov each time, and a re-run of 01
# after a code change skips straight to prep.
#
# The cache key includes the RESOLVED WINDOW, not just the name: an in-season
# fishery re-captured by 00b gets a different window, and a stale cache keyed
# on name alone would silently serve the old span. Keyed on data source too,
# since internal and external need not return identical tables.
#
# The cache is disposable -- delete outputs/cache/dwg/ to force a refetch.

DWG_CACHE_DIR <- here::here("analysis", "bss_bias", "outputs", "cache", "dwg")
USE_DWG_CACHE <- TRUE

fetch_fishery_dwg <- function(fishery_name, est_dates, use_cache = USE_DWG_CACHE) {
  key <- paste0(safe_name(fishery_name), "__", DATA_SOURCE, "__",
                est_dates$est_date_start, "_", est_dates$est_date_end, ".rds")
  path <- file.path(DWG_CACHE_DIR, key)

  if (use_cache && file.exists(path)) {
    cli::cli_alert_info("  Using cached fetch ({.file {basename(path)}}).")
    return(readRDS(path))
  }

  dwg <- creelutils::fetch_data(conn = DB_CONN, fishery_name = fishery_name,
                                data_source = DATA_SOURCE)
  # Normalise BEFORE caching, so every consumer gets Date-typed event_date and
  # the cache cannot hand back the character columns the external path returns.
  dwg <- normalize_dwg_dates(dwg)

  if (use_cache) {
    dir.create(DWG_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
    saveRDS(dwg, path)
  }
  dwg
}

disconnect_fishery_data <- function() {
  if (!is.null(DB_CONN)) try(DBI::dbDisconnect(DB_CONN), silent = TRUE)
}
