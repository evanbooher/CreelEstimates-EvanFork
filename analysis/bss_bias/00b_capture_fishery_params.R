# ==============================================================================
# 00b_capture_fishery_params.R
#
# Purpose:
#   Capture every pipeline parameter that currently requires VPN / internal-DB
#   access into a committed lookup file, so that a subsequent run with
#   DATA_SOURCE <- "external" is genuinely reproducible start to finish by
#   anyone with the repo and plain internet access -- no VPN, no WDFW database
#   credentials.
#
#   RUN THIS ONCE, WITH VPN, and commit the resulting lookup file. Meeting
#   participants never run this script; they run 00 -> 01 -> 02 ... against
#   the captured lookup.
#
# ------------------------------------------------------------------------------
# WHAT IS ACTUALLY VPN-DEPENDENT
#
#   Only one thing, after the switch of DATA_SOURCE to "external":
#
#     resolve_dates() -> creelutils::fishery_lut()  -- the estimation window
#       (fishery_start_date / fishery_end_date). fetch_data() has a public
#       data.wa.gov path; fishery_lut() does not, and the public Socrata
#       fishery registry (vkjc-s5u8, what 00_discover_fisheries.R reads)
#       carries fishery_name but NO start/end dates.
#
#   Everything else downstream -- prep_days(), the prep_dwg_* family,
#   prep_inputs_bss(), fit_bss() -- operates on data already fetched, or on
#   hard-coded constants in 01_fit_bss_bias.R, and needs no database.
#
# ------------------------------------------------------------------------------
# WHY CAPTURE THE *RESOLVED* WINDOW, NOT THE RAW LOOKUP-TABLE COLUMNS
#
#   resolve_dates() does not pass fishery_lut's dates through unchanged:
#
#     resolved_end <- min(as.Date(fishery_row$fishery_end_date), Sys.Date() - 1)
#
#   For an IN-SEASON fishery (end date still in the future -- e.g. the 2026
#   fishery-years) the resolved window therefore depends on WHAT DAY THE
#   SCRIPT IS RUN. Run it today and tomorrow and you get different windows,
#   hence different data, hence a different `b`. That is a reproducibility
#   hole independent of VPN, and it silently makes two runs incomparable.
#
#   So this script records resolve_dates()'s OUTPUT, frozen at capture time,
#   rather than the raw lut columns. Everyone re-running the pipeline --
#   including future you -- then gets exactly the window this capture used.
#   `captured_at` is written alongside so the freeze date is never in doubt.
#
# Usage:
#   Rscript analysis/bss_bias/00b_capture_fishery_params.R
#   REQUIRES VPN / internal DB access. Run after 00_discover_fisheries.R and
#   after hand-editing include_in_run, since it captures the fisheries that
#   01_fit_bss_bias.R will actually attempt.
#
# Outputs:
#   analysis/bss_bias/lookup/fishery_params.csv  -- TRACKED IN GIT; this file
#     is the whole point of the script. Columns: fishery_name, est_date_start,
#     est_date_end, n_days_window, resolve_status, resolve_note, captured_at.
#
#   Rows for fisheries not in this run's target set are PRESERVED, so
#   re-running after a scope change tops the file up rather than truncating it.
# ==============================================================================

library(tidyverse)
library(cli)
library(here)
library(creelutils)
library(rlang)

walk(list.files(here("R_functions"), full.names = TRUE), source)

OUT_DIR     <- here::here("analysis", "bss_bias", "outputs")
LOOKUP_DIR  <- here::here("analysis", "bss_bias", "lookup")
LOOKUP_PATH <- file.path(LOOKUP_DIR, "fishery_params.csv")
dir.create(LOOKUP_DIR, recursive = TRUE, showWarnings = FALSE)

discovery_path <- file.path(OUT_DIR, "fishery_discovery_target.csv")
if (!file.exists(discovery_path)) {
  cli::cli_abort("Run 00_discover_fisheries.R first -- {.file {discovery_path}} not found.")
}

target_fisheries <- read_csv(discovery_path, show_col_types = FALSE) |>
  filter(basin_match == "target", include_in_run) |>
  pull(fishery_name_raw) |>
  unique()

cli::cli_alert_info("Capturing estimation windows for {length(target_fisheries)} fishery-year(s).")

cli::cli_alert_info("Connecting to internal DB...")
conn <- creelutils::connect_creel_db()
on.exit(try(DBI::dbDisconnect(conn), silent = TRUE), add = TRUE)

captured_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

captured <- map(target_fisheries, function(fn) {
  res <- tryCatch(
    {
      d <- resolve_dates(fn, "", "", conn = conn)
      tibble(
        fishery_name   = fn,
        est_date_start = as.character(d$est_date_start),
        est_date_end   = as.character(d$est_date_end),
        resolve_status = "ok",
        resolve_note   = NA_character_
      )
    },
    error = function(e) tibble(
      fishery_name   = fn,
      est_date_start = NA_character_,
      est_date_end   = NA_character_,
      resolve_status = "error",
      # Recorded rather than dropped: a fishery missing from the lookup and a
      # fishery whose lookup FAILED are different problems, and the external
      # run needs to be able to tell them apart.
      resolve_note   = conditionMessage(e)
    )
  )
  if (res$resolve_status == "ok") {
    cli::cli_alert_success("  {.val {fn}}: {res$est_date_start} to {res$est_date_end}")
  } else {
    cli::cli_alert_danger("  {.val {fn}}: {res$resolve_note}")
  }
  res
}) |> bind_rows() |>
  mutate(
    n_days_window = as.integer(as.Date(est_date_end) - as.Date(est_date_start) + 1),
    captured_at   = captured_at
  ) |>
  relocate(n_days_window, .after = est_date_end)

# Preserve any previously-captured fisheries that are not in this run's target
# set, so a scope change tops the file up instead of truncating it.
if (file.exists(LOOKUP_PATH)) {
  prior <- read_csv(LOOKUP_PATH, show_col_types = FALSE) |>
    filter(!fishery_name %in% captured$fishery_name)
  n_kept <- nrow(prior)
  if (n_kept > 0) cli::cli_alert_info("Preserving {n_kept} previously-captured fishery-year(s) not in this run.")
  captured <- bind_rows(prior, captured)
}

captured <- captured |> arrange(fishery_name)
write_csv(captured, LOOKUP_PATH)

cli::cli_h2("Capture summary")
captured |> count(resolve_status) |> print()
if (any(captured$resolve_status == "error")) {
  cli::cli_h3("Failed to resolve")
  captured |> filter(resolve_status == "error") |> select(fishery_name, resolve_note) |> print(n = 50)
}

cli::cli_alert_success("Wrote {nrow(captured)} row(s) to {.file {LOOKUP_PATH}}")
cli::cli_alert_warning(
  "COMMIT THIS FILE. Without it, a DATA_SOURCE = \"external\" run cannot \\
   resolve estimation windows and every fishery-year will skip."
)
