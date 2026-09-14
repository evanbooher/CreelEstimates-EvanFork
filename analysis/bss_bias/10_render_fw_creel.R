# ==============================================================================
# 10_render_fw_creel.R
#
# Purpose:
#   Render template_scripts/fw_creel.Rmd once per fishery-year, with the params
#   you would otherwise type into the YAML header by hand.
#
#   Each render is a full production run: it writes its own analysis folder
#   under fishery_analyses/<project>/<fishery_name>/, fits every catch group,
#   and saves estimates_bss.rds. 09_read_production_estimates.R then collects
#   them for 07.
#
#   ONE RENDER COVERS EVERY CATCH GROUP. fw_creel builds an inputs_bss entry
#   per unique est_cg and loops the fit across all of them, so passing a
#   two-row est_catch_groups fits Chinook and Coho in the same render rather
#   than requiring two.
#
# Usage (R Console):
#   source("analysis/bss_bias/10_render_fw_creel.R")
#
#   Override before sourcing:
#     RENDER_FISHERY_RE <- "Snohomish|Stillaguamish"   # regex over fishery_name
#     RENDER_GROUPS     <- c("chinook_all", "coho_harvest")
#     RENDER_PROJECT    <- "bss_bias"
#     RENDER_SKIP_DONE  <- TRUE     # skip fishery-years that already have output
#     RENDER_DRY_RUN    <- TRUE     # print the plan, render nothing
#
#   Run it serially and watch it. Each render is a full BSS fit; budget
#   accordingly and check the first one before walking away.
#
# Outputs:
#   fishery_analyses/<project>/<fishery_name>/<run folder>/   per the Rmd
#   analysis/bss_bias/outputs/reports/<fishery>.html          the rendered report
#   analysis/bss_bias/outputs/bss_render_ledger.csv           status per run
# ==============================================================================

library(tidyverse)
library(cli)
library(here)

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  cli::cli_abort("Package {.pkg rmarkdown} is required.")
}

OUT_DIR    <- here::here("analysis", "bss_bias", "outputs")
REPORT_DIR <- file.path(OUT_DIR, "reports")
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)

source(here::here("analysis", "bss_bias", "catch_groups.R"))

if (!exists("RENDER_FISHERY_RE", inherits = FALSE)) RENDER_FISHERY_RE <- "Snohomish|Stillaguamish"
if (!exists("RENDER_GROUPS",     inherits = FALSE)) RENDER_GROUPS     <- names(CATCH_GROUPS)
if (!exists("RENDER_PROJECT",    inherits = FALSE)) RENDER_PROJECT    <- "bss_bias"
if (!exists("RENDER_SKIP_DONE",  inherits = FALSE)) RENDER_SKIP_DONE  <- TRUE
if (!exists("RENDER_DRY_RUN",    inherits = FALSE)) RENDER_DRY_RUN    <- FALSE

RMD <- here::here("template_scripts", "fw_creel.Rmd")
if (!file.exists(RMD)) cli::cli_abort("{.file {RMD}} not found.")

safe_name <- function(x) stringr::str_replace_all(x, "[^[:alnum:]]", "_")

# ------------------------------------------------------------------------------
# Params. Everything not named here keeps the Rmd's own YAML default, so this
# stays a thin override rather than a second copy of the header that can drift.
#
# est_date_start / est_date_end are left EMPTY on purpose: fw_creel resolves the
# window itself via resolve_dates(), which is what you get typing the header by
# hand. Set them per fishery in RUN_PARAM_OVERRIDES to pin a window.
# ------------------------------------------------------------------------------

BASE_PARAMS <- list(
  project_name                  = RENDER_PROJECT,
  est_date_start                = "",
  est_date_end                  = "",
  est_catch_groups              = catch_groups_df(RENDER_GROUPS),
  study_design                  = "Standard",
  boat_type_collapse            = "Yes",
  fish_location_determines_type = "No",
  angler_type_kayak_pontoon     = "bank",
  person_count_type             = "group",
  period_pe                     = "week",
  period_bss                    = "day",
  day_length_expansion          = "night closure",
  min_fishing_time              = 0.5,
  bss_model_file_name           = "BSS_creel_model_02_2021-01-22_ppc.stan",
  model_used                    = "Both models",
  data_grade                    = "provisional",
  export                        = "local",
  export_tables                 = "both",
  enable_cache                  = FALSE,
  # REQUIRED. With save_draws FALSE, fw_creel deletes
  # estimates_bss[[ecg]]$draws and $season_results before writing the file, so
  # 09_read_production_estimates.R would find no C_sum, E_sum or b to read.
  save_draws                    = TRUE,
  # Skips every plot and table chunk: fetch, prep, sample, save, nothing else.
  # Faster, and it removes the rendering steps that can fail AFTER a fit has
  # completed -- the trace plot in particular aborts on the NaN generated
  # quantities a low-catch group produces.
  fit_only                      = TRUE
)

# Per-fishery overrides, keyed on exact fishery_name. Anything here wins over
# BASE_PARAMS for that run only -- e.g. a fishery that needs a different
# study_design or a pinned window.
RUN_PARAM_OVERRIDES <- list(
  # "Stillaguamish salmon and gamefish 2025-26" = list(est_date_start = "2025-09-01")
)

# ------------------------------------------------------------------------------
# Which fishery-years
# ------------------------------------------------------------------------------

disc_path <- file.path(OUT_DIR, "fishery_discovery_target.csv")
if (!file.exists(disc_path)) {
  cli::cli_abort("{.file {disc_path}} not found -- run 00_discover_fisheries.R first.")
}

targets <- read_csv(disc_path, show_col_types = FALSE) |>
  filter(basin_match == "target", include_in_run,
         str_detect(fishery_name, RENDER_FISHERY_RE)) |>
  pull(fishery_name) |> unique() |> sort()

if (length(targets) == 0) cli::cli_abort("No fishery matched {.val {RENDER_FISHERY_RE}}.")

already_done <- function(fn) {
  d <- here::here("fishery_analyses", RENDER_PROJECT, fn)
  if (!dir.exists(d)) return(FALSE)
  length(list.files(d, pattern = "^estimates_bss\\.rds$", recursive = TRUE)) > 0
}

skipped <- character(0)
if (RENDER_SKIP_DONE) {
  done <- targets[map_lgl(targets, already_done)]
  skipped <- done
  targets <- setdiff(targets, done)
}

cli::cli_h1("10 -- render fw_creel.Rmd per fishery-year")
cli::cli_alert_info("Project:  {.val {RENDER_PROJECT}}")
cli::cli_alert_info("Filter:   {.val {RENDER_FISHERY_RE}}")
cli::cli_alert_info("Groups:   {.val {RENDER_GROUPS}} (all fitted in ONE render each)")
cli::cli_alert_info("To render: {length(targets)}")
if (length(skipped) > 0) {
  cli::cli_alert_info("Skipping {length(skipped)} already carrying estimates_bss.rds (RENDER_SKIP_DONE):")
  print(skipped)
}
if (length(targets) == 0) {
  cli::cli_alert_success("Nothing to do. Next: 09_read_production_estimates.R")
} else {
  print(targets)
}

if (RENDER_DRY_RUN) {
  cli::cli_alert_info("RENDER_DRY_RUN is TRUE -- stopping before any render.")
} else if (length(targets) > 0) {

  # --------------------------------------------------------------------------
  # Render. Serial on purpose -- see 01b's history with concurrent runs sharing
  # output files.
  # --------------------------------------------------------------------------

  render_one <- function(fn) {
    cli::cli_h2("{fn}")

    # generate_analysis_lut() does assign("analysis_lut", ..., envir = .GlobalEnv)
    # and REUSES an existing one. Left in place, the second render of a session
    # inherits the first fishery's analysis_id and writes into its folder --
    # silently, with both fisheries' outputs landing in one place. Clearing it
    # is what makes one render per fishery actually mean one folder per fishery.
    if (exists("analysis_lut", envir = .GlobalEnv)) {
      rm("analysis_lut", envir = .GlobalEnv)
    }

    p <- utils::modifyList(BASE_PARAMS, list(fishery_name = fn))
    if (!is.null(RUN_PARAM_OVERRIDES[[fn]])) {
      p <- utils::modifyList(p, RUN_PARAM_OVERRIDES[[fn]])
    }

    t0 <- Sys.time()
    out <- tryCatch({
      rmarkdown::render(
        input        = RMD,
        params       = p,
        output_file  = paste0(safe_name(fn), ".html"),
        output_dir   = REPORT_DIR,
        knit_root_dir = here::here(),   # the Rmd uses here() throughout
        envir        = new.env(parent = globalenv()),
        quiet        = FALSE
      )
      list(status = "ok", error = NA_character_)
    }, error = function(e) {
      cli::cli_alert_danger("FAILED: {conditionMessage(e)}")
      list(status = "error", error = conditionMessage(e))
    })

    tibble(
      fishery_name = fn,
      status       = out$status,
      error        = out$error,
      runtime_min  = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
      report       = file.path(REPORT_DIR, paste0(safe_name(fn), ".html")),
      rendered_at  = Sys.time()
    )
  }

  ledger <- map_dfr(targets, render_one)

  ledger_path <- file.path(OUT_DIR, "bss_render_ledger.csv")
  if (file.exists(ledger_path)) {
    prev <- read_csv(ledger_path, show_col_types = FALSE)
    ledger <- bind_rows(anti_join(prev, ledger, by = "fishery_name"), ledger)
  }
  write_csv(ledger, ledger_path)

  cli::cli_h2("Result")
  ledger |> select(fishery_name, status, runtime_min) |> print(n = Inf)

  n_bad <- sum(ledger$status == "error")
  if (n_bad > 0) {
    cli::cli_alert_danger("{n_bad} render{?s} failed -- see {.file bss_render_ledger.csv} for the messages.")
  }
  cli::cli_alert_info("Next: {.code source(\"analysis/bss_bias/09_read_production_estimates.R\")}")
}
