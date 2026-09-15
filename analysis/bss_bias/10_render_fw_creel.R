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
#     RENDER_FISHERY_RE <- "^Stillaguamish salmon and gamefish 2022-23$"  # one fishery-year
#     RENDER_FISHERY_RE <- "Stillaguamish"                                # a whole basin
#     RENDER_GROUPS     <- c("chinook_all", "coho_harvest")
#     RENDER_PROJECT    <- "bss_bias"
#     RENDER_SKIP_DONE  <- TRUE     # skip fishery-years that already have output
#     RENDER_FIT_ONLY   <- TRUE     # TRUE drops the plots/tables; faster, empty HTML
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

# The DEFAULT is the working set, not whatever was last debugged. A one-off
# target belongs in the variable before sourcing, not here -- a narrow default
# means a plain source() silently does less than it appears to.
if (!exists("RENDER_FISHERY_RE", inherits = FALSE)) RENDER_FISHERY_RE <- "Snohomish|Stillaguamish"
if (!exists("RENDER_GROUPS",     inherits = FALSE)) RENDER_GROUPS     <- names(CATCH_GROUPS)
if (!exists("RENDER_PROJECT",    inherits = FALSE)) RENDER_PROJECT    <- "bss_bias"
if (!exists("RENDER_SKIP_DONE",  inherits = FALSE)) RENDER_SKIP_DONE  <- TRUE
# FALSE (the default) renders the full report -- plots and tables -- so the HTML
# per fishery-year is worth opening. TRUE skips all 25 plot/table chunks, which
# is faster but leaves the HTML essentially empty.
if (!exists("RENDER_FIT_ONLY",   inherits = FALSE)) RENDER_FIT_ONLY   <- FALSE
# Render order. A sweep is hours long and the first thing it is asked to do is
# tell you whether it still crashes -- so the cheapest fishery-years go first.
# Basin priority is explicit rather than inferred, with the smallest dataset
# within a basin (fewest interviews, per 00d's inventory) running first.
# Anything not named here sorts after everything that is.
if (!exists("RENDER_BASIN_ORDER", inherits = FALSE)) {
  RENDER_BASIN_ORDER <- c("Stillaguamish", "Snohomish", "Skagit")
}

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
  fit_only                      = RENDER_FIT_ONLY
)

# Per-fishery overrides, keyed on exact fishery_name. Anything here wins over
# BASE_PARAMS for that run only -- e.g. a fishery that needs a different
# study_design or a pinned window.
RUN_PARAM_OVERRIDES <- list(
  # Stillaguamish 2022-23: mainstem sections only, and a window truncated before
  # the late-November tail.
  #
  # SECTIONS -- this year carries sections beyond the mainstem that the rest of
  # its series does not. A b series only means something if every year in it
  # describes the same fishery, and this is also the source of that year's
  # unmatched closure rows.
  #
  # DATES -- the lookup table gives 2022-09-01 to 2022-11-30 (91 days), the
  # longest window in the series. Both ends are pinned here rather than just the
  # end: resolve_dates() queries the database whenever EITHER is blank, so
  # pinning both also makes this run reproducible without a connection.
  #
  # 2022-11-15 IS AN ASSUMPTION, not a value read off the data -- change this one
  # date if the tail starts somewhere else.
  "Stillaguamish salmon and gamefish 2022-23" = list(
    section_filter = c(1, 2, 3),
    est_date_start = "2022-09-01",
    est_date_end   = "2022-11-15"
  )
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
# Ordered further down, once the inventory is available to size each one.

if (length(targets) == 0) cli::cli_abort("No fishery matched {.val {RENDER_FISHERY_RE}}.")

# ------------------------------------------------------------------------------
# Which catch groups are worth fitting, per fishery-year
#
# fw_creel builds inputs_bss from unique(est_cg) on the interview table, and
# prep_dwg_interview_catch() puts est_cg on EVERY interview regardless of
# whether a fish of that group was caught -- it replicates the interviews per
# group and fills fish_count = 0. So the Rmd will happily fit a catch group
# with zero recorded fish, which is a full MCMC run to learn that the posterior
# is the prior truncated by having seen nothing.
#
# 00d's inventory already knows which groups have records. Narrowing
# est_catch_groups per fishery here means those runs never start.
# ------------------------------------------------------------------------------

inv_path <- file.path(OUT_DIR, "bss_catch_inventory.csv")
inventory <- if (file.exists(inv_path)) read_csv(inv_path, show_col_types = FALSE) else NULL
if (is.null(inventory)) {
  cli::cli_alert_warning(
    "No {.file bss_catch_inventory.csv} -- every group will be fitted, including any with \
     zero fish. Run {.file 00d_catch_inventory.R} first to skip those."
  )
}

groups_for <- function(fn) {
  if (is.null(inventory)) return(RENDER_GROUPS)
  rows <- inventory |> filter(fishery_name == fn, catch_group %in% RENDER_GROUPS)
  if (nrow(rows) == 0) {
    cli::cli_alert_warning("{fn}: not in the inventory -- fitting all groups.")
    return(RENDER_GROUPS)
  }
  keep <- rows |> filter(!is.na(n_fish), n_fish > 0) |> pull(catch_group)
  intersect(RENDER_GROUPS, keep)
}

already_done <- function(fn) {
  d <- here::here("fishery_analyses", RENDER_PROJECT, fn)
  if (!dir.exists(d)) return(FALSE)
  length(list.files(d, pattern = "^estimates_bss\\.rds$", recursive = TRUE)) > 0
}

group_plan <- setNames(lapply(targets, groups_for), targets)
empty <- names(group_plan)[lengths(group_plan) == 0]
if (length(empty) > 0) {
  cli::cli_alert_info("Skipping {length(empty)} fishery-year{?s} with no catch group carrying fish:")
  print(empty)
  targets <- setdiff(targets, empty)
  group_plan <- group_plan[targets]
}
narrowed <- names(group_plan)[lengths(group_plan) < length(RENDER_GROUPS)]
if (length(narrowed) > 0) {
  cli::cli_alert_info("Fitting a reduced group set (zero-fish groups dropped) for:")
  for (fn in narrowed) {
    cli::cli_li("{fn}: {paste(group_plan[[fn]], collapse = ', ')}")
  }
}

# ------------------------------------------------------------------------------
# Order: cheapest first, so a crash surfaces in minutes rather than hours
# ------------------------------------------------------------------------------

# Size proxy: the interviews the largest planned group will carry. n_fish is
# the wrong measure -- a group with one Chinook still fits the whole season of
# effort and interview data, and it is the interview and count volume that sets
# sampling time.
size_of <- function(fn) {
  if (is.null(inventory) || !"n_interviews" %in% names(inventory)) return(NA_real_)
  rows <- inventory |> filter(fishery_name == fn, catch_group %in% group_plan[[fn]])
  if (nrow(rows) == 0) return(NA_real_)
  # max() of an all-NA column returns -Inf, which would sort this fishery to
  # the FRONT -- the opposite of what an unknown size should do.
  v <- suppressWarnings(max(rows$n_interviews, na.rm = TRUE))
  if (!is.finite(v)) NA_real_ else as.numeric(v)
}

basin_rank <- function(fn) {
  hit <- which(vapply(RENDER_BASIN_ORDER, function(b) grepl(b, fn, fixed = TRUE), logical(1)))
  if (length(hit) == 0) length(RENDER_BASIN_ORDER) + 1L else min(hit)
}

order_tbl <- tibble(
  fishery_name = targets,
  basin_rank   = map_int(targets, basin_rank),
  n_interviews = map_dbl(targets, size_of)
) |>
  # NA size sorts last within its basin rather than first: an unknown is not a
  # reason to spend the first hour of the sweep on it.
  arrange(basin_rank, is.na(n_interviews), n_interviews, fishery_name)

targets    <- order_tbl$fishery_name
group_plan <- group_plan[targets]

skipped <- character(0)
if (RENDER_SKIP_DONE) {
  done <- targets[map_lgl(targets, already_done)]
  skipped <- done
  targets <- setdiff(targets, done)
}

cli::cli_h1("10 -- render fw_creel.Rmd per fishery-year")
cli::cli_alert_info("Project:  {.val {RENDER_PROJECT}}")
cli::cli_alert_info("Filter:   {.val {RENDER_FISHERY_RE}}")
cli::cli_alert_info("Report:   {if (RENDER_FIT_ONLY) 'fit only -- HTML will be near-empty' else 'full report with plots and tables'}")
cli::cli_alert_info("Groups:   {.val {RENDER_GROUPS}} (fitted in ONE render each; zero-fish groups dropped per fishery)")
cli::cli_alert_info("Order:    cheapest first -- {.val {RENDER_BASIN_ORDER}}, smallest dataset within each")
cli::cli_alert_info("To render: {length(targets)}")
order_tbl |>
  filter(fishery_name %in% targets) |>
  select(fishery_name, n_interviews) |>
  print(n = Inf)
if (length(skipped) > 0) {
  cli::cli_alert_info("Skipping {length(skipped)} already carrying estimates_bss.rds (RENDER_SKIP_DONE):")
  print(skipped)
}
if (length(targets) == 0) {
  cli::cli_alert_success("Nothing to do. Next: 09_read_production_estimates.R")
} else {
  print(targets)
}

if (length(targets) > 0) {

  # --------------------------------------------------------------------------
  # Render. Serial on purpose -- see 01b's history with concurrent runs sharing
  # output files.
  # --------------------------------------------------------------------------

  ledger_file <- function() file.path(OUT_DIR, "bss_render_ledger.csv")

  # Append-or-replace one row, keyed on fishery_name, so a re-run of one
  # fishery updates its row rather than adding a second one.
  write_ledger_row <- function(row) {
    path <- ledger_file()
    all <- if (file.exists(path)) {
      prev <- read_csv(path, show_col_types = FALSE)
      bind_rows(anti_join(prev, row, by = "fishery_name"), row)
    } else {
      row
    }
    write_csv(all, path)
  }

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

    grps <- group_plan[[fn]]
    # NOT folded into the modifyList() below: a data.frame IS a list in R, so
    # modifyList() recurses into it column-by-column instead of replacing it
    # wholesale. Assigning a narrowed (fewer-row) est_catch_groups into a
    # column of BASE_PARAMS's (wider) one via `x[[v]] <- val[[v]]` recycles
    # the shorter column to the original row count -- e.g. one real group
    # narrowed down from the two-group default silently becomes TWO identical
    # rows of that one group. That duplicated est_catch_groups row is what
    # made prep_dwg_interview_catch() replicate every interview twice for a
    # single-group fishery, which is exactly the mismatch ("dims declared=785,
    # dims found=1788") that killed Stillaguamish 2025-26's coho_harvest fit.
    p <- utils::modifyList(BASE_PARAMS, list(fishery_name = fn))
    p$est_catch_groups <- catch_groups_df(grps)
    cli::cli_alert_info("Catch group{?s} for this render: {.val {grps}}")
    if (!is.null(RUN_PARAM_OVERRIDES[[fn]])) {
      # Same hazard applies here if an override ever supplies est_catch_groups
      # with a different row count than `p$est_catch_groups` above.
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

    row <- tibble(
      fishery_name = fn,
      status       = out$status,
      error        = out$error,
      runtime_min  = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
      report       = file.path(REPORT_DIR, paste0(safe_name(fn), ".html")),
      rendered_at  = Sys.time()
    )

    # Write the ledger after EVERY render, not once at the end of map_dfr.
    # A sweep is hours long and gets interrupted; with a single write at the
    # bottom, Ctrl-C loses the error messages for every render that already
    # ran -- which is exactly why the last failure had to be diagnosed from a
    # truncated console paste.
    write_ledger_row(row)
    row
  }

  ledger <- map_dfr(targets, render_one)

  cli::cli_h2("Result")
  ledger |> select(fishery_name, status, runtime_min) |> print(n = Inf)

  n_bad <- sum(ledger$status == "error")
  if (n_bad > 0) {
    cli::cli_alert_danger("{n_bad} render{?s} failed -- see {.file bss_render_ledger.csv} for the messages.")
  }
  cli::cli_alert_info("Next: {.code source(\"analysis/bss_bias/09_read_production_estimates.R\")}")
}
