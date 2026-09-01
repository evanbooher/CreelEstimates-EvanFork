# ==============================================================================
# 01_fit_bss_bias.R
#
# Purpose:
#   For every fishery-year in analysis/bss_bias/outputs/fishery_discovery_target.csv
#   (produced by 00_discover_fisheries.R and hand-reviewed -- only rows with
#   include_in_run == TRUE are attempted), prep BSS inputs, fit the BSS Stan
#   model, and extract the effort-index bias term `b` (vehicle-count bias
#   b[1], trailer-count bias b[2]) into a tidy summary table.
#
#   Also assembles, PER FISHERY-YEAR, a "comparability row" -- date window,
#   CRC areas/sections, p_TI, count-type availability, sample sizes -- and
#   writes it to disk INCREMENTALLY, before any Stan fitting is attempted.
#   This means 02_build_comparability_table.R's deliverable (the standalone
#   comparability table) is available even if MCMC fitting is still running,
#   partially failed, or never gets attempted at all -- see README.md's
#   "go/no-go" section. That table alone is a legitimate meeting deliverable.
#
# Requires (NOT available in a VPN-less / no-internal-DB environment):
#   creelutils::connect_creel_db(), creelutils::fishery_lut(),
#   creelutils::fetch_data(..., data_source = "internal")
#   -- UNLESS the public-data path below (0.A) turns out to carry BSS-grade
#   data, in which case data_source = "internal" can be swapped for whatever
#   value creelutils uses for its data.wa.gov path and this becomes VPN-free.
#
# ------------------------------------------------------------------------------
# BEFORE RUNNING -- two things to check first, in order:
#
# [0.A] CONFIRMED (checked locally, VPN off): creelutils::fetch_data(...,
#       data_source = "external") is the public/data.wa.gov path, and it DOES
#       carry BSS-grade data -- dwg_test$effort (5247x25) and dwg_test$catch
#       (964x11) passed the full column checklist below; dwg_test$interview
#       (1811x39) has every required raw column except fishing_time_total /
#       person_count_final, which are computed downstream from columns that
#       ARE present (fishing_start_time/fishing_end_time, angler_count/
#       total_group_count) -- not a gap. dwg$ll's centroid_lat/centroid_lon
#       and dwg$fishery_manager's p_census_bank/p_census_boat were not
#       independently re-verified in this check, but the latter is moot for
#       `b`: p_census_bank/boat feed the p_TI/census-tie-in spatial-coverage
#       correction, a defunct-in-practice DIFFERENT parameter from `b` (see
#       README.md's "What 'the bias term' is") -- irrelevant to this
#       analysis's identification of `b` either way.
#
#       So DATA_SOURCE <- "external" below is a viable VPN-free option. This
#       run used DATA_SOURCE <- "internal" anyway (VPN/DB already available
#       locally) -- purely a convenience choice, not a fallback from a failed
#       "external" check.
#
# [0.B] Time ONE fishery-year end to end at the "smoke" config (see
#       FIT_CONFIGS below) before committing to a full run. If it takes materially
#       longer than a few minutes, multiply by the number of target fishery-years
#       and compare to the time actually available before the meeting -- see
#       README.md's triage sequence (T0-T4) for what to do if the answer is
#       "not enough time."
# ------------------------------------------------------------------------------
#
# Outputs (analysis/bss_bias/outputs/), all written incrementally (as-you-go,
# not once at the end -- a crash partway through must not lose prior
# results) and UPSERTED by fishery_name (append_csv_row() below) -- re-
# running a fishery-year (e.g. re-attempting at "quick" after "smoke", or
# just retrying after a fix) replaces its prior row(s) rather than
# accumulating duplicates alongside them:
#   bss_b_comparability_raw.csv  -- one row per attempted fishery-year, written
#                                    BEFORE fitting; feeds 02_build_comparability_table.R
#   bss_b_summary.csv            -- one/two rows (b[1], b[2]) per fishery-year that fit
#   bss_b_run_ledger.csv         -- one row per attempted fishery-year: status/stage/reason
#   bss_b_stan_dims.csv          -- D,G,S,H,V_n,T_n,A_n,E_n,IntC,IntA per fishery-year
#   bss_b_na_drops.csv           -- per fishery-year x group, NA rows dropped before Stan
#                                    (n_before/n_dropped/n_after; see R_functions/drop_na_bss_inputs.R)
#   b_draws/<safe_name>.rds      -- raw b[1]/b[2] posterior draws (small; kept for Phase 5)
#   fits/<safe_name>.rds         -- full stanfit object, ONLY if SAVE_FITS == TRUE (large)
# ==============================================================================

library(tidyverse)
library(cli)
library(here)
library(rstan)
# Without this, stan() recompiles the model from C++ source on EVERY call --
# not just the first -- since there is no compiled-model cache to read from
# or write to. With it, stan_model() writes a hashed .rds cache next to the
# .stan file after the first compile and reuses it on every subsequent call
# (same model file = same hash), including across separate Rscript runs.
# Matches template_scripts/fw_creel.Rmd's setup; missing here was turning
# the README's expected "one-time ~1-3 min recompile" into a recompile on
# every one of the ~29 queued fishery-years.
rstan_options(auto_write = TRUE)
library(posterior)
library(creelutils)
library(timeDate)
library(suncalc)
library(lubridate)
library(rlang)

walk(list.files(here("R_functions"), full.names = TRUE), source)

OUT_DIR   <- here::here("analysis", "bss_bias", "outputs")
FITS_DIR  <- file.path(OUT_DIR, "fits")
DRAWS_DIR <- file.path(OUT_DIR, "b_draws")
dir.create(OUT_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(FITS_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(DRAWS_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Configuration -- read this before running
# ------------------------------------------------------------------------------

DATA_SOURCE <- "internal"  # <-- change per the 0.A check above if a public path works

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

# Standardize on ONE Stan model version for every fishery-year. `b`'s
# declaration/prior/likelihood placement is identical across all four model
# files in stan_models/ (verified by inspection), so this choice does not
# change what `b` means -- but see the dimension-mismatch note below.
#
# BSS_creel_model_02_2021-01-22_ppc.stan is used here because it is the file
# actually in current production use (per direct confirmation, overriding an
# earlier draft of this script that defaulted to 2024-04-03 on speed grounds
# alone -- don't re-litigate this without checking with the user first). It
# declares `O` as `matrix[D,S]` (2-D), matching what prep_inputs_bss()
# (unmodified, per the reuse decision below) actually builds -- compatible.
#
# Two things this choice costs, both handled: (1) it's the *_ppc variant, so
# it carries a `*_rep` generated-quantities block for posterior predictive
# checks -- those quantities still get COMPUTED during sampling regardless of
# `pars=`, so per-iteration cost is somewhat higher than the leaner
# non-_ppc/2024-04-03 files, even though MONITOR_PARS below keeps them out of
# what's returned/stored. (2) no compiled .rds cache exists for this exact
# file in stan_models/ (only BSS_creel_model_02_2021-01-22.rds, the non-_ppc
# base, and BSS_creel_model_02_2024-04-03.rds) -- expect a ~1-3 min recompile
# on the FIRST fit only, same as any other first-time rstan compile.
#
# 2024-07-24 remains excluded regardless of which of the other three is
# chosen: it declares `O` as `real O[D,S,G]` (3-D), a hard mismatch against
# prep_inputs_bss()'s 2-D output, not a subtle bug.
BSS_MODEL_FILE <- "BSS_creel_model_02_2021-01-22_ppc.stan"

# Copied verbatim from template_scripts/fw_creel.Rmd's prep_inputs_bss() call.
# Held CONSTANT across every fishery-year -- changing the b prior between
# years would make the whole temporal-stability comparison meaningless.
BSS_PRIORS <- c(
  value_cauchyDF_sigma_eps_C   = 0.5,
  value_cauchyDF_sigma_eps_E   = 0.5,
  value_cauchyDF_sigma_r_E     = 0.5,
  value_cauchyDF_sigma_r_C     = 0.5,
  value_cauchyDF_sigma_mu_C    = 0.5,
  value_cauchyDF_sigma_mu_E    = 0.5,
  value_normal_sigma_omega_C_0 = 1,
  value_normal_sigma_omega_E_0 = 3,
  value_lognormal_sigma_b      = 1,     # <-- this is the number get_bss_bias() needs as prior_sigma_b
  value_normal_sigma_B1        = 5,
  value_normal_mu_mu_C         = log(0.02),
  value_normal_sigma_mu_C      = 1.5,
  value_normal_mu_mu_E         = log(5),
  value_normal_sigma_mu_E      = 2,
  value_betashape_phi_E_scaled = 1,
  value_betashape_phi_C_scaled = 1
)

# Only parameters actually needed downstream are monitored -- excludes
# eps_E_H / lambda_E_S_I / lambda_E_S / lambda_C_S / omega_*, which dominate
# fit size/runtime for a model this shape (see README.md). This is the
# single biggest speed/memory lever available without changing chains/iter.
MONITOR_PARS <- c(
  "b", "R_V", "R_T", "p_I", "mu_E", "mu_C", "B1",
  "phi_E_scaled", "phi_C_scaled", "sigma_eps_E", "sigma_mu_E",
  "E_sum", "C_sum"
)

FIT_CONFIGS <- list(
  smoke = list(n_chain = 1, n_cores = 1, n_iter = 60,   n_warmup = 30,   n_thin = 1, adapt_delta = 0.70, max_treedepth = 10),
  quick = list(n_chain = 2, n_cores = 2, n_iter = 600,  n_warmup = 300,  n_thin = 1, adapt_delta = 0.80, max_treedepth = 11),
  prod  = list(n_chain = 4, n_cores = 4, n_iter = 2000, n_warmup = 1000, n_thin = 1, adapt_delta = 0.95, max_treedepth = 13)
)
FIT_CONFIG_NAME <- "smoke"   # <-- default is the fast/throwaway config; change to "quick" once smoke passes, "prod" for backfill later

SAVE_FITS <- FALSE   # TRUE keeps the full stanfit per fishery-year (large!); the small
                      # b-summary + draws are the actual deliverable and are always saved.

MIN_INTA_INFORMATIVE <- 30   # below this, flag b as weakly informed (see get_bss_bias() / prior_contraction)

BOAT_TYPE_COLLAPSE        <- "Yes"
FISH_LOC_DETERMINES_TYPE  <- "No"
ANGLER_TYPE_KAYAK_PONTOON <- "bank"
MIN_FISHING_TIME          <- 0.5
DAY_LENGTH                <- "night closure"
PERIOD_PE                 <- "week"   # for prep_days()'s PE-style period column (unused by BSS itself)
PERIOD_BSS                <- "day"    # what actually indexes the BSS model; see prep_inputs_bss(period=)

SALMON_SPECIES     <- c("Chinook", "Coho", "Chum", "Pink", "Sockeye")
HARVEST_FATE        <- "Kept|Released"  # BSS catch (`c`) is total encounters, not harvest-only -- broader than the PST harvest scope on purpose
EXCLUDE_LIFE_STAGES <- c("Smolt")

# ------------------------------------------------------------------------------
# Skip/error machinery -- transcribed from analysis/pst/02_ingest/
# multi_fishery_creel_summary.R on chore/multi-fishery-trip-summary (read-only
# reference; NOT imported/merged -- this branch does not depend on that one).
# ------------------------------------------------------------------------------

skip_fishery <- function(reason, stage = "preflight") {
  rlang::abort(message = reason, class = "fishery_skip", stage = stage, reason = reason)
}

run_stage <- function(stage, code) {
  tryCatch(
    code,
    fishery_skip = function(cnd) stop(cnd),
    error = function(e) {
      rlang::abort(
        message = paste0("[", stage, "] ", conditionMessage(e)),
        class = "fishery_error", stage = stage, reason = conditionMessage(e)
      )
    }
  )
}

DESIGN_RULES <- tibble::tribble(
  ~pattern,      ~study_design,
  "Drano Lake",  "Drano"
)
DEFAULT_STUDY_DESIGN <- "Standard"
SUPPORTED_DESIGNS    <- c("Standard", "Drano")

resolve_study_design <- function(fishery_name) {
  hits <- DESIGN_RULES[
    stringr::str_detect(fishery_name, stringr::regex(DESIGN_RULES$pattern, ignore_case = TRUE)), ,
    drop = FALSE
  ]
  design <- if (nrow(hits) == 0) DEFAULT_STUDY_DESIGN else hits$study_design[1]
  if (!design %in% SUPPORTED_DESIGNS) {
    cli::cli_abort("Unsupported study design {.val {design}} for {.val {fishery_name}}.")
  }
  # All three target basins are expected to resolve to "Standard". A hit here
  # (e.g. against a Drano-style name) is a loud signal to re-check scope, not
  # something this script is set up to model correctly for BSS.
  design
}

REQUIRED_INTERVIEW_COLS <- list(
  common   = c("interview_id", "event_date", "section_num", "crc_area",
               "trip_status", "previously_interviewed", "fishing_start_time",
               "interview_time", "vehicle_count", "trailer_count", "boat_used", "boat_type"),
  Standard = c("total_group_count", "fish_from_boat"),
  Drano    = c("angler_count")
)

preflight_fishery <- function(dwg, fishery_name, date_start, date_end, study_design) {
  need <- c(REQUIRED_INTERVIEW_COLS$common, REQUIRED_INTERVIEW_COLS[[study_design]])
  missing_cols <- setdiff(need, names(dwg$interview))
  if (length(missing_cols) > 0) {
    skip_fishery(paste0("Interview table missing column(s) required by '", study_design,
                         "' design: ", paste(missing_cols, collapse = ", "), "."))
  }
  if (is.na(date_start) || is.na(date_end)) skip_fishery("Unresolvable estimation dates.")
  if (date_end < date_start) skip_fishery("Estimation end date precedes start date.")

  int_n <- dwg$interview |> filter(between(event_date, date_start, date_end)) |> nrow()
  if (int_n == 0) skip_fishery("No interviews within the estimation window.")
  eff_n <- dwg$effort |> filter(between(event_date, date_start, date_end)) |> nrow()
  if (eff_n == 0) skip_fishery("No effort count records within the estimation window.")

  sections <- sort(unique(na.omit(
    dwg$interview |> filter(between(event_date, date_start, date_end)) |> pull(section_num)
  )))
  if (length(sections) == 0) skip_fishery("No non-missing section_num values on interviews.")

  lat  <- suppressWarnings(mean(dwg$ll$centroid_lat, na.rm = TRUE))
  long <- suppressWarnings(mean(dwg$ll$centroid_lon, na.rm = TRUE))
  if (!is.finite(lat) || !is.finite(long)) skip_fishery("No usable centroid coordinates in dwg$ll.")

  if (is.null(dwg$catch) || nrow(dwg$catch) == 0) skip_fishery("No catch records returned for this fishery.")

  closures <- dwg$closures
  if (is.null(closures)) closures <- tibble(section_num = double(), event_date = as.Date(character()))
  if (nrow(closures) > 0) {
    closures <- closures |>
      mutate(event_date = as.Date(event_date, format = "%Y-%m-%d"), section_num = as.double(section_num)) |>
      # drop orphan closures (in-window, section not among those being estimated);
      # leave anything outside the window untouched
      filter(!between(event_date, date_start, date_end) | section_num %in% sections) |>
      distinct(section_num, event_date, .keep_all = TRUE)
  }

  list(sections = sections, lat = lat, long = long, closures = closures)
}

# MODIFIED from the PE version: BSS's period is params$period_bss = "day"
# (days$day_index, a running 1..D index), which cannot collide across
# calendar years the way a bare %W week number can. The PE-only
# period-x-year collision check is DOWNGRADED to a warning here rather than
# a skip, because it does not apply to what BSS actually indexes on -- and
# would otherwise wrongly skip every Stillaguamish multi-season fishery
# (e.g. "2022-23"), three of which span calendar years.
validate_days <- function(days, fishery_name) {
  if (!any(grepl("^open_section_", names(days)))) {
    skip_fishery("prep_days() produced no open_section_* columns.", stage = "prep_days")
  }
  if (all(is.na(days$day_length))) skip_fishery("All day_length values are NA.", stage = "prep_days")

  collisions <- days |> distinct(period, year) |> count(period) |> filter(n > 1)
  if (nrow(collisions) > 0) {
    cli::cli_alert_warning(
      "  {.val {fishery_name}}: PE-style week/month `period` values collide \\
       across years ({paste(collisions$period, collapse=', ')}). Not a skip here \\
       -- BSS indexes on day_index (period_bss='day'), not on this column."
    )
  }
  invisible(TRUE)
}

safe_name <- function(x) stringr::str_replace_all(x, "[^[:alnum:]]", "_")

# Upserts by fishery_name rather than blindly appending: a re-run of the same
# fishery-year (routine during troubleshooting, and whenever a fishery gets
# re-attempted at a higher fit_config after "smoke") replaces its prior
# row(s) instead of accumulating duplicates alongside them. Still safe for
# the incremental/crash-resilient design goal -- an interrupted run leaves
# every OTHER fishery's rows untouched; only the fishery-year being written
# right now is ever removed-then-re-added. Every row_df passed in here
# carries a fishery_name column (comparability_row, stan_dims_row,
# bias_summary, na_drop_log all do).
append_csv_row <- function(row_df, path) {
  if (file.exists(path)) {
    existing <- readr::read_csv(path, show_col_types = FALSE) |>
      dplyr::filter(!(fishery_name %in% row_df$fishery_name))
    readr::write_csv(dplyr::bind_rows(existing, row_df), path)
  } else {
    readr::write_csv(row_df, path)
  }
}

# ------------------------------------------------------------------------------
# Catch-group selection rule (see README.md): `b` is fit ONCE per fishery-year
# against the pooled TOTAL SALMON catch group -- all SALMON_SPECIES combined,
# not any single species and not a "most distinct interviewed groups"
# heuristic. This mirrors build_est_catch_groups()'s total_group /
# TOTAL_LABEL convention on chore/multi-fishery-trip-summary's
# multi_fishery_creel_summary.R, and for the same reason given there ([N1]):
# species co-occur within interviews, so pooling must happen at the
# per-interview level BEFORE the CPUE/likelihood calculation, not by summing
# species-level results afterward -- summing would ignore covariance. Species/
# life-stage/fin-mark alternations are built from OBSERVED values per fishery.
# This need not match PST harvest scope exactly since `c`/`h` here just need
# to be a reasonable, consistently-applied catch definition, not a harvest
# total. The chosen est_cg string is recorded in every output row.
# ------------------------------------------------------------------------------

build_catch_groups <- function(dwg_catch, fishery_name) {
  cat_std <- dwg_catch |>
    mutate(across(c(species, life_stage, fin_mark, fate), ~ replace_na(as.character(.), "NA")))

  spp_present <- sort(intersect(unique(cat_std$species), SALMON_SPECIES))
  if (length(spp_present) == 0) {
    skip_fishery(paste0("No PST salmon species in catch records. Species present: ",
                         paste(sort(unique(cat_std$species)), collapse = ", ")), stage = "catch_groups")
  }

  harvest_rows <- cat_std |> filter(species %in% spp_present, str_detect(fate, HARVEST_FATE))
  if (nrow(harvest_rows) == 0) {
    skip_fishery("Salmon present but no records matching the fate scope.", stage = "catch_groups")
  }

  ls_vals <- setdiff(sort(unique(harvest_rows$life_stage)), EXCLUDE_LIFE_STAGES)
  if (length(ls_vals) == 0) skip_fishery("Every record is an excluded life stage.", stage = "catch_groups")
  fm_vals <- sort(unique(harvest_rows$fin_mark[harvest_rows$life_stage %in% ls_vals]))

  bad <- c(spp_present, ls_vals, fm_vals) |> keep(~ str_detect(.x, "[.\\\\+*?\\[\\]^$(){}=!<>|:-]"))
  if (length(bad) > 0) {
    cli::cli_abort("Category value(s) contain regex metacharacters, cannot safely build a catch-group pattern: {.val {bad}}")
  }

  species_groups <- tibble(species = spp_present, life_stage = paste(ls_vals, collapse = "|"),
                            fin_mark = paste(fm_vals, collapse = "|"), fate = HARVEST_FATE)
  total_group <- tibble(species = paste(spp_present, collapse = "|"), life_stage = paste(ls_vals, collapse = "|"),
                         fin_mark = paste(fm_vals, collapse = "|"), fate = HARVEST_FATE)

  # When only one species is present (e.g. "Skagit spring Chinook", "Skagit
  # summer sockeye"), paste(spp_present, collapse="|") on a single element
  # just returns that element -- total_group becomes BYTE-IDENTICAL to
  # species_groups' sole row, same reconstructed est_cg string. Binding both
  # would give prep_dwg_interview_catch() two replicate sets tagged with the
  # same est_cg label; filtering to it then returns both combined, doubling
  # the interview count (the "dims declared vs found" Stan crash). For a
  # single-species fishery, that lone species group already IS the total.
  if (length(spp_present) > 1) {
    bind_rows(species_groups, total_group) |> as.data.frame(stringsAsFactors = FALSE)
  } else {
    as.data.frame(species_groups, stringsAsFactors = FALSE)
  }
}

# Reconstructs the est_cg string EXACTLY as prep_dwg_interview_catch() builds
# it (paste0(unlist(one row of species/life_stage/fin_mark/fate), collapse =
# "_")) for the pooled total-salmon row specifically -- always the LAST row
# of build_catch_groups()'s output, since total_group is bound after
# species_groups. Deterministic: no dependence on which group happens to have
# the most interviews.
total_salmon_est_cg <- function(catch_groups) {
  total_row <- catch_groups[nrow(catch_groups), , drop = FALSE]
  paste0(c(total_row$species, total_row$life_stage, total_row$fin_mark, total_row$fate), collapse = "_")
}

# ------------------------------------------------------------------------------
# Per-fishery-year pipeline
# ------------------------------------------------------------------------------

fit_one_fishery <- function(fishery_name, fit_config_name = FIT_CONFIG_NAME) {

  study_design <- resolve_study_design(fishery_name)
  cli::cli_h2("Processing: {.val {fishery_name}} [design: {.val {study_design}}]")

  est_dates  <- run_stage("resolve_dates", resolve_dates(fishery_name, "", ""))
  date_start <- suppressWarnings(as.Date(est_dates$est_date_start))
  date_end   <- suppressWarnings(as.Date(est_dates$est_date_end))

  dwg <- run_stage("fetch_data", {
    creelutils::fetch_data(conn = DB_CONN, fishery_name = fishery_name, data_source = DATA_SOURCE)
  })

  pf <- run_stage("preflight", preflight_fishery(dwg, fishery_name, date_start, date_end, study_design))

  dwg$effort <- run_stage("patch_p_census", {
    dwg$effort |>
      select(-any_of(c("p_census_bank", "p_census_boat"))) |>
      left_join(
        dwg$fishery_manager |>
          filter(!is.na(p_census_bank) | !is.na(p_census_boat)) |>
          distinct(section_num, p_census_bank, p_census_boat),
        by = "section_num"
      )
  })

  params <- list(fishery_name = fishery_name, project_name = "bss_bias", study_design = study_design)

  dwg$days <- run_stage("prep_days", {
    prep_days(params = params, date_begin = est_dates$est_date_start, date_end = est_dates$est_date_end,
              weekends = c("Saturday", "Sunday"), lat = pf$lat, long = pf$long, period_pe = PERIOD_PE,
              sections = pf$sections, closures = pf$closures, day_length = DAY_LENGTH, day_length_inputs = list())
  })
  run_stage("prep_days", validate_days(dwg$days, fishery_name))

  eff_filt <- dwg$effort |> filter(between(event_date, date_start, date_end))
  int_filt <- dwg$interview |> filter(between(event_date, date_start, date_end))

  interview_fishing_time <- run_stage("interview_fishing_time", {
    prep_dwg_interview_fishing_time(params = params, dwg_interview = int_filt,
                                     min_fishing_time = MIN_FISHING_TIME, study_design = study_design)
  })
  interview_angler_types <- run_stage("interview_angler_types", {
    prep_dwg_interview_angler_types(params = params, interview_fishing_time = interview_fishing_time,
                                     study_design = study_design, boat_type_collapse = BOAT_TYPE_COLLAPSE,
                                     fish_location_determines_type = FISH_LOC_DETERMINES_TYPE,
                                     angler_type_kayak_pontoon = ANGLER_TYPE_KAYAK_PONTOON)
  })

  catch_groups <- run_stage("catch_groups", build_catch_groups(dwg$catch, fishery_name))

  interview_plus_catch <- run_stage("interview_catch", {
    prep_dwg_interview_catch(params = params, interview_plus_angler_types = interview_angler_types,
                              dwg_catch = dwg$catch, study_design = study_design, est_catch_groups = catch_groups)
  })

  chosen_ecg <- run_stage("choose_ecg", total_salmon_est_cg(catch_groups))
  if (!chosen_ecg %in% interview_plus_catch$est_cg) {
    skip_fishery(paste0("Reconstructed total-salmon est_cg ('", chosen_ecg, "') does not match ",
                         "any est_cg actually produced by prep_dwg_interview_catch() -- catch-group ",
                         "string construction mismatch, not a data problem."), stage = "choose_ecg")
  }

  effort_index_summ <- run_stage("effort_index", {
    prep_dwg_effort_index(params = params, eff = eff_filt, study_design = study_design,
                           boat_type_collapse = BOAT_TYPE_COLLAPSE,
                           fish_location_determines_type = FISH_LOC_DETERMINES_TYPE,
                           angler_type_kayak_pontoon = ANGLER_TYPE_KAYAK_PONTOON)
  })
  effort_census_summ <- run_stage("effort_census", {
    prep_dwg_effort_census(params = params, eff = eff_filt, study_design = study_design,
                            boat_type_collapse = BOAT_TYPE_COLLAPSE,
                            fish_location_determines_type = FISH_LOC_DETERMINES_TYPE,
                            angler_type_kayak_pontoon = ANGLER_TYPE_KAYAK_PONTOON)
  })
  if (nrow(effort_index_summ$index_angler_final) == 0) {
    skip_fishery("No index effort counts after angler-type assignment.", stage = "effort_index")
  }

  dwg_summ <- list(
    interview     = interview_plus_catch,
    effort_index  = effort_index_summ$index_angler_final,
    effort_census = effort_census_summ$census_angler_final,
    census_expan  = prep_dwg_census_expan(eff = dwg$effort, days = dwg$days)
  )

  # --- Comparability row: assembled and WRITTEN NOW, before fitting, so it
  #     survives even if the Stan fit fails, is skipped, or never gets run. ---
  section_crc_area <- interview_angler_types |>
    distinct(section_num, crc_area) |>
    group_by(section_num) |>
    filter(!(is.na(crc_area) & any(!is.na(crc_area)))) |>
    ungroup()

  comparability_row <- run_stage("comparability_row", {
    tibble(
      fishery_name           = fishery_name,
      study_design            = study_design,
      date_start               = date_start,
      date_end                 = date_end,
      season_month_span        = paste0(format(date_start, "%b"), "–", format(date_end, "%b")),
      n_days_in_window          = as.integer(date_end - date_start + 1),
      n_days_open                = nrow(dwg$days |> filter(if_any(starts_with("open_section"), ~.))),
      section_nums                = paste(sort(pf$sections), collapse = "|"),
      n_sections                   = length(pf$sections),
      crc_areas                     = paste(sort(unique(na.omit(section_crc_area$crc_area))), collapse = "|"),
      count_types_present            = paste(sort(unique(eff_filt$count_type)), collapse = "|"),
      has_vehicle_counts               = "vehicle" %in% tolower(eff_filt$count_type) || any(str_detect(tolower(eff_filt$count_type), "vehicle")),
      has_trailer_counts                = any(str_detect(tolower(eff_filt$count_type), "trailer")),
      p_TI_bank                          = dwg_summ$census_expan |> filter(angler_final == "bank") |> summarise(m = mean(p_census, na.rm = TRUE)) |> pull(m),
      p_TI_boat                           = dwg_summ$census_expan |> filter(angler_final == "boat") |> summarise(m = mean(p_census, na.rm = TRUE)) |> pull(m),
      p_TI_varies_by_section                = dwg_summ$census_expan |> group_by(angler_final) |> summarise(v = n_distinct(p_census) > 1, .groups="drop") |> pull(v) |> any(),
      chosen_est_cg                          = chosen_ecg,
      angler_types                            = paste(sort(unique(interview_angler_types$angler_final)), collapse = "|")
    )
  })
  append_csv_row(comparability_row, file.path(OUT_DIR, "bss_b_comparability_raw.csv"))
  cli::cli_alert_success("  Comparability row written.")

  # --- BSS input prep + BSS-specific preflight -----------------------------

  fail_counts <- effort_index_summ$index_angler_groups |> filter(angler_final == "fail")
  if (nrow(fail_counts) > 0) {
    frac_fail <- nrow(fail_counts) / nrow(effort_index_summ$index_angler_groups)
    if (frac_fail > 0.5) {
      skip_fishery(paste0(round(100*frac_fail), "% of index effort counts map to 'fail' under design '",
                           study_design, "' -- likely wrong study design."), stage = "effort_index_design")
    }
  }

  inputs_bss <- run_stage("prep_inputs_bss", {
    prep_inputs_bss(
      est_catch_group = chosen_ecg,
      period          = PERIOD_BSS,
      days            = dwg$days,
      dwg_summarized  = dwg_summ,
      census_expan    = dwg_summ$census_expan,
      study_design    = study_design,
      priors          = BSS_PRIORS
    )
  })

  # Stan hard-errors on any NA anywhere in `data` -- drop NA observation rows
  # (for expediency; see R_functions/drop_na_bss_inputs.R for what is and
  # isn't safe to drop) BEFORE the preflight checks below, so V_n/IntA/etc.
  # reflect the post-drop counts.
  inputs_bss <- run_stage("drop_na_bss_inputs", drop_na_bss_inputs(inputs_bss, fishery_name = fishery_name))
  na_drop_log <- attr(inputs_bss, "na_drop_log")
  if (!is.null(na_drop_log) && nrow(na_drop_log) > 0) {
    append_csv_row(na_drop_log, file.path(OUT_DIR, "bss_b_na_drops.csv"))
  }

  # BSS-specific preflight: converts a cryptic Stan crash into a ledger row.
  if (inputs_bss$G < 2) skip_fishery("Only one angler type; b[2]/lambda[...,2] out of bounds in the BSS likelihood.", stage = "bss_preflight")
  if (inputs_bss$V_n == 0) skip_fishery("No vehicle index counts; b[1] is not informed by data (would sample its prior).", stage = "bss_preflight")
  if (inputs_bss$IntA == 0 || sum(inputs_bss$V_A) == 0) skip_fishery("No interview vehicle/angler-group records; R_V unidentified, b[1] confounded with it.", stage = "bss_preflight")
  if (any(is.na(unlist(inputs_bss[c("V_I","T_I","A_A")])))) {
    cli::cli_abort("[bss_preflight] NA values in core count vectors for {.val {fishery_name}} -- data problem, not a skip.")
  }
  if (nrow(inputs_bss$p_TI) != inputs_bss$G || ncol(inputs_bss$p_TI) != inputs_bss$S) {
    cli::cli_abort("[bss_preflight] p_TI dims ({nrow(inputs_bss$p_TI)}x{ncol(inputs_bss$p_TI)}) don't match G x S ({inputs_bss$G}x{inputs_bss$S}) for {.val {fishery_name}}.")
  }
  b_trailer_informed <- inputs_bss$T_n > 0
  if (!b_trailer_informed) cli::cli_alert_warning("  No trailer index counts; b[2] will be prior-dominated.")
  b_weakly_informed <- inputs_bss$IntA < MIN_INTA_INFORMATIVE
  if (b_weakly_informed) cli::cli_alert_warning("  IntA = {inputs_bss$IntA} < {MIN_INTA_INFORMATIVE}; b likely weakly informed.")

  stan_dims_row <- tibble(
    fishery_name = fishery_name, D = inputs_bss$D, G = inputs_bss$G, S = inputs_bss$S, H = inputs_bss$H,
    V_n = inputs_bss$V_n, T_n = inputs_bss$T_n, A_n = inputs_bss$A_n, E_n = inputs_bss$E_n,
    IntC = inputs_bss$IntC, IntA = inputs_bss$IntA,
    b_trailer_informed = b_trailer_informed, b_weakly_informed = b_weakly_informed
  )
  append_csv_row(stan_dims_row, file.path(OUT_DIR, "bss_b_stan_dims.csv"))

  # --- Fit ------------------------------------------------------------------

  cfg <- FIT_CONFIGS[[fit_config_name]]
  t0 <- Sys.time()
  bss_fit <- run_stage("fit_bss", {
    fit_bss(
      model_file_name = BSS_MODEL_FILE, bss_inputs_list = inputs_bss,
      n_chain = cfg$n_chain, n_cores = cfg$n_cores, n_iter = cfg$n_iter, n_warmup = cfg$n_warmup,
      n_thin = cfg$n_thin, adapt_delta = cfg$adapt_delta, max_treedepth = cfg$max_treedepth,
      init = "0", pars = MONITOR_PARS, include = TRUE
    )
  })
  runtime_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  sp <- rstan::get_sampler_params(bss_fit, inc_warmup = FALSE)
  n_div <- sum(vapply(sp, function(x) sum(x[, "divergent__"]), numeric(1)))

  bias_summary <- run_stage("get_bss_bias", {
    get_bss_bias(bss_fit, fishery_name = fishery_name, ecg = chosen_ecg,
                 prior_sigma_b = BSS_PRIORS[["value_lognormal_sigma_b"]])
  }) |>
    mutate(
      fit_config = fit_config_name, model_file = BSS_MODEL_FILE, n_div = n_div, runtime_sec = runtime_sec,
      b_trailer_informed = b_trailer_informed, b_weakly_informed = b_weakly_informed,
      informed_flag = case_when(
        rhat > 1.05 | ess_bulk < 100 ~ "unconverged",
        prior_contraction < 0.10     ~ "prior-dominated",
        b_weakly_informed            ~ "weak",
        TRUE                          ~ "informed"
      )
    )
  append_csv_row(bias_summary, file.path(OUT_DIR, "bss_b_summary.csv"))

  draws <- posterior::as_draws_df(bss_fit) |> posterior::subset_draws(variable = "b")
  saveRDS(draws, file.path(DRAWS_DIR, paste0(safe_name(fishery_name), ".rds")))

  if (SAVE_FITS) saveRDS(bss_fit, file.path(FITS_DIR, paste0(safe_name(fishery_name), ".rds")))

  list(status = "ok", bias_summary = bias_summary, runtime_sec = runtime_sec)
}

# ------------------------------------------------------------------------------
# Driver: read discovery output, filter to include_in_run, run with ledger
# ------------------------------------------------------------------------------

discovery_path <- file.path(OUT_DIR, "fishery_discovery_target.csv")
if (!file.exists(discovery_path)) {
  cli::cli_abort("Run 00_discover_fisheries.R first -- {.file {discovery_path}} not found.")
}

target_fisheries <- read_csv(discovery_path, show_col_types = FALSE) |>
  filter(basin_match == "target", include_in_run) |>
  pull(fishery_name_raw) |>
  unique()

cli::cli_alert_info("{length(target_fisheries)} fishery-year(s) queued at fit_config = {.val {FIT_CONFIG_NAME}}.")

run_ledger <- map(target_fisheries, function(fn) {
  withCallingHandlers(
    tryCatch(
      {
        res <- fit_one_fishery(fn)
        tibble(fishery_name = fn, status = "ok", stage = NA_character_, reason = NA_character_,
               runtime_sec = res$runtime_sec)
      },
      fishery_skip = function(cnd) {
        cli::cli_alert_warning("Skipped [{.val {fn}}]: {conditionMessage(cnd)}")
        tibble(fishery_name = fn, status = "skipped", stage = cnd$stage %||% "unknown",
               reason = cnd$reason %||% conditionMessage(cnd), runtime_sec = NA_real_)
      },
      fishery_error = function(cnd) {
        cli::cli_alert_danger("Failed [{.val {fn}}] at stage {.val {cnd$stage}}: {cnd$reason}")
        tibble(fishery_name = fn, status = "error", stage = cnd$stage %||% "unknown",
               reason = cnd$reason %||% conditionMessage(cnd), runtime_sec = NA_real_)
      },
      error = function(e) {
        cli::cli_alert_danger("Failed [{.val {fn}}] (unstaged): {conditionMessage(e)}")
        tibble(fishery_name = fn, status = "error", stage = "unstaged",
               reason = conditionMessage(e), runtime_sec = NA_real_)
      }
    ),
    warning = function(w) { cli::cli_alert_info("  warning [{.val {fn}}]: {conditionMessage(w)}"); invokeRestart("muffleWarning") }
  )
}) |> bind_rows()

write_csv(run_ledger, file.path(OUT_DIR, "bss_b_run_ledger.csv"))

cli::cli_h2("Run outcomes")
run_ledger |> count(status) |> print()
if (any(run_ledger$status == "error"))   { cli::cli_h3("Errors");  run_ledger |> filter(status=="error")   |> print(n=50) }
if (any(run_ledger$status == "skipped")) { cli::cli_h3("Skipped"); run_ledger |> filter(status=="skipped") |> print(n=50) }

cli::cli_alert_success("Done. See analysis/bss_bias/outputs/ for bss_b_summary.csv, bss_b_comparability_raw.csv, bss_b_run_ledger.csv.")
cli::cli_alert_info("Next: Rscript analysis/bss_bias/02_build_comparability_table.R (fast, no fits needed) and 03_plot_b_series.R / 04_candidate_options.R once bss_b_summary.csv has rows.")

if (!is.null(DB_CONN)) try(DBI::dbDisconnect(DB_CONN), silent = TRUE)
