# ==============================================================================
# 13_snohomish_scope_compare.R -- Snohomish mainstem vs Skykomish, 2023 vs 2024
#
# Same question, same shape, as 12_fork_scope_compare.R's Stillaguamish
# comparison -- whether a reach-specific `b` can be estimated from a common
# window well enough to set 2026's prior -- adapted for a different basin
# with a genuinely different section layout. Kept as a SEPARATE file rather
# than generalising 12 in place: 12's FITS/CENSUS_BOAT_FLOOR/labels are
# Stillaguamish-specific, and this basin's own scopes (SN_MAIN, SN_SKY) and
# window (2023, 2024; Sep 18 - Oct 31) are defined in scope_rules.R
# alongside Stillaguamish's MS/NF, not duplicated here.
#
# WHY SEPARATE OUTPUT FILES, not bss_b_scope_effort.csv /
# bss_b_scope_compare.csv: 12's write_csv() calls OVERWRITE those files
# whole, not append. Running this script after 12 (or vice versa) would
# silently wipe out the other basin's audit rows -- and the Stillaguamish
# memo already reads bss_b_scope_effort.csv by that exact name. Writing to
# *_snohomish.csv avoids that collision entirely rather than requiring the
# two scripts to run in a fixed order forever.
#
#   AUDIT   -- what was actually surveyed under each scope. Reads the cached
#              DWG if present, live otherwise (VPN required for a fresh
#              basin/year not yet cached -- this script has not been run
#              from this environment, which has neither).
#   COMPARE -- the four fitted `b` values (SN_MAIN x 2023/2024, SN_SKY x
#              2023/2024) with those diagnostics beside them. Skipped, with
#              a note, until 01 has been run under both scopes.
#
#   Rscript analysis/bss_bias/13_snohomish_scope_compare.R
#
# Prerequisite fits (01_fit_bss_bias.R, once per scope, on "Snohomish fall
# salmon 2023" and "Snohomish fall salmon 2024"):
#   RUN_SCOPE <- SCOPE_PRESETS$SN_MAIN   # then source 01 for each year
#   RUN_SCOPE <- SCOPE_PRESETS$SN_SKY    # then source 01 for each year
#
# Outputs (analysis/bss_bias/outputs/):
#   bss_b_scope_sites_snohomish.csv     index sites per section x fishery-year
#   bss_b_scope_effort_snohomish.csv    anchors, census anglers, index objects per fit
#   bss_b_scope_compare_snohomish.csv   the four b values + diagnostics   [after fits]
#   figures/fig21_snohomish_scope_b     b with 95% CrI, Mainstem vs Skykomish x year   [after fits]
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(cli); library(here)
})

OUT_DIR <- here::here("analysis", "bss_bias", "outputs")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

source(here::here("analysis", "bss_bias", "common.R"))
source(here::here("analysis", "bss_bias", "fishery_data.R"))
source(here::here("analysis", "bss_bias", "scope_rules.R"))

# Boat anglers in census, below which the trailer term is flagged as
# data-limited. Same floor and same reasoning as 12_fork_scope_compare.R --
# "enough to be a measurement at all", not a tuned threshold.
CENSUS_BOAT_FLOOR <- 20

# The four fits. Sections are resolved per scope from scope_rules.R (SN_MAIN,
# SN_SKY), same as 12 does for MS/NF -- this cannot disagree with what 01
# actually fits.
FITS <- tibble::tribble(
  ~scope,    ~fishery_name,
  "SN_MAIN", "Snohomish fall salmon 2023",
  "SN_MAIN", "Snohomish fall salmon 2024",
  "SN_SKY",  "Snohomish fall salmon 2023",
  "SN_SKY",  "Snohomish fall salmon 2024"
) |>
  mutate(year = as.integer(str_extract(fishery_name, "\\d{4}")),
         fit  = paste(scope, year))

# Census counts are by angler type; index counts are by object type. bank/boat
# comes from angler_final where the fetch provides it, and from the raw label
# otherwise -- "Boat - Motor" is a boat angler, "Shore - Motor Boat" is not.
bank_or_boat <- function(d) {
  if ("angler_final" %in% names(d)) return(as.character(d$angler_final))
  lbl <- as.character(d$count_type %||% rep(NA_character_, nrow(d)))
  dplyr::if_else(str_starts(lbl, "Boat"), "boat", "bank")
}

# ------------------------------------------------------------------------------
# AUDIT -- one scope-resolved pull per fit
# ------------------------------------------------------------------------------

audit_one <- function(scope_tag, fishery_name) {
  RUN_SCOPE <<- SCOPE_PRESETS[[scope_tag]]

  est_dates <- resolve_window(fishery_name)
  win <- fishery_window_limit(fishery_name)
  if (!is.null(win)) {
    est_dates$est_date_start <- win$est_date_start
    est_dates$est_date_end   <- win$est_date_end
  }
  keep <- fishery_section_limit(fishery_name)

  dwg <- fetch_fishery_dwg(fishery_name, est_dates)
  eff <- dwg$effort |>
    mutate(event_date = as.Date(event_date),
           section_num = suppressWarnings(as.double(section_num))) |>
    filter(between(event_date,
                   as.Date(est_dates$est_date_start),
                   as.Date(est_dates$est_date_end)))
  if (!is.null(keep)) eff <- eff |> filter(section_num %in% keep)
  eff$gear <- bank_or_boat(eff)

  is_census <- eff$tie_in_indicator %in% c(1, TRUE, "TRUE", "true")
  cen <- eff[is_census, , drop = FALSE]
  idx <- eff[!is_census, , drop = FALSE]

  # A section's census anchors are the days it has BOTH a census and an index
  # count. Census alone is dropped by the model; index alone carries no anchor.
  anchors <- inner_join(
    cen |> distinct(section_num, event_date),
    idx |> distinct(section_num, event_date),
    by = c("section_num", "event_date")
  )

  list(
    sites = idx |>
      distinct(section_num, location) |>
      arrange(section_num, location) |>
      mutate(scope = scope_tag, fishery_name = fishery_name, .before = 1),
    effort = tibble(
      scope = scope_tag, fishery_name = fishery_name,
      window_start = est_dates$est_date_start,
      window_end   = est_dates$est_date_end,
      sections = paste(sort(unique(eff$section_num)), collapse = ", "),
      n_index_sites = n_distinct(idx$location),
      n_index_days  = n_distinct(idx$event_date),
      n_census_days = n_distinct(cen$event_date),
      n_paired_anchors = n_distinct(anchors$event_date),
      anchor_dates = paste(sort(unique(as.character(anchors$event_date))), collapse = ", "),
      census_anglers_bank = sum(cen$count_quantity[cen$gear == "bank"], na.rm = TRUE),
      census_anglers_boat = sum(cen$count_quantity[cen$gear == "boat"], na.rm = TRUE),
      index_vehicles = sum(idx$count_quantity[idx$count_type == "Vehicle Only"], na.rm = TRUE),
      index_trailers = sum(idx$count_quantity[idx$count_type == "Trailers Only"], na.rm = TRUE)
    )
  )
}

cli_h1("13 -- Snohomish mainstem vs Skykomish, 2023 vs 2024")
.orig_scope <- RUN_SCOPE
audits <- pmap(list(FITS$scope, FITS$fishery_name), audit_one)
RUN_SCOPE <- .orig_scope

sites  <- map_dfr(audits, "sites")
effort <- map_dfr(audits, "effort") |>
  mutate(year = as.integer(str_extract(fishery_name, "\\d{4}")),
         fit  = paste(scope, year), .before = 1)

write_csv(sites,  file.path(OUT_DIR, "bss_b_scope_sites_snohomish.csv"))
write_csv(effort, file.path(OUT_DIR, "bss_b_scope_effort_snohomish.csv"))

cli_h2("Index sites by section")
sites |> select(fit = scope, fishery_name, section_num, location) |> print(n = Inf)

cli_h2("What each fit has to work with")
effort |>
  select(fit, sections, n_index_days, n_paired_anchors,
         census_anglers_bank, census_anglers_boat, index_vehicles, index_trailers) |>
  print(n = Inf)

# The trailer channel is measured against boat anglers in the CENSUS, not
# against trailers in the index: however many trailers were counted, there is
# nothing to compare them to without boat anglers in the census.
thin <- effort |> filter(census_anglers_boat < CENSUS_BOAT_FLOOR)
if (nrow(thin) > 0) {
  n_boat_floor <- CENSUS_BOAT_FLOOR
  cli_alert_warning(
    "Trailer term is DATA-LIMITED for {.val {thin$fit}}: fewer than \\
     {n_boat_floor} boat anglers counted in census across the whole window."
  )
}

# ------------------------------------------------------------------------------
# COMPARE -- the four b values, with the audit beside them
# ------------------------------------------------------------------------------

# IN A FUNCTION for the same reason as 12: quit() ends the whole RStudio
# session when this is sourced from the console rather than run via Rscript,
# not just the script. return() exits the same early-out without that
# difference in behaviour between the two ways of running it.
summary_path <- file.path(OUT_DIR, "bss_b_summary.csv")

run_comparison <- function() {

  if (!file.exists(summary_path)) {
    cli_alert_info(
      "No {.file bss_b_summary.csv} yet -- audit tables written, comparison skipped."
    )
    cli_alert_info("Fit first: {.code RUN_SCOPE <- SCOPE_PRESETS$SN_MAIN} then source 01, and again for $SN_SKY.")
    return(invisible(NULL))
  }

  b_all <- read_csv(summary_path, show_col_types = FALSE)

  # 01 files scoped results under "<fishery_name> [TAG]" -- see output_name()
  # in scope_rules.R. Matching on that is what keeps the SN_MAIN/SN_SKY fits
  # from being confused with any whole-basin Snohomish fit of the same two
  # fishery-years (which bss_b_summary.csv already has, from the earlier
  # variability analysis -- unscoped, so a different key).
  scoped <- FITS |>
    mutate(out_name = paste0(fishery_name, " [", scope, "]")) |>
    left_join(b_all, by = c("out_name" = "fishery_name"))

  fitted_rows <- scoped |> filter(!is.na(median))
  if (nrow(fitted_rows) == 0) {
    cli_alert_warning("No fitted b for any of the four scopes yet -- comparison skipped.")
    cli_alert_info("Audit tables are written. Fit with {.code RUN_SCOPE <- SCOPE_PRESETS$SN_MAIN} then 01, and again for {.code $SN_SKY}.")
    return(invisible(NULL))
  }
  missing <- scoped |> filter(is.na(median)) |> distinct(fit)
  if (nrow(missing) > 0) {
    cli_alert_warning("No fitted b yet for: {.val {missing$fit}}")
  }

  cmp <- fitted_rows |>
    left_join(effort |> select(fit, n_paired_anchors, census_anglers_bank,
                               census_anglers_boat, index_vehicles, index_trailers),
              by = "fit") |>
    mutate(
      limitation = case_when(
        bias_type == "trailer" & census_anglers_boat < CENSUS_BOAT_FLOOR ~
          "data-limited: few boat anglers in census",
        informed_flag == "unconverged"    ~ "data-limited: did not converge",
        informed_flag == "prior-dominated" ~ "data-limited: posterior tracks the prior",
        informed_flag == "weak"            ~ "data-limited: few angler interviews",
        TRUE                               ~ ""
      ),
      data_limited = limitation != "",
      b_ci = sprintf("%.2f (%.2f-%.2f)", median, q2.5, q97.5)
    ) |>
    select(fit, scope, year, bias_type, b = median, b_lo = q2.5, b_hi = q97.5,
           b_ci, data_limited, limitation, informed_flag, prior_contraction,
           n_paired_anchors, census_anglers_bank, census_anglers_boat,
           index_vehicles, index_trailers) |>
    arrange(bias_type, scope, year)

  write_csv(cmp, file.path(OUT_DIR, "bss_b_scope_compare_snohomish.csv"))
  cli_alert_success("Comparison written ({nrow(cmp)} row{?s}).")

  cli_h2("b by scope and year")
  cmp |> select(fit, bias_type, b_ci, limitation, n_paired_anchors,
                census_anglers_boat) |> print(n = Inf)

  # --- The figure ---------------------------------------------------------------
  # Vehicle only, same reasoning as 12: a trailer panel the table already
  # flags as data-limited would be read as a result by anyone who saw the
  # figure alone.
  fig_df <- cmp |> filter(bias_type == "vehicle")

  if (nrow(fig_df) > 0) {
    fig21 <- fig_df |>
      mutate(scope_label = if_else(scope == "SN_MAIN", "Mainstem", "Skykomish"),
             year_label  = factor(year)) |>
      ggplot(aes(x = year_label, y = b, colour = scope_label)) +
      geom_hline(yintercept = 1, colour = BASELINE_COL, linewidth = 0.4) +
      geom_linerange(aes(ymin = b_lo, ymax = b_hi),
                     position = position_dodge(width = 0.4), linewidth = 1.6, alpha = 0.9) +
      geom_point(position = position_dodge(width = 0.4), size = 2.8) +
      geom_text(aes(label = paste0(n_paired_anchors, " anchors")),
                position = position_dodge(width = 0.4),
                vjust = -1.4, size = 3, show.legend = FALSE) +
      scale_colour_manual(values = c(Mainstem = CAT[["blue"]], Skykomish = CAT[["aqua"]]),
                          name = NULL) +
      labs(title = "Vehicle-index bias term by reach, Sep 18 - Oct 31",
           x = NULL, y = "b (vehicle)") +
      theme_bss()

    save_fig(fig21, "fig21_snohomish_scope_b", width = 8, height = 6)
    cli_alert_success("{.file fig21_snohomish_scope_b}")
  } else {
    cli_alert_info("No vehicle b available yet -- figure skipped.")
  }

  cli_rule()
  cli_alert_info("Tables: {.file bss_b_scope_sites_snohomish.csv}, {.file bss_b_scope_effort_snohomish.csv}, {.file bss_b_scope_compare_snohomish.csv}")
}

run_comparison()
