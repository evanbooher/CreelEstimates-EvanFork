# ==============================================================================
# 10a_render_stilly.R -- render fw_creel for Stillaguamish 2024 and 2025 only
#
# A thin loop over 10_render_fw_creel.R, not a second driver. The reason it
# exists: 10 reads its RENDER_* settings with `exists(..., inherits = FALSE)`,
# so a value left in the global environment by one render silently carries into
# the next -- RENDER_PROJECT in particular, which 10 derives from the scope tag
# only when it is absent. Sourcing 10 twice by hand puts the second render's
# output in the first render's folder. This clears them between iterations.
#
# What it produces, per fishery-year: the full fw_creel report, which includes
# plot_census_index_ratio_daily -- the daily paired census-to-index ratio, the
# per-day version of the PE bias term (census anglers over index counts
# expanded to anglers, divided by p_census), faceted by section.
#
#   Rscript analysis/bss_bias/10a_render_stilly.R
#
# Each render is a full BSS fit. Two of them is not a quick job -- start it and
# check the first report before walking away.
#
# Scope: whole Stillaguamish by default, so the ratio plot shows every paired
# census day in the season rather than the three inside the matched comparison
# window. For the fork-scoped version instead, set
#
#   STILLY_SCOPE_TAGS <- c("MS", "NF")
#
# before sourcing, which renders four (fork x year) and applies the
# September 16 - October 31 window from scope_rules.R.
# ==============================================================================

library(cli)
library(here)

if (!exists("STILLY_YEARS", inherits = FALSE)) {
  STILLY_YEARS <- c("Stillaguamish salmon and gamefish 2024-25",
                    "Stillaguamish salmon and gamefish 2025-26")
}
# NULL means no fork scope: the whole fishery, and the lookup's own window.
if (!exists("STILLY_SCOPE_TAGS", inherits = FALSE)) STILLY_SCOPE_TAGS <- NULL

DRIVER <- here::here("analysis", "bss_bias", "10_render_fw_creel.R")
if (!file.exists(DRIVER)) cli::cli_abort("{.file {DRIVER}} not found.")

# Everything 10 sets or reads that must not survive into the next iteration.
CARRIERS <- c("RENDER_FISHERY_RE", "RENDER_SCOPE_TAG", "RENDER_PROJECT",
              "RENDER_GROUPS", "RENDER_SKIP_DONE", "RENDER_FIT_ONLY",
              "RENDER_BASIN_ORDER", "RUN_SCOPE")

clear_carriers <- function() {
  rm(list = intersect(ls(.GlobalEnv), CARRIERS), envir = .GlobalEnv)
}

render_one <- function(fishery_name, scope_tag) {
  clear_carriers()
  # 10 matches on a regex, so anchor it -- "2024-25" would otherwise be a
  # substring of nothing today but of a renamed fishery-year tomorrow.
  assign("RENDER_FISHERY_RE", paste0("^", fishery_name, "$"), envir = .GlobalEnv)
  if (!is.null(scope_tag)) assign("RENDER_SCOPE_TAG", scope_tag, envir = .GlobalEnv)

  cli::cli_h1("Rendering {.val {fishery_name}}{if (is.null(scope_tag)) '' else paste0(' [', scope_tag, ']')}")
  source(DRIVER)
  invisible(NULL)
}

combos <- if (is.null(STILLY_SCOPE_TAGS)) {
  lapply(STILLY_YEARS, function(f) list(fishery_name = f, scope_tag = NULL))
} else {
  unlist(lapply(STILLY_SCOPE_TAGS, function(s)
    lapply(STILLY_YEARS, function(f) list(fishery_name = f, scope_tag = s))),
    recursive = FALSE)
}

cli::cli_alert_info("{length(combos)} render{?s} queued.")
for (cb in combos) render_one(cb$fishery_name, cb$scope_tag)
clear_carriers()

cli::cli_alert_success("Done. Reports in analysis/bss_bias/outputs/reports/; \\
                        the daily ratio figure is in each report.")
