# ==============================================================================
# scope_rules.R
#
# WHICH WATER a fishery-year's `b` is fit to. Sourced by 01_fit_bss_bias.R,
# which applies the rules when building model inputs, and by
# 02a_location_lut_changes.R, which reports year-over-year continuity both
# as-defined and as-fitted.
#
# One file so the two cannot disagree. A continuity measure computed on the
# full lookup describes water `b` never saw -- Skagit spring Chinook upper
# looks like it gained five sites in 2024 when in fact it gained the Cascade,
# which the fits exclude.
#
# Requires, from the sourcing script: dplyr, stringr, readr, cli, here.
# ==============================================================================

# ------------------------------------------------------------------------------
# Section restrictions
# ------------------------------------------------------------------------------
# Two independent reasons to hold a fishery-year to a subset of its sections.
#
# 1. COMPARABILITY (WATER_BODY_RESTRICTIONS). `b` is a single pooled scalar --
#    `vector[G] b` in the Stan model, no section index, applied outside the
#    section-indexed terms of the V_I/T_I likelihoods. So it averages over
#    every block the fishery-year happens to contain. A water body that appears
#    in one year and not the others therefore enters that year's `b` and no
#    other, and the series is no longer comparing like with like.
#
#    Skagit fall salmon is the case: the Cascade River appears only in 2025
#    (one block, section 7). Dropping it holds every year of that series to the
#    Skagit mainstem.
#
#    Keyed on WATER BODY, not section number, because section numbers are
#    labels on blocks and the blocks get redrawn -- the Cascade is section 7 in
#    2025 and does not exist in any other year. The section numbers are
#    resolved per fishery-year from the committed location lookup, so this stays
#    reproducible on the no-VPN path.
#
# 2. INDEXING (SECTION_RESTRICTIONS). Kept as a mechanism, currently EMPTY.
#    prep_inputs_bss() sizes the Stan arrays with length(unique(section_num))
#    but indexes them with the RAW section_num, so any gap reads out of range.
#    Stillaguamish was held to 1-6 for that reason -- a stop-gap to get the
#    model running, not a scope decision with a basis. It has been removed:
#    align_bss_sections() now closes numbering gaps properly by renumbering to
#    a dense 1..S, and the Stillaguamish comparability question is answered by
#    a fork rule under (1), which 1-6 did not do (it dropped the South Fork and
#    NF section 7 in 2022, only the South Fork in 2024, and nothing at all in
#    2023 or 2025).
#
# Both are SCOPE DECISIONS: the excluded water is gone from the fit entirely, so
# the resulting `b` describes the retained reach only. The ledger records the
# resolved sections per fishery-year in `sections_limited_to`. Neither is the
# fix for the indexing bug -- that is align_bss_sections(), which runs
# afterwards and handles the gaps these restrictions do not.

# Empty on purpose -- see (2) above. Add an entry only for a numbering problem
# align_bss_sections() genuinely cannot resolve, never for scope.
SECTION_RESTRICTIONS <- list(
  # NF 2024-25 held to section 4 ONLY, and only under the NF scope.
  #
  # The water-body rule resolves NF 2024-25 to sections 4, 5 and 6. Sections 5
  # and 6 cover rm 9.5-37.5, which has no counterpart in 2025-26 at all -- so
  # keeping them would put 28 river miles into 2024's b and nothing comparable
  # into 2025's. Section 4 (Mouth of the NF to Cicero Bridge, rm 0-9.5) is the
  # same census block as 2025-26's section 3, with the same three index sites.
  #
  # Two further reasons, both from the effort data rather than the lookup:
  #   * s5 and s6 have ONE paired census day each in the comparison window
  #     (Oct 22), against s4's three. They add almost no anchor.
  #   * s6 carries TWO census blocks under one section_num (rm 24.5-30 and
  #     30-37.5). prep_inputs_bss() sizes S, O and p_TI off section_num, so that
  #     is not a shape the model is built for.
  #
  # Gated on the NF tag because it composes by INTERSECTION with whatever else
  # applies: ungated, it would meet the MS scope's {2, 3} and resolve to
  # nothing, aborting the MS run.
  list(pattern   = regex("^Stillaguamish salmon and gamefish 2024-25$"),
       scope_tag = "NF",
       sections  = 4)
)

# Each entry holds a fishery's series to the water bodies present in EVERY year
# of it, so `b` averages over a constant set of water:
#
#   Skagit fall salmon          -- Cascade appears in 2025 only.
#   Skagit spring Chinook upper -- Cascade appears in 2024 and 2025 only.
#
# STILLAGUAMISH IS DELIBERATELY UNRESTRICTED, though its South Fork is present
# in 2022, 2023 and 2024 and absent in 2025. An MS + NF rule was tried and
# removed: it could not be applied evenly. The rule resolves to section
# numbers, and 2023-24 merges the North and South Forks into a single census
# block (section 5), so that year kept five South Fork index sites while 2022
# and 2024 lost theirs. A rule that excludes the South Fork from two years and
# retains it in a third creates the comparability break it exists to prevent.
#
# Nor can it be fixed by dropping those five sites: section 5's census count
# totals both forks, so index counts covering only the North Fork against a
# census covering both would be absorbed by `b` as a low bias. That is a
# corrupted estimate, not a restricted one.
#
# So all four years run whole-basin, and the 2025 South Fork absence is a
# STATED LIMITATION rather than an unevenly-applied rule. See
# bss_b_lut_scope_effect.csv and bss_b_lut_water_body_year.csv for what each
# year actually covers.
#
# FOLLOW-UP, not yet built: fork-specific `b` -- a mainstem-only and a
# North-Fork-only estimate -- is the scope that matches the proposed fishery
# (North Fork gamefish, mainstem coho). It needs more than a rule here, because
# every output is keyed on fishery_name and a second scope for the same
# fishery-year would overwrite the first. The contained way in is a run-level
# scope tag that suffixes the output key, so "Stillaguamish ... 2024-25" and
# "Stillaguamish ... 2024-25 [MS]" coexist in the ledger and the b summary.
WATER_BODY_RESTRICTIONS <- list(
  list(pattern = regex("Skagit fall salmon", ignore_case = TRUE), keep = "Skagit"),
  list(pattern = regex("Skagit spring Chinook.*upper", ignore_case = TRUE), keep = "Skagit")
  # STILLAGUAMISH IS UNRESTRICTED, 2022-23 included -- per direction, the North
  # and South Forks stay in for every year. That year's problem is its date
  # window, handled below, not its water.
  #
  # Recorded because it settles a question the note above leaves open. The note
  # rejects an MS + NF rule, correctly: it cannot be applied evenly, because
  # 2023-24's section 5 is on both the North and South Forks. A MAINSTEM-ONLY
  # rule has no such problem -- per the location lookup the mainstem shares no
  # section with either fork in any year:
  #
  #   2022-23  MS 1,2,3   NF 4-7     SF 8-9
  #   2023-24  MS 1,2     NF 3,4,5   SF 5
  #   2024-25  MS 2,3     NF 4,5,6   SF 7
  #   2025-26  MS 1,2     NF 3
  #
  # So `list(pattern = regex("Stillaguamish"), keep = "Stillaguamish - MS")`
  # would hold all four years to the mainstem evenly, which is the scope that
  # matches the proposed fishery (mainstem coho). That is the shape the
  # fork-specific follow-up should take; it changes every Stillaguamish `b`, so
  # it belongs to that decision rather than to this one.
)

# ------------------------------------------------------------------------------
# Window restrictions -- WHICH DAYS a fishery-year's `b` is fit to
# ------------------------------------------------------------------------------
# Same idea as the section rules, on the other axis. The lookup window is what
# the fishery was open for, which is not always what was surveyed: a long tail
# of unsampled days still enters prep_days()'s day grid, and the BSS estimates
# effort for them from the priors alone. Those days carry no index counts, so
# they contribute nothing to `b` directly -- but they inflate the season totals
# that `b` is then used to rescale.
#
# `trim_to_last_sampled` pulls the end back to the last day actually surveyed
# within the window, so the end date is observed rather than chosen.
#
# One definition, three consumers: 01_fit_bss_bias.R (the `b` fit),
# 10_render_fw_creel.R (the production render) and anything else that resolves a
# window. Without this the b in the brief comes from one window and the season
# totals from another.
WINDOW_RESTRICTIONS <- list(
  # Stillaguamish 2022-23 runs 2022-09-01 to 2022-11-30, 91 days -- the longest
  # in the series, against 76 for 2024-25 and 61 for 2025-26 -- with a tail of
  # late-season dates the other years do not have. Held to September, then
  # trimmed to the last day sampled in the retained (mainstem) sections.
  list(pattern = regex("^Stillaguamish salmon and gamefish 2022-23$"),
       est_date_start = "2022-09-01",
       est_date_end   = "2022-09-30",
       trim_to_last_sampled = TRUE),

  # THE FORK COMPARISON WINDOW -- Sep 16 to Oct 31, both years, scope-gated.
  #
  # The common period of the two years' own survey seasons: 2024-25's lookup
  # window starts 09-16 and 2025-26's ends 10-31, so this is their intersection.
  # Inside it the two years match almost exactly -- three paired census days per
  # section in each (Sep 25 and Oct 5 are the same calendar dates in both), the
  # same census blocks, and the same index sites.
  #
  # scope_tag keeps this OFF the whole-basin fits of the same two fishery-years,
  # whose b feeds the collaborator brief.
  list(pattern = regex("^Stillaguamish salmon and gamefish 2024-25$"),
       scope_tag = c("MS", "NF"),
       est_date_start = "2024-09-16",
       est_date_end   = "2024-10-31",
       trim_to_last_sampled = TRUE),
  list(pattern = regex("^Stillaguamish salmon and gamefish 2025-26$"),
       scope_tag = c("MS", "NF"),
       est_date_start = "2025-09-16",
       est_date_end   = "2025-10-31",
       trim_to_last_sampled = TRUE)
)

# Returns list(est_date_start, est_date_end, trim_to_last_sampled) or NULL.
# First match wins; a second matching rule is an authoring error, not something
# to silently compose, because two windows have no sensible intersection here.
fishery_window_limit <- function(fishery_name) {
  hits <- Filter(function(r) rule_applies(r, fishery_name), WINDOW_RESTRICTIONS)
  if (length(hits) == 0) return(NULL)
  if (length(hits) > 1) {
    cli::cli_abort(
      "{.val {fishery_name}} matches {length(hits)} window restrictions -- \
       make the patterns disjoint."
    )
  }
  hits[[1]]
}

# ------------------------------------------------------------------------------
# Catch groups a fishery-year should not be fitted for
# ------------------------------------------------------------------------------
# 00d's inventory already drops a group with zero fish, but it counts over the
# fishery's FULL window. Narrow the window and a group that had fish can have
# none left -- and fitting it is then a full MCMC run to learn that the
# posterior is the prior truncated by having seen nothing.
#
# `b` is invariant to the catch group (the effort and catch sub-models share no
# parameters), so excluding a group changes nothing about `b` -- only what catch
# is reported alongside it.
CATCH_GROUP_EXCLUSIONS <- list(
  # No Chinook were caught in the September part of the Stillaguamish 2022-23
  # fishery, which is all that survives the window restriction above.
  list(pattern = regex("^Stillaguamish salmon and gamefish 2022-23$"),
       groups = "chinook_all"),
  # No Chinook reported in either fork-comparison year. Stated as a rule rather
  # than left to 00d's inventory, which counts over the full window and would
  # not see a group emptied by the window restriction above.
  list(pattern = regex("^Stillaguamish salmon and gamefish 202[45]-"),
       groups = "chinook_all")
)

fishery_excluded_groups <- function(fishery_name) {
  out <- character(0)
  for (r in CATCH_GROUP_EXCLUSIONS) {
    if (rule_applies(r, fishery_name)) out <- union(out, r$groups)
  }
  out
}

# ------------------------------------------------------------------------------
# RUN SCOPE -- fit `b` to part of a fishery rather than all of it
# ------------------------------------------------------------------------------
# The standing rules above hold a SERIES comparable across its own years. This
# is different: a whole run narrowed to particular water, with its outputs
# TAGGED so they sit alongside the whole-fishery ones instead of overwriting
# them. Every output in 01 is keyed on the fishery name, and the tag becomes
# part of that key -- so "Stillaguamish ... 2024-25" and
# "Stillaguamish ... 2024-25 [MS]" coexist in the ledger, the b summary, the
# draws directory and everything downstream of them.
#
# The motivating case: the proposed Stillaguamish fishery is North Fork
# gamefish and mainstem coho, so a mainstem-only or North-Fork-only `b` is the
# scope that matches it. `b` is a single pooled scalar, so a whole-basin fit
# averages the mainstem together with water that fishery will not touch.
#
# To run one, set this in 01_fit_bss_bias.R BEFORE it sources this file, or
# edit it here, and re-run:
#
#   RUN_SCOPE <- list(tag = "MS",
#                     pattern = regex("Stillaguamish", ignore_case = TRUE),
#                     keep    = "Stillaguamish - MS")
#
# EXPECT FEWER ANCHORS. Census counts set the scale of the latent effort that
# index counts are measured against, and are recorded only alongside a same-day
# index count in the same section -- 2 to 6 such days per year for the whole
# Stillaguamish basin, and a fork-only scope keeps a fraction of those. Read
# prior_contraction in bss_b_summary.csv before drawing anything from the
# result: a fork-scoped `b` may be mostly prior.
#
# Left alone by whatever sources this file first, so 01 can set it beforehand
# and 02a (which has no run scope) still gets a definition.
if (!exists("RUN_SCOPE", inherits = FALSE)) RUN_SCOPE <- NULL

# Named scopes, so a run is `RUN_SCOPE <- SCOPE_PRESETS$MS` rather than a
# hand-edited regex each time -- the fork comparison needs two runs that differ
# in exactly one field, and hand-editing is how they drift.
SCOPE_PRESETS <- list(
  MS = list(tag = "MS",
            pattern = regex("Stillaguamish salmon and gamefish 202[45]-"),
            keep    = "Stillaguamish - MS"),
  NF = list(tag = "NF",
            pattern = regex("Stillaguamish salmon and gamefish 202[45]-"),
            keep    = "Stillaguamish - NF")
)

# Does a rule apply to this fishery-year, under the scope currently in force?
#
# `scope_tag` on a rule means "only when RUN_SCOPE has one of these tags". It is
# what keeps the fork comparison's narrower window off the whole-basin series:
# the same fishery_name is fitted both ways, and re-cutting the whole-basin b as
# a side effect of a different question would silently move the numbers already
# in the collaborator brief.
rule_applies <- function(r, fishery_name) {
  if (!str_detect(fishery_name, r$pattern)) return(FALSE)
  if (is.null(r$scope_tag)) return(TRUE)
  !is.null(RUN_SCOPE) && RUN_SCOPE$tag %in% r$scope_tag
}

# The name an output is filed under. Identical to fishery_name when RUN_SCOPE
# is NULL or does not match, so a default run is exactly what it always was.
output_name <- function(fishery_name) {
  if (is.null(RUN_SCOPE) || !str_detect(fishery_name, RUN_SCOPE$pattern)) return(fishery_name)
  paste0(fishery_name, " [", RUN_SCOPE$tag, "]")
}

# Committed by 00c_probe_location_lut.R. Read lazily and cached: a fishery-name
# run with no water-body rule never needs it, and a run of 27 fisheries should
# not read the same file 27 times.
LOCATION_LUT_PATH <- here::here("analysis", "bss_bias", "lookup", "fishery_location_lut.csv")
.location_lut <- NULL
location_lut <- function() {
  if (is.null(.location_lut)) {
    if (!file.exists(LOCATION_LUT_PATH)) {
      cli::cli_abort(c(
        "A water-body restriction needs {.file {LOCATION_LUT_PATH}}, which is missing.",
        "i" = "Run {.file analysis/bss_bias/00c_probe_location_lut.R} (VPN) and commit the result,",
        "i" = "or clear {.code WATER_BODY_RESTRICTIONS} to run without it."
      ))
    }
    .location_lut <<- readr::read_csv(LOCATION_LUT_PATH, show_col_types = FALSE)
  }
  .location_lut
}

# Sections of `fishery_name` that sit in one of `keep`, per the location lookup.
sections_in_water_bodies <- function(fishery_name, keep) {
  lut <- location_lut()
  rows <- lut |> dplyr::filter(.data$fishery_name == !!fishery_name)
  if (nrow(rows) == 0) {
    cli::cli_abort(
      "{.val {fishery_name}} has no rows in the location lookup, so its \\
       water-body restriction cannot be resolved. Re-run 00c, or exempt this \\
       fishery from WATER_BODY_RESTRICTIONS."
    )
  }
  kept <- rows |> dplyr::filter(.data$water_body_code %in% keep)
  if (nrow(kept) == 0) {
    cli::cli_abort(
      "{.val {fishery_name}} has no sections in water bodies {.val {keep}}; \\
       restricting would leave nothing to fit."
    )
  }
  dropped <- setdiff(unique(rows$water_body_code), keep)
  if (length(dropped) > 0) {
    cli::cli_alert_info("  Water-body scope: keeping {.val {keep}}, dropping {.val {dropped}}.")
  }

  keep_sections <- sort(unique(as.double(kept$section_num)))

  # A restriction expressed in water bodies has to be APPLIED in section
  # numbers, and a section can straddle two water bodies -- Stillaguamish
  # 2023-24 section 5 is on both the North and South Forks. Keeping that
  # section therefore keeps some of the water the rule meant to drop. Say so:
  # the restriction is partial for that fishery-year, and a `b` from it is not
  # quite the clean comparison the rule was written to produce.
  bleed <- rows |>
    dplyr::filter(as.double(.data$section_num) %in% keep_sections,
                  !.data$water_body_code %in% keep) |>
    dplyr::distinct(section_num, water_body_code)
  if (nrow(bleed) > 0) {
    # Stated without cli pluralization: two {?} markers over two different
    # quantities in one string is ambiguous, and cli errors rather than guessing.
    bleed_sections <- sort(unique(bleed$section_num))
    bleed_water    <- sort(unique(bleed$water_body_code))
    cli::cli_alert_warning(
      "  Restriction is PARTIAL. Kept sections {.val {bleed_sections}} also carry \\
       {.val {bleed_water}}, which cannot be separated by section number."
    )
  }

  keep_sections
}

# The two rule sets COMPOSE: a fishery matching both is held to the
# intersection. Returns NULL when neither applies, which restrict_dwg_sections()
# reads as "use everything".
fishery_section_limit <- function(fishery_name) {
  limits <- list()
  for (r in SECTION_RESTRICTIONS) {
    if (rule_applies(r, fishery_name)) limits <- c(limits, list(as.double(r$sections)))
  }
  for (r in WATER_BODY_RESTRICTIONS) {
    if (rule_applies(r, fishery_name)) {
      limits <- c(limits, list(sections_in_water_bodies(fishery_name, r$keep)))
    }
  }
  # A run-level scope composes with the standing rules like any other limit.
  if (!is.null(RUN_SCOPE) && str_detect(fishery_name, RUN_SCOPE$pattern)) {
    limits <- c(limits, list(sections_in_water_bodies(fishery_name, RUN_SCOPE$keep)))
  }
  if (length(limits) == 0) return(NULL)
  keep <- Reduce(intersect, limits)
  if (length(keep) == 0) {
    cli::cli_abort(
      "Section restrictions for {.val {fishery_name}} intersect to nothing -- \\
       check SECTION_RESTRICTIONS against WATER_BODY_RESTRICTIONS."
    )
  }
  sort(keep)
}

