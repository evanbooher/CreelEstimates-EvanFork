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
SECTION_RESTRICTIONS <- list()

# Each entry holds a fishery's series to the water bodies present in EVERY year
# of it, so `b` averages over a constant set of water:
#
#   Skagit fall salmon       -- Cascade appears in 2025 only.
#   Skagit spring Chinook upper -- Cascade appears in 2024 and 2025 only.
#   Stillaguamish            -- mainstem and North Fork are in all four years;
#                               the South Fork is in 2022, 2023 and 2024 but
#                               not 2025. MS + NF is the largest constant set.
WATER_BODY_RESTRICTIONS <- list(
  list(pattern = regex("Skagit fall salmon", ignore_case = TRUE), keep = "Skagit"),
  list(pattern = regex("Skagit spring Chinook.*upper", ignore_case = TRUE), keep = "Skagit"),
  list(pattern = regex("Stillaguamish", ignore_case = TRUE),
       keep = c("Stillaguamish - MS", "Stillaguamish - NF"))
)

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
      "{.val {fishery_name}} has no sections in water bod{?y/ies} {.val {keep}}; \\
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
    cli::cli_alert_warning(
      "  Restriction is PARTIAL: kept section{?s} {.val {sort(unique(bleed$section_num))}} \\
       also carr{?ies/y} {.val {sort(unique(bleed$water_body_code))}}, which cannot be \\
       separated by section number."
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
    if (str_detect(fishery_name, r$pattern)) limits <- c(limits, list(as.double(r$sections)))
  }
  for (r in WATER_BODY_RESTRICTIONS) {
    if (str_detect(fishery_name, r$pattern)) {
      limits <- c(limits, list(sections_in_water_bodies(fishery_name, r$keep)))
    }
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

