# ==============================================================================
# 02a_location_lut_changes.R
#
# Purpose:
#   Answer, precisely, what changed about a fishery's location definition
#   across years -- and what did not. Reads only the committed capture from
#   00c_probe_location_lut.R, so it needs no VPN, no database and no network.
#
# ------------------------------------------------------------------------------
# THE DATA MODEL, AND WHY IT DECIDES WHAT TO MEASURE
#
#   SITE     -- a discrete location where INDEX counts occur.
#   SECTION  -- the block sites are nested in. CENSUS (tie-in) counts total
#               EVERYTHING in the block, paired with a scheduled index count.
#               Interviews can occur anywhere in the block, not only at sites.
#
#   Sections are the spatial unit the data aggregate to for estimation. They
#   often follow a real boundary -- a sport-fishing pamphlet break, a CRC
#   break, a closure line -- but not necessarily.
#
#   That structure is the whole reason this script matters to `b`. `b` is the
#   ratio linking an index count to the census count it is paired with, so:
#
#     * the SITES in a block determine what the index count sees;
#     * the SECTION boundary determines what the census count totals;
#     * `b` is what absorbs the difference between the two.
#
#   So re-drawing a section is NOT neutral for `b`, even when every site is
#   kept. Split one block into two and each new block has its own
#   sites-to-block coverage, which need not equal the original's -- and `b` is
#   a single scalar fit across all of them. sites_per_section is therefore
#   tracked as a first-class quantity here, not as a by-product.
#
# ------------------------------------------------------------------------------
# WHICH SIGNAL TO TRUST, AND IN WHAT ORDER
#
#   1. SURVEY SITES (location_id) -- the most reliable. A site is a physical,
#      named place. The same sites in two years means index counts were looking
#      at the same places.
#
#   2. SECTIONS -- reliable as recorded, but they change meaning between years:
#      a section number is a label on a block, and blocks get redrawn. Track
#      what is IN a section, never the number itself.
#
#   3. RIVER MILES -- reported here, but treat as unvalidated. They need
#      checking in their own right, so nothing in this script's conclusions
#      rests on them. bss_b_lut_rm_consistency.csv is the check: the same site
#      should carry the same river mile every year, and where it does not, the
#      river miles are not yet trustworthy for that fishery.
#
#   Water body is part of the grain throughout, because a section number only
#   identifies a block within one, and every tributary restarts at river mile 0.
#
# ------------------------------------------------------------------------------
# DEFINED SCOPE vs FITTED SCOPE
#
#   The continuity measures here run on the sections 01_fit_bss_bias.R actually
#   fits (scope_rules.R), not on everything the lookup defines. Otherwise they
#   answer the wrong question: Skagit spring Chinook upper reads as gaining
#   five sites in 2024, and Skagit fall salmon a seventh section in 2025, when
#   both gained the CASCADE -- a separate water body the fits exclude precisely
#   so the series stays comparable.
#
#   bss_b_lut_scope_effect.csv keeps both views side by side. "What did the
#   program change" and "what changed in the data behind b" are different
#   questions, and the gap between them is itself worth seeing.
#
# Usage:
#   Rscript analysis/bss_bias/02a_location_lut_changes.R
#
# Inputs:
#   analysis/bss_bias/lookup/fishery_location_lut.csv  (from 00c; tracked)
#
# Outputs (analysis/bss_bias/outputs/):
#   bss_b_lut_site_stability.csv    -- PRIMARY: per year, sites retained/added/dropped vs prior and first year
#   bss_b_lut_index_density.csv     -- index sites per section-block, the coverage ratio `b` absorbs
#   bss_b_lut_index_density_summary.csv -- per fishery-year: blocks, sites per block, spread
#   bss_b_lut_site_moves.csv        -- sites kept but re-blocked, so paired with a different census count
#   bss_b_lut_year_over_year.csv    -- NAMED diff between consecutive years: sites added, dropped, re-blocked
#   bss_b_lut_core_sites.csv        -- sites present in every year, and whether their block held
#   bss_b_lut_block_year_over_year.csv -- census blocks added/dropped between consecutive years
#   bss_b_lut_scope_effect.csv      -- defined vs fitted scope per fishery-year, and what was excluded
#   bss_b_lut_site_series_summary.csv -- core vs union sites, permanent departures vs skipped years
#   bss_b_lut_location_changes.csv  -- per site: years present, persistent/added/dropped/intermittent
#   bss_b_lut_section_year.csv      -- fishery_type x year x water body x section
#   bss_b_lut_water_body_year.csv   -- fishery_type x year x water body
#   bss_b_lut_changes.csv           -- per fishery_type: which attributes varied, with per-year values
#   bss_b_lut_stability.csv         -- per fishery_type: which attributes held constant
#   bss_b_lut_rm_consistency.csv    -- sites whose river mile moved between years (RM validation)
#   bss_b_lut_anomalies.csv         -- case-variant labels, split sections, zero-length spans
#   bss_b_lut_changes.html          -- gt rendering of the change table
#   figures/fig9_lut_site_persistence.{png,pdf}
#   figures/fig10_lut_site_turnover.{png,pdf}
#   figures/fig11_lut_section_year.{png,pdf}
#   figures/fig12_lut_index_density.{png,pdf}
# ==============================================================================

library(tidyverse)
library(gt)
library(cli)
library(here)

source(here::here("analysis", "bss_bias", "common.R"))
# The same scope rules 01_fit_bss_bias.R applies when building model inputs, so
# the continuity reported here can be read against the water b is actually fit
# to -- not against water the fits exclude.
source(here::here("analysis", "bss_bias", "scope_rules.R"))

OUT_DIR <- here::here("analysis", "bss_bias", "outputs")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

LUT_PATH <- here::here("analysis", "bss_bias", "lookup", "fishery_location_lut.csv")
if (!file.exists(LUT_PATH)) {
  cli::cli_abort(c(
    "{.file {LUT_PATH}} not found.",
    "i" = "Run {.file analysis/bss_bias/00c_probe_location_lut.R} (VPN required) and commit the result."
  ))
}

lut_raw <- read_csv(LUT_PATH, show_col_types = FALSE)
cli::cli_alert_info("Read {nrow(lut_raw)} row(s) covering {n_distinct(lut_raw$fishery_name)} fishery-year(s).")

# ------------------------------------------------------------------------------
# Normalisation
# ------------------------------------------------------------------------------
# location_type and survey_type carry case-variant duplicates in the source
# ("Site"/"site", "Census"/"census"). Left alone they split a group in two, and
# the pipeline's exact-match tests -- prep_dwg_census_expan()'s
# `location_type == "Site"`, for one -- silently drop the odd-cased rows.

title_case <- function(x) str_replace(str_to_lower(x), "^(.)", toupper)

case_anomalies <- lut_raw |>
  filter(location_type != title_case(location_type) | survey_type != title_case(survey_type)) |>
  transmute(anomaly = "case-variant type label", fishery_name,
            section_num = as.integer(section_num), water_body_code, location_code,
            detail = paste0("location_type=", location_type, ", survey_type=", survey_type))

lut <- lut_raw |>
  mutate(
    location_type = title_case(location_type),
    survey_type   = title_case(survey_type),
    fishery_type  = fishery_type_from_name(fishery_name),
    # fishery_start_year, not a regex over the name: the database's own answer,
    # and it resolves the multi-season names ("2022-23" -> 2022) the same way
    # the b-series parity anchor does.
    year          = as.integer(fishery_start_year),
    section_num   = as.integer(section_num),
    across(c(upper_rm, lower_rm), as.numeric)
  )

basin_of <- function(x) {
  case_when(
    str_detect(x, regex("skagit",        ignore_case = TRUE)) ~ "Skagit",
    str_detect(x, regex("snohomish",     ignore_case = TRUE)) ~ "Snohomish",
    str_detect(x, regex("stillaguamish", ignore_case = TRUE)) ~ "Stillaguamish",
    TRUE ~ "Other"
  )
}

# Water bodies take a fixed colour each, assigned once and used everywhere --
# colour follows the entity, so the Skykomish is the same hue in every figure.
# The eight CAT hues are used in their documented order; that palette passes the
# categorical checks at all eight slots.
WATER_BODIES <- sort(unique(lut$water_body_code))
WB_COLORS <- setNames(unname(CAT)[seq_along(WATER_BODIES)], WATER_BODIES)

finite_or_na <- function(x) if_else(is.finite(x), x, NA_real_)

# location_id is the key, not location_code: 141 ids carry only 138 distinct
# codes, so codes are not unique. Codes are for display.
#
# INDEX SITES AND CENSUS BLOCKS ARE DIFFERENT THINGS and must not be pooled.
# The lookup holds both: location_type "Site" rows are the discrete places
# index counts happen, "Section" rows define the block a census count totals.
# Counting them together makes a block being split look like sites arriving --
# Skagit spring Chinook lower keeps the same eight index sites throughout while
# its 2024 block split adds a census-block row, which pooled would read as
# turnover in the water surveyed. It is not.
#
# So `sites` below is INDEX SITES ONLY, and every continuity measure built on
# it is about the places index counts actually happen. Census-block changes are
# reported separately, further down.
all_locations <- lut |>
  distinct(basin = basin_of(fishery_type), fishery_type, fishery_name, year, location_type,
           location_id, location_code, section_num, water_body_code, survey_type)

# IN SCOPE = the sections 01_fit_bss_bias.R actually fits, per scope_rules.R.
# Without this the continuity measures describe water `b` never sees: Skagit
# spring Chinook upper reads as gaining five sites in 2024 when what it gained
# was the Cascade, which the fits exclude, and Skagit fall salmon reads as
# gaining a seventh section in 2025 for the same reason.
#
# Reported BOTH ways throughout. As-defined answers "what did the program
# change"; as-fitted answers "what changed in the data behind b". They are
# different questions and the gap between them is itself informative.
# Resolved once per fishery-year, not per row: fishery_section_limit() reads the
# lookup and can abort, and there are 28 names against 600-odd rows.
scope_by_fishery <- tibble(fishery_name = sort(unique(all_locations$fishery_name))) |>
  mutate(section_limit = map(fishery_name, fishery_section_limit))

all_locations <- all_locations |>
  left_join(scope_by_fishery, by = "fishery_name") |>
  mutate(in_scope = map2_lgl(section_limit, section_num,
                             ~ is.null(.x) || as.double(.y) %in% .x)) |>
  select(-section_limit)

n_out <- sum(!all_locations$in_scope)
if (n_out > 0) {
  cli::cli_alert_info("{n_out} location row(s) fall outside the fitted scope:")
  all_locations |>
    filter(!in_scope) |>
    count(fishery_type, year, water_body_code, name = "n_rows") |>
    print(n = 30)
}

# Everything downstream runs on the FITTED scope, because the question the
# whole script serves is what a `b` comparison across years is comparing.
# sites_defined keeps the unrestricted view for the scope-effect table below.
sites_defined  <- all_locations |> filter(location_type == "Site")
sites          <- sites_defined |> filter(in_scope)
census_blocks  <- all_locations |> filter(location_type == "Section", in_scope)

cli::cli_alert_info(
  "{nrow(sites)} index-site rows and {nrow(census_blocks)} census-block rows \\
   ({sum(sites$in_scope)} sites in fitted scope)."
)

# ------------------------------------------------------------------------------
# What the scope rules removed
# ------------------------------------------------------------------------------
# Side by side, per fishery-year: what the program defined, and what the fits
# actually use. A row where the two differ is a year whose raw construction
# change is partly or wholly an artifact of water the model never sees.

scope_effect <- sites_defined |>
  group_by(basin, fishery_type, year) |>
  summarise(
    n_sites_defined = n_distinct(location_id),
    n_sites_fitted  = n_distinct(location_id[in_scope]),
    water_defined   = paste(sort(unique(water_body_code)), collapse = ", "),
    water_excluded  = paste(sort(unique(water_body_code[!in_scope])), collapse = ", "),
    .groups = "drop"
  ) |>
  mutate(n_sites_excluded = n_sites_defined - n_sites_fitted) |>
  arrange(basin, fishery_type, year)

write_csv(scope_effect, file.path(OUT_DIR, "bss_b_lut_scope_effect.csv"))
cli::cli_alert_success("Wrote bss_b_lut_scope_effect.csv.")

cli::cli_h2("Defined vs fitted scope")
scope_effect |>
  filter(n_sites_excluded > 0) |>
  select(fishery_type, year, n_sites_defined, n_sites_fitted, n_sites_excluded, water_excluded) |>
  print(n = 40)

# ==============================================================================
# PRIMARY: is it the same set of sites?
# ==============================================================================

jaccard <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

site_sets <- sites |>
  group_by(basin, fishery_type, year) |>
  summarise(ids = list(sort(unique(location_id))), .groups = "drop") |>
  arrange(basin, fishery_type, year)

site_stability <- site_sets |>
  group_by(fishery_type) |>
  mutate(
    # Explicit prior-row shift rather than lag(): lag() on a list column fills
    # the first element with NA, and intersect(x, NA) is character(0) -- a
    # first year would report 0 retained instead of "no prior year".
    prev_ids  = c(list(NULL), ids[-n()]),
    first_ids = rep(ids[1], n()),
    n_sites   = map_int(ids, length),
    # vs the previous year: the year-on-year turnover
    n_retained_prev = map2_int(ids, prev_ids, ~ if (is.null(.y)) NA_integer_ else length(intersect(.x, .y))),
    n_added_prev    = map2_int(ids, prev_ids, ~ if (is.null(.y)) NA_integer_ else length(setdiff(.x, .y))),
    n_dropped_prev  = map2_int(ids, prev_ids, ~ if (is.null(.y)) NA_integer_ else length(setdiff(.y, .x))),
    jaccard_prev    = map2_dbl(ids, prev_ids, ~ if (is.null(.y)) NA_real_ else jaccard(.x, .y)),
    # vs the fishery's FIRST year: cumulative drift, which year-on-year hides
    jaccard_first   = map2_dbl(ids, first_ids, jaccard),
    n_shared_first  = map2_int(ids, first_ids, ~ length(intersect(.x, .y)))
  ) |>
  ungroup() |>
  mutate(
    # The judgement this table exists to support. "Same sites" is the aggregate
    # footprint test; it deliberately says nothing about how they were grouped.
    site_verdict = case_when(
      is.na(jaccard_prev)   ~ "first year",
      jaccard_prev == 1     ~ "same sites as prior year",
      jaccard_prev >= 0.9   ~ "near-identical (>=90% shared)",
      jaccard_prev >= 0.7   ~ "mostly shared (70-90%)",
      TRUE                  ~ "substantially different (<70%)"
    )
  ) |>
  select(basin, fishery_type, year, n_sites, n_retained_prev, n_added_prev, n_dropped_prev,
         jaccard_prev, n_shared_first, jaccard_first, site_verdict)

write_csv(site_stability, file.path(OUT_DIR, "bss_b_lut_site_stability.csv"))
cli::cli_alert_success("Wrote bss_b_lut_site_stability.csv ({nrow(site_stability)} fishery-years).")

cli::cli_h2("Site-set stability (the aggregate-footprint test)")
site_stability |>
  select(fishery_type, year, n_sites, n_added_prev, n_dropped_prev, jaccard_prev, site_verdict) |>
  print(n = 100)

# ------------------------------------------------------------------------------
# Sites that stayed but were re-grouped
# ------------------------------------------------------------------------------
# The Hwy 9 signature: the site is still surveyed, the block it sits in changed.
# Separating this from "the site went away" matters, but NEITHER is benign for
# `b`: a re-blocked site keeps feeding an index count, and now that count is
# paired with a different census total. Whatever site-to-block coverage ratio
# `b` had absorbed for the old block no longer describes either new one.
#
# It is also the list to take to the closure record. Sections are often drawn
# on a closure line, and a closure change would move fishing dynamics on top of
# the coverage change.

site_moves <- sites |>
  group_by(basin, fishery_type, location_id, location_code) |>
  filter(n_distinct(paste(water_body_code, section_num)) > 1) |>
  summarise(
    n_years     = n_distinct(year),
    assignments = paste0(year, ": ", water_body_code, " s", section_num, collapse = "; "),
    n_sections  = n_distinct(section_num),
    changed_water_body = n_distinct(water_body_code) > 1,
    .groups = "drop"
  ) |>
  arrange(basin, fishery_type, location_code)

write_csv(site_moves, file.path(OUT_DIR, "bss_b_lut_site_moves.csv"))
cli::cli_alert_success("Wrote bss_b_lut_site_moves.csv ({nrow(site_moves)} sites re-grouped between years).")
if (nrow(site_moves) > 0) {
  cli::cli_h3("Sites kept but re-assigned -- check the closures for these years")
  site_moves |> count(fishery_type, name = "n_sites_regrouped") |> print()
}

# ------------------------------------------------------------------------------
# Index-site density per block -- the ratio `b` absorbs
# ------------------------------------------------------------------------------
# An index count sees the sites; the census count it is paired with totals the
# whole block. So how many index sites sit in a block is a direct proxy for the
# coverage fraction `b` has to carry -- and `b` is one scalar fit across every
# block in the fishery. A year whose blocks hold a consistent number of sites
# is one where a single `b` is a reasonable description; a year with 1 site in
# one block and 9 in another is asking `b` to average over a wide spread.
#
# LUT convention, confirmed against the capture: location_type "Site" rows are
# the index sites, "Section" rows are the census blocks.

index_density <- lut |>
  filter(location_type == "Site") |>
  group_by(basin = basin_of(fishery_type), fishery_type, year, water_body_code, section_num) |>
  summarise(n_index_sites = n_distinct(location_id), .groups = "drop") |>
  left_join(
    lut |>
      filter(location_type == "Section") |>
      distinct(fishery_type, year, water_body_code, section_num) |>
      mutate(has_census_block = TRUE),
    by = c("fishery_type", "year", "water_body_code", "section_num")
  ) |>
  mutate(has_census_block = coalesce(has_census_block, FALSE)) |>
  arrange(basin, fishery_type, year, water_body_code, section_num)

write_csv(index_density, file.path(OUT_DIR, "bss_b_lut_index_density.csv"))

# Per fishery-year: how evenly the sites are spread across blocks. The SPREAD
# is the point, not the mean -- `b` is a single scalar, so an uneven year is
# one where no single value fits every block equally well.
density_summary <- index_density |>
  group_by(basin, fishery_type, year) |>
  summarise(
    n_blocks              = n(),
    n_index_sites         = sum(n_index_sites),
    sites_per_block_mean  = round(mean(n_index_sites), 2),
    sites_per_block_min   = min(n_index_sites),
    sites_per_block_max   = max(n_index_sites),
    blocks_without_census = sum(!has_census_block),
    .groups = "drop"
  ) |>
  group_by(fishery_type) |>
  mutate(mean_change_vs_prior = round(sites_per_block_mean - lag(sites_per_block_mean), 2)) |>
  ungroup()

write_csv(density_summary, file.path(OUT_DIR, "bss_b_lut_index_density_summary.csv"))
cli::cli_alert_success("Wrote bss_b_lut_index_density.csv and bss_b_lut_index_density_summary.csv.")

cli::cli_h2("Index sites per census block")
density_summary |>
  select(fishery_type, year, n_blocks, n_index_sites, sites_per_block_mean,
         sites_per_block_min, sites_per_block_max, mean_change_vs_prior) |>
  print(n = 100)

# ------------------------------------------------------------------------------
# Per-site presence
# ------------------------------------------------------------------------------

fishery_year_span <- sites |>
  distinct(fishery_type, year) |>
  group_by(fishery_type) |>
  summarise(all_years = list(sort(unique(year))), .groups = "drop")

location_changes <- sites |>
  group_by(basin, fishery_type, location_id, location_code) |>
  summarise(
    years_present = list(sort(unique(year))),
    sections      = paste(sort(unique(section_num)), collapse = "|"),
    water_bodies  = paste(sort(unique(water_body_code)), collapse = "|"),
    survey_types  = paste(sort(unique(survey_type)), collapse = "|"),
    .groups = "drop"
  ) |>
  left_join(fishery_year_span, by = "fishery_type") |>
  mutate(
    n_years_present = map_int(years_present, length),
    n_years_fishery = map_int(all_years, length),
    first_year      = map_int(years_present, ~ .x[1]),
    last_year       = map_int(years_present, ~ tail(.x, 1)),
    years           = map_chr(years_present, ~ paste(.x, collapse = ",")),
    status = case_when(
      n_years_present == n_years_fishery                                       ~ "persistent",
      map2_lgl(years_present, all_years, ~ identical(.x, tail(.y, length(.x)))) ~ "added",
      map2_lgl(years_present, all_years, ~ identical(.x, head(.y, length(.x)))) ~ "dropped",
      TRUE                                                                     ~ "intermittent"
    )
  ) |>
  select(basin, fishery_type, location_code, sections, water_bodies, survey_types,
         status, years, n_years_present, n_years_fishery, first_year, last_year) |>
  arrange(basin, fishery_type, status, location_code)

write_csv(location_changes, file.path(OUT_DIR, "bss_b_lut_location_changes.csv"))
cli::cli_alert_success("Wrote bss_b_lut_location_changes.csv ({nrow(location_changes)} sites).")

cli::cli_h3("Site turnover by fishery")
location_changes |>
  count(fishery_type, status) |>
  pivot_wider(names_from = status, values_from = n, values_fill = 0) |>
  print()

# ------------------------------------------------------------------------------
# Year over year: what was built the same, and what was built differently
# ------------------------------------------------------------------------------
# The named diff between consecutive years of one fishery. The set-overlap
# numbers above say HOW MUCH changed; this says WHAT -- which sites arrived,
# which left, and which stayed put but were re-blocked. Those are three
# different things to follow up and only the first two change the footprint.

site_year_section <- sites |>
  group_by(basin, fishery_type, year, location_id, location_code) |>
  summarise(assigned = paste0(water_body_code, " s", sort(unique(section_num)), collapse = " + "),
            .groups = "drop")

year_pairs <- site_year_section |>
  distinct(basin, fishery_type, year) |>
  arrange(basin, fishery_type, year) |>
  group_by(fishery_type) |>
  mutate(prev_year = lag(year)) |>
  ungroup() |>
  filter(!is.na(prev_year))

# A consecutive-year diff alone OVER-COUNTS a temporary contraction. Skagit
# fall salmon dropped 8 upper-river sites in 2024 and restored 7 of them in
# 2025: read pair by pair that is "8 dropped" then "7 added", which sounds like
# a program in flux, when it is one year's contraction reversed. Only ONE site
# actually left for good.
#
# So each arrival and departure is classified against the fishery's WHOLE
# history, not just the adjacent year:
#   new       -- never present before
#   returning -- present in an earlier year, absent, now back
#   permanent -- gone and never returns
#   temporary -- absent this year, present again later
describe_pair <- function(ftype, y_prev, y_now) {
  hist <- site_year_section |> filter(fishery_type == ftype)
  prev <- hist |> filter(year == y_prev)
  now  <- hist |> filter(year == y_now)

  seen_before <- hist |> filter(year < y_prev)  |> pull(location_id) |> unique()
  seen_after  <- hist |> filter(year > y_now)   |> pull(location_id) |> unique()

  added   <- now  |> filter(!location_id %in% prev$location_id) |>
    mutate(kind = if_else(location_id %in% seen_before, "returning", "new"))
  dropped <- prev |> filter(!location_id %in% now$location_id) |>
    mutate(kind = if_else(location_id %in% seen_after, "temporary", "permanent"))
  # A site in both years whose block assignment differs: kept, re-blocked.
  moved <- now |>
    inner_join(prev, by = c("location_id", "location_code"), suffix = c("_now", "_prev")) |>
    filter(assigned_now != assigned_prev)

  # From `sites` rather than section_year, which is built further down -- and
  # `sites` is the same source, so the block set is identical either way.
  sec <- function(y) {
    sites |>
      filter(fishery_type == ftype, year == y) |>
      distinct(water_body_code, section_num) |>
      mutate(lbl = paste0(water_body_code, " s", section_num)) |>
      pull(lbl) |> sort()
  }
  sec_prev <- sec(y_prev); sec_now <- sec(y_now)

  tibble(
    fishery_type   = ftype,
    year_prev      = y_prev,
    year           = y_now,
    n_retained     = nrow(now) - nrow(added),
    n_added        = nrow(added),
    n_added_new       = sum(added$kind == "new"),
    n_added_returning = sum(added$kind == "returning"),
    n_dropped      = nrow(dropped),
    n_dropped_permanent = sum(dropped$kind == "permanent"),
    n_dropped_temporary = sum(dropped$kind == "temporary"),
    n_moved_block  = nrow(moved),
    sites_added    = paste(sort(paste0(added$location_code, " [", added$kind, "]")), collapse = "; "),
    sites_dropped  = paste(sort(paste0(dropped$location_code, " [", dropped$kind, "]")), collapse = "; "),
    sites_moved    = paste(sort(paste0(moved$location_code, " (", moved$assigned_prev,
                                       " -> ", moved$assigned_now, ")")), collapse = "; "),
    blocks_added   = paste(setdiff(sec_now, sec_prev), collapse = "; "),
    blocks_dropped = paste(setdiff(sec_prev, sec_now), collapse = "; "),
    # The verdict reads NET-OF-REVERSAL change: a site that leaves and comes
    # back has not changed the fishery, it has interrupted it for a year.
    verdict = dplyr::case_when(
      nrow(added) == 0 & nrow(dropped) == 0 & nrow(moved) == 0 ~ "identical construction",
      nrow(added) == 0 & nrow(dropped) == 0                    ~ "same sites, re-blocked",
      sum(added$kind == "new") == 0 &
        sum(dropped$kind == "permanent") == 0                  ~ "temporary contraction or restoration",
      nrow(moved) == 0                                         ~ "sites changed, blocking held",
      TRUE                                                     ~ "sites changed and re-blocked"
    )
  )
}

year_over_year <- pmap(
  list(year_pairs$fishery_type, year_pairs$prev_year, year_pairs$year),
  describe_pair
) |> bind_rows()

write_csv(year_over_year, file.path(OUT_DIR, "bss_b_lut_year_over_year.csv"))
cli::cli_alert_success("Wrote bss_b_lut_year_over_year.csv ({nrow(year_over_year)} year pairs).")

cli::cli_h2("Year-over-year construction")
year_over_year |>
  select(fishery_type, year_prev, year, n_retained, n_added, n_dropped, n_moved_block, verdict) |>
  print(n = 60)

# ------------------------------------------------------------------------------
# Net site stability over the whole series
# ------------------------------------------------------------------------------
# The answer that survives temporary contractions. Core = present in every
# year. Union = every site the fishery ever used. Sites that left for good and
# sites that ever skipped a year are counted separately, because "the program
# stopped surveying here" and "one lean season" are different facts.

# Reuses location_changes' classification rather than re-deriving it, so the
# two cannot disagree about what "dropped" means:
#   persistent   -- every year
#   dropped      -- present from the start, then gone and never back
#   added        -- arrived partway and stayed
#   intermittent -- skipped at least one year and returned
site_series_summary <- location_changes |>
  group_by(basin, fishery_type) |>
  summarise(
    n_years         = first(n_years_fishery),
    core_all_years  = sum(status == "persistent"),
    union_ever_used = n(),
    left_for_good   = sum(status == "dropped"),
    arrived_later   = sum(status == "added"),
    skipped_a_year  = sum(status == "intermittent"),
    core_pct_of_union = round(100 * sum(status == "persistent") / n(), 1),
    .groups = "drop"
  )

write_csv(site_series_summary, file.path(OUT_DIR, "bss_b_lut_site_series_summary.csv"))
cli::cli_alert_success("Wrote bss_b_lut_site_series_summary.csv.")

cli::cli_h2("Net site stability across the whole series")
print(site_series_summary |> select(-basin), n = 20)

# ------------------------------------------------------------------------------
# Census blocks, year over year
# ------------------------------------------------------------------------------
# The other half of the construction. A block appearing or disappearing changes
# what a census count totals, and therefore the coverage ratio `b` absorbs --
# independently of whether any index site moved.

block_year_over_year <- census_blocks |>
  group_by(basin, fishery_type, year) |>
  summarise(blocks = list(sort(unique(paste0(water_body_code, " s", section_num)))),
            n_blocks = n_distinct(paste(water_body_code, section_num)), .groups = "drop") |>
  arrange(basin, fishery_type, year) |>
  group_by(fishery_type) |>
  mutate(prev_blocks = c(list(NULL), blocks[-n()]), prev_year = lag(year)) |>
  ungroup() |>
  filter(!is.na(prev_year)) |>
  mutate(
    n_added   = map2_int(blocks, prev_blocks, ~ length(setdiff(.x, .y))),
    n_dropped = map2_int(blocks, prev_blocks, ~ length(setdiff(.y, .x))),
    added     = map2_chr(blocks, prev_blocks, ~ paste(setdiff(.x, .y), collapse = "; ")),
    dropped   = map2_chr(blocks, prev_blocks, ~ paste(setdiff(.y, .x), collapse = "; "))
  ) |>
  select(basin, fishery_type, prev_year, year, n_blocks, n_added, n_dropped, added, dropped)

write_csv(block_year_over_year, file.path(OUT_DIR, "bss_b_lut_block_year_over_year.csv"))
cli::cli_alert_success("Wrote bss_b_lut_block_year_over_year.csv.")

cli::cli_h2("Census blocks, year over year")
block_year_over_year |> select(-basin, -added, -dropped) |> print(n = 40)

# ------------------------------------------------------------------------------
# The constant core
# ------------------------------------------------------------------------------
# Sites present in EVERY year of a fishery, and whether their block assignment
# also held. This is the "what stayed the same" answer at its most concrete: a
# site in the core with one assignment throughout is a piece of the fishery that
# is genuinely comparable across the whole series.

core_sites <- site_year_section |>
  group_by(basin, fishery_type, location_id, location_code) |>
  summarise(n_years = n_distinct(year),
            assignments = paste0(year, ": ", assigned, collapse = "; "),
            n_assignments = n_distinct(assigned),
            .groups = "drop") |>
  left_join(fishery_year_span, by = "fishery_type") |>
  mutate(n_years_fishery = map_int(all_years, length),
         in_core = n_years == n_years_fishery,
         core_status = case_when(
           !in_core             ~ "not in every year",
           n_assignments == 1   ~ "core, same block throughout",
           TRUE                 ~ "core, re-blocked"
         )) |>
  select(basin, fishery_type, location_code, core_status, n_years, n_years_fishery,
         n_assignments, assignments) |>
  arrange(basin, fishery_type, core_status, location_code)

write_csv(core_sites, file.path(OUT_DIR, "bss_b_lut_core_sites.csv"))
cli::cli_alert_success("Wrote bss_b_lut_core_sites.csv.")

cli::cli_h2("Constant core: sites in every year, and whether the block held")
core_sites |>
  count(fishery_type, core_status) |>
  pivot_wider(names_from = core_status, values_from = n, values_fill = 0) |>
  print()

# ==============================================================================
# SUPPORTING: how the sites were grouped
# ==============================================================================
# Water body is part of the grain, not a label hung off it: Stillaguamish 2023
# carries section 5 on both the North and South Forks, so (year, section_num)
# alone does not identify a row.

section_year <- lut |>
  group_by(basin = basin_of(fishery_type), fishery_type, year, water_body_code, section_num) |>
  summarise(
    n_locations = n_distinct(location_id),
    n_index     = n_distinct(location_id[survey_type == "Index"]),
    n_census    = n_distinct(location_id[survey_type == "Census"]),
    rm_lower    = suppressWarnings(min(lower_rm, na.rm = TRUE)),
    rm_upper    = suppressWarnings(max(upper_rm, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  # min/max over an all-NA river-mile group returns +/-Inf, not NA.
  mutate(across(c(rm_lower, rm_upper), finite_or_na), rm_span = rm_upper - rm_lower) |>
  arrange(basin, fishery_type, year, water_body_code, section_num)

write_csv(section_year, file.path(OUT_DIR, "bss_b_lut_section_year.csv"))

water_body_year <- lut |>
  group_by(basin = basin_of(fishery_type), fishery_type, year, water_body_code, water_body_desc) |>
  summarise(
    n_sections  = n_distinct(section_num),
    n_locations = n_distinct(location_id),
    n_index     = n_distinct(location_id[survey_type == "Index"]),
    n_census    = n_distinct(location_id[survey_type == "Census"]),
    rm_lower    = suppressWarnings(min(lower_rm, na.rm = TRUE)),
    rm_upper    = suppressWarnings(max(upper_rm, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  mutate(across(c(rm_lower, rm_upper), finite_or_na)) |>
  arrange(basin, fishery_type, year, water_body_code)

write_csv(water_body_year, file.path(OUT_DIR, "bss_b_lut_water_body_year.csv"))
cli::cli_alert_success("Wrote bss_b_lut_section_year.csv and bss_b_lut_water_body_year.csv.")

# ------------------------------------------------------------------------------
# Attribute-level diff
# ------------------------------------------------------------------------------

fishery_year_attrs <- lut |>
  group_by(basin = basin_of(fishery_type), fishery_type, year) |>
  summarise(
    n_sites          = n_distinct(location_id),
    n_index_sites    = n_distinct(location_id[survey_type == "Index"]),
    n_census_units   = n_distinct(location_id[survey_type == "Census"]),
    n_sections       = n_distinct(section_num),
    sections         = paste(sort(unique(section_num)), collapse = ","),
    n_water_bodies   = n_distinct(water_body_code),
    water_bodies     = paste(sort(unique(water_body_code)), collapse = ","),
    n_surveyors      = n_distinct(na.omit(surveyor_num)),
    fishery_start    = min(fishery_start_date),
    fishery_end      = max(fishery_end_date),
    .groups = "drop"
  )

attr_long <- fishery_year_attrs |>
  mutate(across(everything(), as.character)) |>
  pivot_longer(-c(basin, fishery_type, year), names_to = "attribute", values_to = "value")

attr_summary <- attr_long |>
  # `values` pastes year=value in row order, so the per-year list has to be
  # chronological rather than whatever order the pivot left.
  arrange(basin, fishery_type, attribute, year) |>
  group_by(basin, fishery_type, attribute) |>
  summarise(
    n_years  = n_distinct(year),
    n_values = n_distinct(value),
    values   = paste0(year, "=", value, collapse = "; "),
    .groups = "drop"
  ) |>
  mutate(changed = n_values > 1)

lut_changes <- attr_summary |> filter(changed)  |> arrange(basin, fishery_type, desc(n_values), attribute)
lut_stable  <- attr_summary |> filter(!changed) |> arrange(basin, fishery_type, attribute)

write_csv(lut_changes, file.path(OUT_DIR, "bss_b_lut_changes.csv"))
write_csv(lut_stable,  file.path(OUT_DIR, "bss_b_lut_stability.csv"))
cli::cli_alert_success(
  "Wrote bss_b_lut_changes.csv ({nrow(lut_changes)} varying) and \\
   bss_b_lut_stability.csv ({nrow(lut_stable)} constant)."
)

# ==============================================================================
# River-mile validation
# ==============================================================================
# Not a result -- a check on whether the river miles can be used at all. The
# same site should carry the same river mile every year. Where it does not,
# either the site moved or the river miles were re-measured, and until that is
# resolved no distance-based comparison for that fishery means anything.

rm_consistency <- lut |>
  filter(!is.na(lower_rm) | !is.na(upper_rm)) |>
  group_by(basin = basin_of(fishery_type), fishery_type, location_id, location_code) |>
  summarise(
    n_years      = n_distinct(year),
    n_rm_values  = n_distinct(paste(lower_rm, upper_rm)),
    rm_by_year   = paste0(year, ": ", lower_rm, "-", upper_rm, collapse = "; "),
    .groups = "drop"
  ) |>
  filter(n_years > 1, n_rm_values > 1) |>
  arrange(basin, fishery_type, location_code)

write_csv(rm_consistency, file.path(OUT_DIR, "bss_b_lut_rm_consistency.csv"))
if (nrow(rm_consistency) > 0) {
  cli::cli_alert_warning(
    "{nrow(rm_consistency)} site(s) carry more than one river mile across years -- \\
     river miles are NOT yet reliable for {.val {sort(unique(rm_consistency$fishery_type))}}."
  )
} else {
  cli::cli_alert_success("Every site carries a consistent river mile across years.")
}

# ------------------------------------------------------------------------------
# Anomalies
# ------------------------------------------------------------------------------

split_sections <- section_year |>
  group_by(fishery_type, year, section_num) |>
  filter(n_distinct(water_body_code) > 1) |>
  summarise(wbs = paste(sort(unique(water_body_code)), collapse = " + "), .groups = "drop") |>
  transmute(anomaly = "section spans >1 water body",
            fishery_name = paste(fishery_type, year), section_num,
            water_body_code = wbs, location_code = NA_character_,
            detail = "(year, section_num) is not a unique key for this fishery-year")

zero_spans <- section_year |>
  filter(!is.na(rm_span), rm_span == 0) |>
  transmute(anomaly = "zero-length river-mile span",
            fishery_name = paste(fishery_type, year), section_num, water_body_code,
            location_code = NA_character_,
            detail = paste0("rm ", rm_lower, "-", rm_upper))

# The capture is otherwise consistent: "Site" rows are Index, "Section" rows
# are Census. A row breaking that pairing means a block is being treated as an
# index location or vice versa, which crosses the two count types `b` relates.
type_mismatches <- lut |>
  filter((location_type == "Site"    & survey_type != "Index") |
         (location_type == "Section" & survey_type != "Census")) |>
  transmute(anomaly = "location_type / survey_type mismatch", fishery_name,
            section_num, water_body_code, location_code,
            detail = paste0(location_type, " + ", survey_type))

anomalies <- bind_rows(case_anomalies, type_mismatches, split_sections, zero_spans)
if (nrow(anomalies) > 0) {
  write_csv(anomalies, file.path(OUT_DIR, "bss_b_lut_anomalies.csv"))
  cli::cli_alert_warning("{nrow(anomalies)} anomal{?y/ies} in the lookup:")
  print(anomalies, width = Inf)
} else {
  cli::cli_alert_success("No lookup anomalies found.")
}

# ==============================================================================
# Figures
# ==============================================================================

# ------------------------------------------------------------------------------
# Fig 9 -- site persistence, coloured by section
# ------------------------------------------------------------------------------
# One chart, two questions. A row that is filled in every year is a site the
# fishery kept; a gap is a site added or dropped. A row that stays filled but
# CHANGES COLOUR is a site kept and re-grouped -- the Hwy 9 case, where the
# footprint held and only the paperwork moved.
#
# Section number takes the sequential ramp, not categorical hues: sections are
# ordinal (numbered downstream to upstream), and there are up to nine of them,
# past where a categorical palette should be asked to go.

site_grid <- sites |>
  group_by(fishery_type, location_id) |>
  mutate(order_key = min(section_num)) |>
  ungroup() |>
  arrange(fishery_type, order_key, location_code) |>
  mutate(site_label = fct_rev(fct_inorder(paste0(location_code, "  [", water_body_code, "]"))))

fig9 <- site_grid |>
  ggplot(aes(x = factor(year), y = site_label, fill = section_num)) +
  geom_tile(colour = SURFACE, linewidth = 0.8) +
  facet_wrap(~fishery_type, scales = "free_y", ncol = 2) +
  scale_fill_gradientn(colours = SEQ_RAMP, name = "Section", breaks = scales::breaks_width(2)) +
  labs(
    title = "Which survey sites each fishery used, by year",
    subtitle = "A gap is a site added or dropped. A row that stays filled but changes colour is a site kept and re-blocked -- its index count is now paired with a different census total.",
    x = NULL, y = NULL
  ) +
  theme_bss() +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 6),
        legend.position = "right")

save_fig(fig9, "fig9_lut_site_persistence", width = 13, height = 14)
cli::cli_alert_success("Wrote fig9_lut_site_persistence.")

# ------------------------------------------------------------------------------
# Fig 10 -- year-on-year site turnover
# ------------------------------------------------------------------------------
# The summary of fig 9: how many sites carried over, how many are new, how many
# went away. Retained/added/dropped is a status-like distinction, so it takes
# the reserved status hues with a legend rather than arbitrary series colours.

turnover <- site_stability |>
  filter(!is.na(n_retained_prev)) |>
  select(fishery_type, year, Retained = n_retained_prev, Added = n_added_prev, Dropped = n_dropped_prev) |>
  pivot_longer(c(Retained, Added, Dropped), names_to = "change", values_to = "n") |>
  mutate(change = factor(change, levels = c("Retained", "Added", "Dropped")))

fig10 <- turnover |>
  ggplot(aes(x = factor(year), y = n, fill = change)) +
  geom_col(width = 0.7, colour = SURFACE, linewidth = 1) +
  facet_wrap(~fishery_type, scales = "free_y", ncol = 3) +
  scale_fill_manual(
    values = c(Retained = BASELINE_COL, Added = STATUS[["good"]], Dropped = STATUS[["serious"]]),
    name = NULL
  ) +
  labs(
    title = "Site turnover against the prior year",
    subtitle = "Retained sites are the fishery's continuity. A year that is mostly grey looked at the same places as the one before it -- read it with the block density below, which can move even when every site is kept.",
    x = NULL, y = "Sites"
  ) +
  theme_bss() +
  theme(legend.position = "top")

save_fig(fig10, "fig10_lut_site_turnover", width = 12, height = 7)
cli::cli_alert_success("Wrote fig10_lut_site_turnover.")

# ------------------------------------------------------------------------------
# Fig 11 -- section x year, by water body
# ------------------------------------------------------------------------------
# The supporting view: how the sites were partitioned. Read this AFTER fig 9 --
# a change here with no change there is a re-grouping, not a footprint change.

fig11 <- section_year |>
  mutate(row_label = fct_rev(fct_inorder(paste0(water_body_code, "  s", section_num)))) |>
  ggplot(aes(x = factor(year), y = row_label, fill = n_locations)) +
  geom_tile(colour = SURFACE, linewidth = 1.2) +
  geom_text(aes(label = n_locations), colour = INK, size = 3) +
  facet_wrap(~fishery_type, scales = "free_y", ncol = 2) +
  scale_fill_gradientn(colours = SEQ_RAMP, name = "Sites") +
  labs(
    title = "How sites were blocked into sections, by year",
    subtitle = "Sections are the estimation unit: a census count totals the whole block and is paired with an index count of its sites. A change here moves the coverage ratio b absorbs, even when every site is kept.",
    x = NULL, y = NULL
  ) +
  theme_bss() +
  theme(panel.grid = element_blank(), legend.position = "right")

save_fig(fig11, "fig11_lut_section_year", width = 12, height = 10)
cli::cli_alert_success("Wrote fig11_lut_section_year.")

# ------------------------------------------------------------------------------
# Fig 12 -- index sites per census block
# ------------------------------------------------------------------------------
# The `b`-relevant view. Each point is one block; the vertical spread within a
# year is how unevenly the index sites are distributed across the blocks a
# single `b` has to cover. A year where the spread widens is a year where one
# scalar fits the fishery less well, regardless of what happened to the sites.

fig12 <- index_density |>
  ggplot(aes(x = factor(year), y = n_index_sites)) +
  geom_line(aes(group = 1), data = ~ .x |> group_by(fishery_type, year) |>
              summarise(n_index_sites = mean(n_index_sites), .groups = "drop"),
            colour = INK_MUTED, linewidth = 0.6) +
  geom_point(aes(colour = water_body_code), size = 2.4,
             position = position_jitter(width = 0.12, height = 0, seed = 1)) +
  facet_wrap(~fishery_type, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = WB_COLORS, name = NULL) +
  labs(
    title = "Index sites per census block",
    subtitle = "One point per block; the grey line is the fishery's mean. A census count totals the whole block, an index count sees only its sites -- so the spread here is the range of coverage a single b has to average over.",
    x = NULL, y = "Index sites in block"
  ) +
  theme_bss() +
  theme(legend.position = "top")

save_fig(fig12, "fig12_lut_index_density", width = 12, height = 7)
cli::cli_alert_success("Wrote fig12_lut_index_density.")

# ------------------------------------------------------------------------------
# gt rendering
# ------------------------------------------------------------------------------

gt_tbl <- site_stability |>
  select(basin, fishery_type, year, n_sites, n_added_prev, n_dropped_prev,
         jaccard_prev, jaccard_first, site_verdict) |>
  mutate(across(c(jaccard_prev, jaccard_first), ~ round(.x, 3))) |>
  gt(groupname_col = "fishery_type") |>
  tab_header(
    title = "Is it the same fishery, year to year?",
    subtitle = "Measured on the set of survey sites, which are physical places -- not on sections, which are administrative and get redrawn."
  ) |>
  cols_label(n_sites = "Sites", n_added_prev = "Added", n_dropped_prev = "Dropped",
             jaccard_prev = "Overlap vs prior", jaccard_first = "Overlap vs first yr",
             site_verdict = "Verdict") |>
  tab_footnote(
    footnote = "Overlap is |shared| / |union| of site sets. 1.00 means exactly the same sites.",
    locations = cells_column_labels(columns = jaccard_prev)
  ) |>
  opt_row_striping()

gtsave(gt_tbl, file.path(OUT_DIR, "bss_b_lut_changes.html"))
cli::cli_alert_success("Wrote bss_b_lut_changes.html")
