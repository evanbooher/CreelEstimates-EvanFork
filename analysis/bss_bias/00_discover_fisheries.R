# ==============================================================================
# 00_discover_fisheries.R
#
# Purpose:
#   Fast, VPN-free, DB-free discovery of every historical Stillaguamish /
#   Snohomish / Skagit creel fishery-name x year, using the PUBLIC "WDFW Creel
#   Fishery Manager" dataset on data.wa.gov (Socrata dataset id vkjc-s5u8):
#     https://data.wa.gov/Natural-Resources-Environment/WDFW-Creel-Fishery-Manager/vkjc-s5u8/about_data
#
#   This determines the SIZE of the whole downstream job before anything slow
#   (VPN, DB, Stan) is touched. Run this first and eyeball the output.
#
# IMPORTANT -- schema is NOT verified. This script was written without the
#   ability to query data.wa.gov directly (network-restricted sandbox), so it
#   is deliberately defensive: it fetches a small sample first, prints the
#   field names, and STOPS if it cannot find a plausible fishery-name column
#   rather than guessing further. Read the printed sample before trusting the
#   rest of the run.
#
# ALSO WORTH 5 MINUTES: browse https://data.wa.gov/browse?q=creel by hand.
#   vkjc-s5u8 is a fishery *registry* (names/metadata). There may be sibling
#   published datasets carrying actual effort counts / interview-level data --
#   if the "public" fetch_data() path in analysis/bss_bias/01_fit_bss_bias.R
#   Phase 0.A works, one of those siblings is probably what it's hitting.
#
# ALSO SEE Phase 0.A in the header of 01_fit_bss_bias.R: before assuming a
#   VPN-only path is required for the actual BSS fits, try
#   creelutils::fetch_data(fishery_name = ..., data_source = "public") (or
#   whatever value creelutils::fetch_data's data_source arg actually accepts --
#   check `?creelutils::fetch_data` / `formals(creelutils::fetch_data)` locally)
#   against ONE target fishery first. If it returns interview/effort-count-level
#   data (not just published estimates), the VPN dependency may disappear
#   entirely for this analysis.
#
# Usage:
#   Rscript analysis/bss_bias/00_discover_fisheries.R
#   No VPN, no DB connection, no rstan required. Needs internet access only.
#
# Outputs (analysis/bss_bias/outputs/):
#   fishery_discovery_raw.csv     -- everything fetched from the registry, unfiltered
#   fishery_discovery_target.csv  -- the 3-basin (+ adjacent) subset, with an
#                                     `include_in_run` column for hand review
# ==============================================================================

library(tidyverse)
library(httr2)
library(cli)
library(here)

OUT_DIR <- here::here("analysis", "bss_bias", "outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SOCRATA_BASE <- "https://data.wa.gov/resource/vkjc-s5u8.json"
SOCRATA_META <- "https://data.wa.gov/api/views/vkjc-s5u8.json"

# Primary basins the meeting is about.
TARGET_BASIN_PATTERN <- regex("skagit|snohomish|stillaguamish", ignore_case = TRUE)

# Adjacent/tributary systems worth surfacing for the group to rule in or out
# (e.g. Cascade folds into Skagit spring Chinook upper in some years per the
# PST crosswalk -- see analysis/bss_bias/README.md).
ADJACENT_BASIN_PATTERN <- regex(
  "sauk|cascade|skykomish|snoqualmie|wallace|pilchuck|canyon creek",
  ignore_case = TRUE
)

# ------------------------------------------------------------------------------
# 1. Probe schema first. Do NOT assume column names.
# ------------------------------------------------------------------------------

cli::cli_h1("Step 1: probing dataset schema")

sample_resp <- tryCatch(
  request(SOCRATA_BASE) |>
    req_url_query(`$limit` = 5) |>
    req_perform(),
  error = function(e) {
    cli::cli_abort(c(
      "Could not reach {.url {SOCRATA_BASE}}.",
      "x" = "{conditionMessage(e)}",
      "i" = "This needs plain internet access (not VPN). If this machine also",
      "i" = "blocks data.wa.gov, open the URL in a browser instead and adapt",
      "i" = "this script's column names by hand from what you see."
    ))
  }
)

sample_df <- resp_body_json(sample_resp, simplifyVector = TRUE)

if (!is.data.frame(sample_df) || nrow(sample_df) == 0) {
  cli::cli_abort(
    "Sample request to {.url {SOCRATA_BASE}} returned no rows. Check the URL \\
     by hand in a browser -- the dataset id (vkjc-s5u8) may have changed."
  )
}

cli::cli_alert_success("Sample fetch OK. Fields present:")
print(names(sample_df))
cli::cli_h3("First few rows (inspect before trusting the filters below)")
print(sample_df)

# Also pull the dataset's field metadata/descriptions -- helps disambiguate a
# column like "fishery" vs "fishery_name" vs "project_name" without guessing.
meta <- tryCatch(
  request(SOCRATA_META) |> req_perform() |> resp_body_json(simplifyVector = TRUE),
  error = function(e) NULL
)
if (!is.null(meta) && !is.null(meta$columns)) {
  cli::cli_h3("Column metadata (name / description) from the dataset's own schema")
  meta$columns |>
    as_tibble() |>
    select(any_of(c("name", "fieldName", "description"))) |>
    print(n = 100)
}

# ------------------------------------------------------------------------------
# 2. Resolve the fishery-name and year/date columns defensively.
# ------------------------------------------------------------------------------

cli::cli_h1("Step 2: resolving key columns")

candidate_name_cols <- c(
  "fishery_name", "fishery", "name", "fishery_title", "project_name", "waterbody"
)
name_col <- intersect(candidate_name_cols, names(sample_df))[1]

if (is.na(name_col)) {
  cli::cli_abort(c(
    "Could not find a plausible fishery-name column automatically.",
    "x" = "Checked for: {.val {candidate_name_cols}}",
    "i" = "Fields actually present: {.val {names(sample_df)}}",
    "i" = "Set `name_col` manually below to the correct field and re-run from",
    "i" = "Step 3, or inspect https://data.wa.gov/d/vkjc-s5u8 in a browser."
  ))
}
cli::cli_alert_success("Using {.val {name_col}} as the fishery-name column.")

candidate_date_cols <- c("start_date", "date_start", "season_start", "year", "season")
date_col <- intersect(candidate_date_cols, names(sample_df))[1]
if (!is.na(date_col)) {
  cli::cli_alert_info("Using {.val {date_col}} as a supplementary date/year field (year is also parsed from the name).")
} else {
  cli::cli_alert_warning("No obvious date/year column found; year will be parsed from {.val {name_col}} only (regex \\d{{4}}).")
}

# ------------------------------------------------------------------------------
# 3. Paginate the FULL registry (safer than a server-side $where on an
#    unverified column name -- filter locally instead).
# ------------------------------------------------------------------------------

cli::cli_h1("Step 3: paginating full registry")

fetch_page <- function(offset, limit = 5000) {
  request(SOCRATA_BASE) |>
    req_url_query(`$limit` = limit, `$offset` = offset, `$order` = ":id") |>
    req_perform() |>
    resp_body_json(simplifyVector = TRUE)
}

all_pages <- list()
offset <- 0L
limit <- 5000L
repeat {
  page <- fetch_page(offset, limit)
  if (!is.data.frame(page) || nrow(page) == 0) break
  all_pages[[length(all_pages) + 1]] <- page
  cli::cli_alert_info("  fetched {nrow(page)} rows at offset {offset}")
  if (nrow(page) < limit) break
  offset <- offset + limit
}

if (length(all_pages) == 0) {
  cli::cli_abort("Pagination returned zero rows total -- something is wrong upstream of the filters.")
}

registry_raw <- bind_rows(all_pages) |> as_tibble()
cli::cli_alert_success("Total rows fetched: {nrow(registry_raw)}")

write_csv(registry_raw, file.path(OUT_DIR, "fishery_discovery_raw.csv"))

# ------------------------------------------------------------------------------
# 4. Filter to target + adjacent basins. Match on BASIN TOKEN ONLY.
#
#    Do NOT filter on "salmon" in the name -- analysis/pst/02_ingest/
#    multi_fishery_creel_summary.R (chore/multi-fishery-trip-summary branch)
#    documents that exact mistake dropping "Skagit spring Chinook 2024 upper"
#    and "Skagit summer sockeye 2023" because they don't contain "salmon".
# ------------------------------------------------------------------------------

cli::cli_h1("Step 4: filtering to target basins")

registry <- registry_raw |>
  mutate(fishery_name_raw = .data[[name_col]])

target <- registry |>
  filter(str_detect(fishery_name_raw, TARGET_BASIN_PATTERN)) |>
  mutate(basin_match = "target")

adjacent <- registry |>
  filter(
    str_detect(fishery_name_raw, ADJACENT_BASIN_PATTERN),
    !str_detect(fishery_name_raw, TARGET_BASIN_PATTERN)
  ) |>
  mutate(basin_match = "adjacent")

combined <- bind_rows(target, adjacent) |>
  mutate(
    basin = case_when(
      str_detect(fishery_name_raw, regex("skagit", ignore_case = TRUE)) ~ "Skagit",
      str_detect(fishery_name_raw, regex("snohomish", ignore_case = TRUE)) ~ "Snohomish",
      str_detect(fishery_name_raw, regex("skykomish", ignore_case = TRUE)) ~ "Skykomish",
      str_detect(fishery_name_raw, regex("stillaguamish", ignore_case = TRUE)) ~ "Stillaguamish",
      TRUE ~ str_extract(fishery_name_raw, ADJACENT_BASIN_PATTERN) |> str_to_title()
    ),
    year_extracted = str_extract(fishery_name_raw, "\\d{4}"),
    # multi-season fisheries (e.g. "Stillaguamish salmon and gamefish 2024-25")
    season_label = str_extract(fishery_name_raw, "\\d{4}(-\\d{2,4})?"),
    spans_calendar_years = str_detect(coalesce(season_label, ""), "-"),
    year_start = suppressWarnings(as.integer(year_extracted)),
    include_in_run = TRUE  # hand-edit this column to prune before Phase 1
  ) |>
  arrange(basin_match, basin, year_start, fishery_name_raw) |>
  select(basin_match, basin, fishery_name_raw, season_label, year_start,
         spans_calendar_years, include_in_run, everything())

write_csv(combined, file.path(OUT_DIR, "fishery_discovery_target.csv"))

cli::cli_h1("Result: basin x year FISHERY-YEAR count")
cli::cli_alert_info(
  "This dataset is one row per (fishery, location/section), not one row per \\
   fishery-year -- distinct()'d on fishery_name_raw below so this count \\
   matches what 01_fit_bss_bias.R will actually iterate over."
)
combined |>
  filter(basin_match == "target") |>
  distinct(basin, year_start, fishery_name_raw) |>
  count(basin, year_start) |>
  print(n = 100)

cli::cli_h2("Distinct target fishery-years (review THIS list by hand)")
combined |>
  filter(basin_match == "target") |>
  distinct(basin, year_start, fishery_name_raw) |>
  arrange(basin, year_start, fishery_name_raw) |>
  print(n = 300)

cli::cli_alert_success(
  "Wrote {nrow(combined)} target/adjacent rows to \\
   analysis/bss_bias/outputs/fishery_discovery_target.csv"
)
cli::cli_alert_warning(
  "REVIEW THIS FILE BY HAND before running 01_fit_bss_bias.R -- edit \\
   `include_in_run` to prune anything that shouldn't be fit (e.g. a fishery \\
   name that matched by coincidence, or a year with known no-data)."
)
