# ==============================================================================
# 00c_probe_location_lut.R
#
# Purpose:
#   Find the fishery/location lookup table in the internal creel database,
#   report its schema, and capture the rows for this analysis's target
#   fishery-years into a COMMITTED lookup file -- so that comparing how a
#   fishery's location definition changed across years is reproducible by
#   anyone with the repo, no VPN and no database credentials.
#
#   Same contract as 00b_capture_fishery_params.R: RUN ONCE, WITH VPN, and
#   commit the resulting lookup. Meeting participants never run this script.
#
# ------------------------------------------------------------------------------
# WHY IT PROBES INSTEAD OF JUST QUERYING
#
#   `fishery_location_lut` is referenced nowhere in this repo -- not in
#   R_functions/, not in template_scripts/, not in any analysis script -- and
#   creelutils exposes no accessor for it (only fishery_lut(), which carries
#   the season dates 00b captures). So the table's real name, schema and
#   fishery key are unknown here.
#
#   Rather than guess a name and fail with a confusing SQL error, this script
#   asks the database what it has: it looks for tables whose COLUMN SIGNATURE
#   fits a fishery-by-location lookup (something fishery-shaped plus something
#   section- or location-shaped), scores them, and stops with the evidence
#   printed if the answer is not obvious. Name matching alone is too brittle --
#   the table may not have "lut" in its name at all.
#
#   If you already know the table, set LUT_TABLE below and the search is
#   skipped entirely.
#
# Usage:
#   Rscript analysis/bss_bias/00c_probe_location_lut.R
#   REQUIRES VPN / internal DB access. Run after 00_discover_fisheries.R,
#   since it captures the fisheries 01_fit_bss_bias.R will actually attempt.
#
# Outputs:
#   analysis/bss_bias/outputs/location_lut_candidates.csv -- every table whose
#     columns fit the pattern, with its score and matched columns. Diagnostic.
#   analysis/bss_bias/outputs/location_lut_schema.csv     -- columns + types of
#     the chosen table. Diagnostic.
#
#   analysis/bss_bias/lookup/fishery_location_lut.csv     -- TRACKED IN GIT.
#     The captured rows. This file is the point of the script.
# ==============================================================================

library(tidyverse)
library(cli)
library(here)
library(creelutils)
library(DBI)

OUT_DIR     <- here::here("analysis", "bss_bias", "outputs")
LOOKUP_DIR  <- here::here("analysis", "bss_bias", "lookup")
LOOKUP_PATH <- file.path(LOOKUP_DIR, "fishery_location_lut.csv")
dir.create(OUT_DIR,    recursive = TRUE, showWarnings = FALSE)
dir.create(LOOKUP_DIR, recursive = TRUE, showWarnings = FALSE)

# Set this to skip the search, e.g. "dbo.fishery_location_lut". Schema-qualify
# it if the database needs that to resolve the name.
LUT_TABLE <- NULL

# Cap on rows pulled when the table cannot be filtered by fishery name, so a
# probe against an unexpectedly large table cannot hang the session.
MAX_UNFILTERED_ROWS <- 50000

# ------------------------------------------------------------------------------
# Target fisheries (same source of truth as 00b)
# ------------------------------------------------------------------------------

discovery_path <- file.path(OUT_DIR, "fishery_discovery_target.csv")
if (!file.exists(discovery_path)) {
  cli::cli_abort("Run 00_discover_fisheries.R first -- {.file {discovery_path}} not found.")
}

target_fisheries <- read_csv(discovery_path, show_col_types = FALSE) |>
  filter(basin_match == "target", include_in_run) |>
  pull(fishery_name_raw) |>
  unique()

cli::cli_alert_info("Target fishery-years: {length(target_fisheries)}")

cli::cli_alert_info("Connecting to internal DB...")
conn <- creelutils::connect_creel_db()
on.exit(try(DBI::dbDisconnect(conn), silent = TRUE), add = TRUE)

# ------------------------------------------------------------------------------
# Schema discovery
# ------------------------------------------------------------------------------
# information_schema is the portable route and gives table + column in one
# query, which is what the scoring needs. dbListTables()/dbListFields() is the
# fallback for a driver that does not expose it -- slower (one round trip per
# table) but it works anywhere.

all_columns <- tryCatch(
  DBI::dbGetQuery(conn, "
    SELECT table_schema, table_name, column_name, data_type
    FROM information_schema.columns
  ") |> as_tibble() |> rename_with(tolower),
  error = function(e) {
    cli::cli_alert_warning("information_schema unavailable ({conditionMessage(e)}); falling back to dbListTables().")
    tbls <- DBI::dbListTables(conn)
    map(tbls, function(t) {
      flds <- tryCatch(DBI::dbListFields(conn, t), error = function(e) character(0))
      if (length(flds) == 0) return(NULL)
      tibble(table_schema = NA_character_, table_name = t, column_name = flds, data_type = NA_character_)
    }) |> bind_rows()
  }
)

if (nrow(all_columns) == 0) {
  cli::cli_abort("Could not read any table metadata from the connection.")
}
cli::cli_alert_success("Read metadata for {n_distinct(all_columns$table_name)} table(s).")

# What a fishery-by-location lookup looks like, in columns rather than in its
# name. Deliberately loose -- the point is to surface candidates for a human to
# read, not to auto-pick.
PAT_FISHERY  <- regex("fishery", ignore_case = TRUE)
PAT_SECTION  <- regex("section", ignore_case = TRUE)
PAT_LOCATION <- regex("location|site|reach|water_body|waterbody|river", ignore_case = TRUE)

candidates <- all_columns |>
  group_by(table_schema, table_name) |>
  summarise(
    n_cols          = n(),
    fishery_cols    = paste(column_name[str_detect(column_name, PAT_FISHERY)],  collapse = "|"),
    section_cols    = paste(column_name[str_detect(column_name, PAT_SECTION)],  collapse = "|"),
    location_cols   = paste(column_name[str_detect(column_name, PAT_LOCATION)], collapse = "|"),
    all_cols        = paste(column_name, collapse = "|"),
    .groups = "drop"
  ) |>
  mutate(
    has_fishery  = nzchar(fishery_cols),
    has_section  = nzchar(section_cols),
    has_location = nzchar(location_cols),
    # Name is a tiebreaker, not the test: a table can be the right one without
    # "lut" anywhere in its name, and a table named "..._lut" can be unrelated.
    name_hint    = str_detect(table_name, regex("lut|lookup|location", ignore_case = TRUE)),
    score        = 2L * has_fishery + 2L * has_section + has_location + name_hint
  ) |>
  filter(has_fishery, has_section | has_location) |>
  arrange(desc(score), table_name)

write_csv(candidates, file.path(OUT_DIR, "location_lut_candidates.csv"))

cli::cli_h2("Candidate tables")
if (nrow(candidates) == 0) {
  cli::cli_abort(c(
    "No table has both a fishery-shaped and a section/location-shaped column.",
    "i" = "Read {.file {file.path(OUT_DIR, 'location_lut_candidates.csv')}} and set LUT_TABLE by hand."
  ))
}
candidates |> select(table_schema, table_name, score, fishery_cols, section_cols, location_cols) |> print(n = 30)

chosen <- if (!is.null(LUT_TABLE)) {
  LUT_TABLE
} else if (nrow(candidates) == 1 || candidates$score[1] > candidates$score[2]) {
  paste(na.omit(c(candidates$table_schema[1], candidates$table_name[1])), collapse = ".")
} else {
  # A tie means the evidence does not single one out. Stopping here with the
  # table printed is more useful than picking one and capturing the wrong rows.
  cli::cli_abort(c(
    "Top {sum(candidates$score == candidates$score[1])} candidates tie on score -- not guessing.",
    "i" = "Pick one from the table above, set {.code LUT_TABLE} at the top of this script, and re-run."
  ))
}
cli::cli_alert_success("Using table {.val {chosen}}")

# ------------------------------------------------------------------------------
# Schema of the chosen table
# ------------------------------------------------------------------------------

chosen_parts  <- str_split(chosen, fixed("."))[[1]]
chosen_schema <- if (length(chosen_parts) == 2) chosen_parts[1] else NA_character_
chosen_name   <- tail(chosen_parts, 1)

lut_schema <- all_columns |>
  filter(table_name == chosen_name,
         is.na(chosen_schema) | is.na(table_schema) | table_schema == chosen_schema)
write_csv(lut_schema, file.path(OUT_DIR, "location_lut_schema.csv"))

cli::cli_h2("Schema of {chosen}")
lut_schema |> select(column_name, data_type) |> print(n = 100)

# ------------------------------------------------------------------------------
# Capture
# ------------------------------------------------------------------------------
# Filter server-side on the fishery-name column when there is one, so the
# capture stays small and the committed lookup contains only this analysis's
# scope. Quoting goes through the driver (dbQuoteIdentifier / dbQuoteString) --
# fishery names carry spaces and can carry apostrophes.

# as.character() on the DBI SQL objects: they are S4, and paste()'s coercion of
# them is driver-dependent.
quote_id <- function(x) as.character(DBI::dbQuoteIdentifier(conn, x))
qualified <- if (is.na(chosen_schema)) {
  quote_id(chosen_name)
} else {
  paste0(quote_id(chosen_schema), ".", quote_id(chosen_name))
}

# Prefer an exact "fishery_name"; otherwise the narrowest fishery+name column.
name_col <- {
  fcols <- lut_schema$column_name[str_detect(lut_schema$column_name, PAT_FISHERY)]
  exact <- fcols[tolower(fcols) == "fishery_name"]
  if (length(exact)) exact[1]
  else {
    named <- fcols[str_detect(fcols, regex("name", ignore_case = TRUE))]
    if (length(named)) named[1] else NA_character_
  }
}

lut <- if (!is.na(name_col)) {
  cli::cli_alert_info("Filtering on {.field {name_col}} for {length(target_fisheries)} fishery-year(s).")
  sql <- paste0(
    "SELECT * FROM ", qualified,
    " WHERE ", quote_id(name_col),
    " IN (", paste(as.character(DBI::dbQuoteString(conn, target_fisheries)), collapse = ", "), ")"
  )
  DBI::dbGetQuery(conn, sql) |> as_tibble()
} else {
  cli::cli_alert_warning(
    "No fishery-name column found; pulling the whole table (capped at \\
     {MAX_UNFILTERED_ROWS} rows) so you can see how it keys to a fishery."
  )
  DBI::dbGetQuery(conn, paste0("SELECT * FROM ", qualified)) |>
    as_tibble() |>
    head(MAX_UNFILTERED_ROWS)
}

if (nrow(lut) == 0) {
  cli::cli_abort(c(
    "{.val {chosen}} returned no rows for the target fishery-years.",
    "i" = "The table may key on a fishery id rather than a name, or use different name strings.",
    "i" = "Check {.file {file.path(OUT_DIR, 'location_lut_schema.csv')}} and set LUT_TABLE / name_col accordingly."
  ))
}

lut <- lut |>
  mutate(
    source_table = chosen,
    captured_at  = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )

write_csv(lut, LOOKUP_PATH)

cli::cli_h2("Capture summary")
cli::cli_alert_success("{nrow(lut)} row(s), {ncol(lut)} column(s) from {.val {chosen}}")
if (!is.na(name_col)) {
  missing <- setdiff(target_fisheries, unique(lut[[name_col]]))
  if (length(missing) > 0) {
    # Recorded rather than silently tolerated: a fishery-year absent from the
    # location lookup is itself a finding for the across-year comparison.
    cli::cli_alert_warning("{length(missing)} target fishery-year(s) have NO row in this table:")
    cli::cli_ul(missing)
  }
  lut |> count(.data[[name_col]], name = "n_rows") |> print(n = 60)
}

cli::cli_h3("First rows")
print(head(lut, 10), width = Inf)

cli::cli_alert_success("Wrote {.file {LOOKUP_PATH}}")
cli::cli_alert_warning(
  "COMMIT THIS FILE, plus outputs/location_lut_schema.csv if you want the \\
   schema visible in review. Without the lookup, the across-year location \\
   comparison cannot run on the external (VPN-free) path."
)
