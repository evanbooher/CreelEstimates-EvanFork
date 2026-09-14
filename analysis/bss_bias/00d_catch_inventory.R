# ==============================================================================
# 00d_catch_inventory.R
#
# Purpose:
#   Before spending MCMC on a catch group, find out whether that group HAS any
#   data in each fishery-year. Totals of Chinook encounters and Coho harvest
#   per fishery-year, straight from the catch records.
#
#   The prompting case: Stillaguamish is a "salmon and gamefish" fishery whose
#   current proposal is NF gamefish and mainstem coho. Chinook encounters there
#   may be zero or near-zero, in which case fitting that group is wasted time
#   and any C_sum it produced would be a number nobody should quote.
#
# What these numbers ARE:
#   RAW INTERVIEW CATCH -- the fish actually reported by interviewed anglers,
#   summed. NOT an expanded season estimate. They are a sample, so they answer
#   "is there anything here, and roughly how much" and nothing more. The
#   expanded estimate is C_sum, which only the BSS fit produces.
#
# Usage:
#   Rscript analysis/bss_bias/00d_catch_inventory.R
#   or, from the Console:  source("analysis/bss_bias/00d_catch_inventory.R")
#
#   Reads the same cached DWG fetches 01 uses, so no DB and no VPN for any
#   fishery-year already fetched. One that has never been fetched is reported
#   as uncached rather than triggering a fetch, so this stays a fast read-only
#   pass. Run 01 (or clear USE_DWG_CACHE) to populate one.
#
# Outputs:
#   bss_catch_inventory.csv         -- fishery-year x catch group totals
#   bss_catch_inventory_species.csv -- the raw species x life_stage x fin_mark x
#                                      fate breakdown, for seeing what the
#                                      group patterns actually matched
# ==============================================================================

library(tidyverse)
library(cli)
library(here)

OUT_DIR <- here::here("analysis", "bss_bias", "outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

source(here::here("analysis", "bss_bias", "fishery_data.R"))
source(here::here("analysis", "bss_bias", "catch_groups.R"))

# Default: every target fishery-year in the discovery CSV. Narrow with a regex.
if (!exists("FISHERY_RE", inherits = FALSE)) FISHERY_RE <- "."

disc_path <- file.path(OUT_DIR, "fishery_discovery_target.csv")
if (!file.exists(disc_path)) {
  cli::cli_abort("{.file {disc_path}} not found -- run 00_discover_fisheries.R first.")
}

targets <- read_csv(disc_path, show_col_types = FALSE) |>
  filter(basin_match == "target", include_in_run, str_detect(fishery_name, FISHERY_RE)) |>
  pull(fishery_name) |>
  unique() |>
  sort()

cli::cli_h1("00d -- catch inventory")
# Echo the effective filter. FISHERY_RE uses the exists() idiom so it can be
# set before sourcing -- which also means a value left over from an earlier
# 01b run silently narrows this scan. Printing it is the difference between
# noticing that and quietly concluding a basin has no data.
cli::cli_alert_info("Fishery filter: {.val {FISHERY_RE}}{if (FISHERY_RE == '.') ' (all)' else ' -- set FISHERY_RE or rm() it to widen'}")
cli::cli_alert_info("Scanning {length(targets)} fishery-year{?s}.")

# Read-only: report an uncached fishery-year rather than fetching it, so this
# pass cannot silently turn into a long DB run.
cached_dwg <- function(fishery_name) {
  est_dates <- resolve_window(fishery_name)
  if (is.null(est_dates)) return(NULL)
  key <- paste0(safe_name(fishery_name), "__", DATA_SOURCE, "__",
                est_dates$est_date_start, "_", est_dates$est_date_end, ".rds")
  path <- file.path(DWG_CACHE_DIR, key)
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

one_fishery <- function(fishery_name) {
  dwg <- cached_dwg(fishery_name)
  if (is.null(dwg) || is.null(dwg$catch) || !is.data.frame(dwg$catch)) {
    return(tibble(fishery_name = fishery_name, catch_group = NA_character_,
                  est_cg = NA_character_, n_records = NA_integer_,
                  n_fish = NA_real_, n_interviews = NA_integer_,
                  status = "not cached -- run 01 for this fishery-year first"))
  }
  catch <- dwg$catch
  needed <- c("species", "life_stage", "fin_mark", "fate", "fish_count")
  if (!all(needed %in% names(catch))) {
    return(tibble(fishery_name = fishery_name, catch_group = NA_character_,
                  est_cg = NA_character_, n_records = NA_integer_,
                  n_fish = NA_real_, n_interviews = NA_integer_,
                  status = paste0("catch table missing: ",
                                  paste(setdiff(needed, names(catch)), collapse = ", "))))
  }

  imap_dfr(CATCH_GROUPS, function(g, key) {
    hits <- match_catch_group(catch, g)
    tibble(
      fishery_name = fishery_name,
      catch_group  = key,
      est_cg       = catch_group_label(g),
      n_records    = nrow(hits),
      n_fish       = sum(hits$fish_count, na.rm = TRUE),
      n_interviews = if ("interview_id" %in% names(hits)) n_distinct(hits$interview_id) else NA_integer_,
      status       = if (nrow(hits) == 0) "NO RECORDS -- do not fit this group" else "ok"
    )
  })
}

inventory <- map_dfr(targets, function(fn) {
  out <- try(one_fishery(fn), silent = TRUE)
  if (inherits(out, "try-error")) {
    cli::cli_alert_warning("{fn}: {conditionMessage(attr(out, 'condition'))}")
    return(tibble(fishery_name = fn, catch_group = NA_character_, est_cg = NA_character_,
                  n_records = NA_integer_, n_fish = NA_real_, n_interviews = NA_integer_,
                  status = "error reading cache"))
  }
  out
})

inventory <- inventory |>
  mutate(
    year_start = as.integer(str_extract(fishery_name, "\\d{4}")),
    fishery_type = str_squish(str_replace(fishery_name, "\\d{4}(-\\d{2,4})?", " "))
  ) |>
  arrange(fishery_type, year_start, catch_group) |>
  relocate(fishery_type, year_start, .after = fishery_name)

write_csv(inventory, file.path(OUT_DIR, "bss_catch_inventory.csv"))

# The raw breakdown, so the group patterns can be checked against what is
# actually coded in the data rather than trusted.
species_detail <- map_dfr(targets, function(fn) {
  dwg <- cached_dwg(fn)
  if (is.null(dwg) || is.null(dwg$catch) || !is.data.frame(dwg$catch)) return(NULL)
  if (!all(c("species", "fish_count") %in% names(dwg$catch))) return(NULL)
  dwg$catch |>
    mutate(across(any_of(c("species", "life_stage", "fin_mark", "fate")),
                  ~replace_na(as.character(.), "NA"))) |>
    count(across(any_of(c("species", "life_stage", "fin_mark", "fate"))),
          wt = fish_count, name = "n_fish") |>
    mutate(fishery_name = fn, .before = 1)
})

if (!is.null(species_detail) && nrow(species_detail) > 0) {
  write_csv(species_detail, file.path(OUT_DIR, "bss_catch_inventory_species.csv"))
}

# ------------------------------------------------------------------------------
# Console view -- the wide table to actually look at
# ------------------------------------------------------------------------------

cli::cli_h2("Fish counted in interviews, by fishery-year and catch group")

wide <- inventory |>
  filter(!is.na(catch_group)) |>
  select(fishery_type, year_start, catch_group, n_fish) |>
  pivot_wider(names_from = catch_group, values_from = n_fish)

print(wide, n = Inf)

empty <- inventory |> filter(!is.na(catch_group), n_records == 0)
if (nrow(empty) > 0) {
  cli::cli_h2("Groups with NO records -- skip these in 01b")
  empty |> select(fishery_name, catch_group) |> print(n = Inf)
}

uncached <- inventory |> filter(str_detect(status, "not cached"))
if (nrow(uncached) > 0) {
  cli::cli_h2("Not cached -- fetched by 01, not by this script")
  uncached |> distinct(fishery_name) |> print(n = Inf)
}

cli::cli_alert_info(
  "These are RAW INTERVIEW CATCH totals, not expanded season estimates. \\
   They say whether a group is worth fitting, not how many fish were caught."
)
