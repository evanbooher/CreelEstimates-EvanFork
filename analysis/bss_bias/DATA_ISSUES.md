# Data issues found while building the `b` analysis

Findings about the underlying creel data, not about this analysis's code. Kept
here so they survive the session they were found in. Each one names how it was
found so it can be re-checked.

Nothing here is fixed. Some are worked around in `01_fit_bss_bias.R` /
`scope_rules.R` and say so; the rest are open.

---

## 1. Missing closure records — affects the ESTIMATES, not just this analysis

**Status: open. The most consequential item here.**

**Confirmed against the raw tables.** Stillaguamish 2022-23 has **no effort
row and no interview row anywhere in the basin from 1 October to 15 November
2022** — 46 consecutive days out of a 91-day estimation window, **51% of it**.
Sampling is two blocks: 23 dates in September, 10 dates from 16 to 28 November.
October is entirely absent. The fishery lead's reading is that this stretch was
**closed**, and that the closure dates are absent from the database.

Neither of the ways a record could have hidden from `bss_b_survey_days.csv`
applies here: `tie_in_indicator` is clean (2772 `FALSE`, 77 `TRUE`, no `NA`)
and neither table has a missing `section_num`. Nothing was dropped — the data
are simply not there.

Both census dates that year, 1 and 17 September, fall in the opening fortnight,
so that fishery-year's `b` is anchored entirely in the first 17 days of the
window.

### Why it matters

`prep_days()` assumes a day is OPEN unless a closure record exists:

```r
rows_update(
  expand_grid(event_date = days$event_date, section_num = sections, open = TRUE),
  closures |> ... |> mutate(open = FALSE),
  by = c("section_num", "event_date")
)
```

So a missing closure record is indistinguishable from a genuinely open day. The
resulting `O[s,d]` is 1 for every day of that gap, and `O` multiplies latent
effort in the Stan model:

```stan
lambda_E_S[s][d,g] = exp(mu_E[g,s] + omega_E[period[d]][g,s] + B1*w[d]) * O[d,s];
```

With `O = 1` the model fits low-but-nonzero effort across a stretch the
fishery was shut, instead of switching it off. Consequences, in rough order of
severity:

- **the fishery-level effort and catch expansion includes ~5 weeks of
  non-fishery**, which is a bias in the estimates themselves, not a
  presentation problem;
- `n_days_open` in `bss_b_comparability.csv` overstates the season;
- the AR(1) effort process spends a third of a 91-day window on a period that
  should not be in it.

It does **not** bias `b` directly — `b` is a ratio and there are no counts in
the gap to contribute to it — but it does move the effort estimate that the
expansion rests on.

### How to confirm

Re-run of the check that confirmed it, kept for the next fishery-year that
needs the same treatment:

```r
library(tidyverse); library(here)
source(here("analysis/bss_bias/fishery_data.R"))
fn  <- "Stillaguamish salmon and gamefish 2022-23"
d   <- resolve_window(fn); dwg <- fetch_fishery_dwg(fn, d)
win <- function(x) filter(x, between(event_date, as.Date(d$est_date_start), as.Date(d$est_date_end)))

bind_rows(
  win(dwg$effort)    |> transmute(event_date, src = "effort"),
  win(dwg$interview) |> transmute(event_date, src = "interview")
) |> count(event_date, src) |>
  pivot_wider(names_from = src, values_from = n, values_fill = 0) |>
  arrange(event_date) |> print(n = 100)

win(dwg$effort) |> count(tie_in_indicator)          # NA -> counted as neither index nor census
win(dwg$effort) |> summarise(na_section = sum(is.na(section_num)))
```

### Options

1. **Add the closure dates to the database.** The real fix; everything
   downstream picks it up with no code change.
2. **Trim the estimation window** for that fishery-year in
   `lookup/fishery_params.csv`. Faster, but it changes what the window means
   and has to be recorded as a deviation.
3. **Check every other fishery-year for the same pattern.** Now automatic:
   `02b_survey_coverage.R` writes `bss_b_no_data_gaps.csv`, every run of 7+
   consecutive days with no record in any section, with the run's share of the
   estimation window. It warns on the console when any are found.

`bss_b_no_data_gaps.csv` is a FLAG, not a finding — a run may be a real
closure, a genuine mid-season break, or missing records. It says where to look.
Worth clearing before trusting `n_days_open` anywhere.

---

## 2. Case-variant type labels in the location lookup

**Status: worked around in `02a`, not fixed at source.**

`creel.vw_fishery_location` carries `location_type` as both `Site` and `site`,
and `survey_type` as both `Census` and `census` — two rows, both in Skagit fall
salmon 2023. The pipeline does exact-match string tests on these, e.g.
`prep_dwg_census_expan()` filters `location_type == "Site"`, so an odd-cased row
is **silently dropped**.

`02a_location_lut_changes.R` folds the case and writes the offending rows to
`bss_b_lut_anomalies.csv`.

---

## 3. One section on two water bodies

**Status: open; blocks a clean scope rule.**

Stillaguamish 2023-24 section 5 sits on **both the North and South Forks**, so
`(year, section_num)` is not a unique key for that fishery-year, and its census
count totals both forks together.

This is why Stillaguamish runs whole-basin: a rule excluding the South Fork
resolves to section numbers, so it drops the SF in 2022 and 2024 but keeps five
SF index sites in 2023 — creating the comparability break it was meant to
prevent. Dropping those five sites instead would be worse, since the census
count still covers both forks and `b` would absorb the mismatch. See
`scope_rules.R`.

---

## 4. Section numbers are not stable between years

**Status: handled in analysis; a property of the data.**

A section number is a label on a block, and inserting a block renumbers
everything above it. Skagit fall salmon 2021→2022 shifts `s3→s4, s4→s5, s5→s6`
when the Hwy 9 split adds a block low down; Stillaguamish renumbers in every
year-pair.

`bss_b_lut_block_lineage.csv` matches blocks between years on **site
membership** rather than number, and separates a renumbering from a genuine
redraw. Any join or comparison keyed on `(year, section_num)` across years is
wrong without it.

---

## 5. `p_TI` is 1 everywhere

**Status: fact, not a defect — but worth knowing.**

Every `p_census_bank` / `p_census_boat` value in the captured lookup is 1 or
missing, and `prep_dwg_census_expan()` fills missing with 1. So the census
spatial-coverage expansion is inert in these fisheries, and `b` carries the
whole index-vs-census discrepancy.

---

## 6. A fishery-year with interviews but no effort counts

**Status: fails cleanly; cause not investigated.**

Skagit spring Chinook 2021 lower has 8 surveyed days in a 13-day window, all
interviews, and **zero effort count records**. `01` skips it at preflight ("No
effort count records within the estimation window") and it has no `b`. Whether
the counts were never done or never entered is unknown.

---

## 7. Angler types that do not classify

**Status: worked around in `01`, recorded per fishery-year.**

`prep_dwg_interview_angler_types()` ends its `case_when` with `TRUE ~ "fail"`
and never filters those rows, so an interview the rules cannot classify — a
missing `boat_used`, or an `NA` `boat_type`, since `str_detect(NA, .)` is `NA`
rather than `FALSE` — becomes a third angler type and breaks the BSS
likelihood, which is written for two.

`01_fit_bss_bias.R` recodes against fixed `bank`/`boat` levels and writes the
counts to `bss_b_angler_recode.csv`. A high drop rate there is a data-entry
signal worth following up.

---

## 8. `fish_count` is stored as text in some fishery-years

**Found:** 2026-09-14, while building `00d_catch_inventory.R`.

**What happens:** `dwg$catch$fish_count` comes back as `character` rather than
numeric for at least Skagit spring Chinook 2021 lower. Any arithmetic on it
aborts — `sum()` raises `invalid 'type' (character) of argument`.

**Why it matters:** a column whose type varies by fishery-year is a trap for
every consumer, not just this one. Code that happens to be written against a
numeric year works, and then fails the first time it meets a character year.
The failure is at least loud; the dangerous version is a silent coercion that
turns unparseable entries into `NA` and quietly drops fish from a total.

**Handled here:** `00d_catch_inventory.R` coerces explicitly, counts the
fishery-years affected and the number of values that would not parse, and
reports both at the end of its run rather than absorbing them.

**Not yet investigated:** whether the character values are all clean integers
(so the type is cosmetic) or whether some carry text that represents real
information — a range, a qualifier, a sentinel. Until someone looks, treat
totals for the affected fishery-years as provisional. The BSS fits reach
`fish_count` through a different path and are not known to be affected, but
that has not been verified either.

## 9. Stillaguamish 2022-23 covers sections the rest of the series does not

**What happens:** the 2022-23 fetch returns sections beyond the mainstem
(1–3) that the other years of `Stillaguamish salmon and gamefish` do not
carry. It is also the year whose closure table produced unmatched
`(section_num, event_date)` keys — 122 of them — which aborted `prep_days()`
until `rows_update(unmatched = "ignore")` went in.

**Why it matters:** a `b` series is only meaningful if every year in it
describes the same fishery. A year that includes water the others exclude has
a different mix of bank and boat effort and a different relationship between
index counts and anglers, so its `b` is not exchangeable with the rest — which
is precisely the assumption the whole comparability argument rests on.

It also runs the longest window in the series -- 2022-09-01 to 2022-11-30, 91
days, against 76 for 2024-25 and 61 for 2025-26 -- with a sparse tail of late
November dates that the other years do not have.

**Handled here:** a window restriction in `scope_rules.R` holds that
fishery-year to 2022-09-01 / 2022-09-30 with `trim_to_last_sampled = TRUE`, and
a catch-group exclusion drops `chinook_all` (no Chinook were caught in the
September part of the fishery). One definition, read by both
`01_fit_bss_bias.R` (the `b` fit) and `10_render_fw_creel.R` (the production
render) — without that, the `b` in the brief would come from one window and the
season totals from another.

**Its sections are NOT restricted.** Per direction the North and South Forks
stay in for every Stillaguamish year; the problem with 2022-23 is its window,
not its water. The section-restriction mechanism exists and `scope_rules.R`
records that a mainstem-only rule *could* be applied evenly across all four
years (the mainstem shares no section with either fork in any year), which the
rejected MS + NF rule could not — but that changes every Stillaguamish `b` and
belongs to the fork-specific follow-up.

Trailing unsampled days are not harmless: `prep_days()` puts every one of them in
the day grid, and the BSS estimates effort and catch for days nobody visited on
the strength of the priors alone.

**Not yet investigated:** whether the extra sections are tributaries, a
different reach naming convention, or a data-entry artifact; and whether any
*other* fishery-year in the three basins has the same problem. The section
sets per fishery-year have not been compared systematically — if they differ
elsewhere, the affected `b` estimates carry the same objection and the fix is
the same override.
