# BSS effort-bias (`b`) temporal-stability analysis

Branch: `explore/bss-effort-bias-stability` only. Unrelated to, and does not
depend on, `chore/multi-fishery-trip-summary` (the PST/consultant deliverable
branch) -- some code patterns below were transcribed from that branch's
`analysis/pst/02_ingest/multi_fishery_creel_summary.R` for reference, but
nothing here imports or merges from it.

Built for a DFW/tribal technical-staff meeting deciding how to use historical
years of data for the bias-correction component of the freshwater creel
Bayesian State-Space (BSS) estimation framework, for the Stillaguamish,
Snohomish, and (regional-context) Skagit basins.

## What "the bias term" is

In the BSS Stan model (`vector<lower=0>[G] b`, prior `lognormal(0,
value_lognormal_sigma_b)`, prior median 1 = "no bias"), `b` is used
POSITIONALLY, not by angler type:

- `b[1]` -- bias in the **vehicle** index count (road counts of cars vs. the
  true angler-vehicle relationship)
- `b[2]` -- bias in the **trailer** index count

This is the "effort bias term" the meeting is about. It is a **different**
parameter from `p_TI` / `TI_expan`, the census:index tie-in spatial-coverage
ratio (`prep_dwg_census_expan.R` calls that "the bias parameter/expansion
factor" in its own comments) -- the two sit adjacent in the code and in the
same likelihood line, but this analysis is about `b` specifically, per an
explicit scoping decision. `p_TI` is tracked in the comparability table
(Phase 3) because a `p_TI` change year-to-year would distort `b`, not because
`p_TI` itself is the thing being compared.

`b` is confounded with `R_V` (vehicles-per-angler) in the vehicle-count
likelihood and is only pinned down by the interview-side binomial on
`V_A`/`A_A`. A fishery-year with few/no vehicle-count interviews will
produce a `b` posterior that mostly reflects its prior, not its data --
every output table carries a `prior_contraction` / `informed_flag` column
for exactly this reason. Read it before trusting a point estimate.

## Stan model used

**`BSS_creel_model_02_2021-01-22_ppc.stan`**, standardized across every
fishery-year -- this is the file confirmed to be in current production use.
`b`'s declaration/prior/likelihood placement is identical across all four
model files in `stan_models/` (verified directly against the `_ppc` file's
source), so this choice doesn't change what `b` means. It does declare `O`
as `matrix[D,S]` (2-D), matching `prep_inputs_bss()`'s unmodified output --
compatible. `BSS_creel_model_02_2024-07-24.stan` is excluded from this
analysis entirely: it declares `O` as `real O[D,S,G]` (3-D), a hard
dimension mismatch against `prep_inputs_bss()`, not a subtle bug -- don't
point this pipeline at it without first updating `prep_inputs_bss()` itself.

The `_ppc` variant carries a `*_rep` generated-quantities block (still
computed during sampling regardless of `pars=`, so somewhat slower
per-iteration than the non-`_ppc` files) and has no compiled `.rds` cache in
`stan_models/` yet -- expect a one-time ~1-3 min recompile on the first fit.

The full prior vector is copied verbatim from `template_scripts/fw_creel.Rmd`
and held constant across every fishery-year (see `BSS_PRIORS` in
`01_fit_bss_bias.R`) -- changing `value_lognormal_sigma_b` between years
would make the whole comparison meaningless.

## Catch-group selection

`b` is fit once per fishery-year, against a **fixed target catch group keyed
off the fishery name** -- not pooled total salmon, not derived from observed
species in the data. Per explicit direction, `fishery_target_catch_group()`
in `01_fit_bss_bias.R` applies, in order:

| Fishery name contains | Target catch group (`est_cg`) |
|---|---|
| `Chinook` | `Chinook_Adult_AD_Kept` |
| `fall salmon`, or basin = Stillaguamish (its "salmon and gamefish" naming never contains "fall salmon" literally) | `Coho_Adult_AD\|UM_Kept` |
| `sockeye` | `Sockeye_Adult_AD\|UM_Kept` |

Fate is **harvest only** (`Kept`) for every rule -- a deliberate narrowing
from an earlier pooled-total design that used `Kept|Released` on the
reasoning that BSS catch is total encounters, not harvest. Chinook is
`AD`-only (hatchery-marked) while Coho/Sockeye are `AD|UM` (both marks) --
intentional given wild/unmarked Chinook are typically release-only under
regulation in these fisheries. A fishery name matching none of these rules
hard-stops the whole run (`cli_abort`, not a per-fishery skip) rather than
silently proceeding with no catch group -- surfacing a naming gap
immediately. This matters for `b` specifically (not just for catch) because
the interview rows feeding `IntA`/`V_A`/`T_A` are filtered by `est_cg`
upstream of `prep_inputs_bss()`. The chosen `est_cg` string is recorded in
every output row (`chosen_est_cg` / `est_cg` columns).

## Pipeline

Run in order:

1. **`00_discover_fisheries.R`** -- no VPN, no DB, no Stan. Pulls the full
   historical fishery-name registry from the public "WDFW Creel Fishery
   Manager" dataset (`data.wa.gov`, Socrata id `vkjc-s5u8`), filters to
   Skagit/Snohomish/Stillaguamish (+ adjacent systems), writes
   `outputs/fishery_discovery_target.csv`. **Review this file by hand** and
   edit its `include_in_run` column before step 2.

   Also worth 5 minutes: browse `https://data.wa.gov/browse?q=creel` --
   `vkjc-s5u8` is a fishery *registry*; a sibling dataset with actual
   effort/interview-level data may exist and could feed step 2 directly
   without a DB pull at all. See the "public-data path" note below.

2. **`01_fit_bss_bias.R`** -- the main driver. For each included
   fishery-year: fetch data, prep BSS inputs, fit the Stan model, extract
   `b`. Writes a **comparability row per fishery-year BEFORE fitting**, so
   step 3 doesn't need any fit to have succeeded. Everything is written
   incrementally (append-as-you-go) so an interrupted overnight run doesn't
   lose completed fishery-years.

3. **`02_build_comparability_table.R`** -- no MCMC required, only needs step
   2 to have run its data-prep stage for at least one fishery. Produces the
   standalone comparability table (date windows, CRC areas, `p_TI`, count
   availability, an explicit `comparability_tier` per fishery-year). **This
   table alone is a legitimate meeting deliverable** even with zero
   completed fits -- see "Go/no-go" below.

4. **`03_plot_b_series.R`** -- the top-priority deliverable. Six figures;
   Figure 1 + 2 composite is the actual slide. Needs `dataviz`-skill-guided
   design choices already baked in (see the script's header for how a
   web/CSS-oriented skill was adapted to static R/ggplot output).

5. **`04_candidate_options.R`** -- the four candidate bias-correction options
   (mean / most-recent / precision-weighted / hierarchical-shrinkage) per
   fishery-series. Re-run `03_plot_b_series.R` after this to get Figure 6
   (options overlaid on the series).

## Public-data path -- confirmed working, VPN not required

**Checked locally (VPN off):** `creelutils::fetch_data(fishery_name = "Skagit
fall salmon 2024", data_source = "external")` is the public/`data.wa.gov`
path, and it returns BSS-grade data. `dwg$effort` (5247x25) and `dwg$catch`
(964x11) pass the full column checklist in `01_fit_bss_bias.R`'s header
(including `species`/`life_stage`/`fin_mark`/`fate` on `catch`). `dwg$interview`
(1811x39) has every required raw column except `fishing_time_total` /
`person_count_final`, which are computed downstream from columns that ARE
present (`fishing_start_time`/`fishing_end_time`, `angler_count`/
`total_group_count`) -- not a gap.

`dwg$ll`'s `centroid_lat`/`centroid_lon` and `dwg$fishery_manager`'s
`p_census_bank`/`p_census_boat` were not independently re-verified against
the public path, but per direct confirmation the latter is moot for `b`
either way: `p_census_bank`/`boat` feed the `p_TI`/census-tie-in spatial-
coverage correction (a defunct-in-practice, DIFFERENT parameter -- see "What
'the bias term' is" above), not `b`'s identification.

So `DATA_SOURCE <- "external"` in `01_fit_bss_bias.R` is a viable VPN-free
option for the whole pipeline. The default remains `DATA_SOURCE <-
"internal"` since VPN/DB access is available locally anyway -- switch to
`"external"` only if running somewhere without VPN.

## Same-day feasibility -- read this before starting a full historical run

**No cached BSS posterior fits exist anywhere in this repo for these three
basins** (only rstan's compiled-model `.rds` caches in `stan_models/`, not
posterior draws). A full historical run -- 3 basins x N years x up to 4
concurrent named Skagit fisheries, each a fresh MCMC fit -- is unlikely to
finish before the meeting. Triage in this order:

- **T0 (15 min, do this first):** search your own machine / network drives /
  prior project folders for any saved `stanfit` `.rds` from a past BSS run of
  these fisheries. If any turn up, `get_bss_bias()` extracts `b` from them in
  seconds and steps 3-5 complete immediately for those years.
- **T1 (30 min):** the public-data-path check above -- DONE, confirmed
  working (`data_source = "external"`); not on this run's critical path
  since VPN is available, but ruled out as a blocker either way.
- **T2 (30-60 min):** run ONE fishery-year (suggest `"Skagit fall salmon
  2024"`) end-to-end at `FIT_CONFIG_NAME <- "smoke"` in `01_fit_bss_bias.R`
  to prove the pipeline, then re-run the SAME fishery at `"quick"` and time
  it. That number x the fishery-year count is your actual budget. **If a
  single `quick` fit exceeds ~45 minutes, stop the full historical run and
  go straight to T3 step 1 only** -- put remaining time into the
  comparability table and plots instead.
- **T3:** prioritized subset, in order, stopping when time runs out:
  1. Most recent 3 years x 3 basins, primary fishery only (Skagit fall
     salmon, Snohomish fall salmon, Stillaguamish salmon-and-gamefish) -- 9
     fits. Minimum viable story.
  2. Skagit's other named fisheries (spring Chinook upper/lower, summer
     sockeye), most recent 2 years -- 6 fits. This is what makes the
     "separate series per named fishery" point concrete.
  3. Backfill earlier years, most-recent-first, until out of time.
- **T4:** kick off remaining fits overnight in a detached R session; steps
  3-5 re-run cleanly against whatever's in `outputs/` at 7am, no edits needed.

### Go/no-go -- what ships even in the worst case

| Tier | Needs | Contents |
|---|---|---|
| Floor | Steps 1 + 3 only, **no fits** | The comparability table alone -- a legitimate deliverable: it defines which fishery-years are even comparable, which the group has to agree on before any `b` number matters |
| Minimum viable | T3 step 1 (~9 fits) | + Figure 1 (3-yr series, 3 basins) + Figure 3 + options (a)/(b) |
| Target | T3 steps 1-2 (~15 fits) | + Skagit's concurrent tracks + Figures 2/4/5 + all four options with I²/tau² |
| Stretch | Full backfill | + complete series + a scoped (d2) in-model hierarchical `b` proposal |

Frame the meeting outcome accordingly: **the defensible ask for tomorrow is
agreement on the *rule* (which option) and the *comparability criteria*,
with numbers backfilled at production fit quality afterward.**

## Net-new vs. reused

| | |
|---|---|
| Reused unchanged | `resolve_dates.R`, `prep_days.R`, `prep_dwg_interview_fishing_time.R`, `prep_dwg_interview_angler_types.R`, `prep_dwg_interview_catch.R`, `prep_dwg_effort_index.R`, `prep_dwg_effort_census.R`, `prep_dwg_census_expan.R`, `prep_inputs_bss.R` |
| Reused, one additive change | `fit_bss.R` -- added `pars`/`include` args, forwarded to `stan()` (restricts monitored parameters; the main speed/memory lever available without touching chains/iter) |
| Net-new function | `R_functions/get_bss_bias.R` -- extracts `b[1]`/`b[2]` posterior summaries from a stanfit |
| Net-new function | `R_functions/drop_na_bss_inputs.R` -- drops NA rows from prep_inputs_bss()'s observation-level Stan inputs before fitting (Stan hard-errors on any NA in `data`); logs what it dropped per fishery-year to `bss_b_na_drops.csv` |
| Net-new scripts | everything in `analysis/bss_bias/` |
| Transcribed patterns (not imports) | `skip_fishery()`/`run_stage()` condition-class error handling, the run ledger, `resolve_study_design()`, `preflight_fishery()` -- from `chore/multi-fishery-trip-summary`'s `multi_fishery_creel_summary.R`, read-only reference |
| Explicitly deferred | (d2), an in-model hierarchical re-specification of `b` itself across years -- flagged as a principled follow-on, not attempted here |

## Outputs (`analysis/bss_bias/outputs/`)

Tracked in git (small, and they ARE the deliverable): all `.csv`, `.html`,
`figures/`, `b_draws/`. Gitignored (large): `fits/` (full stanfit objects,
only written if `SAVE_FITS <- TRUE`).
