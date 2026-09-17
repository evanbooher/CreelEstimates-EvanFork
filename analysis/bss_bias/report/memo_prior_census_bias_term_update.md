# Using prior census data to set the effort-index bias term

Draft — 2026-09-16 — Evan Booher, WDFW Freshwater Fisheries Policy

*Technical basis for the effort bias-correction approach described in the 2026 Snohomish and Stillaguamish basin monitoring briefs.*

## Purpose and recommendation

Following the elimination of District 13 freshwater salmon harvest monitoring funding (MI 5M413, a $278,007 annual reduction), aerial census counts are not funded for the Snohomish mainstem, Skykomish, Wallace, and Snoqualmie fisheries or for the Stillaguamish mainstem and North Fork Stillaguamish gamefish fishery in 2026. Where census (tie-in) counts are reduced or eliminated this season, WDFW proposes to carry the effort-index bias term forward from prior years as an informative prior rather than leave it unconstrained.

The bias term, written b, converts index counts of vehicles and boat trailers at access sites into anglers. Measuring it requires census counts run alongside index counts on the same day and in the same section. Without either census counts or a prior, the model has nothing to pin b against, and the effort and catch estimates inherit that.

The recommendation has three parts:

1. Set an informative prior on the vehicle bias term from the fishery's own prior years, wide enough to span the year-to-year variation actually observed rather than the precision of any single year's estimate.
2. Fix the trailer bias term at 1 in fisheries with negligible boat effort, rather than estimate it from data that cannot support it.
3. Direct the tie-in effort that is available to the reaches and gear types where it changes the estimate most.

This refines, and should replace, the "multi-year average" description of the effort bias correction currently in Table 1 of the Snohomish and Stillaguamish monitoring briefs. An informative prior, combined with whatever census data is collected this season, is a materially different and more defensible statement than a flat historical average, and the briefs' language should be updated to match before they go to co-managers.

This commits WDFW to reporting, for each fishery-year, whether b was measured or imported, and to publishing the prior used alongside the estimate.

## Background: what the bias term does

Creel effort is estimated from index counts, which are counts of objects at a standard set of access sites rather than counts of anglers. The bias term b is the correction between what those counts see and how many anglers are actually present. There are two: one on the vehicle index and one on the trailer index.

Census counts are what measure it. A census count totals a whole section, and is scheduled against an index count on the same day, so the pair gives a direct reading of the ratio. A census day without a matching index count, or the reverse, contributes nothing to b.

One consequence governs everything below. In a year without census counts, estimated catch is proportional to 1/b: double b and the catch estimate is cut in half, halve it and the estimate doubles. Estimated effort moves by exactly the same factor, so the catch rate is unchanged. The bias term rescales the size of the fishery without changing the rate at which fish are caught.

That also sets what an error in b costs. A b that is 10 percent too high pulls estimated catch and effort 10 percent too low, for every species and catch group in that fishery, in the same direction and by the same proportion.

## Evidence: how close a borrowed bias term gets

The vehicle term carries over between years well enough to import. The trailer term does not.

Each fishery-year with a measured b was predicted from the other years of its own series and nothing else, then compared against what was actually measured. Leave-one-out, so no year is scored against an average containing itself. Forty-three such tests across Skagit, Snohomish and Stillaguamish.

| Index channel | Inside the 95 percent interval | Typical error | Worst error |
| --- | --- | --- | --- |
| Vehicle | 20 of 21 | 3.7 to 18.2 percent | 36 percent |
| Trailer | 17 of 22 | 9.5 to 44.5 percent | 168 percent |

Because catch and effort track b one-for-one, those are also the errors carried straight into the estimate.

Two qualifications belong with these figures. Every historical year in the test had census counts, so the comparison measures how well the prediction reproduces a measured value, while the effect on catch assumes no census at all; where some census is collected this season, it will pull the estimate back toward the data and the error will be smaller. And the fisheries in the series differ in reach, season length and access, which is the assumption the prediction interval rests on and the reason it is wide.

## Evidence: Stillaguamish mainstem and North Fork, 2024 and 2025

TO BE COMPLETED when the current fits finish. Values below are placeholders.

The two most recent Stillaguamish seasons were compared on a common window, 16 September to 31 October, restricted to the water each fork actually covers. The two years match closely enough that a difference in b can be read as a difference rather than an artifact of coverage: the same census blocks, the same index sites, and three paired census days in each.

| Fit | Sections | Index days | Paired census days | Census anglers, bank | Census anglers, boat |
| --- | --- | --- | --- | --- | --- |
| Mainstem 2024 | 2, 3 | 24 | 3 | 145 | 0 |
| Mainstem 2025 | 1, 2 | 24 | 3 | 227 | 15 |
| North Fork 2024 | 4 | 25 | 3 | 7 | 2 |
| North Fork 2025 | 3 | 28 | 3 | 12 | 0 |

The survey totals are the constraint, and they are known independently of any model fit. Each estimate rests on three days that carry both a census and an index count. The North Fork rests on 7 and 12 anglers counted in census across the whole window, and no fit can add information the survey did not collect. Boat anglers in census are what the trailer term is measured against, and there are 0 to 15 of them in each fishery-year.

| Fit | Vehicle b | 95 percent interval | Data limitation |
| --- | --- | --- | --- |
| Mainstem 2024 | TBD | TBD | TBD |
| Mainstem 2025 | TBD | TBD | TBD |
| North Fork 2024 | TBD | TBD | TBD |
| North Fork 2025 | TBD | TBD | TBD |

Trailer estimates are not reported for any of the four. With fewer than 20 boat anglers counted in census, there is nothing for that term to be measured against, whatever the trailer counts show.

The North Fork Stillaguamish gamefish fishery carries a standing co-manager commitment, under LOAF section 1.16, to monitor through November 30 whenever salmon fishing is open in that river system, unlike the discretionary Snohomish tributary seasons. That commitment is what makes the North Fork's thin census base a live issue rather than a background limitation: it is the fishery where reduced monitoring is least optional and the bias term is least well measured.

## Proposed approach for this season

Vehicle term. Place a lognormal prior on b, centred on the pooled estimate from that fishery's own prior years and given a spread taken from the leave-one-out prediction interval, not from the posterior width of any single year. The distinction matters: a single year's posterior can be narrow while the year-to-year variation is wide, and it is the year-to-year variation that describes what an unmeasured year might be.

Trailer term. Fix at 1 where census boat-angler coverage falls below the threshold used in the comparison above. This is a stated assumption, not an estimate, and it is reported as such. In a fishery that is 2 to 3 percent boat effort, the choice moves the season total very little, and the alternative is a number with no data behind it.

Reaches without their own census. Apply the fishery's bias term rather than fitting a reach-specific one, and state that the reach estimate inherits a correction measured elsewhere. Do not set a reach-specific prior from a reach-specific fit that the survey could not support; that would import a number and read it back as a measurement.

Where census is collected this season. The prior and the census data combine in the usual way, and the census will dominate wherever it is adequate. The prior is what keeps the estimate stable where it is not.

Specific values follow once the current fits complete.

## Where tie-in effort should go

A census count only informs b if an index count was taken in the same section on the same day. Scheduling census days against index days is therefore worth more than adding census days, and an unpaired census count is effort that does not reach the bias term at all.

Three points follow for this season's schedule:

1. Pair every census day with an index count in the same section. The Stillaguamish comparison above turned on three paired days per fit, and the count of paired days, not the count of census surveys, is what sets how well b is known.
2. Spread paired days across the season rather than clustering them. A single scalar is fitted per fishery, so paired days at one point in the season describe that point and are assumed to hold for the rest.
3. Cover boat anglers in census wherever the trailer term is meant to be estimated. Trailer counts alone cannot do it. A reach with boat effort and no boat-angler census will return a trailer term that reflects its prior, not the fishery.

Where census capacity is limited, concentrating it on one reach with adequate pairing is more useful than spreading it thinly across several, because a reach without its own b can borrow one while a reach with a poorly measured b cannot be distinguished from one that is genuinely different.

This is in tension with the front-loaded monitoring priority proposed in the Snohomish brief, which concentrates survey days in the window of highest Chinook catchability (the 9/26 opener through approximately mid-October) and proposes tapering after an October co-manager check-in. If tapering means paired census stops after that point, any bias term applied to the tapered period rests entirely on pairs drawn from the high-catchability window and carried into a period that may have different angler behavior and access conditions. Whether that gap matters depends on how much effort composition actually changes between the two windows, and it is worth raising explicitly with co-managers alongside the tapering decision, rather than leaving it as an implicit assumption of the front-loaded schedule.

## Limitations and open items

What the approach assumes. That a fishery's bias term is exchangeable across years, so an unmeasured year can be predicted from measured ones. The leave-one-out test is the check on that assumption rather than an argument for it, and it passes for the vehicle term and fails for the trailer term. It also assumes the fishery is the same fishery across those years; a series whose reach, access sites or season length changed materially is not comparable, and the Stillaguamish comparison above was restricted to a common window and common sections for that reason.

What would change the recommendation. A vehicle term that moves substantially between two adjacent, near-identically surveyed years would argue for a wider prior than the pooled series implies, or against importing a point estimate at all. A fishery with adequate paired census this season needs no prior and should not be given a tight one.

Decisions needed:

- Confirm which fisheries and reaches will carry reduced census coverage this season, since the approach only applies where that is true.
- Confirm the fixed value for the trailer term where boat effort is negligible, and whether that is acceptable to co-managers as a stated assumption.
- Confirm how imported bias terms are disclosed in the estimate tables and in any pre-season or post-season reporting.
- Reconcile this document's bias-term proposal (informative prior, not a flat multi-year average) with the "Effort bias correction" row in Table 1 of the Snohomish and Stillaguamish monitoring briefs before those briefs are finalized.
- Decide whether paired census days should be spread across the full season for bias-term stability, or concentrated in the front-loaded high-catchability window as currently proposed, and whether co-managers need to weigh in on that tradeoff directly.