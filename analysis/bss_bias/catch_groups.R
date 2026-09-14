# ==============================================================================
# catch_groups.R -- the named catch groups, shared by 00d and 01b
#
# Sourced by 00d_catch_inventory.R (what is in the data) and
# 01b_fit_catch_groups.R (what gets fitted), so the inventory always describes
# exactly the groups the fits will use. Previously these lived only in 01b,
# which made the inventory a separate transcription waiting to drift.
#
# Each field is a str_detect PATTERN, matched against the catch table's
# species / life_stage / fin_mark / fate columns. Alternation works, and "NA"
# matches the literal string that prep_dwg_interview_catch() coerces NA to.
#
# NOTE these are SUBSTRING matches, not exact. "AD" matches any fin_mark
# containing "AD"; that is the same behaviour prep_dwg_interview_catch() has
# always had, and 00d prints the distinct raw values it matched so an
# unintended catch is visible rather than assumed away.
# ==============================================================================

CATCH_GROUPS <- list(
  # Any encounter of Chinook -- adults and jacks, every mark status including
  # unknown and unrecorded, kept and released.
  #
  # DO NOT DROP THIS GROUP ON SAMPLE SIZE, IN ANY BASIN. 00d shows single-digit
  # Chinook encounters in most fishery-years and zero in three, which makes this
  # look like a group not worth fitting. That reasoning is a category error.
  #
  # Chinook encounters are an IMPACT-LIMITED quantity: the management question
  # is whether an incidental-impact estimate exists and what its upper bound is,
  # against a threshold that decides whether a fishery opens or continues. That
  # is not a precision problem, and a small n is not a reason to skip it -- in
  # Snohomish and Stillaguamish alike, the low-count years are exactly the ones
  # the decision turns on. A coho series with 400 fish is statistically
  # comfortable and decides much less.
  #
  # Two consequences for how these get reported:
  #   - Always with an interval, never as a point. A C_sum built on one
  #     encounter has a very wide posterior and the UPPER bound is the
  #     decision-relevant end. A median quoted alone reads as a precise small
  #     number when it means "possibly near zero, possibly a good deal more".
  #   - A fishery-year with zero encounters is a real and reportable result,
  #     not a missing one. 01 skips it at the "no matching records" guard, so
  #     the absence must be carried over from 00d by hand.
  chinook_all = list(
    species    = "Chinook",
    life_stage = "Adult|Jack",
    fin_mark   = "UM|AD|UNK|NA",
    fate       = "Released|Kept"
  ),
  # Coho harvest. NOT the pipeline default (Coho_Adult_AD|UM_Kept) -- this adds
  # jacks, so it is a genuinely different group.
  coho_harvest = list(
    species    = "Coho",
    life_stage = "Adult|Jack",
    fin_mark   = "UM|AD",
    fate       = "Kept"
  )
)

# The est_cg string 01/01b build for a group, and the key 07 joins on.
catch_group_label <- function(g) paste0(unlist(g), collapse = "_")

# Rows of a catch table matching one group. Kept here so the inventory counts
# exactly what a fit would count.
match_catch_group <- function(catch_df, g) {
  catch_df |>
    dplyr::mutate(dplyr::across(c(species, life_stage, fin_mark, fate),
                                ~tidyr::replace_na(as.character(.), "NA"))) |>
    dplyr::filter(
      stringr::str_detect(species,    g$species),
      stringr::str_detect(life_stage, g$life_stage),
      stringr::str_detect(fin_mark,   g$fin_mark),
      stringr::str_detect(fate,       g$fate)
    )
}

# The catch groups as the data.frame that fw_creel.Rmd's `est_catch_groups`
# param expects -- one row per group, columns species / life_stage / fin_mark /
# fate. Passing several rows in ONE render is deliberate: fw_creel builds an
# inputs_bss entry per unique est_cg and loops the fit over all of them, so a
# single render produces every group for that fishery-year.
catch_groups_df <- function(keys = names(CATCH_GROUPS)) {
  bad <- setdiff(keys, names(CATCH_GROUPS))
  if (length(bad) > 0) stop("Unknown catch group(s): ", paste(bad, collapse = ", "))
  do.call(rbind, lapply(keys, function(k) {
    g <- CATCH_GROUPS[[k]]
    data.frame(species = g$species, life_stage = g$life_stage,
               fin_mark = g$fin_mark, fate = g$fate, stringsAsFactors = FALSE)
  }))
}
