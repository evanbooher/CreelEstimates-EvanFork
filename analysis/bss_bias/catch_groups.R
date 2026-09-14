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
  # DO NOT DROP THIS GROUP ON SAMPLE SIZE. 00d shows single-digit Chinook
  # encounters in most fishery-years and zero in three, which looks like a
  # group not worth fitting. It is not: in Stillaguamish the low-count years
  # are the ones that decide whether the fishery continues at all. The question
  # there is whether an incidental-impact estimate EXISTS and what its upper
  # bound is -- a management threshold, not a precision problem. An estimate
  # built on one encounter is wide and still decision-relevant; report it with
  # its interval, never as a point.
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
