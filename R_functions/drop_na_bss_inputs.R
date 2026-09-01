# Drop NA-containing observation rows from a prep_inputs_bss() output list
# before handing it to fit_bss()/stan(), which hard-errors on any NA
# anywhere in `data` (e.g. "Stan does not support NA (in countnum_E) in
# data"). For expediency, not as a substitute for fixing the raw data.
#
# Only touches the six "one row per observation" groups in the Stan data
# block -- vehicle/trailer/angler index counts (V/T/A), census tie-in counts
# (E), and the two interview tables (IntC/IntA). Each observation in these
# groups stands alone: nothing else in the model indexes INTO them by
# position, so dropping a row and decrementing its `_n` count is safe.
#
# Deliberately does NOT touch the day-indexed arrays (w, period, L) or the
# O (open/closed) matrix: every group above references THOSE by position
# via day_V/day_T/.../day_IntA, so silently dropping a day would corrupt
# every other group's day index. NA there is a different, deeper data
# problem -- this skips the fishery-year loudly instead of masking it.
drop_na_bss_inputs <- function(inputs_bss, fishery_name = NA_character_) {

  groups <- list(
    list(n = "V_n",  vars = c("day_V", "section_V", "countnum_V", "V_I")),
    list(n = "T_n",  vars = c("day_T", "section_T", "countnum_T", "T_I")),
    list(n = "A_n",  vars = c("day_A", "gear_A", "section_A", "countnum_A", "A_I")),
    list(n = "E_n",  vars = c("day_E", "gear_E", "section_E", "countnum_E", "E_s")),
    list(n = "IntC", vars = c("day_IntC", "gear_IntC", "section_IntC", "c", "h")),
    list(n = "IntA", vars = c("day_IntA", "gear_IntA", "section_IntA", "V_A", "T_A", "A_A"))
  )

  for (g in groups) {
    present_vars <- intersect(g$vars, names(inputs_bss))
    if (length(present_vars) == 0) next

    na_mask <- Reduce(`|`, lapply(present_vars, function(v) is.na(inputs_bss[[v]])))
    n_bad <- sum(na_mask)
    if (n_bad == 0) next

    cli::cli_alert_warning(
      "  {.val {fishery_name}}: dropping {n_bad} NA row(s) from {.val {g$n}} \\
       ({paste(present_vars, collapse=', ')}) before Stan -- {inputs_bss[[g$n]]} -> \\
       {inputs_bss[[g$n]] - n_bad}."
    )
    for (v in present_vars) inputs_bss[[v]] <- inputs_bss[[v]][!na_mask]
    inputs_bss[[g$n]] <- inputs_bss[[g$n]] - n_bad
  }

  day_level <- intersect(c("w", "period", "L"), names(inputs_bss))
  bad_day <- day_level[vapply(day_level, function(v) any(is.na(inputs_bss[[v]])), logical(1))]
  if (length(bad_day) > 0 || (!is.null(inputs_bss$O) && any(is.na(inputs_bss$O)))) {
    skip_fishery(
      paste0("NA found in day-indexed/status input(s) (",
             paste(c(bad_day, if (any(is.na(inputs_bss$O))) "O"), collapse = ", "),
             ") -- these are referenced by position from every count/interview group and ",
             "can't be safely dropped. Data problem upstream of Stan (prep_days()/",
             "prep_inputs_bss()), not something to paper over here."),
      stage = "drop_na_bss_inputs"
    )
  }

  inputs_bss
}
