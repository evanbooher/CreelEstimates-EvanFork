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
#
# The returned list carries a `na_drop_log` attribute -- a tibble, one row
# per group present in this fishery-year's inputs (n_dropped may be 0) --
# so callers can write a persistent audit record rather than relying on the
# console warnings alone. See its use in 01_fit_bss_bias.R.
drop_na_bss_inputs <- function(inputs_bss, fishery_name = NA_character_) {

  groups <- list(
    list(n = "V_n",  vars = c("day_V", "section_V", "countnum_V", "V_I")),
    list(n = "T_n",  vars = c("day_T", "section_T", "countnum_T", "T_I")),
    list(n = "A_n",  vars = c("day_A", "gear_A", "section_A", "countnum_A", "A_I")),
    list(n = "E_n",  vars = c("day_E", "gear_E", "section_E", "countnum_E", "E_s")),
    list(n = "IntC", vars = c("day_IntC", "gear_IntC", "section_IntC", "c", "h")),
    list(n = "IntA", vars = c("day_IntA", "gear_IntA", "section_IntA", "V_A", "T_A", "A_A"))
  )

  drop_log <- tibble::tibble(
    fishery_name = character(), group = character(), vars = character(),
    n_before = integer(), n_dropped = integer(), n_after = integer()
  )

  for (g in groups) {
    present_vars <- intersect(g$vars, names(inputs_bss))
    if (length(present_vars) == 0) next

    n_before <- inputs_bss[[g$n]]
    na_mask <- Reduce(`|`, lapply(present_vars, function(v) is.na(inputs_bss[[v]])))
    n_bad <- sum(na_mask)

    if (n_bad > 0) {
      cli::cli_alert_warning(
        "  {.val {fishery_name}}: dropping {n_bad} NA row(s) from {.val {g$n}} \\
         ({paste(present_vars, collapse=', ')}) before Stan -- {n_before} -> {n_before - n_bad}."
      )
      for (v in present_vars) inputs_bss[[v]] <- inputs_bss[[v]][!na_mask]
      inputs_bss[[g$n]] <- inputs_bss[[g$n]] - n_bad
    }

    drop_log <- dplyr::bind_rows(drop_log, tibble::tibble(
      fishery_name = fishery_name, group = g$n, vars = paste(present_vars, collapse = "|"),
      n_before = n_before, n_dropped = n_bad, n_after = n_before - n_bad
    ))
  }

  attr(inputs_bss, "na_drop_log") <- drop_log

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
