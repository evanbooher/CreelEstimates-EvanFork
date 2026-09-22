# ==============================================================================
# bss_input_fixes.R -- make BSS inputs safe to hand to Stan
#
# Two classes of defect turn up in real creel data and are invisible until
# Stan aborts with an out-of-range index, or -- worse -- does not abort and
# silently fits the wrong thing:
#
#   * SECTION NUMBERING. prep_inputs_bss() takes DIMENSIONS as counts
#     (S = length(unique(...))) but INDICES as raw values (section_V =
#     section_num). Those agree only when sections are a dense 1..S.
#   * ANGLER TYPE. angler_final_int is computed positionally, so a stray
#     "fail" label becomes a third level, or bank/boat land on the wrong
#     numbers entirely.
#
# Both are fixed here rather than in either caller, so
# analysis/bss_bias/01_fit_bss_bias.R and template_scripts/fw_creel.Rmd apply
# exactly the same corrections. A second copy in the Rmd would drift from this
# one the first time either was touched.
#
# EVERY FUNCTION HERE IS AN IDENTITY TRANSFORM when the data is already well
# formed -- which is every fishery-year that fits today, since they fit
# precisely because it is. Adding these to a working analysis does not change
# its results.
#
# Failure handling: where a defect cannot be repaired, these call
# skip_fishery() if the caller defines one (01 does, and skips that
# fishery-year while the run continues) and otherwise stop() with the same
# message (fw_creel, where one fishery-year is the whole render).
# ==============================================================================

bss_fix_fail <- function(msg, stage = "bss_input_fixes") {
  if (exists("skip_fishery", mode = "function")) {
    skip_fishery(msg, stage = stage)
  } else {
    stop("[", stage, "] ", msg, call. = FALSE)
  }
}

# BUG FIXED 2026-09-22: align_bss_sections() (below) called skip_fishery()
# directly in two places instead of through bss_fix_fail() above -- bypassing
# the fallback this file exists to provide. Every fishery-year fit until now
# had well-formed sections, so the branch was never exercised; rendering
# fw_creel.Rmd standalone (no 01_fit_bss_bias.R sourced first, so
# skip_fishery() is undefined) turned a real, reportable data problem into
# "could not find function" instead of the intended stop() with the actual
# message. Same fix applied to drop_na_bss_inputs.R's one direct call.

# ------------------------------------------------------------------------------
# Angler-type coding
# ------------------------------------------------------------------------------
# The BSS likelihood is written for exactly two angler types in a fixed order:
# p_TI[1,]/R_V[1]/R_T[1] are bank, p_TI[2,]/R_V[2]/R_T[2] are boat. Two of the
# shared prep functions assign that index with a bare
#
#     angler_final_int = as.integer(factor(angler_final))
#
# which numbers whatever levels happen to be present, alphabetically. Two ways
# that goes wrong:
#
#   * prep_dwg_interview_angler_types() ends its case_when with TRUE ~ "fail"
#     and never filters those rows out. Any interview the rules do not classify
#     -- a missing boat_used, or an NA boat_type, since str_detect(NA, .) is NA
#     rather than FALSE -- becomes a third level, and G = 3. That is the "three
#     angler types" in the ledger. There is no third gear.
#   * prep_dwg_effort_census() computes the code BEFORE dropping "fail", so
#     bank/boat land on 1/2 whenever a fail row happens to exist. But a fishery
#     with only boat and fail rows gives boat = 1 -- boat counts handed to the
#     model as bank. Silent, and invisible to any dimension check.
#
# Recoding against fixed levels removes both. Rows that are neither bank nor
# boat are dropped: the model has nowhere to put them.
ANGLER_LEVELS <- c("bank", "boat")

recode_angler_final_int <- function(d, table_name) {
  if (!is.data.frame(d) || !"angler_final" %in% names(d)) {
    return(list(data = d, n_before = 0L, n_dropped = 0L, n_recoded = 0L, dropped_labels = ""))
  }
  n_before <- nrow(d)
  dropped_labels <- setdiff(unique(d$angler_final), ANGLER_LEVELS)
  kept <- d |> dplyr::filter(angler_final %in% ANGLER_LEVELS)
  n_dropped <- n_before - nrow(kept)

  new_int <- as.integer(factor(kept$angler_final, levels = ANGLER_LEVELS))
  n_recoded <- if ("angler_final_int" %in% names(kept)) {
    sum(kept$angler_final_int != new_int, na.rm = TRUE)
  } else 0L
  kept$angler_final_int <- new_int

  if (n_dropped > 0) {
    cli::cli_alert_warning(
      "  {.field {table_name}}: dropped {n_dropped} of {n_before} rows labelled \\
       {.val {dropped_labels}} -- neither bank nor boat."
    )
  }
  if (n_recoded > 0) {
    cli::cli_alert_warning(
      "  {.field {table_name}}: recoded angler_final_int on {n_recoded} rows -- the \\
       original numbering did not put bank at 1 and boat at 2."
    )
  }
  list(data = kept, n_before = n_before, n_dropped = n_dropped, n_recoded = n_recoded,
       dropped_labels = paste(sort(dropped_labels), collapse = "|"))
}
# ------------------------------------------------------------------------------
# Section alignment
# ------------------------------------------------------------------------------
# prep_inputs_bss() sizes the Stan arrays from the CENSUS sections --
#   S      = length(unique(effort_census$section_num))
#   O cols = any_of(paste0("open_section_", unique(effort_census$section_num)))
#   p_TI   = pivot_wider(census_expan, names_from = section_num)
# -- but indexes them with the RAW section_num carried on every observation row
# (section_V, section_T, section_E, section_IntC, ...). Three conditions must
# therefore hold, and nothing in the pipeline enforces any of them:
#
#   1. every observed section is also a census section, or the index exceeds S;
#   2. the census sections are exactly 1..S with no gaps -- Stillaguamish is
#      1-6 and 8;
#   3. census_expan covers the same sections as effort_census, in the same
#      order, or p_TI's columns are misaligned with the section they price.
#
# A water-body restriction can leave the surviving numbers with gaps too --
# Stillaguamish 2024-25 starts at section 2, and dropping the South Fork leaves
# 2-6 -- and it does nothing at all when the census counts themselves skip a
# section, or when census_expan (built from location_type == "Site" rows)
# covers a different set than effort_census. Those are this function's job, not
# the restriction's.
#
# align_bss_sections() enforces all three: it keeps the sections that have BOTH
# census effort counts and a p_census entry, drops the rest, and renumbers what
# is left to a dense 1..S. It is an IDENTITY TRANSFORM whenever the three
# conditions already hold -- which is every fishery-year that fits today, since
# they fit precisely because they hold. It also sorts effort_census by section,
# so unique() hands O's columns back in the same ascending order that
# census_expan's arrange() gives p_TI's.
#
# Renumbering is internal to the model inputs. `b` is one scalar per index count
# type, so it does not depend on a section's label; the mapping is written to
# bss_b_section_map.csv so any section-level quantity can be traced back.
align_bss_sections <- function(dwg_summ, days, fishery_name) {
  secs_of <- function(d) sort(unique(as.double(na.omit(d$section_num))))
  census_secs <- secs_of(dwg_summ$effort_census)
  expan_secs  <- secs_of(dwg_summ$census_expan)
  index_secs  <- secs_of(dwg_summ$effort_index)

  # CENSUS-FREE YEAR (2026-09-22): census_secs is legitimately empty when no
  # tie-in counts were collected at all -- not a data defect. p_census/p_TI is
  # NOT gated on a live census event: prep_dwg_census_expan() builds it from
  # dwg$effort at large (index rows included), with p_census_bank/p_census_boat
  # joined from the static fishery_manager table and defaulting to 1 (full
  # coverage) where unset -- see fw_creel.Rmd's own comment, "so this data can
  # be used in situations where census counts have not ocurred yet". V_I/T_I's
  # likelihood still needs p_TI by section regardless of whether census data
  # exists, so the right intersection here is against whichever sections
  # actually have data THIS year -- index sections when there is no census,
  # census sections otherwise, unchanged from before.
  if (length(census_secs) == 0) {
    cli::cli_alert_warning(
      "  Section alignment: no census effort counts this fishery-year -- \
       aligning to index sections against the p_census lookup instead of \
       census sections. p_census still applies (default 1 where unset); this \
       is not skipping the census-index bias correction, just its section \
       reconciliation step, which has nothing to reconcile against with zero \
       census events."
    )
    primary_secs <- index_secs
  } else {
    primary_secs <- census_secs
  }
  usable <- intersect(primary_secs, expan_secs)

  if (length(usable) == 0) {
    bss_fix_fail(
      paste0("No section has both usable effort counts and a p_census entry. ",
             "Census sections: ", paste(census_secs, collapse = ", "),
             "; index sections: ", paste(index_secs, collapse = ", "),
             "; census_expan sections: ", paste(expan_secs, collapse = ", "), "."),
      stage = "align_sections"
    )
  }

  open_cols <- grep("^open_section_", names(days), value = TRUE)
  want_open <- paste0("open_section_", as.character(usable))
  missing_open <- setdiff(want_open, open_cols)
  if (length(missing_open) > 0) {
    bss_fix_fail(
      paste0("`days` has no open/closed column for section(s) ",
             paste(sub("^open_section_", "", missing_open), collapse = ", "),
             ", which carry ", if (length(census_secs) == 0) "index" else "census",
             " counts. prep_days() was given sections: ",
             paste(sort(unique(na.omit(c(primary_secs, expan_secs)))), collapse = ", "), "."),
      stage = "align_sections"
    )
  }

  dropped <- setdiff(union(primary_secs, expan_secs), usable)
  if (length(dropped) > 0) {
    cli::cli_alert_warning(
      "  Section alignment: dropping {.val {dropped}} -- present in the \\
       {if (length(census_secs) == 0) 'index counts' else 'census counts'} or \\
       the p_census lookup, but not both."
    )
  }
  if (!identical(usable, as.double(seq_along(usable)))) {
    cli::cli_alert_info(
      "  Section alignment: renumbering {.val {usable}} -> {.val {seq_along(usable)}} \\
       so the Stan indices match S = {length(usable)}."
    )
  }

  remap <- function(d) {
    if (!is.data.frame(d) || !"section_num" %in% names(d)) return(d)
    d |>
      dplyr::filter(as.double(section_num) %in% usable) |>
      dplyr::mutate(section_num = as.integer(match(as.double(section_num), usable)))
  }
  for (nm in c("interview", "effort_index", "effort_census", "census_expan")) {
    dwg_summ[[nm]] <- remap(dwg_summ[[nm]])
  }
  # unique() on effort_census is what orders O's columns; arrange() is what
  # orders p_TI's. Sort both so the two agree.
  dwg_summ$effort_census <- dwg_summ$effort_census |> dplyr::arrange(section_num)
  dwg_summ$census_expan  <- dwg_summ$census_expan  |> dplyr::arrange(angler_final, section_num)

  keep <- setdiff(names(days), open_cols)
  open_new <- days[, want_open, drop = FALSE]
  names(open_new) <- paste0("open_section_", seq_along(usable))
  days <- dplyr::bind_cols(days[, keep, drop = FALSE], open_new)

  list(
    dwg_summ = dwg_summ,
    days     = days,
    map      = tibble(
      fishery_name    = fishery_name,
      section_num_src = usable,
      section_num_bss = seq_along(usable),
      sections_dropped = paste(dropped, collapse = "|")
    )
  )
}
# ------------------------------------------------------------------------------
# Apply everything, in the order the defects have to be resolved: recode angler
# types first (align works on the recoded tables), then sections.
#
# Returns the corrected dwg_summ and days plus the section map, so a caller can
# record what was changed instead of it happening invisibly.
# ------------------------------------------------------------------------------

apply_bss_input_fixes <- function(dwg_summ, days, fishery_name) {
  # interview and effort_census ONLY -- deliberately NOT effort_index.
  #
  # effort_index legitimately carries `total` rows: an index count of total
  # anglers, not split into bank and boat. Those are not unclassifiable
  # angler types, they are a different KIND of count, and prep_inputs_bss()
  # routes them to A_I rather than through angler_final_int. Recoding that
  # table drops about half its rows and trips the >50% guard below -- which is
  # the guard doing its job on a table that should never have been passed to
  # it. 01_fit_bss_bias.R has always recoded just these two.
  recode_log <- list()
  for (nm in c("interview", "effort_census")) {
    if (is.null(dwg_summ[[nm]])) next
    r <- recode_angler_final_int(dwg_summ[[nm]], nm)
    dwg_summ[[nm]] <- r$data
    recode_log[[nm]] <- r

    # Over half the rows unclassifiable is not a recode, it is a different
    # data problem wearing a recode's clothes.
    if (r$n_before > 0 && r$n_dropped / r$n_before > 0.5) {
      bss_fix_fail(
        paste0(nm, ": ", r$n_dropped, " of ", r$n_before,
               " rows are neither bank nor boat (", r$dropped_labels,
               "). Check prep_dwg_interview_angler_types() for this fishery."),
        stage = "recode_angler"
      )
    }
  }

  aligned <- align_bss_sections(dwg_summ, days, fishery_name)
  aligned$recode_log <- recode_log
  aligned
}

# ------------------------------------------------------------------------------
# Preflight on the assembled Stan input list. Catches what survived the fixes
# above -- and would otherwise surface as "index N out of range" from inside
# the sampler, with nothing naming the fishery or the cause.
# ------------------------------------------------------------------------------

preflight_bss_inputs <- function(inputs_bss, fishery_name) {
  # Declared-count vs. actual-length checks. Stan does not error when these
  # disagree until sampling starts, and then only on the FIRST chain to touch
  # the mismatched variable -- reported as "mismatch in dimension declared and
  # found", with every chain otherwise returning silently with no draws and no
  # R-level error. Catching it here names the variable and the fishery instead.
  dim_checks <- list(
    V_n   = c("day_V", "section_V", "countnum_V", "V_I"),
    T_n   = c("day_T", "section_T", "countnum_T", "T_I"),
    A_n   = c("day_A", "gear_A", "section_A", "countnum_A", "A_I"),
    B_n   = c("day_B", "gear_B", "section_B", "countnum_B", "B_s"),
    E_n   = c("day_E", "gear_E", "section_E", "countnum_E", "E_s"),
    IntC  = c("day_IntC", "gear_IntC", "section_IntC", "c", "h"),
    IntA  = c("day_IntA", "gear_IntA", "section_IntA", "V_A", "T_A", "B_A", "A_A")
  )
  for (dim_name in names(dim_checks)) {
    declared <- inputs_bss[[dim_name]]
    if (is.null(declared)) next
    for (vec_name in dim_checks[[dim_name]]) {
      actual <- length(inputs_bss[[vec_name]])
      if (actual != declared) {
        bss_fix_fail(
          paste0(dim_name, " = ", declared, " but ", vec_name, " has length ", actual,
                 ". A catch group or interview table upstream is duplicating or ",
                 "dropping rows before this fishery's Stan data was assembled."),
          stage = "bss_preflight"
        )
      }
    }
  }

  sec_idx <- unlist(inputs_bss[c("section_V", "section_T", "section_A", "section_E",
                                 "section_IntC", "section_IntA")], use.names = FALSE)
  sec_idx <- sec_idx[!is.na(sec_idx)]
  gear_idx <- unlist(inputs_bss[c("gear_A", "gear_E", "gear_IntC", "gear_IntA")],
                     use.names = FALSE)
  gear_idx <- gear_idx[!is.na(gear_idx)]

  if (length(sec_idx) > 0 && max(sec_idx) > inputs_bss$S) {
    bss_fix_fail(
      paste0("Section index ", max(sec_idx), " exceeds S = ", inputs_bss$S,
             " (sections present: ", paste(sort(unique(sec_idx)), collapse = ", "),
             "). Stan indexes p_TI/O positionally, so this reads out of range."),
      stage = "bss_preflight"
    )
  }
  if (length(gear_idx) > 0 && max(gear_idx) > inputs_bss$G) {
    bss_fix_fail(
      paste0("Angler-type index ", max(gear_idx), " exceeds G = ", inputs_bss$G, "."),
      stage = "bss_preflight"
    )
  }
  if (inputs_bss$G > 2) {
    bss_fix_fail(
      paste0("G = ", inputs_bss$G, " angler types, but the BSS likelihood hard-codes ",
             "exactly two. recode_angler_final_int() should have made this impossible."),
      stage = "bss_preflight"
    )
  }
  p_dim <- dim(inputs_bss$p_TI)
  if (!is.null(p_dim) && !identical(as.integer(p_dim), c(inputs_bss$G, inputs_bss$S))) {
    bss_fix_fail(
      paste0("p_TI is ", p_dim[1], "x", p_dim[2], " but G x S is ",
             inputs_bss$G, "x", inputs_bss$S, "."),
      stage = "bss_preflight"
    )
  }
  invisible(TRUE)
}
