# ==============================================================================
# render_briefs.R -- one per-basin collaborator packet, rendered from one source
#
# Part 1 (11_share_brief.qmd) and part 2 (12_estimates_update.qmd) are the same
# document for every basin; only the `focus` param changes. Rendering them from
# a loop rather than by hand is what keeps the Stillaguamish and Snohomish
# copies from quietly drifting apart as the text is edited.
#
# Part 1 needs only 07_catch_sensitivity.R to have run. Part 2 additionally
# needs the production fits (09_read_production_estimates.R, then 07 again), so
# it is skipped -- not failed -- while the sweep is still going.
#
#   Rscript analysis/bss_bias/report/render_briefs.R
#
# Narrow or widen the set first if wanted:
#   BRIEF_FOCI <- "Stillaguamish"; source(".../render_briefs.R")
# ==============================================================================

suppressPackageStartupMessages({
  library(cli)
  library(here)
})

# Works both as `Rscript .../render_briefs.R` (--file= is present) and as
# `source(...)` from the Console (it is not, so fall back to here()).
.file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
REPORT_DIR <- if (length(.file_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", .file_arg), mustWork = FALSE))
} else {
  here::here("analysis", "bss_bias", "report")
}
if (!dir.exists(REPORT_DIR)) {
  REPORT_DIR <- here::here("analysis", "bss_bias", "report")
}
OUT_DIR <- file.path(dirname(REPORT_DIR), "outputs")

# `focus` is matched against fishery_type with str_detect, so a basin name is
# enough. Single words only -- the value crosses a command line, and a label
# with spaces is a quoting hazard on Windows. The document derives its own
# display label from the data it kept.
if (!exists("BRIEF_FOCI", inherits = FALSE)) {
  BRIEF_FOCI <- c("Stillaguamish", "Snohomish")
}

quarto_bin <- Sys.which("quarto")
if (!nzchar(quarto_bin)) {
  cli_abort("{.code quarto} is not on PATH. Install Quarto, or render by hand:
             {.code quarto render 11_share_brief.qmd -P focus:Stillaguamish}")
}

have_part2 <- file.exists(file.path(OUT_DIR, "bss_b_T8_gear_split.csv"))
if (!have_part2) {
  cli_alert_info(
    "No {.file bss_b_T8_gear_split.csv} -- rendering part 1 only. Part 2 needs the \\
     production fits: run {.file 09_read_production_estimates.R}, then \\
     {.file 07_catch_sensitivity.R}, then re-run this."
  )
}

jobs <- expand.grid(
  focus = BRIEF_FOCI,
  qmd   = c("11_share_brief.qmd", if (have_part2) "12_estimates_update.qmd"),
  stringsAsFactors = FALSE
)

failures <- character(0)

for (i in seq_len(nrow(jobs))) {
  focus <- jobs$focus[i]
  qmd   <- jobs$qmd[i]
  out   <- sub("\\.qmd$", paste0("_", tolower(focus), ".html"), qmd)

  cli_h2("{qmd} -- {focus}")
  status <- system2(
    quarto_bin,
    c("render", shQuote(file.path(REPORT_DIR, qmd)),
      "-P", paste0("focus:", focus),
      "--output", shQuote(out))
  )

  if (identical(status, 0L)) {
    cli_alert_success("{.file {file.path(REPORT_DIR, out)}}")
  } else {
    cli_alert_danger("{qmd} failed for {focus} (exit {status}).")
    failures <- c(failures, paste0(qmd, " [", focus, "]"))
  }
}

cli_rule()
if (length(failures) == 0) {
  cli_alert_success("{nrow(jobs)} brief{?s} rendered.")
} else {
  cli_alert_danger("{length(failures)} of {nrow(jobs)} failed: {.val {failures}}")
}
