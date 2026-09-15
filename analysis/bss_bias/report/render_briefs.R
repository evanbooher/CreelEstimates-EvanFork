# ==============================================================================
# render_briefs.R -- one per-basin collaborator brief, rendered from one source
#
# 11_share_brief.qmd is the same document for every basin; only the `focus`
# param changes. Rendering from a loop rather than by hand is what keeps the
# Stillaguamish and Snohomish copies from quietly drifting apart as the text is
# edited.
#
# The brief needs 07_catch_sensitivity.R to have run. It GROWS an effort-and-
# catch section once the production fits exist (09_read_production_estimates.R,
# then 07 again) -- the document gates on the data itself, so this script does
# not have to know or care which version it is producing.
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
# Absolute, because the render loop setwd()s into REPORT_DIR and a relative
# OUT_DIR would then point somewhere else.
REPORT_DIR <- normalizePath(REPORT_DIR, winslash = "/", mustWork = TRUE)
OUT_DIR    <- file.path(dirname(REPORT_DIR), "outputs")

# `focus` is matched against fishery_type with str_detect, so a basin name is
# enough. Single words only -- the value crosses a command line, and a label
# with spaces is a quoting hazard on Windows. The document derives its own
# display label from the data it kept.
if (!exists("BRIEF_FOCI", inherits = FALSE)) {
  BRIEF_FOCI <- c("Stillaguamish", "Snohomish")
}

quarto_bin <- Sys.which("quarto")
if (!nzchar(quarto_bin)) {
  cli_abort("{.code quarto} is not on PATH. Install Quarto, or render by hand \
             FROM {.file {REPORT_DIR}}:
             {.code quarto render 11_share_brief.qmd -P focus:Stillaguamish}")
}

# Reported, not acted on: the document decides for itself. This line just says
# which of the two versions is about to come out, so a brief that quietly lacks
# the fish section is not a surprise.
if (file.exists(file.path(OUT_DIR, "bss_b_T8_gear_split.csv"))) {
  cli_alert_info("Gear split found -- briefs will include effort and catch in fish.")
} else {
  cli_alert_info(
    "No {.file bss_b_T8_gear_split.csv} -- briefs will end at the percentages, with \\
     the follow-up note. For the in-fish section run \\
     {.file 09_read_production_estimates.R}, then {.file 07_catch_sensitivity.R}, \\
     then re-run this."
  )
}

jobs <- expand.grid(
  focus = BRIEF_FOCI,
  qmd   = "11_share_brief.qmd",
  stringsAsFactors = FALSE
)

# RUN FROM THE REPORT DIRECTORY, always.
#
# Quarto resolves --output relative to the INPUT file's directory, but writes
# the _files support directory relative to the invocation directory. Call it
# from the project root with a path into report/ and those two disagree: pandoc
# gets `--output ..\..\..\11_share_brief_stillaguamish.html`, the support
# directory lands somewhere else, and embed-resources then dies looking for
# quarto-html/quarto.js. Passing bare basenames from inside report/ keeps both
# in the same place.
#
# The documents use here::here() throughout, which walks up to the project root
# from wherever it starts, so nothing inside them depends on the cwd.
#
# In a function, so on.exit() actually restores the directory -- at top level it
# would not, and a failed run sourced from the Console would leave the session
# sitting in report/.
render_all <- function() {
  owd <- setwd(REPORT_DIR)
  on.exit(setwd(owd), add = TRUE)
  failures <- character(0)

  for (i in seq_len(nrow(jobs))) {
    focus <- jobs$focus[i]
    qmd   <- jobs$qmd[i]
    out   <- sub("\\.qmd$", paste0("_", tolower(focus), ".html"), qmd)

    cli_h2("{qmd} -- {focus}")

    # No --output. Quarto's default name is the input's, which lands
    # unambiguously next to the qmd; renaming afterwards avoids the output-path
    # resolution entirely rather than trying to get it right on two platforms.
    default_out <- sub("\\.qmd$", ".html", qmd)
    status <- system2(
      quarto_bin,
      c("render", qmd, "-P", paste0("focus:", focus))
    )

    if (!identical(status, 0L)) {
      cli_alert_danger("{qmd} failed for {focus} (exit {status}).")
      failures <- c(failures, paste0(qmd, " [", focus, "]"))
    } else if (!file.exists(default_out)) {
      cli_alert_danger("{qmd} reported success for {focus} but {.file {default_out}} is not there.")
      failures <- c(failures, paste0(qmd, " [", focus, "]"))
    } else {
      # Renders are sequential, so the default name is free each time.
      if (file.exists(out)) file.remove(out)
      if (file.rename(default_out, out)) {
        cli_alert_success("{.file {file.path(REPORT_DIR, out)}}")
      } else {
        cli_alert_danger("Rendered {focus} but could not rename to {.file {out}}.")
        failures <- c(failures, paste0(qmd, " [", focus, "]"))
      }
    }
  }

  failures
}

failures <- render_all()

cli_rule()
if (length(failures) == 0) {
  cli_alert_success("{nrow(jobs)} brief{?s} rendered.")
} else {
  cli_alert_danger("{length(failures)} of {nrow(jobs)} failed: {.val {failures}}")
}
