# ==============================================================================
# make_reference_docx.R
#
# Builds the Word reference document that the census-bias memos
# (memo_prior_census_bias_term_stillaguamish.qmd,
# memo_prior_census_bias_term_snohomish.qmd) render against: Calibri
# throughout, black text, one body size, headings distinguished by weight
# rather than colour or size jumps. Both memos point reference-doc at the
# same file, so restyling once restyles both.
#
# Run this ONCE and commit the result. Quarto's stock reference doc puts
# headings in blue Calibri Light at four different sizes, which is why the
# rendered memo arrives looking like three documents stapled together.
#
# Why a script and not Word: restyling by hand is not reproducible, and the
# next person to touch the memo would have no way to tell which of the two
# dozen pandoc styles had been changed. A .docx is a zip of XML -- the styles
# are editable directly, so the whole restyle is the twenty lines below.
#
#   Rscript analysis/bss_bias/report/make_reference_docx.R
#
# Requires quarto (or pandoc) on PATH to dump the stock file, and the xml2 and
# zip packages to rewrite it.
# ==============================================================================

suppressPackageStartupMessages({
  library(xml2)
})

stopifnot_pkg <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("Package '", p, "' is required: install.packages(\"", p, "\")", call. = FALSE)
  }
}
stopifnot_pkg("xml2")
stopifnot_pkg("zip")

REPORT_DIR <- here::here("analysis", "bss_bias", "report")
OUT        <- file.path(REPORT_DIR, "reference-plain.docx")

FONT      <- "Calibri"
BODY_HALF <- 22L   # half-points, so 11 pt

# ONE SIZE, title included. Headings and the title differ from body text by
# weight alone. Leaving the title at 14 pt was still two sizes on the page.
HEADING_HALF <- c(Title = BODY_HALF, Heading1 = BODY_HALF,
                  Heading2 = BODY_HALF, Heading3 = BODY_HALF)

# ------------------------------------------------------------------------------
# 1. Dump the stock reference doc
# ------------------------------------------------------------------------------
dump_stock <- function(dest) {
  for (cmd in list(c("quarto", "pandoc"), "pandoc")) {
    exe <- cmd[1]
    if (nzchar(Sys.which(exe))) {
      args <- c(cmd[-1], "-o", dest, "--print-default-data-file", "reference.docx")
      status <- suppressWarnings(
        system2(exe, shQuote(args), stdout = FALSE, stderr = FALSE)
      )
      # A .docx is a zip, so it starts "PK". Anything else means the command
      # printed a diagnostic where the file should be.
      ok <- status == 0 && file.exists(dest) && file.size(dest) > 1000 &&
        identical(readBin(dest, "raw", 2L), as.raw(c(0x50, 0x4b)))
      if (ok) return(invisible(exe))
    }
  }
  stop("Neither quarto nor pandoc is on PATH -- cannot dump the stock ",
       "reference.docx. Open a terminal where `quarto --version` works.",
       call. = FALSE)
}

work <- file.path(tempdir(), "refdocx")
unlink(work, recursive = TRUE)
dir.create(work, recursive = TRUE)

stock <- file.path(work, "stock.docx")
via <- dump_stock(stock)
message("Stock reference doc dumped via ", via, ".")

# ------------------------------------------------------------------------------
# 2. Rewrite word/styles.xml
# ------------------------------------------------------------------------------
xdir <- file.path(work, "x")
utils::unzip(stock, exdir = xdir)

styles_path <- file.path(xdir, "word", "styles.xml")
doc <- read_xml(styles_path)
ns  <- xml_ns(doc)

# Font. The theme attributes outrank the explicit ones where both are set, so
# setting w:ascii alone changes nothing -- the headings kept coming out in the
# theme's major font. The theme attributes are stripped after serialisation
# (see below) because xml2 will not remove a namespace-prefixed attribute:
# `xml_attr(n, "w:asciiTheme") <- NULL` is a silent no-op.
for (n in xml_find_all(doc, "//w:rFonts", ns)) {
  xml_attr(n, "w:ascii")    <- FONT
  xml_attr(n, "w:hAnsi")    <- FONT
  xml_attr(n, "w:cs")       <- FONT
  xml_attr(n, "w:eastAsia") <- FONT
}

# Colour. Black everywhere, including the theme-coloured headings.
for (n in xml_find_all(doc, "//w:color", ns)) {
  for (a in c("w:themeColor", "w:themeTint", "w:themeShade")) xml_attr(n, a) <- NULL
  xml_attr(n, "w:val") <- "000000"
}

# Size. Flatten everything to the body size first, then raise only the title.
for (n in xml_find_all(doc, "//w:sz | //w:szCs", ns)) {
  xml_attr(n, "w:val") <- as.character(BODY_HALF)
}

set_style_size <- function(style_id, half) {
  st <- xml_find_first(doc, sprintf("//w:style[@w:styleId='%s']", style_id), ns)
  if (inherits(st, "xml_missing")) return(invisible(FALSE))
  rpr <- xml_find_first(st, "./w:rPr", ns)
  if (inherits(rpr, "xml_missing")) {
    xml_add_child(st, "w:rPr")
    rpr <- xml_find_first(st, "./w:rPr", ns)
  }
  for (tag in c("w:sz", "w:szCs")) {
    node <- xml_find_first(rpr, paste0("./", tag), ns)
    if (inherits(node, "xml_missing")) {
      xml_add_child(rpr, tag)
      node <- xml_find_first(rpr, paste0("./", tag), ns)
    }
    xml_attr(node, "w:val") <- as.character(half)
  }
  invisible(TRUE)
}

for (nm in names(HEADING_HALF)) set_style_size(nm, HEADING_HALF[[nm]])

write_xml(doc, styles_path)

# The theme-attribute strip, on the serialised text where the prefix is just
# characters. 64 of them in the stock file.
txt <- readLines(styles_path, warn = FALSE)
txt <- gsub("[[:space:]]+w:(ascii|hAnsi|cs|eastAsia)Theme=\"[^\"]*\"", "", txt)
writeLines(txt, styles_path)

# Belt and braces: point the theme itself at Calibri, so anything still
# resolving through majorHAnsi/minorHAnsi lands in the same place. The stock
# theme is Aptos.
theme_path <- file.path(xdir, "word", "theme", "theme1.xml")
if (file.exists(theme_path)) {
  # One string, not one per line: the element and its typeface attribute can
  # sit on different lines, and a per-line gsub would never see the pair.
  th <- paste(readLines(theme_path, warn = FALSE), collapse = "\n")
  th <- gsub("(<a:(major|minor)Font>[[:space:]]*<a:latin typeface=\")[^\"]*(\")",
             paste0("\\1", FONT, "\\3"), th)
  writeLines(th, theme_path)
}

# ------------------------------------------------------------------------------
# 3. Rezip
# ------------------------------------------------------------------------------
# Paths relative to the extraction root, or Word sees a nested folder and
# refuses the file.
# mode = "mirror", NOT "cherry-pick". cherry-pick adds each file at its
# basename, which flattens word/styles.xml to styles.xml and produces a .docx
# that Word and pandoc both reject. mirror keeps the paths relative to root.
files <- list.files(xdir, recursive = TRUE, all.files = TRUE, no.. = TRUE)
if (file.exists(OUT)) unlink(OUT)
zip::zip(zipfile = OUT, files = files, root = xdir, mode = "mirror")

# A flattened archive looks fine until Word opens it, so check the one path
# that proves the structure survived.
if (!"word/styles.xml" %in% zip::zip_list(OUT)$filename) {
  stop("Built ", OUT, " but word/styles.xml is not in it -- the archive was ",
       "flattened and Word will refuse the file.", call. = FALSE)
}

message("Wrote ", OUT)
message("Commit it: the memo's YAML points at reference-plain.docx and will ",
        "fail to render without it.")
