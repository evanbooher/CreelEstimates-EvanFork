# Background job script - runs .Rmd without rendering final document
rmarkdown::render(
  input = "fw_creel_HWS24_CD72_20260318.Rmd",  # Change to your .Rmd filename
  run_pandoc = TRUE,            # Skip final rendering
  quiet = FALSE,                 # Set TRUE to suppress messages
  envir = new.env()              # Run in clean environment
)