get_bss_overview <- function(bss_fit, ecg, ...){
  # Divergent transitions, summed across chains.
  #
  # This was:
  #   get_sampler_params(inc_warmup = FALSE) |>
  #     purrr::set_names(paste0("n_div_", 1:length(bss_fit@stan_args))) |>
  #     purrr::map_dbl(~.x[, "divergent__"] |> sum()) |> sum()
  #
  # which aborted with "`x` must be a vector" out of purrr::set_names() at four
  # chains, taking the whole render with it. The names it assigned were never
  # used -- map_dbl() summed the list and sum() collapsed it -- so the naming
  # step existed only to be fragile.
  #
  # vapply() over the list is what analysis/bss_bias/01_fit_bss_bias.R has
  # always used for the same quantity, and it does not care how many chains
  # there are. Same number, no purrr.
  sp <- rstan::get_sampler_params(bss_fit, inc_warmup = FALSE)
  n_div <- sum(vapply(sp, function(x) sum(x[, "divergent__"]), numeric(1)))

  bss_fit |>
    summary(pars = c("E_sum", "C_sum")) |>
    pluck("summary") |>
    as.data.frame() |>
    rownames_to_column("estimate") |>
    as_tibble() |>
    mutate(n_div = n_div, est_cg = ecg) |>
    relocate(estimate, est_cg)
}
