get_bss_overview <- function(bss_fit, ecg, ...){
  # Divergent transitions, summed across chains.
  #
  # This was:
  #   get_sampler_params(inc_warmup = FALSE) |>
  #     purrr::set_names(paste0("n_div_", 1:length(bss_fit@stan_args))) |>
  #     purrr::map_dbl(~.x[, "divergent__"] |> sum()) |> sum()
  #
  # which aborted with "`x` must be a vector" out of purrr::set_names(), taking
  # the whole render with it.
  #
  # set_names() raises that when its input is not a vector -- NULL being the
  # case that matters here, since get_sampler_params() returns NULL for a fit
  # with no usable post-warmup samples. So the message was never about the
  # naming; it was a degenerate fit reported in the least informative way
  # available. The names were not used either: map_dbl() summed the list and
  # sum() collapsed it, so the naming step existed only to be fragile.
  #
  # vapply() over the list is what analysis/bss_bias/01_fit_bss_bias.R has
  # always used for the same quantity. It returns 0 for an empty list instead
  # of aborting, which lets the caller see the fit and judge it.
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
