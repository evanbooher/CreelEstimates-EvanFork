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

  # n_eff/Rhat for E_sum/C_sum specifically (2026-09-22): rstan's summary()
  # computes these two columns with a variance-based split-Rhat calculation
  # that is NOT NaN-tolerant, unlike the mean/sd/quantile columns next to
  # them -- ONE non-finite draw anywhere in a chain (poisson_rng rate
  # overflow on low-catch data; see the "Missing or NaN values detected"
  # warning this same render emits) silently poisons Rhat/n_eff for that
  # whole parameter, while the rest of the row still looks like a normal,
  # well-behaved summary. Left as NaN here rather than recomputed on a
  # filtered subset -- a from-scratch NaN-tolerant Rhat is easy to get
  # subtly wrong, and this codebase has no rstan/Stan environment to verify
  # one against. n_finite/n_draws makes the cause visible in the table
  # itself instead of requiring a trip back through console output: a low
  # n_div with n_finite == n_draws and Rhat still NaN is a genuine mystery
  # worth escalating; NaN Rhat alongside even one non-finite draw is this
  # known, explained case, not evidence sampling failed.
  finite_frac <- function(par) {
    draws <- unlist(rstan::extract(bss_fit, pars = par), use.names = FALSE)
    tibble(estimate = par, n_finite = sum(is.finite(draws)), n_draws = length(draws))
  }
  finite_counts <- bind_rows(finite_frac("E_sum"), finite_frac("C_sum"))

  bss_fit |>
    summary(pars = c("E_sum", "C_sum")) |>
    pluck("summary") |>
    as.data.frame() |>
    rownames_to_column("estimate") |>
    as_tibble() |>
    left_join(finite_counts, by = "estimate") |>
    mutate(n_div = n_div, est_cg = ecg) |>
    relocate(estimate, est_cg)
}
