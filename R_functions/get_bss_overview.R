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
  # well-behaved summary. n_finite/n_draws makes the cause visible in the
  # table itself instead of requiring a trip back through console output.
  finite_frac <- function(par) {
    draws <- unlist(rstan::extract(bss_fit, pars = par), use.names = FALSE)
    tibble(estimate = par, n_finite = sum(is.finite(draws)), n_draws = length(draws))
  }

  # Per-chain draws as an iteration x chain matrix, the input shape
  # posterior::rhat()/ess_bulk()/ess_tail() expect for a single variable
  # (same functions this codebase already calls, via summarize_draws(),
  # for params_keep -- not a fresh guess). permuted = FALSE keeps chains
  # separate (required for split-Rhat); [, , 1] drops the trailing
  # parameter dimension, which extract() gives length 1 for a scalar
  # generated quantity like E_sum/C_sum.
  draws_matrix <- function(par) {
    arr <- rstan::extract(bss_fit, pars = par, permuted = FALSE, inc_warmup = FALSE)
    arr[, , 1]
  }

  # Rhat/ESS recomputed on the same non-finite draws that poison rstan's
  # summary(), after dropping them -- rather than a from-scratch Rhat
  # implementation (easy to get subtly wrong with no rstan/Stan environment
  # here to check one against), this drops whole ITERATIONS (same iteration
  # index across every chain, so the per-chain matrix stays rectangular --
  # split-Rhat needs equal chain lengths) that contain a non-finite draw in
  # ANY chain, then hands the reduced matrix to posterior's own algorithm.
  # With n_finite within a few of n_draws (the poisson_rng overflow is rare),
  # this drops a handful of iterations out of thousands -- not a meaningful
  # loss of effective sample, and a real Rhat/ESS instead of NaN.
  finite_rhat_ess <- function(par) {
    m <- draws_matrix(par)
    bad_iter <- apply(m, 1, function(row) any(!is.finite(row)))
    m_finite <- m[!bad_iter, , drop = FALSE]
    enough <- nrow(m_finite) > 1
    tibble(
      estimate = par,
      n_dropped_iter = sum(bad_iter),
      rhat_finite = if (enough) posterior::rhat(m_finite) else NA_real_,
      ess_bulk_finite = if (enough) posterior::ess_bulk(m_finite) else NA_real_,
      ess_tail_finite = if (enough) posterior::ess_tail(m_finite) else NA_real_
    )
  }

  # Pareto k-hat (Vehtari, Gelman, Simpson, Yao & Gabry) -- the standard
  # answer, in the same posterior/loo ecosystem this codebase already uses
  # for rhat/ess_bulk, to "is the sample mean of these draws even reliable."
  # It fits a generalised Pareto distribution to the upper tail and returns
  # its shape parameter; k_hat < 0.5 means the tail is thin enough for the
  # mean to converge at the usual sqrt(n) rate (trust mean/sd); 0.5-0.7 means
  # it is estimable but converges slowly (treat with caution); >= 0.7 means
  # the tail is heavy enough that the sample mean/variance may not even be
  # finite -- exactly what "mean 185,578 next to a median of 48" is a
  # symptom of. Median/quantiles, being rank-based rather than moment-based,
  # stay trustworthy regardless of k_hat -- this diagnoses WHICH columns to
  # trust, it does not fix the mean/sd columns themselves.
  #
  # First render returned khat = NA / khat_flag = "not computed" with no
  # error, for every row -- i.e. pareto_khat() exists and ran, it just
  # returns NA silently rather than throwing, on the same rows where
  # n_finite < n_draws. Same root cause as the Rhat/n_eff NaNs above:
  # non-finite draws poisoning a diagnostic that isn't itself NaN-tolerant.
  # Filtering to is.finite() draws before the call (pareto_khat pools all
  # chains into one vector already, so no rectangular-matrix constraint
  # here, unlike rhat/ess above) is the same fix in the same spirit.
  pareto_k <- function(par) {
    draws <- unlist(rstan::extract(bss_fit, pars = par), use.names = FALSE)
    draws_finite <- draws[is.finite(draws)]
    err_msg <- NULL
    k <- tryCatch(
      {
        val <- posterior::pareto_khat(draws_finite)
        if (is.list(val)) val$khat else as.numeric(val)[1]
      },
      error = function(e) {
        err_msg <<- conditionMessage(e)
        NA_real_
      }
    )
    tibble(
      estimate = par,
      khat = k,
      # err_msg surfaced in khat_flag itself (not a separate silently-dropped
      # column) so a call signature/version mismatch shows up in the
      # rendered table directly.
      khat_flag = dplyr::case_when(
        !is.null(err_msg) ~ paste0("error: ", err_msg),
        is.na(k)  ~ "not computed",
        k < 0.5   ~ "mean reliable",
        k < 0.7   ~ "mean reliable, converges slowly",
        TRUE      ~ "mean UNRELIABLE -- use median/CI"
      )
    )
  }

  finite_counts <- bind_rows(finite_frac("E_sum"), finite_frac("C_sum")) |>
    left_join(bind_rows(finite_rhat_ess("E_sum"), finite_rhat_ess("C_sum")), by = "estimate") |>
    left_join(bind_rows(pareto_k("E_sum"), pareto_k("C_sum")), by = "estimate")

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
