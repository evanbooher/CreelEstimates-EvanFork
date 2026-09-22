# Extract the BSS effort-index bias term ("b") from a fitted stanfit object.
#
# NOTE ON WHAT "b" IS: in the BSS Stan model (see e.g.
# stan_models/BSS_creel_model_02_2024-04-03.stan), `b` is declared
# `vector<lower=0>[G] b;` and priored `b[g] ~ lognormal(0, value_lognormal_sigma_b)`
# (prior median = 1, "no bias"). Despite being indexed 1:G, it is NOT an
# angler-type-indexed quantity -- it is used positionally as:
#   b[1] multiplies the expected VEHICLE index count (V_I likelihood)
#   b[2] multiplies the expected TRAILER index count (T_I likelihood)
# This function labels them "vehicle" / "trailer" accordingly, never "bank" /
# "boat". Do not relabel by angler_final without re-reading the Stan block.
#
# `b` is confounded with R_V (vehicles-per-angler) in the vehicle likelihood
# and is only separated by the interview-side binomial on V_A/A_A -- so a
# fishery-year with few/no vehicle-count interviews (IntA, V_A) will produce
# a `b` posterior that mostly reflects its prior rather than data. See
# `prior_contraction` below, and the IntA/V_n preflight checks in
# analysis/bss_bias/01_fit_bss_bias.R.

get_bss_bias <- function(
    bss_fit,          # a stanfit object returned by fit_bss()
    fishery_name,
    ecg,              # the est_cg string this fit was run on (see build_est_catch_groups-style helper)
    prior_sigma_b = 1, # value_lognormal_sigma_b used for this fit; must match what was passed into prep_inputs_bss()
    probs = c(0.025, 0.10, 0.50, 0.90, 0.975),
    ...) {

  draws <- posterior::as_draws_df(bss_fit) |>
    posterior::subset_draws(variable = "b")

  if (ncol(dplyr::select(draws, dplyr::starts_with("b["))) == 0) {
    rlang::abort(
      paste0("No 'b' parameter found in stanfit draws for fishery '", fishery_name,
             "'. Was fit_bss() called with pars = ... that excluded \"b\"?")
    )
  }

  # lognormal(0, sigma) prior variance, used to gauge how much the posterior
  # has moved away from the prior (see `prior_contraction` below).
  var_prior <- (exp(prior_sigma_b^2) - 1) * exp(prior_sigma_b^2)

  # NOTE: summarise_draws() ignores the `...` argument names below for
  # formula-wrapped stats::quantile() calls -- it names the resulting column
  # after quantile()'s own output name ("2.5%", not "q2.5"), confirmed
  # against a live posterior::summarise_draws() call. Rename explicitly
  # afterward rather than relying on the argument names, since q2.5/q10/
  # q90/q97.5 are load-bearing column names downstream (03_plot_b_series.R,
  # 04_candidate_options.R both read them directly from bss_b_summary.csv).
  summ <- posterior::summarise_draws(
    draws,
    mean, sd, median = ~stats::median(.x),
    q2.5  = ~stats::quantile(.x, probs = 0.025),
    q10   = ~stats::quantile(.x, probs = 0.10),
    q90   = ~stats::quantile(.x, probs = 0.90),
    q97.5 = ~stats::quantile(.x, probs = 0.975),
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  ) |>
    dplyr::rename(q2.5 = `2.5%`, q10 = `10%`, q90 = `90%`, q97.5 = `97.5%`)

  summ |>
    dplyr::mutate(
      fishery_name = fishery_name,
      est_cg       = ecg,
      bias_type    = dplyr::case_when(
        variable == "b[1]" ~ "vehicle",
        variable == "b[2]" ~ "trailer",
        TRUE ~ variable  # future-proofing if G ever has a b[3]+; should not occur under current model
      ),
      post_var          = sd^2,
      prior_sigma       = prior_sigma_b,
      prior_var         = var_prior,
      prior_contraction = 1 - (post_var / prior_var)
    ) |>
    dplyr::rename(param = variable) |>
    dplyr::select(
      fishery_name, est_cg, param, bias_type,
      mean, sd, median, q2.5, q10, q90, q97.5,
      rhat, ess_bulk, ess_tail,
      post_var, prior_sigma, prior_var, prior_contraction
    )
}
