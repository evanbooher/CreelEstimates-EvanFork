#thin wrapper on stan()
fit_bss <- function(
  # model_file = here::here("stan_models/BSS_creel_model_02_2021-01-22.stan"),
  #model_file_name = here::here(paste0("stan_models/", model_file_name)), #BSS_creel_model_02_2021-01-22.stan"),
  # model_file_name = here::here("stan_models/BSS_creel_model_02_2021-01-22_ppc.stan"),
  model_file_name,
  bss_inputs_list,
  n_chain = n_chain,
  n_cores = n_cores,
  n_iter = n_iter,
  n_warmup = n_warmup,
  n_thin = n_thin,
  adapt_delta = adapt_delta,
  max_treedepth = max_treedepth,
  init = "0",
  pars = NA,     # optional character vector restricting which parameters are monitored/returned.
                 # NA (default) preserves prior behavior (every parameter kept). Restricting this
                 # is the main lever for shrinking fit size/runtime when only a few parameters
                 # (e.g. "b") are actually needed downstream -- see R_functions/get_bss_bias.R.
  include = TRUE, # TRUE keeps exactly `pars`; FALSE keeps everything EXCEPT `pars` (rstan::stan() semantics)
  ...){

  model_path <- here::here("stan_models", model_file_name)

  if(!file.exists(model_path)){
    stop(
      "Stan model file not found: ", model_path,
      "\nCheck params$bss_model_file_name matches a file in stan_models/"
    )
  }

  # rstan reads any length-1 R vector as a SCALAR, but these data are declared
  # vector[D] or int x[n] in the models. A one-day window (D = 1) -- or any
  # single-row table (V_n = 1, IntA = 1, ...) -- then fails at data
  # initialization: "dims declared=(1); dims found=()". as.array() keeps the
  # dimension. Scalars are not in this list, and neither is
  # value_lognormal_*_b: a real in the original model and a length-2 vector in
  # the b-prior models, never length 1. Checked against the data blocks of
  # BSS_creel_model_02_2021-01-22_ppc, _2026-09-22_bprior and
  # _2026-09-22_univariate_bprior.
  vector_data <- c(
    "w", "period", "L",
    "day_V", "section_V", "countnum_V", "V_I",
    "day_T", "section_T", "countnum_T", "T_I",
    "day_A", "gear_A", "section_A", "countnum_A", "A_I",
    "day_E", "gear_E", "section_E", "countnum_E", "E_s",
    "day_IntC", "gear_IntC", "section_IntC", "c", "h",
    "day_IntA", "gear_IntA", "section_IntA", "V_A", "T_A", "A_A"
  )
  for (nm in intersect(vector_data, names(bss_inputs_list))) {
    if (length(bss_inputs_list[[nm]]) == 1) {
      bss_inputs_list[[nm]] <- as.array(bss_inputs_list[[nm]])
    }
  }

  stan_args <- list(
    file = model_path,
    data = bss_inputs_list,
    chains = n_chain,
    cores = n_cores,
    iter = n_iter,
    warmup = n_warmup,
    thin = n_thin, init = init, include = include,
    control = list(
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    )
  )
  # rstan::stan() errors if `pars` is passed as NA rather than omitted, so only
  # attach it when the caller actually wants to restrict the monitored set.
  if (!identical(pars, NA)) stan_args$pars <- pars

  do.call(stan, stan_args)

}