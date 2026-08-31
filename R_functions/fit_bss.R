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