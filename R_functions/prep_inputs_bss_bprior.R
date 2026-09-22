# ==============================================================================
# prep_inputs_bss_bprior.R -- fork of prep_inputs_bss.R, 2026-09-22
#
# ONE CHANGE: b's prior is assembled as two length-2 vectors (mu and sigma,
# vehicle then trailer) from four scalar prior names, instead of one
# value_lognormal_sigma_b passed straight through. Pairs with
# stan_models/BSS_creel_model_02_2026-09-22_univariate_bprior.stan, which declares
# value_lognormal_mu_b / value_lognormal_sigma_b as vector[G] instead of a
# shared real -- the original prep_inputs_bss() would build Stan data this
# model can't read (no value_lognormal_mu_b field, and a scalar where a
# length-2 vector is declared).
#
# A DIFFERENT FUNCTION NAME, not an edit to prep_inputs_bss() in place:
# fw_creel.Rmd sources every file in R_functions/ (walk(list.files(...),
# source)), so this loads alongside the original with no collision, and
# every other basin's report keeps calling the unmodified prep_inputs_bss()
# against the unmodified model. This fork exists for one purpose: a
# census-free fishery-year that needs a history-informed vehicle prior with
# the trailer channel left at the uninformative default.
#
# Everything below prep_inputs_bss() is identical -- diff against it to
# confirm -- except the study_design branches are omitted here since this
# fork is Standard-only (see the source function for Drano handling), and
# the two `#priors` lines marked below.
# ==============================================================================

# create object/list of data inputs for BSS model

prep_inputs_bss_bprior <- function(
    period,
    days,            # tibble with time strata and closure fields
    dwg_summarized,  # list with shared interview, index and census tibbles
    est_catch_group, # data.frame passed from params of aggregated catch groups of interest to estimate
    census_expan,    # tibble summarizing p_census by angler_final and section_num where p_census is a hard coded value in the database specifying the proportion of a a section that is covered during a census count; values less than 1 with result in census counts being expanded (e.g., census count divided by p_census)
    priors,          # list of hyperpriors specified by user
    study_design,    # string passed from params denoting which study design was followed during data collection
    ...){

if(str_detect(study_design, "tandard" )){ 
  #in-function-scope intermediates
  effort_index_vehicle <- dwg_summarized$effort_index |> filter(str_detect(angler_final, "total"))   #angler_final = "total" is the same data as count_type = "vehicle" (see prep_dwg_effort_index function)
  effort_index_trailer <- dwg_summarized$effort_index |> filter(str_detect(angler_final, "boat"))    #angler_final = "boat" is the same data as count_type = "trailer" (see prep_dwg_effort_index function)
  effort_index_angler <- dwg_summarized$effort_index |> filter(str_detect(angler_final, "ngler")) |> #KB: as of April 2024, no projects following the "Standard" study design are counting anglers during index effort counts
    #based on outdated study design/database values that do not include current, more detailed census count levels
    mutate(
      angler_final_int = case_when(
        angler_final == "Bank Anglers"  ~ as.integer(1),
        angler_final == "Boat Anglers"  ~ as.integer(2)
      )
    )
  effort_census_boats <- dwg_summarized$effort_census |> filter(angler_final == "xxx") #KB: there were no boat census counts built into the standard creel study design; filtering "xxx" as a data wrangling trick to create empty dataframe 
  effort_census_anglers<- dwg_summarized$effort_census
  
  interview_cg <- 
    dwg_summarized$interview |> 
    filter(est_cg == est_catch_group)
  
  interview_cg_intA <- 
    interview_cg |> 
    drop_na(vehicle_count, trailer_count, person_count_final) |> 
    mutate(
      boat_count = as.integer(0)
    )
  
  interview_cg_boat<- interview_cg_intA # KB: interview_cg_boat isn't actually used in the "standard" study design analysis but has to be created/duplicated so that the case_when argument works
  
  interview_cg_daily_summ <- #KB NOTE: this object is used to create data objects that are no longer used in the most up-to-date model (perhaps could be omitted via deletion or at least commented out)
    interview_cg |> 
    group_by(event_date, section_num, angler_final_int) |> 
    summarise(across(c(fishing_time_total, fish_count), sum), .groups = "drop")

}else if(study_design == "Drano"){ 
  
  effort_index_vehicle <- dwg_summarized$effort_index |> filter(angler_final == "xxx") # there were no vehicle index counts built into the Drano Lake creel study design; filtering "xxx" as a data wrangling trick to create empty dataframe 
  effort_index_trailer <- dwg_summarized$effort_index |> filter(angler_final == "xxx") # there were no trailer index counts built into the Drano Lake creel study design; filtering "xxx" as a data wrangling trick to create empty dataframe
  effort_index_angler <- dwg_summarized$effort_index |>  filter(angler_final == "xxx") # there were no angler index counts built into the Drano Lake creel study design; technically, angler counts were recorded as index counts but really assumed to be census coutns;  filtering "xxx" as a data wrangling trick to create empty dataframe
  effort_census_boats <- dwg_summarized$effort_census |> filter(str_detect(angler_final, "boat"))
  effort_census_anglers<- dwg_summarized$effort_census |> filter(str_detect(angler_final, "bank"))
  
  interview_cg <- 
    dwg_summarized$interview |> 
    filter(est_cg == est_catch_group)

  interview_cg_intA <- 
    interview_cg |> 
    mutate(
      vehicle_count = as.integer(0), 
      trailer_count = as.integer(0),
      boat_count = ifelse(str_detect(angler_final, "boat"), 1, 0)
      )# KB: interview_cg_intA isn't actually used in the "Drano" study design analysis but has to be created/duplicated so that the case_when argument works

  interview_cg_daily_summ <- #KB NOTE: this object is used to create data objects that are no longer used in the most up-to-date model (perhaps could be omitted via deletion or at least commented out)
    interview_cg |> 
    group_by(event_date, section_num, angler_final_int) |> 
    summarise(across(c(fishing_time_total, fish_count), sum), .groups = "drop")

}  

#returned list object
stan_list <- list(
  est_cg = est_catch_group,
  D = nrow(days), # int; number of fishing days
  G = length(unique(interview_cg$angler_final_int)),  # int; final number of unique gear/angler types 
  # CENSUS-FREE YEAR (2026-09-22): S sourced from census_expan, not
  # effort_census. effort_census is legitimately empty when no tie-in counts
  # were collected this year, which made S = 0 while section_V/section_T
  # (from effort_index) still held real section indices -- caught by
  # preflight_bss_inputs() as "Section index 1 exceeds S = 0", which is what it
  # is for. census_expan (p_census/p_TI) is not gated on a live census event --
  # prep_dwg_census_expan() builds it from dwg$effort at large, defaulting to
  # full coverage where unset -- and after align_bss_sections() runs, it always
  # has exactly one row per (angler_final, section) for every section in the
  # aligned set. Same table p_TI pivots its columns from below, so S and p_TI's
  # column count now agree by construction rather than by coincidence.
  S = as.integer(length(unique(dwg_summarized$census_expan$section_num))),  # int; final number of river sections  
  H = max(dwg_summarized$effort_index$count_sequence), # int; max number of index counts within a sample day
  
  P_n = case_when( #int; total number of periods
    tolower(period) == 'day' ~ max(days$day_index),
    tolower(period) == 'week' ~ max(days$week_index),
    tolower(period) == 'month' ~ max(days$month_index),
    tolower(period) == 'duration' ~ as.integer(1)
  ),
  period = case_when( # int vec; index denoting fishing day/period
    tolower(period) == 'day' ~ days$day_index,
    tolower(period) == 'week' ~ days$week_index,
    tolower(period) == 'month' ~ days$month_index,
    tolower(period) == 'duration' ~ rep(as.integer(1), nrow(days))
  ),

  w = days$day_type_num, # int vec; 0/1 denoting Weekday/end for model offset
  L = days$day_length,  # num vec, daylength (model offset; assumption)
  # num mat; open/closed by section; 0 defined as 1E-6 for model
  O = days |> # only generates closure columns for sections with at least one observation from effort census counts
    select(contains("section_")) |>
    select(any_of(paste0("open_section_", unique(dwg_summarized$effort_census$section_num)))) |>
    mutate(across(everything(), ~if_else(., 1, 0.000001))) |> 
    as.matrix(),
    
  # Vehicle index effort counts 
  V_n = nrow(effort_index_vehicle), # int; total number of individual vehicle index effort counts 
  day_V = left_join(effort_index_vehicle, days, by = "event_date") |> pull(day_index),   # int; index for day/period 
  section_V = as.integer(effort_index_vehicle$section_num), # int; index for section
  countnum_V = as.integer(effort_index_vehicle$count_sequence),     # int; index for count_sequence  
  V_I = effort_index_vehicle$count_index, # num vec; observed # of vehicles 

  # Trailer index effort counts
  T_n = nrow(effort_index_trailer), # int; total number of boat trailer index effort counts 
  day_T = left_join(effort_index_trailer, days, by = "event_date") |> pull(day_index), # int; index for day/period
  section_T = as.integer(effort_index_trailer$section_num), # int vec; index for section
  countnum_T = as.integer(effort_index_trailer$count_sequence), # int vec; index for count_sequence  
  T_I = effort_index_trailer$count_index, # num vec; observed # of boat trailers 
    
  # Angler index effort counts
  A_n = nrow(effort_index_angler), # int; total number of angler index effort counts
  day_A = left_join(effort_index_angler, days, by = "event_date") |> pull(day_index), # int; index for day/period
  gear_A = effort_index_angler$angler_final_int, # int vec; index denoting "gear/angler type"
  section_A = as.integer(effort_index_angler$section_num), # int vec; index for section
  countnum_A = as.integer(effort_index_angler$count_sequence), # int vec; index for count_num
  A_I = effort_index_angler$count_index, #num vec; observed # of anglers

  # Census (tie-in) effort counts for boats 
  B_n = nrow(effort_census_boats), # int; total number of census effort counts for boats
  day_B = left_join(effort_census_boats, days, by = "event_date") |> pull(day_index), # int; index for day/period
  gear_B = effort_census_boats$angler_final_int, # int vec; index denoting "gear/angler type"
  section_B = as.integer(effort_census_boats$section_num), # int vec; index for section
  countnum_B = as.integer(effort_census_boats$count_sequence), # int vec; index for count_num
  B_s = effort_census_boats$count_census, #num vec; observed # of boat vessels

  # Census (tie-in) effort counts for anglers
  E_n = nrow(effort_census_anglers), # int; total number of census effort counts for anglers
  day_E = left_join(effort_census_anglers, days, by = "event_date") |> pull(day_index), # int vec; index denoting day/period
  gear_E = effort_census_anglers$angler_final_int, # int vec; index denoting "gear/angler type"  
  section_E = as.integer(effort_census_anglers$section_num), # int vec; index for section
  countnum_E = as.integer(effort_census_anglers$count_sequence), # int vec; index for count_sequence
  E_s = effort_census_anglers$count_census, # num vec; observed # of anglers
  
  #proportion spatial coverage during census (tie-in; TI) counts  
  # p_TI's ROW COUNT must equal G exactly -- Stan declares matrix[G,S] p_TI,
  # and census_expan ALWAYS carries both "bank" and "boat" rows regardless of
  # this fishery-year's actual gear-type composition (prep_dwg_census_expan()
  # builds it from the static per-section p_census lookup, not from which
  # gear types were actually interviewed). A bank-only fishery-year (G=1) fed
  # the unfiltered 2-row census_expan would hand Stan a matrix[2,S] against a
  # declared matrix[1,S] -- a guaranteed dimension mismatch, not a maybe.
  # Filtered here to the gear types actually present in THIS fishery-year's
  # interviews, bank(1) before boat(2), matching angler_final_int's own
  # convention -- so p_TI's row count tracks G by construction, generally,
  # not just for the bank-only case this fork was built for.
  p_TI = {
    gear_labels <- c("1" = "bank", "2" = "boat")
    needed <- gear_labels[as.character(sort(unique(interview_cg$angler_final_int)))]
    out <- census_expan |>
      filter(angler_final %in% needed) |>
      mutate(angler_final = factor(angler_final, levels = needed)) |>
      arrange(angler_final) |>
      select(angler_final, section_num, p_census) |>
      pivot_wider(names_from = section_num, values_from = p_census) |>
      select(-angler_final) |>
      as.matrix()
    if (nrow(out) != length(needed)) {
      stop("p_TI has ", nrow(out), " row(s) but this fishery-year's interviews ",
           "need ", length(needed), " (", paste(needed, collapse = ", "), "). ",
           "census_expan is missing a p_census entry for a gear type this ",
           "fishery-year's interviews actually have -- check prep_dwg_census_expan() ",
           "and the fishery_manager p_census_bank/p_census_boat lookup.", call. = FALSE)
    }
    out
  },
    
  # interview data - CPUE 
  IntC = nrow(distinct(interview_cg, interview_id)),  # int; total number of angler interviews with c & h data; distinct() here should be redundant
  day_IntC = left_join(interview_cg, days, by = "event_date") |> pull(day_index), # int vec; index denoting day/period   
  gear_IntC = interview_cg$angler_final_int, # int vec; index denoting "gear/angler type"
  section_IntC = interview_cg$section_num, # int vec; index for section
  c = interview_cg$fish_count, # num vec; total catch
  h = interview_cg$fishing_time_total, # num vec; total hours fished as fishing_time * person_count_final
    
  # # interview data - Total Effort & Catch Creeled (#KB: as of April 2024, the following data/parameters are not used in the most up-to-date BSS model so commented out)
  # IntCreel = nrow(interview_cg_daily_summ), # int; totals from interviews aggregated by date-section-anglertype	
  # day_Creel = left_join(interview_cg_daily_summ, days, by = "event_date") |> pull(day_index), # int vec; index denoting day/period
  # gear_Creel = interview_cg_daily_summ$angler_final_int,  # int vec; index denoting "gear/angler type"
  # section_Creel = interview_cg_daily_summ$section_num, # int vec; index for section 
  # C_Creel = interview_cg_daily_summ$fish_count, # num vec; total reported catch by day-section-anglertype
  # E_Creel = interview_cg_daily_summ$fishing_time_total,  #num vec; total hours fished by day-section-anglertype
    
  # interview data - objects per anglers
  IntA = nrow(distinct(interview_cg_intA, interview_id)),     # int; total number of angler interviews where V_A, T_A, A_A were collected
  day_IntA = left_join(interview_cg_intA, days, by = "event_date") |> pull(day_index), # int vec; index denoting day/period
  gear_IntA = interview_cg_intA$angler_final_int, # int vec; index denoting "gear/angler type"
  section_IntA = interview_cg_intA$section_num, # int vec; index for section
  V_A = as.integer(interview_cg_intA$vehicle_count),  # num vec; total number of vehicles an angler group brought
  T_A = as.integer(interview_cg_intA$trailer_count),  # num vec; total number of trailers an angler group brought
  B_A = as.integer(interview_cg_intA$boat_count),     # num vec; total number of boats an angler group brought
  A_A = as.integer(interview_cg_intA$person_count_final),  # num vec; total number of anglers in the groups interviewed

  #priors
  #hyperhyper scale (degrees of freedom) parameters
  value_cauchyDF_sigma_eps_C = priors["value_cauchyDF_sigma_eps_C"] , #for the hyperprior distribution sigma_eps_C; default = 1  
  value_cauchyDF_sigma_eps_E = priors["value_cauchyDF_sigma_eps_E"], #for the hyperprior distribution sigma_eps_E; default = 1  
  value_cauchyDF_sigma_r_E = priors["value_cauchyDF_sigma_r_E"],   #for the hyperprior distribution sigma_r_E; default = 1  
  value_cauchyDF_sigma_r_C = priors["value_cauchyDF_sigma_r_C"],   #for the hyperprior distribution sigma_r_C; default = 1 
  value_cauchyDF_sigma_mu_C = priors["value_cauchyDF_sigma_mu_C"],  #the hyperhyper SD parameter in the hyperprior distribution sigma_mu_C
  value_cauchyDF_sigma_mu_E = priors["value_cauchyDF_sigma_mu_E"],   #the hyperhyper SD parameter in the hyperprior distribution sigma_mu_E
  
  value_normal_sigma_omega_C_0 = priors["value_normal_sigma_omega_C_0"], #the SD hyperparameter in the prior distribution omega_C_0; normal sd (log-space); default = 1   
  value_normal_sigma_omega_E_0 = priors["value_normal_sigma_omega_E_0"], #the SD hyperparameter in the prior distribution omega_E_0; normal sd (log-space);; default = 3  
  # b's prior, PER CHANNEL -- the one change this fork makes. b[1] is the
  # vehicle index, b[2] the trailer index (positional in the Stan likelihood,
  # not angler-type-indexed -- see get_bss_bias.R's header comment). Order
  # matches that positional use, not alphabetical.
  value_lognormal_mu_b = c(
    priors[["value_lognormal_mu_b_vehicle"]],
    priors[["value_lognormal_mu_b_trailer"]]
  ),
  value_lognormal_sigma_b = c(
    priors[["value_lognormal_sigma_b_vehicle"]],
    priors[["value_lognormal_sigma_b_trailer"]]
  ),
  value_normal_sigma_B1 = priors["value_normal_sigma_B1"],        #the SD hyperparameter in the prior distribution B1; default = 5  
  value_normal_mu_mu_C = priors["value_normal_mu_mu_C"], #the mean hyperparameter in the prior distribution mu_C; median (log-space); default = 0.02 (was originally  0.05) 
  value_normal_sigma_mu_C = priors["value_normal_sigma_mu_C"],    #the SD hyperparameter in the prior distribution mu_C; normal sd (log-space); default = 1.5 (was originally 5)
  value_normal_mu_mu_E = priors["value_normal_mu_mu_E"],    #the mean hyperparameter in the prior distribution mu_E; median effort (log-space); default = 15 
  value_normal_sigma_mu_E = priors["value_normal_sigma_mu_E"],      #the SD hyperparameter in the prior distribution mu_E; normal sd (log-space); default = 2 (was originally 5) 
  value_betashape_phi_E_scaled = priors["value_betashape_phi_E_scaled"], #the rate (alpha) and shape (beta) hyperparameters in phi_E_scaled; default = 1 (i.e., beta(1,1) which is uniform), alternative beta(2,2) 
  value_betashape_phi_C_scaled = priors["value_betashape_phi_C_scaled"] #the rate (alpha) and shape (beta) hyperparameters in phi_C_scaled; default = 1 (i.e., beta(1,1) which is uniform), alternative beta(2,2)
  )
  
 return(stan_list)   
}