// ==============================================================================
// BSS_creel_model_02_2026-09-22_univariate_bprior.stan
//
// "univariate": no cross-section/cross-gear correlation structure (no
// Lcorr_E/Lcorr_C) -- each section's effort/CPUE process error is its own
// independent AR(1) series, as in the univariate parent below, not the
// jointly-correlated version *_2021-01-22_ppc.stan uses. "bprior": b's prior
// is vectorised per channel (see change 1 below).
//
// Based on BSS_creel_model_09_2022-09-30_1G_1S.stan (supplied 2026-09-22), not
// on the newer *_2021-01-22_ppc.stan the rest of the pipeline compiles today.
// This older file already carries the property this fork needs: its model-
// block likelihoods for V_I and T_I use only a single (g=1) term, rather than
// summing a bank and a boat term the way the current production model does --
// i.e. it already assumes ONE gear type, which is what a fishery-year with
// zero boat anglers in every interview this season actually has. Started from
// this file and changed only what running it for real required.
//
// TWO CHANGES from the supplied file:
//
// 1. b's prior. b is declared vector<lower=0>[2] b -- fixed size 2 -- but was
//    priored inside `for(g in 1:G){ b[g] ~ lognormal(0,value_lognormal_sigma_b) }`,
//    which for this fishery-year's G=1 only ever touches b[1]: b[2], the
//    trailer channel's bias, got no explicit prior statement at all. Moved to
//    its own for(g in 1:2) loop and vectorised per channel, so the vehicle
//    term can carry a history-derived prior while the trailer term keeps
//    lognormal(0,1) explicitly, rather than an implicit, undocumented default.
//
// 2. log_lik's V_I/T_I terms summed a g=2 (boat) component the model block
//    itself never uses. R_V, R_T, p_TI and lambda_E_S_I are all sized by G
//    (not hardcoded to 2 the way b is), so indexing their second gear-type
//    column is out of bounds whenever G=1 -- a runtime crash, since Stan
//    executes the whole generated quantities block every iteration whether
//    or not log_lik is ever read downstream (it isn't, here). Reduced to the
//    same single g=1 term the model block's likelihood already uses.
//
// UNTESTED: no Stan compiler is available in the environment that wrote this
// fork. Watch the full compiler error on first run, not just "it failed" --
// Stan's messages name the exact line.
//
// T_n/T_I/R_T/T_A (trailer index counts, kept, not removed): this season
// still records occasional nonzero trailer index counts even with zero
// boat anglers interviewed, so the supplied file's structure -- both
// channels tied to the single gear type rather than the trailer channel
// removed outright -- is kept as given, on the reasoning that dropping real
// observations is a bigger, less reversible choice than modelling them.
//
// No PPD (_rep) variables -- the supplied file doesn't have them, and
// nothing downstream in fw_creel_bprior_2026.Rmd reads _rep or log_lik, so
// their absence changes nothing there.
// ==============================================================================
data{
	//Day attributes
	int<lower=0> D; //number of fishing days (sampling frame)
	int<lower=0> G; //number of unique gear/angler types 
	int<lower=0> S; //number of river sections
	int<lower=0> H; //max number of angler effort counts within a sample day across entire sampling frame (max countnum)
	int<lower=0> P_n; //number of periods (number of states for our state variable); P_n can equal D or some other interval (e.g., weekly)
	vector<lower=0,upper=1>[D] w; //index denoting daytype, where 0=weekday, 1= weekend/holiday.  
	int<lower=0> period[D]; //index denoting period    
	vector<lower=0>[D] L; // total amount of available fishing hours per day (e.g., day length - sunrise to sunset)
	matrix<lower=0>[D,S] O; //index denoting fishery status, where 1=open, 0 = closed 
	//Vehicle index effort counts
	int<lower=0> V_n; //total number of individual vehicle index effort counts
	int<lower=0> day_V[V_n]; //index denoting the "day" for an individual vehicle index effort count
	int<lower=0> section_V[V_n]; //index denoting the "section" for an individual vehicle index effort count
	int<lower=0> countnum_V[V_n]; //index denoting the "count number" for an individual vehicle index effort count
	int<lower=0> V_I[V_n]; //number of vehicles enumerated during an individual index effort survey
	//Trailer index effort counts 
	int<lower=0> T_n; //total number of boat trailer index effort counts
	int<lower=0> day_T[T_n]; //index denoting the "day" for an individual boat trailer index effort count
	int<lower=0> section_T[T_n]; //index denoting the "section" for an individual boat trailer index effort count
	int<lower=0> countnum_T[T_n]; //index denoting the "count number" for an individual boat trailer index effort count
	int<lower=0> T_I[T_n]; //number of boat trailers enumerated during an individual index effort survey
	//Angler index effort counts
	int<lower=0> A_n; //total number of angler index effort counts             
	int<lower=0> day_A[A_n]; //index denoting the "day" for an individual angler index effort count     
	int<lower=0> gear_A[A_n]; //index denoting the "gear/angler type" for an individual angler index effort count (e.g., 1=bank and 2=boat)     
	int<lower=0> section_A[A_n]; //index denoting the "section" for an individual angler index effort count  
	int<lower=0> countnum_A[A_n]; //index denoting the "count number" for an individual angler index effort count
	int<lower=0> A_I[A_n]; //number of anglers enumerated during an individual index effort survey         
	//Census (tie-in) effort counts
	int<lower=0> E_n; //total number of angler census effort counts 
	int<lower=0> day_E[E_n]; //index denoting the "day" for an individual angler census effort count
	int<lower=0> gear_E[E_n]; //index denoting the "gear/angler type" for an individual angler census effort count (e.g., 1=bank and 2=boat) 
	int<lower=0> section_E[E_n]; //index denoting the "section" for an individual angler census effort count  
	int<lower=0> countnum_E[E_n]; //index denoting the "count number" for an individual angler census effort count
	int<lower=0> E_s[E_n]; //number of anglers enumerated during an individual census effort survey
	//Proportion tie-in expansion
	matrix<lower=0,upper=1>[G,S]p_TI; //proportion of section covered by tie in counts (serves as an expansion if p_TI != 1)
	//interview data - CPUE
	int<lower=0> IntC; //total number of angler interviews conducted across all surveys dates where CPUE data (c & h) were collected                                       	
	int<lower=0> day_IntC[IntC]; //index denoting the "day" for an individual angler interview                       	
	int<lower=0> gear_IntC[IntC]; //index denoting the "gear/angler type" for an individual angler interview (e.g., 1=bank and 2=boat)      						
	int<lower=0> section_IntC[IntC]; //index denoting the "section" for an individual angler interview    						
	int<lower=0> c[IntC]; // total number of fish caught by an angler (group) collected from an individual angler interview           						
	vector<lower=0>[IntC] h; // total number of hours fish by an angler (group) collected from an individual angler interview        						
	//interview data - angler expansion
	int<lower=0> IntA; //total number of angler interviews conducted across all surveys dates where angler expansion data (V_A, T_A, A_A) were collected                                         	
	int<lower=0> day_IntA[IntA]; //index denoting the "day" for an individual angler interview                       	
	int<lower=0> gear_IntA[IntA]; //index denoting the "gear/angler type" for an individual angler interview (e.g., 1=bank and 2=boat)      						
	int<lower=0> section_IntA[IntA]; //index denoting the "section" for an individual angler interview    
	int<lower=0> V_A[IntA]; // total number of vehicles by an individual angler/group brought to the fishery on a given survey date		
	int<lower=0> T_A[IntA]; // total number of boat trailers by an individual angler/group brought to the fishery on a given survey date			
	int<lower=0> A_A[IntA]; // total number of anglers in the group interviewed 		
	//hyper and hyperhyper parameters
	real value_cauchyDF_sigma_eps_C; //the hyperhyper scale (degrees of freedom) parameter in the hyperprior distribution sigma_eps_C 
	real value_cauchyDF_sigma_eps_E; //the hyperhyper scale (degrees of freedom) parameter in the hyperprior distribution sigma_eps_E
	real value_cauchyDF_sigma_r_E; //the hyperhyper scale (degrees of freedom) parameter in the hyperprior distribution sigma_r_E
	real value_normal_sigma_B1; //the SD hyperparameter in the prior distribution B1
	real value_cauchyDF_sigma_r_C; //the hyper scale (degrees of freedom) parameter in the prior distribution sigma_r_C
	real value_betashape_phi_C_scaled; //the rate (alpha) and shape (beta) hyperparameters in phi_C_scaled 
	real value_betashape_phi_E_scaled; //the rate (alpha) and shape (beta) hyperparameters in phi_E_scaled 
	real value_normal_sigma_omega_C_0; // the SD hyperparameter in the prior distribution omega_C_0
	real value_normal_sigma_omega_E_0; // the SD hyperparameter in the prior distribution omega_E_0
	vector[2] value_lognormal_mu_b; //the mean hyperparameter (on the log scale) in the prior distribution b, per channel (b[1] vehicle, b[2] trailer). Sized 2, matching b's own fixed size -- NOT sized G, since b is priored regardless of how many gear types this fishery-year has.
	vector<lower=0>[2] value_lognormal_sigma_b; //the SD hyperparameter in the prior distribution b, per channel. Was a single real shared by both, mu fixed at 0.
	real value_normal_mu_mu_C; //the hyperhyper mean parameter in the hyperprior distribution mu_mu_C
	real value_normal_sigma_mu_C; //the hyperhyper SD parameter in the hyperprior distribution mu_mu_C
	real value_normal_mu_mu_E; //the hyperhyper mean parameter in the hyperprior distribution mu_mu_E
	real value_normal_sigma_mu_E; //the hyperhyper SD parameter in the hyperprior distribution mu_mu_E
	real value_cauchyDF_sigma_mu_C; //the hyperhyper SD parameter in the hyperprior distribution sigma_mu_C
	real value_cauchyDF_sigma_mu_E; //the hyperhyper SD parameter in the hyperprior distribution sigma_mu_E
}
transformed data{
}
parameters{
	//Effort
	real B1; //fixed effect accounting for the effect of day type on effort   
	real<lower=0> sigma_eps_E; //effort process error standard deviation 
	real<lower=0> sigma_r_E; //prior on r_E
	real<lower=0,upper=1> phi_E_scaled; //prior on a transformation of phi_E								    			
	matrix[P_n-1,G] eps_E[S]; //effort process errors  
	matrix[G,S] omega_E_0; //effort residual for initial time step (p=1)
	vector<lower=0,upper=1>[G] R_V; //true angler vehicles per angler
	vector<lower=0,upper=1>[G] R_T; //true angler trailers per angler
	vector<lower=0>[2] b; //bias in counts of vehicles and trailers per angler group from road counts
	matrix<lower=0>[D,G] eps_E_H[S,H]; //gamma random variate accounting for overdispersion in the census effort counts due to within-day variability in angler pressure
	matrix<lower=0,upper=1>[G,S] p_I; //fixed proportion of angler effort observed in an index area
	real mu_mu_E[G]; //hyper-prior on mean of mu_E //TB 5/3/2019
	real<lower=0>sigma_mu_E; //hyper-prior on SD of mu_E //TB 5/3/2019
	matrix[G,S] eps_mu_E;
	//Catch rates
	real<lower=0> sigma_eps_C; //catch rate (CPUE) process error standard deviation 
	real<lower=0,upper=1> phi_C_scaled; //Prior on a transformation of phi_C									    			
	real<lower=0> sigma_r_C; //Prior on r_C
	matrix[P_n-1,G] eps_C[S]; //CPUE process errors 
	matrix[G,S] omega_C_0; //CPUE residual for initial time step (p=1)
	real mu_mu_C[G] ; //hyper-prior on mean of mu_C //TB 5/3/2019
	real<lower=0>sigma_mu_C; //hyper-prior on SD of mu_C
	matrix[G,S] eps_mu_C;
}
transformed parameters{
	//Effort
	matrix[G,S] mu_E; //season-long effort intercept  
	real<lower=-1,upper=1> phi_E; //auto-regressive (AR), mean-reverting lag-1 coefficient for effort 
	matrix[P_n,G] omega_E[S]; //residual in effort 
	matrix<lower=0>[D,G] lambda_E_S[S]; //mean daily effort
	matrix<lower=0>[D,G] lambda_E_S_I[S,H]; //mean hourly effort
	//Catch rates
	matrix[G,S] mu_C; //season-long catch rate intercept  
	real<lower=-1,upper=1> phi_C; //auto-regressive (AR), mean-reverting lag-1 coefficient for CPUE 
	real<lower=0> r_C; //over-dispersion parameter accounting for among angler (group) variability in CPUE
	matrix[P_n,G] omega_C[S]; //Residual in CPUE
	matrix<lower=0>[D,G] lambda_C_S[S]; //mean daily CPUE
	r_C = 1 / square(sigma_r_C);
	phi_C = (phi_C_scaled * 2)-1;
	phi_E = (phi_E_scaled * 2)-1;
	for(g in 1:G){
		for(s in 1:S){
		  mu_C[g,s] = mu_mu_C[g] + eps_mu_C[g,s] * sigma_mu_C;
		  mu_E[g,s] = mu_mu_E[g] + eps_mu_E[g,s] * sigma_mu_E;
			omega_C[s][1,g] = omega_C_0[g,s];
			omega_E[s][1,g] = omega_E_0[g,s];
		}
		for(p in 2:P_n){ 
			for(s in 1:S){
				omega_C[s][p,g] = phi_C * omega_C[s][p-1,g] + eps_C[s][p-1,g] * sigma_eps_C; 
				omega_E[s][p,g] = phi_E * omega_E[s][p-1,g] + eps_E[s][p-1,g] * sigma_eps_E; 
			}													
		}
		for(d in 1:D){       
			for(s in 1:S){	
				lambda_C_S[s][d,g] = exp(mu_C[g,s] + omega_C[s][period[d],g]) * O[d,s];
				lambda_E_S[s][d,g] = exp(mu_E[g,s] + omega_E[s][period[d],g] + B1 * w[d])* O[d,s];
				for(i in 1:H){
					lambda_E_S_I[s,i][d,g] = lambda_E_S[s][d,g] * eps_E_H[s,i][d,g];									
				}
			}
		}								    
	}	
}
model{
	//Hyperpriors (effort hyperparameters)
	sigma_eps_E ~ cauchy(0,value_cauchyDF_sigma_eps_E);                                						
  phi_E_scaled ~ beta(value_betashape_phi_E_scaled,value_betashape_phi_E_scaled);
	sigma_r_E ~ cauchy(0,value_cauchyDF_sigma_r_E);
	B1 ~ normal(0,value_normal_sigma_B1);
	//Hyperpriors (CPUE hyperparameters)
	sigma_eps_C ~ cauchy(0,value_cauchyDF_sigma_eps_C);                     						
  phi_C_scaled ~ beta(value_betashape_phi_C_scaled,value_betashape_phi_C_scaled);
	sigma_r_C ~ cauchy(0, value_cauchyDF_sigma_r_C);
	sigma_mu_C~cauchy(0,value_cauchyDF_sigma_mu_C);//TB 5/3/2019
	sigma_mu_E~cauchy(0,value_cauchyDF_sigma_mu_E);//TB 5/3/2019
	//Priors 
	for(g in 1:G){
		mu_mu_C[g] ~ normal(value_normal_mu_mu_C,value_normal_sigma_mu_C); //TB 5/3/2019
		mu_mu_E[g] ~ normal(value_normal_mu_mu_E,value_normal_sigma_mu_E); //TB 5/3/2019
		for(p in 2:P_n){ 
			for(s in 1:S){
				eps_C[s][p-1,g] ~ std_normal();
				eps_E[s][p-1,g] ~ std_normal();
			}
		}
		for(d in 1:D){
			for(s in 1:S){  					  
				for(i in 1:H){
					eps_E_H[s,i][d,g] ~ gamma(1/square(sigma_r_E),1/square(sigma_r_E)); 
				}
			}
		}
		for(s in 1:S){
			omega_C_0[g,s] ~ normal(0,value_normal_sigma_omega_C_0); 
			omega_E_0[g,s] ~ normal(0,value_normal_sigma_omega_E_0); 
			eps_mu_C[g,s] ~ std_normal();
			eps_mu_E[g,s] ~ std_normal();
			p_I[g,s] ~ beta(0.5,0.5);
		}
		R_V[g] ~ beta(0.5,0.5); //Note: leaving constant among days AND sections...may need to tweak; can make beta because is "true" angler cars or angler trailers per angler!
		R_T[g] ~ beta(0.5,0.5); //Note: leaving constant among days AND sections...may need to tweak; can make beta because is "true" angler cars or angler trailers per angler!
	}
	// b is priored on its OWN loop, sized 2 (not G): the loop above only runs
	// 1:G, which for a single-gear-type fishery-year (G=1) would leave b[2] --
	// the trailer channel's bias -- with no explicit prior statement at all.
	for(g in 1:2){
		b[g] ~ lognormal(value_lognormal_mu_b[g],value_lognormal_sigma_b[g]); //Note: leaving constant among days AND sections...may need to tweak could go as low as 0.25 for sigma
	}
	//Likelihoods
	//Index effort counts - vehicles
	for(i in 1:V_n){
		V_I[i] ~ poisson((lambda_E_S_I[section_V[i],countnum_V[i]][day_V[i],1] * p_TI[1,section_V[i]] * R_V[1]) * b[1]);
	}
	//Index effort counts - trailers
	for(i in 1:T_n){
		T_I[i] ~ poisson((lambda_E_S_I[section_T[i],countnum_T[i]][day_T[i],1] * p_TI[1,section_T[i]] * R_T[1]) * b[2]); 
	}
	//Index effort counts - anglers
	for(i in 1:A_n){
		A_I[i] ~ poisson(lambda_E_S_I[section_A[i],countnum_A[i]][day_A[i],gear_A[i]] * p_TI[gear_A[i],section_A[i]] * p_I[gear_A[i],section_A[i]]);
	}
	//Census (tie-in) effort counts - anglers
	for(e in 1:E_n){
		E_s[e] ~ poisson(lambda_E_S_I[section_E[e],countnum_E[e]][day_E[e],gear_E[e]] * p_TI[gear_E[e],section_E[e]]);				
	}
	//Angler interviews - CPUE
	for(a in 1:IntC){
		c[a] ~ neg_binomial_2(lambda_C_S[section_IntC[a]][day_IntC[a], gear_IntC[a]] * h[a] , r_C);
	}
	//Angler interviews - Angler expansions
	for(a in 1:IntA){
		//vehicles
		V_A[a] ~ binomial(A_A[a], R_V[gear_IntA[a]]);  //Note: leaving ratio of cars per angler constant among days since was invariant!
		//trailers
		T_A[a] ~ binomial(A_A[a], R_T[gear_IntA[a]]);  //Note: leaving ratio of cars per angler constant among days since was invariant!
	}												
}
generated quantities{ 
 	matrix<lower=0>[D,G] lambda_Ctot_S[S]; //total daily catch
	matrix<lower=0>[D,G] C[S]; //realized total daily catch
	matrix<lower=0>[D,G] E[S]; //realized total daily effort
	real<lower=0> C_sum; //season-total catch
 	real<lower=0> E_sum; //season-total effort
	vector[V_n + T_n + A_n + E_n + IntC + IntA + IntA] log_lik;
	C_sum = 0;
 	E_sum = 0;
	for(g in 1:G){
		for(d in 1:D){
			for(s in 1:S){
				lambda_Ctot_S[s][d,g] = lambda_E_S[s][d,g] * L[d] * lambda_C_S[s][d,g]; 
				C[s][d,g] = poisson_rng(lambda_Ctot_S[s][d,g]); 
 				C_sum = C_sum + C[s][d,g];
 				E[s][d,g] = lambda_E_S[s][d,g] * L[d]; 
 				E_sum = E_sum + E[s][d,g];
			}                                                                                                              
		}
	}
// 	//point-wise log likelihood for LOO-IC
// 	//Index effort counts - vehicles
	for (i in 1:V_n){
		// Matches the model block's actual likelihood (single g=1 term) -- the
		// g=2 term this file's log_lik previously summed here is out of bounds
		// whenever G=1: R_V/p_TI/lambda_E_S_I are all sized by G, not hardcoded
		// to 2 the way b is, so this crashed at runtime.
		log_lik[i] = poisson_lpmf(V_I[i]|(lambda_E_S_I[section_V[i],countnum_V[i]][day_V[i],1] * p_TI[1,section_V[i]] * R_V[1]) * b[1]);
	}
	//Index effort counts - trailers
	for(i in 1:T_n){
		// Same fix, matching the model block's actual T_I likelihood.
		log_lik[V_n +i] = poisson_lpmf(T_I[i]|(lambda_E_S_I[section_T[i],countnum_T[i]][day_T[i],1] * p_TI[1,section_T[i]] * R_T[1]) * b[2]); 
	}
	//Index effort counts - anglers
	for(i in 1:A_n){ 
		log_lik[V_n + T_n + i] = poisson_lpmf(A_I[i]|lambda_E_S_I[section_A[i],countnum_A[i]][day_A[i],gear_A[i]] * p_TI[gear_A[i],section_A[i]] * p_I[gear_A[i],section_A[i]]);
	}
	//Census (tie-in) effort counts - anglers
	for(e in 1:E_n){
		log_lik[V_n + T_n + A_n + e] = poisson_lpmf(E_s[e]|lambda_E_S_I[section_E[e],countnum_E[e]][day_E[e],gear_E[e]] * p_TI[gear_E[e],section_E[e]]);				
	}
	//Angler interviews - catch (number of fish)
	for(a in 1:IntC){
		log_lik[V_n + T_n + A_n + E_n + a] = neg_binomial_2_lpmf(c[a]|lambda_C_S[section_IntC[a]][day_IntC[a],gear_IntC[a]] * h[a] , r_C);
	}
	//Angler interviews - number of vehicles
	for(a in 1:IntA){
		log_lik[V_n + T_n + A_n + E_n + IntC + a] = binomial_lpmf(V_A[a]|A_A[a], R_V[gear_IntA[a]]);
	}
	//Angler interviews - number of trailers
	for(a in 1:IntA){
		log_lik[V_n + T_n + A_n + E_n + IntC + IntA + a] = binomial_lpmf(T_A[a]|A_A[a], R_T[gear_IntA[a]]);
	}												
}

