

source("/Users/danielrud/Desktop/USC/Targeted Learning/tlgxe/tmle_main_function.R")

# function for effect modification
# effect_modifier should be a numeric 0 1 variable 
TMLE_effect_mod = function(Y, A,effect_modifier, W_outcome = NULL, W_exposure = NULL, 
                           family = "binomial",
                           outcome_formula = NULL,
                           propensity_scores = NULL, 
                           case_control_design = F, disease_prevalence = NULL, 
                           obs.weights = rep(1,length(Y)),
                           TMLE_args_list = list(
                             outcome_method = c("glmnet_int", "glmnet", "gesso", "logistf", "SL"), 
                             npv_thresh = (5/sqrt(length(Y)))/log(length(Y)), 
                             near_positivity_method = c("trim", "rebound"), 
                             nfolds_cv_Q_init = 10, 
                             nfolds_cv_glmnet_outcome = 3,
                             alpha_outcome = .5, 
                             clever_cov_propensity_wt = T, 
                             outcome_SL.library = c("SL.glmnet"), 
                             outcome_SL.cvControl = list(V = 10L, 
                                                         stratifyCV = ifelse(family == "gaussian", F, T), 
                                                         shuffle = TRUE, 
                                                         validRows = NULL)), 
                           verbose = F)
{
  
  arg_names = c( 'outcome_method',
                 'npv_thresh',
                 'near_positivity_method',
                 'nfolds_cv_Q_init',
                 'nfolds_cv_glmnet_outcome',
                 'alpha_outcome',
                 'clever_cov_propensity_wt', 
                 'outcome_SL.library', 
                 'outcome_SL.cvControl')
  
  if(!identical(names(TMLE_args_list), arg_names))
  {
    stop(paste0("TMLE_args_list must be specified with elements for all included arguments in list. One argument may have been missing or misspelled.  Check to see if the following arguments 
               are included: ", paste0(arg_names, collapse = ", ")))
  }
  
  # put elements of TMLE_args_list in local environment
  list2env(TMLE_args_list,envir = environment())
  
  # match args 
  outcome_method = match.arg(outcome_method, 
                             choices = c("glmnet_int", "glmnet", "gesso", "logistf", "SL"))
  
  near_positivity_method = match.arg(near_positivity_method, c("trim", "rebound"))
  
  
  # for effect modification, we fit ATE on data where effect_modifier = 0, 
  # then on effect modifier = 1, then on effect modifier = 2, then process ACEs.  
  E0_indices = which(effect_modifier == 0)
  E1_indices = which(effect_modifier == 1)
  E2_indices = which(effect_modifier == 2)
  
  tmle_E0 = TMLE(Y = Y[E0_indices], A = A[E0_indices], 
                 W_outcome = if(is.null(W_outcome)){NULL}else{data.frame(W_outcome)[E0_indices, ]}, 
                 W_exposure = if(is.null(W_exposure)){NULL}else{data.frame(W_exposure)[E0_indices, ]}, 
                 family = family, 
                 outcome_formula = outcome_formula, 
                 case_control_design = case_control_design, 
                 disease_prevalence = disease_prevalence,
                 propensity_scores = propensity_scores[E0_indices],
                 obs.weights = obs.weights[E0_indices],
                 outcome_method = outcome_method, 
                 npv_thresh = npv_thresh, 
                 near_positivity_method = near_positivity_method, 
                 nfolds_cv_Q_init = nfolds_cv_Q_init, 
                 nfolds_cv_glmnet_outcome = nfolds_cv_glmnet_outcome,
                 alpha_outcome = alpha_outcome,
                 clever_cov_propensity_wt = clever_cov_propensity_wt, 
                 outcome_SL.library = outcome_SL.library, 
                 outcome_SL.cvControl = outcome_SL.cvControl)
  
  tmle_E1 = TMLE(Y = Y[E1_indices], A = A[E1_indices], 
                 W_outcome = if(is.null(W_outcome)){NULL}else{data.frame(W_outcome)[E1_indices, ]}, 
                 W_exposure = if(is.null(W_exposure)){NULL}else{data.frame(W_exposure)[E1_indices, ]}, 
                 family = family, 
                 outcome_formula = outcome_formula, 
                 case_control_design = case_control_design, 
                 propensity_scores = propensity_scores[E1_indices],
                 obs.weights = obs.weights[E1_indices],
                 disease_prevalence = disease_prevalence,
                 outcome_method = outcome_method, 
                 npv_thresh = npv_thresh, 
                 near_positivity_method = near_positivity_method, 
                 nfolds_cv_Q_init = nfolds_cv_Q_init, 
                 nfolds_cv_glmnet_outcome = nfolds_cv_glmnet_outcome,
                 alpha_outcome = alpha_outcome,
                 clever_cov_propensity_wt = clever_cov_propensity_wt, 
                 outcome_SL.library = outcome_SL.library, 
                 outcome_SL.cvControl = outcome_SL.cvControl)
  
  tmle_E2 = TMLE(Y = Y[E2_indices], A = A[E2_indices], 
                 W_outcome = if(is.null(W_outcome)){NULL}else{data.frame(W_outcome)[E2_indices, ]}, 
                 W_exposure = if(is.null(W_exposure)){NULL}else{data.frame(W_exposure)[E2_indices, ]}, 
                 family = family, 
                 outcome_formula = outcome_formula, 
                 case_control_design = case_control_design,
                 propensity_scores = propensity_scores[E2_indices],
                 obs.weights = obs.weights[E2_indices], 
                 disease_prevalence = disease_prevalence,
                 outcome_method = outcome_method, 
                 npv_thresh = npv_thresh, 
                 near_positivity_method = near_positivity_method, 
                 nfolds_cv_Q_init = nfolds_cv_Q_init, 
                 nfolds_cv_glmnet_outcome = nfolds_cv_glmnet_outcome,
                 alpha_outcome = alpha_outcome,
                 clever_cov_propensity_wt = clever_cov_propensity_wt, 
                 outcome_SL.library = outcome_SL.library, 
                 outcome_SL.cvControl = outcome_SL.cvControl)
  
  
  
  # Perform anova for difference in ATE by SNP level 
  
  n0 = length(E0_indices)
  n1 = length(E1_indices)
  n2 = length(E2_indices)
  
  var_0 = n0*tmle_E0$ATE_var
  var_1 = n1*tmle_E1$ATE_var
  var_2 = n2*tmle_E2$ATE_var
  
  xbar_0 = tmle_E0$ATE 
  xbar_1 = tmle_E1$ATE 
  xbar_2 = tmle_E2$ATE 
  
  n_vec = c(n0, n1, n2)
  xbar_vec = c(xbar_0, xbar_1, xbar_2)
  var_vec = c(var_0, var_1, var_2) 
  xbar_grand = n_vec %*% xbar_vec / sum(n_vec)
  # Setup for typical ANOVA ##############
  
  SSB = (n_vec) %*% (xbar_vec - rep(xbar_grand, 3))^2
  SSW = (n_vec - 1) %*% var_vec 
  
  MSB = SSB / (3-1)
  MSW = SSW / (sum(n_vec) - 3)
  
  ATE_EM_F_statistic = MSB / MSW
  
  ATE_EM_pvalue = pf(ATE_EM_F_statistic, df1 = 3 - 1,sum(n_vec) - 3 , lower.tail = F)
  
  # Perform anova for MOR results 
  
  var_0 = n0*tmle_E0$MOR_var
  var_1 = n1*tmle_E1$MOR_var
  var_2 = n2*tmle_E2$MOR_var
  
  xbar_0 = tmle_E0$MOR
  xbar_1 = tmle_E1$MOR
  xbar_2 = tmle_E2$MOR
  
  n_vec = c(n0, n1, n2)
  xbar_vec = c(xbar_0, xbar_1, xbar_2)
  var_vec = c(var_0, var_1, var_2)
  
  xbar_grand = n_vec %*% xbar_vec / sum(n_vec)
  
  # Setup for typical ANOVA ##############
  
  SSB = (n_vec) %*% (xbar_vec - rep(xbar_grand, 3))^2
  SSW = (n_vec - 1) %*% var_vec
  
  MSB = SSB / (3-1)
  MSW = SSW / (sum(n_vec) - 3)
  
  MOR_EM_F_statistic = MSB / MSW
  
  MOR_EM_pvalue = pf(MOR_EM_F_statistic, df1 = 3 - 1,sum(n_vec) - 3 , lower.tail = F)
  
  # 1 df linear ATE test 
  m0 = tmle_E0$ATE
  m1 = tmle_E1$ATE
  m2 = tmle_E2$ATE
  
  var_m0 = tmle_E0$ATE_var
  var_m1 = tmle_E1$ATE_var
  var_m2 = tmle_E2$ATE_var
  
  ATE_linear_results = linear_EM_test(m0 = m0, m1 = m1, m2 = m2, 
                                      var_m0 = var_m0, var_m1 = var_m1, var_m2 = var_m2)
  
  
  ATE_EM_lin_baseline_est = ATE_linear_results$EM_lin_baseline
  ATE_EM_lin_est = ATE_linear_results$EM_hat
  ATE_EM_lin_Z_stat = ATE_linear_results$Z_EM
  ATE_EM_lin_pvalue = ATE_linear_results$pvalue
  
  
  # EM for MOR 1 df in log scale  
  
  # apply delta method to variance estimates 
  # var* = var* 1/MOR^2
  
  m0 = tmle_E0$MOR %>% log
  m1 = tmle_E1$MOR %>% log
  m2 = tmle_E2$MOR %>% log
  
  var_0 = tmle_E0$MOR_var / (tmle_E0$MOR)^2
  var_1 = tmle_E1$MOR_var / (tmle_E1$MOR)^2
  var_2 = tmle_E2$MOR_var / (tmle_E2$MOR)^2
  
  MOR_linear_results = linear_EM_test(m0 = m0, m1 = m1, m2 = m2, 
                                      var_m0 = var_0, var_m1 = var_1, var_m2 = var_2)
  
  MOR_EM_mult_baseline_est = MOR_linear_results$EM_lin_baseline %>% exp
  MOR_EM_mult_est = MOR_linear_results$EM_hat %>% exp
  MOR_EM_mult_Z_stat = MOR_linear_results$Z_EM
  MOR_EM_mult_pvalue = MOR_linear_results$pvalue

  result_vector = c(ATE_E0 = tmle_E0$ATE %>% as.numeric, 
                    ATE_E1 = tmle_E1$ATE %>% as.numeric,
                    ATE_E2 = tmle_E2$ATE %>% as.numeric, 
                    var_ATE_E0 = tmle_E0$ATE_var %>% as.numeric, 
                    var_ATE_E1 = tmle_E1$ATE_var %>% as.numeric, 
                    var_ATE_E2 = tmle_E2$ATE_var %>% as.numeric, 
                    MOR_E0 = tmle_E0$MOR %>% as.numeric,
                    MOR_E1 = tmle_E1$MOR %>% as.numeric,
                    MOR_E2 = tmle_E2$MOR %>% as.numeric, 
                    var_MOR_E0 = tmle_E0$MOR_var %>% as.numeric,
                    var_MOR_E1 = tmle_E1$MOR_var %>% as.numeric,
                    var_MOR_E2 = tmle_E2$MOR_var %>% as.numeric, 
                    ATE_EM_F_statistic = ATE_EM_F_statistic,
                    ATE_EM_pvalue = ATE_EM_pvalue, 
                    MOR_EM_F_statistic = MOR_EM_F_statistic,
                    MOR_EM_pvalue = MOR_EM_pvalue, 
                    ATE_EM_lin_baseline_est = ATE_EM_lin_baseline_est,
                    ATE_EM_lin_est = ATE_EM_lin_est,
                    ATE_EM_lin_Z_stat = ATE_EM_lin_Z_stat,
                    ATE_EM_lin_pvalue = ATE_EM_lin_pvalue,
                    MOR_EM_mult_baseline_est = MOR_EM_mult_baseline_est,
                    MOR_EM_mult_est = MOR_EM_mult_est,
                    MOR_EM_mult_Z_stat = MOR_EM_mult_Z_stat,
                    MOR_EM_mult_pvalue =MOR_EM_mult_pvalue)
  
  return(result_vector)
}

# will be called for 1df ATE and 1df MOR tests
linear_EM_test = function(m0, m1, m2, var_m0, var_m1, var_m2)
{
  c02 = 1/var_m0 + 1/var_m2
  c12 = 1/var_m1 + 4/var_m2
  
  u1_hat = ((m1 / var_m1) + (2/var_m2)*(m2 + m0/(var_m0*c02) - m2/(var_m2*c02))) /
    (c12* (1- 4/(c12*c02*var_m2^2) ) )
  u0_hat = ((m0/var_m0) - (m2-2*u1_hat)/var_m2 ) / c02
  
  ## variance calculations:
  # u1_hat
  d1 = c12*(1-4/(c12*c02*(var_m2^2)))
  a10 = 2/(d1*var_m2*var_m0 * c02)
  a11 = 1/(d1*var_m1)
  a12 = (2/(d1*var_m2))*(1-(1/(var_m2 * c02)))
  
  var_u1_hat = a10^2*var_m0 + a11^2*var_m1 + a12^2*var_m2
  
  # u0_hat
  a00 = (1/c02)*(1/var_m0 + 2*a10 / var_m2)
  a01 = (1/c02)*(2*a11/var_m2)
  a02 = (1/c02)*(-1/var_m2 + 2*a12 / var_m2)
  
  var_u0_hat = a00^2*var_m0 +  a01^2*var_m1 + a02^2*var_m2
  
  # Z_u1 = (u1_hat) / sqrt(var_u1_hat)
  # Z_u0 = (u0_hat) / sqrt(var_u0_hat)
  
  # Effect Modification estimate
  
  EM_hat = u1_hat - u0_hat
  
  EM_lin_baseline = u0_hat
  
  cov_u1hat_u0hat = (c(a10,a11, a12)*c(a00, a01, a02)) %*% c(var_m0, var_m1, var_m2)
  
  var_u1hat_minus_u0hat = var_u1_hat + var_u0_hat -2*cov_u1hat_u0hat
  
  Z_EM = (EM_hat) / sqrt(var_u1hat_minus_u0hat) # under null, u1 - u0 = 0
  
  lin_EM_pvalue = 2*pnorm(-1*abs(Z_EM), mean = 0, sd = 1)
  
  return(list(Z_stat = Z_EM, pvalue =  lin_EM_pvalue, 
         EM_hat = EM_hat, EM_lin_baseline = EM_lin_baseline))
}

