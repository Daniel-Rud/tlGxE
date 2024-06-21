

source("/Users/danielrud/Desktop/USC/Targeted Learning/TMLE_function/tmle_main_function.R")

# function for effect modification
# effect_modifier should be a numeric 0 1 variable 
TMLE_effect_mod = function(Y, A,effect_modifier, W_outcome = NULL, W_exposure = NULL, family = "binomial",
                           propensity_scores = NULL, 
                           case_control_design = F, disease_prevalence = NULL, 
                           weights = rep(1,length(Y)),
                           TMLE_args_list = list(
                             outcome_method = c("glmnet_int", "glmnet", "gesso", "logistf"), 
                             npv_thresh = (5/sqrt(length(Y)))/log(length(Y)), 
                             near_positivity_method = c("trim", "rebound"), 
                             nfolds_cv_Q_init = 10, 
                             nfolds_cv_glmnet_propensity = 3, 
                             nfolds_cv_glmnet_outcome = 3,
                             alpha_outcome = .5, 
                             alpha_propensity = .5, 
                             clever_cov_propensity_wt = T
                           ))
{
  
  if(!identical(names(TMLE_args_list), c( 'outcome_method',
                                          'npv_thresh',
                                          'near_positivity_method',
                                          'nfolds_cv_Q_init',
                                          'nfolds_cv_glmnet_propensity',
                                          'nfolds_cv_glmnet_outcome',
                                          'alpha_outcome',
                                          'alpha_propensity',
                                          'clever_cov_propensity_wt')))
  {
    stop("TMLE_args_list must be specified with elements for all included arguments in list.  
         One argument may have been missing or misspelled.")
  }
  
  # put elements of TMLE_args_list in local environment
  list2env(TMLE_args_list,envir = environment())
  
  # match args 
  outcome_method = match.arg(outcome_method, 
                             choices = c("glmnet_int", "glmnet", "gesso", "logistf"))
  
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
                 case_control_design = case_control_design, 
                 disease_prevalence = disease_prevalence,
                 propensity_scores = propensity_scores,
                 weights = weights[E0_indices],
                 outcome_method = outcome_method, 
                 npv_thresh = npv_thresh, 
                 near_positivity_method = near_positivity_method, 
                 nfolds_cv_Q_init = nfolds_cv_Q_init, 
                 nfolds_cv_glmnet_propensity = nfolds_cv_glmnet_propensity, 
                 nfolds_cv_glmnet_outcome = nfolds_cv_glmnet_outcome,
                 alpha_outcome = alpha_outcome,
                 alpha_propensity = alpha_propensity, 
                 clever_cov_propensity_wt = clever_cov_propensity_wt)
  
  tmle_E1 = TMLE(Y = Y[E1_indices], A = A[E1_indices], 
                 W_outcome = if(is.null(W_outcome)){NULL}else{data.frame(W_outcome)[E1_indices, ]}, 
                 W_exposure = if(is.null(W_exposure)){NULL}else{data.frame(W_exposure)[E1_indices, ]}, 
                 family = family, 
                 case_control_design = case_control_design, 
                 propensity_scores = propensity_scores,
                 weights = weights[E1_indices],
                 disease_prevalence = disease_prevalence,
                 outcome_method = outcome_method, 
                 npv_thresh = npv_thresh, 
                 near_positivity_method = near_positivity_method, 
                 nfolds_cv_Q_init = nfolds_cv_Q_init, 
                 nfolds_cv_glmnet_propensity = nfolds_cv_glmnet_propensity, 
                 nfolds_cv_glmnet_outcome = nfolds_cv_glmnet_outcome,
                 alpha_outcome = alpha_outcome,
                 alpha_propensity = alpha_propensity, 
                 clever_cov_propensity_wt = clever_cov_propensity_wt)
  
  tmle_E2 = TMLE(Y = Y[E2_indices], A = A[E2_indices], 
                 W_outcome = if(is.null(W_outcome)){NULL}else{data.frame(W_outcome)[E2_indices, ]}, 
                 W_exposure = if(is.null(W_exposure)){NULL}else{data.frame(W_exposure)[E2_indices, ]}, 
                 family = family, 
                 case_control_design = case_control_design, 
                 disease_prevalence = disease_prevalence,
                 propensity_scores = propensity_scores,
                 weights = weights[E2_indices],
                 outcome_method = outcome_method, 
                 npv_thresh = npv_thresh, 
                 near_positivity_method = near_positivity_method, 
                 nfolds_cv_Q_init = nfolds_cv_Q_init, 
                 nfolds_cv_glmnet_propensity = nfolds_cv_glmnet_propensity, 
                 nfolds_cv_glmnet_outcome = nfolds_cv_glmnet_outcome,
                 alpha_outcome = alpha_outcome,
                 alpha_propensity = alpha_propensity, 
                 clever_cov_propensity_wt = clever_cov_propensity_wt)
  
  
  
  
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
  
  ########################################
  
  # Setup for Welch Unequal Variance ANOVA ######
  
  # w = n_vec / var_vec
  # 
  # w_sum = sum(w)
  # 
  # xbar_grand = as.numeric((w %*% xbar_vec) / w_sum)
  # 
  # F_numerator = (1/(3-1)) * w %*% (xbar_vec - xbar_grand)^2
  # 
  # F_denom = 1 + 2*((3-2) / (3^2 - 1)) * (1 / (n_vec - 1)) %*% (1 - (w / w_sum))^2
  # 
  # df2 = (3^2 - 1) / (3* (1 / (n_vec - 1)) %*% (1 - w / w_sum))
  # 
  # ATE_EM_F_statistic = F_numerator / F_denom
  # 
  # ATE_EM_pvalue = pf(ATE_EM_F_statistic, df1 = 3 - 1, df2, lower.tail = F)
  
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
  
  # First, compute pooled EM estimator
  m0 = tmle_E0$ATE; m1 = tmle_E1$ATE; m2 = tmle_E2$ATE
  var_m0 = tmle_E0$ATE_var; var_m1 = tmle_E1$ATE_var; var_m2 = tmle_E2$ATE_var
  
  
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
  
  # CI_2.5 = EM_hat - qnorm(0.975, mean = 0, sd = 1) * sqrt(var_u1hat_minus_u0hat)
  # CI_97.5 = EM_hat + qnorm(0.975, mean = 0, sd = 1) * sqrt(var_u1hat_minus_u0hat)
  
  #EM Test for MOR 1 df 
  
  var_0 = tmle_E0$MOR_var
  var_1 = tmle_E1$MOR_var
  var_2 = tmle_E2$MOR_var
  
  m0 = tmle_E0$MOR
  m1 = tmle_E1$MOR
  m2 = tmle_E2$MOR
  
  initial_estimates = c(m0, 
                        weighted.mean(x = c(m1/m0, sqrt(m2/m0)),
                                      w = 1/c(var_1, var_2))
  )
  
  u0_u1_mles = nleqslv(x = initial_estimates, 
                       fn = grad_func_MOR, 
                       method = "Newton", 
                       jac = MOR_mle_jacobian,
                       var0 = var_0, 
                       var1 = var_1, 
                       var2 = var_2, 
                       m0 = m0, 
                       m1 = m1, 
                       m2 = m2, 
                       jacobian = F, 
                       control = list(
                         maxit = 10000, 
                         ftol = 1E-14
                       ))
  
  u1_est = u0_u1_mles$x[2]
  
  # compute fisher info for MLEs
  fisher_info = MOR_mle_fisher_info(u =u0_u1_mles$x,
                                    var0 = var_0, 
                                    var1 = var_1, 
                                    var2 = var_2, 
                                    m0 = m0, 
                                    m1 = m1, 
                                    m2 = m2)
  var_est_u1 = fisher_info[2,2]
  
  Z_stat_u1 = (u1_est - 1) / sqrt(var_est_u1) # null is that MOR is 1
  
  pval_u1 = 2*pnorm(abs(Z_stat_u1), lower = F)
  
  if(is.na(pval_u1)) # incase there is an NA for pval_u1 -- indicates fisher info was negative
  {
    save_obj = list(MOR = c(m0, m1, m2), 
                    VAR = c(var_0, var_1, var_2),
                    initial_estimates = initial_estimates,
                    u0_u1_mles = u0_u1_mles, 
                    fisher_info= fisher_info, 
                    var_est_u1 = var_est_u1)
    
    saveRDS(save_obj, file = "/Users/danielrud/Desktop/USC/Targeted Learning/debug_tmle_MOR_mult_list.RDS")
    
  }
  
  
  result_vector = c(ATE_E2 = tmle_E2$ATE %>% as.numeric, 
                    ATE_E1 = tmle_E1$ATE %>% as.numeric, 
                    ATE_E0 = tmle_E0$ATE %>% as.numeric,
                    var_ATE_E2 = tmle_E2$ATE_var %>% as.numeric, 
                    var_ATE_E1 = tmle_E1$ATE_var %>% as.numeric, 
                    var_ATE_E0 = tmle_E0$ATE_var %>% as.numeric,
                    MOR_E2 = tmle_E2$MOR %>% as.numeric, 
                    MOR_E1 = tmle_E1$MOR %>% as.numeric,
                    MOR_E0 = tmle_E0$MOR %>% as.numeric,
                    var_MOR_E2 = tmle_E2$MOR_var %>% as.numeric, 
                    var_MOR_E1 = tmle_E1$MOR_var %>% as.numeric, 
                    var_MOR_E0 = tmle_E0$MOR_var %>% as.numeric,
                    ATE_EM_F_statistic = ATE_EM_F_statistic,
                    ATE_EM_pvalue = ATE_EM_pvalue, 
                    MOR_EM_F_statistic = MOR_EM_F_statistic,
                    MOR_EM_pvalue = MOR_EM_pvalue, 
                    ATE_EM_lin_est = EM_hat, 
                    ATE_EM_lin_baseline_est = EM_lin_baseline,
                    ATE_EM_lin_pvalue = lin_EM_pvalue, 
                    MOR_EM_mult_Z_stat = Z_stat_u1, 
                    MOR_EM_mult_pvalue = pval_u1, 
                    MOR_EM_mult_baseline = u0_u1_mles$x[1], 
                    MOR_EM_mult_factor = u1_est)
  
  return(result_vector)
}