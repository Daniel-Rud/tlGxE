

# function for effect modification
# effect_modifier should be a numeric 0 1 variable
TMLE_effect_mod = function(Y, E,effect_modifier, W_outcome = NULL, W_exposure = NULL,
                           family = "binomial",
                           outcome_formula = NULL,
                           propensity_scores = NULL,
                           case_control_design = F, disease_prevalence = NULL,
                           obs.weights = rep(1,length(Y)),
                           TMLE_args_list = list(
                             outcome_method = c("glmnet", "glmnet_int", "gesso", "SL"),
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

  # put elements of TMLE_args_list in local environment
  list2env(TMLE_args_list,envir = environment())

  # match args
  if(length(outcome_method) >= 1)
  {
    outcome_method = "glmnet"
  }else
  {
    outcome_method = match.arg(outcome_method,
                               choices = c("glmnet_int", "glmnet", "gesso", "SL"))
  }

  near_positivity_method = match.arg(near_positivity_method, c("trim", "rebound"))


  # for effect modification, we fit ATE on data where effect_modifier = 0,
  # then on effect modifier = 1, then on effect modifier = 2, then process ACEs.
  E0_indices = which(effect_modifier == 0)
  E1_indices = which(effect_modifier == 1)
  E2_indices = which(effect_modifier == 2)

  # TMLE for SNP = 0
  tmle_E0 = TMLE(Y = Y[E0_indices], A = E[E0_indices],
                 W_outcome = if(is.null(W_outcome)){NULL}else{data.frame(W_outcome)[E0_indices, ]},
                 W_exposure = if(is.null(W_exposure)){NULL}else{data.frame(W_exposure)[E0_indices, ]},
                 family = family,
                 outcome_formula = outcome_formula,
                 case_control_design = case_control_design,
                 disease_prevalence = disease_prevalence,
                 propensity_scores = if(is.null(propensity_scores)) { NULL } else { propensity_scores[E0_indices] },
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
  # TMLE for SNP = 1
  tmle_E1 = TMLE(Y = Y[E1_indices], A = E[E1_indices],
                 W_outcome = if(is.null(W_outcome)){NULL}else{data.frame(W_outcome)[E1_indices, ]},
                 W_exposure = if(is.null(W_exposure)){NULL}else{data.frame(W_exposure)[E1_indices, ]},
                 family = family,
                 outcome_formula = outcome_formula,
                 case_control_design = case_control_design,
                 propensity_scores = if(is.null(propensity_scores)) { NULL } else { propensity_scores[E1_indices] },
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

  # TMLE for SNP = 2
  tmle_E2 = TMLE(Y = Y[E2_indices], A = E[E2_indices],
                 W_outcome = if(is.null(W_outcome)){NULL}else{data.frame(W_outcome)[E2_indices, ]},
                 W_exposure = if(is.null(W_exposure)){NULL}else{data.frame(W_exposure)[E2_indices, ]},
                 family = family,
                 outcome_formula = outcome_formula,
                 case_control_design = case_control_design,
                 propensity_scores = if(is.null(propensity_scores)) { NULL } else { propensity_scores[E2_indices] },
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

  ATE_welch_aov = welch_anova(means = c(tmle_E0$ATE, tmle_E1$ATE, tmle_E2$ATE),
                              var_means = c(tmle_E0$ATE_var, tmle_E1$ATE_var, tmle_E2$ATE_var),
                              sample_sizes = c(n0, n1, n2))

  # ATE_welch_aov = codominant_LRT(means = c(tmle_E0$ATE, tmle_E1$ATE, tmle_E2$ATE),
  #                             var_means = c(tmle_E0$ATE_var, tmle_E1$ATE_var, tmle_E2$ATE_var))


  ATE_EM_F_statistic = ATE_welch_aov$F_statistic

  ATE_EM_pvalue = ATE_welch_aov$p_value

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

  # only computer MOR stuff if binary outcome
  if(family == "binomial")
  {

    # Perform anova for MOR results

    MOR_welch_aov = welch_anova(means = c(tmle_E0$MOR, tmle_E1$MOR, tmle_E2$MOR),
                                var_means = c(tmle_E0$MOR_var, tmle_E1$MOR_var, tmle_E2$MOR_var),
                                sample_sizes = c(n0, n1, n2))

    # # trying LRT instead
    # MOR_welch_aov = codominant_LRT(means = c(tmle_E0$MOR, tmle_E1$MOR, tmle_E2$MOR),
    #                             var_means = c(tmle_E0$MOR_var, tmle_E1$MOR_var, tmle_E2$MOR_var))

    MOR_EM_F_statistic = MOR_welch_aov$F_statistic

    MOR_EM_pvalue = MOR_welch_aov$p_value


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
  }

  if(family == "binomial")
  {
    result_vector = c(ATE_G0 = tmle_E0$ATE %>% as.numeric,
                      ATE_G1 = tmle_E1$ATE %>% as.numeric,
                      ATE_G2 = tmle_E2$ATE %>% as.numeric,
                      var_ATE_G0 = tmle_E0$ATE_var %>% as.numeric,
                      var_ATE_G1 = tmle_E1$ATE_var %>% as.numeric,
                      var_ATE_G2 = tmle_E2$ATE_var %>% as.numeric,
                      MOR_G0 = tmle_E0$MOR %>% as.numeric,
                      MOR_G1 = tmle_E1$MOR %>% as.numeric,
                      MOR_G2 = tmle_E2$MOR %>% as.numeric,
                      var_MOR_G0 = tmle_E0$MOR_var %>% as.numeric,
                      var_MOR_G1 = tmle_E1$MOR_var %>% as.numeric,
                      var_MOR_G2 = tmle_E2$MOR_var %>% as.numeric,
                      ATE_codominant_F_statistic = ATE_EM_F_statistic,
                      ATE_codominant_pvalue = ATE_EM_pvalue,
                      MOR_codominant_F_statistic = MOR_EM_F_statistic,
                      MOR_codominant_pvalue = MOR_EM_pvalue,
                      ATE_additive_baseline_est = ATE_EM_lin_baseline_est,
                      ATE_additive_lin_est = ATE_EM_lin_est,
                      ATE_additive_Z_stat = ATE_EM_lin_Z_stat,
                      ATE_additive_pvalue = ATE_EM_lin_pvalue,
                      MOR_additive_baseline_est = MOR_EM_mult_baseline_est,
                      MOR_additive_mult_est = MOR_EM_mult_est,
                      MOR_additive_Z_stat = MOR_EM_mult_Z_stat,
                      MOR_additive_pvalue = MOR_EM_mult_pvalue)
  }else
  {
    result_vector = c(ATE_G0 = tmle_E0$ATE %>% as.numeric,
                      ATE_G1 = tmle_E1$ATE %>% as.numeric,
                      ATE_G2 = tmle_E2$ATE %>% as.numeric,
                      var_ATE_G0 = tmle_E0$ATE_var %>% as.numeric,
                      var_ATE_G1 = tmle_E1$ATE_var %>% as.numeric,
                      var_ATE_G2 = tmle_E2$ATE_var %>% as.numeric,
                      ATE_codominant_F_statistic = ATE_EM_F_statistic,
                      ATE_codominant_pvalue = ATE_EM_pvalue,
                      ATE_additive_baseline_est = ATE_EM_lin_baseline_est,
                      ATE_additive_lin_est = ATE_EM_lin_est,
                      ATE_additive_Z_stat = ATE_EM_lin_Z_stat,
                      ATE_additive_pvalue = ATE_EM_lin_pvalue)
  }

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

welch_anova = function(means, var_means, sample_sizes)
{
  # Convert variances of means to sample variances
  sample_variances = var_means * sample_sizes

  # Number of groups
  k = 3

  # Calculate weights and grand mean
  weights = sample_sizes / sample_variances
  grand_mean = sum(weights * means) / sum(weights)

  # Between-group sum of squares (SSB) and mean square (MSB)
  SSB = sum(weights * (means - grand_mean)^2)
  MSB = SSB / (k - 1)

  # Within-group mean square (MSW)

  MSW = 1 + (2*(k-2)/(k^2-1)) * sum((1/(sample_sizes - 1)) * (1 - weights/sum(weights))^2)

  # Calculate Welch's F-statistic
  F_stat = MSB / MSW

  # Degrees of freedom
  df1 = k - 1
  df2 = (k^2 - 1) / (3*sum((1 - weights / sum(weights))^2 / (sample_sizes - 1)))

  # Calculate p-value
  p_value = pf(F_stat, df1, df2, lower.tail = FALSE)

  # Return results
  return(list(F_statistic = F_stat, p_value = p_value, df1 = df1, df2 = df2))
}

# LRT based codominant effect modification test
codominant_LRT = function(means, var_means)
{
  mu0 = means[1]; mu1 = means[2]; mu2 = means[3]
  sigma0 = sqrt(var_means[1]); sigma1 = sqrt(var_means[2]); sigma2 = sqrt(var_means[3])

  # Log-likelihood under the alternative hypothesis
  # I write the value as 0 and the mean as 0 because under H1, we will have dnorm(m_i, m_i, sigma_0)
  logL_alt = dnorm(0, 0, sigma0, log = TRUE) +
    dnorm(0, 0, sigma1, log = TRUE) +
    dnorm(0, 0, sigma2, log = TRUE)

  # MLEs under the null hypothesis (common mean)
  sigma_vec = c(sigma1^2 *sigma2^2,sigma0^2 *sigma2^2, sigma0^2 *sigma1^2)
  mu_null = (sigma_vec %*% c(mu0, mu1, mu2)) / sum(sigma_vec)

  # Log-likelihood under the null hypothesis
  logL_null = dnorm(mu0, mu_null, sigma0, log = TRUE) +
    dnorm(mu1, mu_null, sigma1, log = TRUE) +
    dnorm(mu2, mu_null, sigma2, log = TRUE)

  # Likelihood ratio statistic
  LRT_stat = -2 * (logL_null - logL_alt)

  # Compare with chi-square distribution with 2 degrees of freedom
  p_value = pchisq(LRT_stat, df = 2, lower.tail = FALSE)

  # Return results
  return(list(F_statistic = LRT_stat, p_value = p_value))
}






