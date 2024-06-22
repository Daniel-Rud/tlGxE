
library(glmnet)
library(gesso)
library(logistf)
library(caret)
library(ggpubr)
library(latex2exp)
library(SuperLearner)
library(parallel)
library(future)
library(future.apply)

source("/Users/danielrud/Desktop/USC/Targeted Learning/TMLE_function/tmle_helper_functions.R")
source("/Users/danielrud/Desktop/USC/Targeted Learning/TMLE_function/tmle_propensity_model_functions.R")
source("/Users/danielrud/Desktop/USC/Targeted Learning/TMLE_function/tmle_outcome_model_functions.R")



# npv_thresh stands for near positivity violation threshold for propensity model 
# clever_cov_propensity_wt will put the denominator of the clever covariate as a weight
#                          in glm instead of in the covariate 
# for continuous outcome, uses logistic fluctuation 
# we always use logistic fluctuation, no linear fluctuation 
TMLE = function(Y, A, W_outcome = NULL, W_exposure = NULL, 
                family = c("gaussian", "binomial"),
                case_control_design = F, 
                disease_prevalence = NULL, 
                propensity_scores = NULL,
                weights = NULL,
                outcome_method = c("glmnet_int", "glmnet", "gesso", "logistf"), 
                npv_thresh = (5/sqrt(length(Y)))/log(length(Y)), 
                near_positivity_method = c("trim", "rebound"), 
                nfolds_cv_Q_init = 10,
                nfolds_cv_glmnet_propensity = 3, 
                nfolds_cv_glmnet_outcome = 3, 
                alpha_outcome = .5, alpha_propensity = .5, 
                clever_cov_propensity_wt = T, 
                propensity_SL.library = c("SL.glmnet"))
{
  
  # match args 
  near_positivity_method = match.arg(near_positivity_method, 
                                     choices = c("trim", "rebound"))
  
  # if case control design, check to make sure disease prevalence is specified
  if(case_control_design && is.null(disease_prevalence))
  {
    stop("For case-control study design, must supply estimate of disease prevalence.  Please specify argument `disease_prevalence`.")
  }
  
  #######################################################
  # Data Preprocessing 
  #######################################################
  
  # if Y is continuous, rebound and use logistic fluctuation 
  # Y_star will remain at Y if binomial outcome 
  Y_star = Y
  Y_bounds = NULL
  
  # if continuous outcome, use logistic fluctuation, rebound Y to Y_star in (0,1)
  if(family == "gaussian")
  {
    Y_bounds = range(Y)
    Y_star = rebound(Y, 0, 1)
  }
  
  # create dataframe to be used for outcome model 
  outcome_data = data.frame(Y = Y_star, A = A)
  if(!is.null(W_outcome)) # if there are covariates for propensity model 
  {
    outcome_data = cbind(outcome_data, W_outcome)
  }
  
  # create dataframe to be used for propensity model
  exposure_data = data.frame(A = A)
  if(!is.null(W_exposure)) # if there are covariates for propensity model 
  {
    exposure_data = cbind(exposure_data, W_exposure)
  }
  
  # establish weights for potential case control design 
  
  if(is.null(weights)) # if weights are not predefined
  {
    weights = numeric(length(Y_star))
    
    # for cohort/cross sectional design
    if(!case_control_design)
    {
      weights[1:length(weights)] = 1
    }else # for case control design 
    {
      outcome_table = table(Y_star)
      J = outcome_table[1] / outcome_table[2]
      # weight will be q0 for Y = 1 (cases) and (1-q0)*(1/J) for Y = 0 (controls)
      weights = Y_star*disease_prevalence  + (1-Y_star)*(1 - disease_prevalence)*(1/J)
      
    }
  }
  #######################################################
  # STEP 1: Propensity Score Model -- I change order for dealing with 
  #         near positivity violators
  
  PS = NULL
  if(!is.null(propensity_scores)) # if propensity scores provided to function 
  {
    PS = propensity_scores
  }else # otherwise compute propensity scores 
  {
    PS = generate_propensity_SL(exposure_data = exposure_data, 
                                               weights = weights, 
                                               SL.library = propensity_SL.library, 
                                               parallel = F, 
                                               ncores = NULL) 
  }
  
  # trimming and rebound is only for lower endpoints 
  # this is because propensity is in denominator, low 
  # propensities can explode the EIC computation
  if(near_positivity_method == "rebound")
  {
    # rebound lower and upper enda of interval 
    PS = bound_limits(PS,npv_thresh, 1-npv_thresh)
    
  }else # if near_positivity_method = "trim"
  {
    trim_indices = which(PS < npv_thresh | PS > 1 - npv_thresh)
    if(length(trim_indices) != 0)
    {
      PS = PS[-trim_indices]
      outcome_data = outcome_data[-trim_indices, ]
      A = A[-trim_indices]
      Y = Y[-trim_indices]
      Y_star = Y_star[-trim_indices]
      weights = weights[-trim_indices]
    }
  }
  
  #######################################################
  # STEP 2: G computation for initial estimates
  #######################################################
  # binomial and contionuous outcome
  suppressWarnings(Q0 <- generate_Q0_cv(outcome_data = outcome_data, family = family, 
                                        alpha = alpha_outcome, lambda = "lambda.min",
                                        nfolds_cv_Q_init = nfolds_cv_Q_init, nfolds_cv_glmnet = nfolds_cv_glmnet_outcome, 
                                        outcome_method = outcome_method, weights = weights))
  #######################################################
  # STEP 3: Clever Covariate
  #######################################################
  # note: for each observation, only one of H_1W or 
  # H_0W will be nonzero!
  
  epsilon = NULL
  # note again that we are always using logistic fluctuation, so 
  # we always treat offset like this
  
  offset = qlogis(Q0$Q0AW)
  
  #if we want to use regular clever covariate where propensity is in denominator 
  if(!clever_cov_propensity_wt)
  {
    
    H_1W = A / PS
    H_0W = (1-A) / (1-PS)
    
    # estimate epsilons 
    # again, family is binomial for logistic fluctuation
    suppressWarnings(epsilon <- glm(Y_star ~ -1 + H_1W + H_0W,
                                    offset = offset,
                                    family = binomial, 
                                    weights = weights) %>% coef)
  }else
    # if we choose to put propensity as a weight like IPTW instead of the usual clever covariate
  {
    
    clev_cov_weights = A / PS + (1-A) / (1-PS)
    H_1W = A
    H_0W = (1-A)
    
    # estimate epsilons
    # again, family is binomial for logistic fluctuation
    suppressWarnings(epsilon <- glm(Y_star ~ -1 + H_1W + H_0W,
                                    offset = offset,
                                    weights = clev_cov_weights*weights, # multiply clever covariate weights and original case control weights
                                    family = binomial) %>% coef)
  }
  #######################################################
  # STEP 4: Update Q0AW to Q1AW
  #######################################################
  
  Q1_1W = numeric(length(Y_star))
  Q1_0W = numeric(length(Y_star))
  if(!clever_cov_propensity_wt)
  {
    Q1_1W =  plogis(qlogis(Q0$Q1W) + epsilon[1] / PS)
    Q1_0W = plogis(qlogis(Q0$Q0W) + epsilon[2] / (1-PS))
  }else
  {
    Q1_1W =  plogis(qlogis(Q0$Q1W) + epsilon[1])
    Q1_0W = plogis(qlogis(Q0$Q0W) + epsilon[2])
  }
  
  # remap from (0,1) range to Y_bounds range if continuous outcome
  if(family == "gaussian")
  {
    a = Y_bounds[1]
    b = Y_bounds[2]
    Q1_1W = Q1_1W * (b-a) + a
    Q1_0W = Q1_0W * (b-a) + a
  }
  
  #######################################################
  # Step 5: Compute ATE and MOR 
  #######################################################
  
  # use weighted mean, only different for case control design
  EY1_tmle = weighted.mean(Q1_1W, weights)
  EY0_tmle = weighted.mean(Q1_0W, weights)
  
  ATE_tmle = EY1_tmle - EY0_tmle
  MOR_tmle = (EY1_tmle / (1-EY1_tmle)) / (EY0_tmle / (1-EY0_tmle))
  
  #######################################################
  # Step 6: Inference for Estimates 
  # ATE EIC
  
  ### NEED TO FIGURE OUT WEIGHTS FOR EIC!!!!
  D1 = ((A / PS) * (Y - Q1_1W) + Q1_1W)*weights - EY1_tmle
  D0 = (((1-A) / (1-PS)) * (Y - Q1_0W) + Q1_0W)*weights - EY0_tmle
  
  EIC_ATE = D1 - D0
  
  Var_EIC_ATE = var(EIC_ATE) / length(Y)
  
  qZ = qnorm(.975, mean = 0, sd = 1)
  
  ATE_CI = c("2.5%" = ATE_tmle - qZ*sqrt(Var_EIC_ATE) , 
             "97.5%" = ATE_tmle + qZ*sqrt(Var_EIC_ATE))
  ATE_pvalue = 2*(1-pnorm(abs(ATE_tmle / sqrt(Var_EIC_ATE))))
  
  # MOR EIC 
  
  EIC_MOR = ((1 - EY0_tmle) / (EY0_tmle * (1-EY1_tmle)^2)) * D1 -
    (EY1_tmle / ((1 - EY1_tmle) * EY0_tmle^2))* D0
  Var_EIC_MOR = var(EIC_MOR) / length(Y)
  
  MOR_CI = c("2.5%" = MOR_tmle - qZ*sqrt(Var_EIC_MOR), 
             "97.5%" = MOR_tmle + qZ*sqrt(Var_EIC_MOR))
  MOR_pvalue = 2*(1-pnorm(abs((MOR_tmle - 1) / sqrt(Var_EIC_MOR))))
  
  # include untargeted MOR 
  
  m1 = 0
  m0 = 0
  if(family == "gaussian")
  {
    a = Y_bounds[1]
    b = Y_bounds[2]
    m1 = mean(Q0$Q1W) * (b-a) + a
    m0 = mean(Q0$Q0W) * (b-a) + a
  }else
  {
    m1 = weighted.mean(Q0$Q1W, weights)
    m0 = weighted.mean(Q0$Q0W, weights)
  }
  
  ATE_untargeted = m1 - m0
  MOR_untargeted = (m1 / (1-m1)) / (m0 / (1-m0))
  
  
  result_list = list(
    ATE = ATE_tmle, 
    ATE_var = Var_EIC_ATE, 
    ATE_CI_2.5 = ATE_CI[1], 
    ATE_CI_97.5 = ATE_CI[2],
    ATE_pvalue = ATE_pvalue, 
    ATE_untargeted = ATE_untargeted,
    MOR = MOR_tmle, 
    MOR_var = Var_EIC_MOR, 
    MOR_CI_2.5 = MOR_CI[1], 
    MOR_CI_97.5 = MOR_CI[2],
    MOR_pvalue = MOR_pvalue, 
    MOR_untargeted = MOR_untargeted
  )
  
  return(result_list)
}