
library(progressr)
source("/Users/danielrud/Desktop/USC/Targeted Learning/tlgxe/tmle_effect_mod_function.R")
source("/Users/danielrud/Desktop/USC/Targeted Learning/tlgxe/tmle_gxe_post_process.R")
source("/Users/danielrud/Desktop/USC/Targeted Learning/tlgxe/tmle_helper_functions.R")
source("/Users/danielrud/Desktop/USC/Targeted Learning/tlgxe/tmle_propensity_model_functions.R")
source("/Users/danielrud/Desktop/USC/Targeted Learning/tlgxe/tmle_outcome_model_functions.R")
source("/Users/danielrud/Desktop/USC/Targeted Learning/tlgxe/tmle_main_function.R")

suppress_output <- function(expr) {
  sink(tempfile())
  on.exit(sink())  # Ensure sink is reset even if there's an error
  invisible(capture.output(expr))
}

# propensity formula and outcome formula should be of form:
# A ~ X1 + ...
# Y ~ A + X1
# model matrix will be used to fit them 
# refer to the outcome with Y and the exposure with A 

tmle_gxe = function(Y, A, G, W = NULL, family = "binomial",
                    case_control_design = F, disease_prevalence = NULL,
                    obs.weights = NULL, 
                    propensity_scores = NULL, 
                    TMLE_args_list = list(
                      outcome_method = c("glmnet_int", "glmnet", "gesso", "logistf", "SL"),
                      npv_thresh = (5/sqrt(length(Y)))/log(length(Y)),
                      near_positivity_method = c("trim", "rebound"),
                      nfolds_cv_Q_init = 10,
                      nfolds_cv_glmnet_outcome = 3,
                      alpha_outcome = .5,
                      clever_cov_propensity_wt = T, 
                      outcome_SL.library = c("SL.glmnet", 
                                             "SL.rpart", 
                                             "SL.bartMachine"), 
                      outcome_SL.cvControl = list(V = 10L, 
                                                  stratifyCV = ifelse(family == "gaussian", F, T), 
                                                  shuffle = TRUE, 
                                                  validRows = NULL)), 
                    propensity_SL.library = c("SL.glmnet", 
                                              "SL.rpart", 
                                              "SL.bartMachine"), 
                    propensity_SL.cvControl = list(V = 10L, 
                                                   stratifyCV = T, 
                                                   shuffle = TRUE, 
                                                   validRows = NULL), 
                    include_G_propensity = T,
                    include_W_outcome = T, 
                    propensity_formula = NULL, 
                    outcome_formula = NULL,
                    SNP_results = 1:(ifelse(is.null(ncol(G)), 1, ncol(G))),
                    parallel = T,
                    ncores = ifelse(parallel, availableCores(), 1), 
                    progress = T, 
                    verbose = T)
{
  # need to do all checks for function params
  # fill in later
  ###########################################
  
  # account for when only including a single G
  num_G = length(SNP_results)
  
  
  # set up W_propensity for propensity model 
  W_propensity = W
  if(include_G_propensity)
  {
    if(is.null(W_propensity))
    {
      W_propensity = G
    }else
    {
      W_propensity = cbind(W_propensity, G) 
    }
  }
  
  # initialize observation weights 

  if(case_control_design && is.null(obs.weights)) # if case control design and weights not specified 
  {
    outcome_table = table(Y)
    J = outcome_table[1] / outcome_table[2]
    
    # weight will be q0 for Y = 1 (cases) and (1-q0)*(1/J) for Y = 0 (controls)
    # where q0 is the disease prevalence
    obs.weights = Y*disease_prevalence  + (1-Y)*(1 - disease_prevalence)*(1/J)
    
  }else if (case_control_design && !is.null(obs.weights)) # if case control study AND user weights specified
  {
    if(verbose)
    {
      message("Using prespecified weights for case control design specified.  \nTo use default case control weights, set `obs.weights = NULL`.")
    }
  }else if(!is.null(obs.weights))
  {
    if(verbose)
    {
      message("Using prespecified weights.")
    }
    
  }else # if weights are not prespecified and not case control design
  {
    obs.weights = rep(1, length(Y))
  }
  
  #########################################################################
  # Fit overall propensity model ##########################################
  #########################################################################
  if(is.null(propensity_scores)) # if no propensity scores are supplied 
  {
    
    if(!is.null(propensity_formula))
    {
      
      full_exposure_data = cbind(A = A, G)
      if(!is.null(W))
      {
        full_exposure_data = cbind(full_exposure_data, W)
      }
      full_exposure_data = full_exposure_data %>% data.frame()
      
      # this line is needed because glm was otherwise giving problems 
      # with the scoping of the obs.weight arg.  Lot of headache....
      environment(propensity_formula) = environment()
      
      propensity_scores = glm(formula = propensity_formula, data = full_exposure_data, 
                              weights = obs.weights, family = "binomial") %>% predict(type = "response")
      
    }else
    {
      if(verbose == T)
      {
        message("Fitting Propensity model using Super Learner...")
      }
      suppress_output(propensity_scores <- generate_propensity_SL(exposure_data = data.frame(if(is.null(W_propensity)){A}else{cbind(A, W_propensity)}), 
                                                                  obs.weights = obs.weights, 
                                                                  SL.library = propensity_SL.library, 
                                                                  SL.cvControl = propensity_SL.cvControl, 
                                                                  parallel = parallel, 
                                                                  ncores = ncores))
    }
  }
  
  #########################################################################
  # Perform TMLE EM analysis for each SNP as an effect modifier ###########
  #########################################################################
  tmle_results = NULL
  
  if(verbose == T)
  {
    message("Fitting tmle_GxE models...")
  }
  
  if(progress == T)
  {
    with_progress(tmle_results <- tmle_EM_iterator(Y = Y, 
                                                   A = A,
                                                   G = G, 
                                                   W = W,
                                                   propensity_scores = propensity_scores, 
                                                   family = family,
                                                   case_control_design = case_control_design, 
                                                   disease_prevalence = disease_prevalence,
                                                   obs.weights = obs.weights,  
                                                   include_W_outcome = include_W_outcome, 
                                                   outcome_formula = outcome_formula, 
                                                   TMLE_args_list = TMLE_args_list, 
                                                   SNP_results = SNP_results,
                                                   progress = progress, 
                                                   parallel = parallel, 
                                                   ncores = ncores))
  }else
  {
    tmle_results <- tmle_EM_iterator(Y = Y, 
                                     A = A,
                                     G = G, 
                                     W = W,
                                     propensity_scores = propensity_scores, family = family,
                                     case_control_design = case_control_design, disease_prevalence = disease_prevalence,
                                     obs.weights = obs.weights,  
                                     include_W_outcome = include_W_outcome, 
                                     outcome_formula = outcome_formula,
                                     TMLE_args_list = TMLE_args_list, 
                                     SNP_results = SNP_results,
                                     progress = progress, 
                                     parallel = parallel, 
                                     ncores = ncores)
  }
  
  
  result_frame = do.call(cbind, tmle_results)
  if(num_G > 1)
  {
    colnames(result_frame) = if(is.null(colnames(G))){paste0("SNP_", SNP_results)}else{colnames(G)[SNP_results]} 
  }
  
  class(result_frame) = "tmle_gxe"
  
  return(result_frame)
}

tmle_EM_iterator = function(Y, A, G, W = NULL,propensity_scores = NULL, family = "binomial",
                            case_control_design = F, disease_prevalence = NULL,
                            obs.weights = NULL,  
                            include_W_outcome = T, 
                            outcome_formula = NULL, 
                            TMLE_args_list = list(
                              outcome_method = c("glmnet_int", "glmnet", "gesso", "logistf"),
                              npv_thresh = (5/sqrt(length(Y)))/log(length(Y)),
                              near_positivity_method = c("trim", "rebound"),
                              nfolds_cv_Q_init = 10,
                              nfolds_cv_glmnet_propensity = 3,
                              nfolds_cv_glmnet_outcome = 3,
                              alpha_outcome = .5,
                              clever_cov_propensity_wt = T, 
                              outcome_SL.library = c("SL.glmnet", 
                                                     "SL.rpart", 
                                                     "SL.bartMachine"), 
                              outcome_SL.cvControl = list(V = 10L, 
                                                          stratifyCV = ifelse(family == "gaussian", F, T), 
                                                          shuffle = TRUE, 
                                                          validRows = NULL)), 
                            SNP_results = 1:(ifelse(is.null(ncol(G)), 1, ncol(G))),
                            progress = T, 
                            parallel = T, 
                            ncores = availableCores())
{
  
  p <- NULL
  
  num_G = length(SNP_results)
  if(progress)
  {
    p <- progressor(steps = num_G)
  }
  
  tmle_results = NULL

  if(parallel)
  {
    plan(multisession, workers = ncores)
  }

    tmle_results = future_lapply(SNP_results, FUN = function(i)
    {

      W_outcome_curr = NULL
      # W_exposure_curr is set to NULL since propensity model is done for all data in beginning.
      # We pass the propensities to TMLE
      W_exposure_curr = NULL

      if(include_W_outcome && !is.null(W))
      {
        if(!is.null(ncol(G)))
        {
          W_outcome_curr = cbind(G[,-i], W)
        }else
        {
          W_outcome_curr = W
        }
      }else
      {
        if(!is.null(ncol(G)))
        {
          W_outcome_curr = G[,-i]
        }else
        {
          W_outcome_curr = NULL
        }
      }

      effect_modifier =  if(!is.null(ncol(G))){G[,i]}else{G}
      
      # round the effect modifier SNP in case the data is supplied as imputed dosages
      effect_modifier = round(effect_modifier)
      
      # if the current SNP does not have representation for all of SNP = 0,1,2
      if(length(unique(effect_modifier)) != 3)
      {
        message("SNP ",i, " does not have representation of all levels {0, 1, 2} and is therefore excluded from effect modification testing.")
        return(NULL)
      }else
      {
        tmle_mod = TMLE_effect_mod(Y = Y,
                                   A = A,
                                   effect_modifier = effect_modifier,
                                   W_outcome = W_outcome_curr,
                                   W_exposure = W_exposure_curr,
                                   family = family,
                                   outcome_formula = outcome_formula,
                                   propensity_scores = propensity_scores,
                                   case_control_design = case_control_design,
                                   disease_prevalence = disease_prevalence,
                                   obs.weights = obs.weights,
                                   TMLE_args_list = TMLE_args_list)
        if(progress)
        {
          p()
        }
        
        return(tmle_mod)
      }
    }, future.seed = 2024)
    
    if(parallel)
    {
      plan(sequential)
    }
  
  return(tmle_results)
}











