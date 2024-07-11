
source("/Users/danielrud/Desktop/USC/Targeted Learning/tlgxe/tmle_effect_mod_function.R")
source("/Users/danielrud/Desktop/USC/Targeted Learning/tlgxe/tmle_gxe_post_process.R")
source("/Users/danielrud/Desktop/USC/Targeted Learning/tlgxe/tmle_helper_functions.R")
source("/Users/danielrud/Desktop/USC/Targeted Learning/tlgxe/tmle_propensity_model_functions.R")
source("/Users/danielrud/Desktop/USC/Targeted Learning/tlgxe/tmle_outcome_model_functions.R")
source("/Users/danielrud/Desktop/USC/Targeted Learning/tlgxe/tmle_main_function.R")

tmle_gxe = function(Y, A, G, W = NULL, family = "binomial",
                    case_control_design = F, disease_prevalence = NULL,
                    weights = NULL, 
                    
                    propensity_SL.library = c("SL.gam", 
                                              "SL.glm.interaction", 
                                              "SL.glmnet", 
                                              "SL.ranger", 
                                              "SL.xgboost"), 
                    propensity_SL.cvControl = list(V = 10L, 
                                                   stratifyCV = T, 
                                                   shuffle = TRUE, 
                                                   validRows = NULL), 
                    include_G_propensity = T,
                    include_W_outcome = T, 
                    TMLE_args_list = list(
                      outcome_method = c("glmnet_int", "glmnet", "gesso", "logistf"),
                      npv_thresh = (5/sqrt(length(Y)))/log(length(Y)),
                      near_positivity_method = c("trim", "rebound"),
                      nfolds_cv_Q_init = 1,
                      nfolds_cv_glmnet_propensity = 3,
                      nfolds_cv_glmnet_outcome = 3,
                      alpha_outcome = .5,
                      alpha_propensity = .5,
                      clever_cov_propensity_wt = T),
                    parallel = T,
                    ncores = ifelse(parallel, availableCores(), 1), 
                    progress = T, 
                    verbose = T)
{
  # need to do all checks for function params
  # fill in later
  ###########################################
  
  # account for when only including a single G
  num_G = ncol(G)
  num_G = ifelse(!is.null(num_G),num_G, 1)
  
  
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

  if(case_control_design && is.null(weights)) # if case control design and weights not specified 
  {
    outcome_table = table(Y)
    J = outcome_table[1] / outcome_table[2]
    
    # weight will be q0 for Y = 1 (cases) and (1-q0)*(1/J) for Y = 0 (controls)
    # where q0 is the disease prevalence
    weights = Y*disease_prevalence  + (1-Y)*(1 - disease_prevalence)*(1/J)
    
  }else if (case_control_design && !is.null(weights)) # if case control study AND user weights specified
  {
    if(verbose)
    {
      message("Using prespecified weights for case control design specified.  \nTo use default case control weights, set `weights = NULL`.")
    }
  }else if(!is.null(weights))
  {
    if(verbose)
    {
      message("Using prespecified weights.")
    }
    
  }else # if weights are not prespecified and not case control design
  {
    weights = rep(1, length(Y))
  }
  #########################################################################
  # Fit overall propensity model ##########################################
  #########################################################################
  if(verbose == T)
  {
    message("Fitting Propensity model using Super Learner...")
  }
  
  propensity_scores = generate_propensity_SL(exposure_data = data.frame(if(is.null(W_propensity)){A}else{cbind(A, W_propensity)}), 
                                             weights = weights, 
                                             SL.library = propensity_SL.library, 
                                             SL.cvControl = propensity_SL.cvControl, 
                                             parallel = parallel, 
                                             ncores = ncores) 
  
  #########################################################################
  # Perform TMLE EM analysis for each SNP as an effect modifier ###########
  #########################################################################
  tmle_results = NULL
  
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
                                                   weights = weights,  
                                                   include_W_outcome = include_W_outcome, 
                                                   TMLE_args_list = TMLE_args_list, 
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
                                     weights = weights,  
                                     include_W_outcome = include_W_outcome, 
                                     TMLE_args_list = TMLE_args_list, 
                                     progress = progress, 
                                     parallel = parallel, 
                                     ncores = ncores)
  }
  
  
  result_frame = do.call(cbind, tmle_results)
  if(num_G > 1)
  {
    colnames(result_frame) = if(is.null(colnames(G))){paste0("SNP_", 1:G)}else{colnames(G)} 
  }
  
  class(result_frame) = "tmle_gxe"
  
  return(result_frame)
}

tmle_EM_iterator = function(Y, A, G, W = NULL,propensity_scores = NULL, family = "binomial",
                            case_control_design = F, disease_prevalence = NULL,
                            weights = NULL,  
                            include_W_outcome = T, 
                            TMLE_args_list = list(
                              outcome_method = c("glmnet_int", "glmnet", "gesso", "logistf"),
                              npv_thresh = (5/sqrt(length(Y)))/log(length(Y)),
                              near_positivity_method = c("trim", "rebound"),
                              nfolds_cv_Q_init = 1,
                              nfolds_cv_glmnet_propensity = 3,
                              nfolds_cv_glmnet_outcome = 3,
                              alpha_outcome = .5,
                              alpha_propensity = .5,
                              clever_cov_propensity_wt = T), 
                            progress = T, 
                            parallel = T, 
                            ncores = availableCores())
{
  
  p <- NULL
  # if G is a single column, we need to set num_G to 1
  num_G = ifelse(!is.null(ncol(G)),ncol(G), 1)
  if(progress)
  {
    p <- progressor(steps = num_G)
  }
  
  tmle_results = NULL

  if(parallel)
  {
    plan(multisession, workers = ncores)

    tmle_results = future_lapply(1:num_G, FUN = function(i)
    {

      W_outcome_curr = NULL
      # W_exposure_curr is set to NULL since propensity model is done for all data in beginning.
      # We pass the propensities to TMLE
      W_exposure_curr = NULL

      if(include_W_outcome && !is.null(W))
      {
        if(num_G>1)
        {
          W_outcome_curr = cbind(G[,-i], W)
        }else
        {
          W_outcome_curr = W
        }
      }else
      {
        if(num_G>1)
        {
          W_outcome_curr = G[,-i]
        }else
        {
          W_outcome_curr = NULL
        }
      }

      effect_modifier =  if(num_G>1){G[,i]}else{G}

      tmle_mod = TMLE_effect_mod(Y = Y,
                                 A = A,
                                 effect_modifier = effect_modifier,
                                 W_outcome = W_outcome_curr,
                                 W_exposure = W_exposure_curr,
                                 family = family,
                                 propensity_scores = propensity_scores,
                                 case_control_design = case_control_design,
                                 disease_prevalence = disease_prevalence,
                                 weights = weights,
                                 TMLE_args_list = TMLE_args_list)
      if(progress)
      {
        p()
      }

      return(tmle_mod)
    }, future.seed = 2024)

    plan(sequential)
  }else
  {
    tmle_results = lapply(1:num_G, FUN = function(i)
    {

      W_outcome_curr = NULL
      # W_exposure_curr is set to NULL since propensity model is done for all data in beginning.
      # We pass the propensities to TMLE
      W_exposure_curr = NULL
      
      if(include_W_outcome && !is.null(W))
      {
        if(num_G>1)
        {
          W_outcome_curr = cbind(G[,-i], W)
        }else
        {
          W_outcome_curr = W
        }
      }else
      {
        if(num_G>1)
        {
          W_outcome_curr = G[,-i]
        }else
        {
          W_outcome_curr = NULL
        }
      }
      
      effect_modifier =  if(num_G>1){G[,i]}else{G}

      tmle_mod = TMLE_effect_mod(Y = Y,
                                 A = A,
                                 effect_modifier = effect_modifier,
                                 W_outcome = W_outcome_curr,
                                 W_exposure = W_exposure_curr,
                                 family = family,
                                 propensity_scores = propensity_scores,
                                 case_control_design = case_control_design,
                                 disease_prevalence = disease_prevalence,
                                 weights = weights,
                                 TMLE_args_list = TMLE_args_list)
      if(progress)
      {
        p()
      }

      return(tmle_mod)
    })
  }

  
  
  return(tmle_results)
}











