
source("/Users/danielrud/Desktop/USC/Targeted Learning/TMLE_function/tmle_effect_mod_function.R")
source("/Users/danielrud/Desktop/USC/Targeted Learning/TMLE_function/tmle_gxe_post_process.R")

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
  num_G = ifelse(!is.null(num_G),num_G, length(G))
  
  # set up confounding frames W_outcome and W_exposure
  # only need to worry if W is specified
  
  if(include_G_propensity)
  {
    if(is.null(W))
    {
      W = G
    }else
    {
      W = cbind(W, G) 
    }
  }
  
  # initialize progressor
  p <- NULL
  if(progress == T)
  {
    p <- progressor(steps = num_G)
  }
  
  
  if(verbose == T)
  {
    message("Fitting Propensity model using Super Learner...")
  }
  
  # initialize observation weights 
  weights = numeric(length(Y))
  if(case_control_design && is.null(weights)) # if case control design and weights not specified 
  {
    outcome_table = table(Y)
    J = outcome_table[1] / outcome_table[2]
    # weight will be q0 for Y = 1 (cases) and (1-q0)*(1/J) for Y = 0 (controls)
    weights = Y*disease_prevalence  + (1-Y)*(1 - disease_prevalence)*(1/J)
    
  }else if (case_control_design && !is.null(weights)) # if case control study AND user weights specified
  {
    if(verbose)
    {
      message("Using prescpecified weights for case control design specified.  \nTo use default case control weights, set `weights = NULL`.")
    }
  }else
  {
    weights[1:length(weights)] = 1
  }
  
  # fit overall propensity model 
  
  propensity_scores = generate_propensity_SL(exposure_data = data.frame(if(is.null(W)){A}else{cbind(A, W)}), 
                                             weights = weights, 
                                             SL.library = SL.library, 
                                             SL.cvControl = SL.cvControl, 
                                             parallel = parallel, 
                                             ncores = ncores) 
  
  # include propensity scored in TMLE effect mod 
  
  plan(multisession, workers = cores)
  
  tmle_results = future_lapply(1:num_G, FUN = function(i)
  {
    
    W_outcome_curr = if(!is.null(W_outcome)){cbind(G[,-i], W_outcome)}else{G[,-i]}
    W_exposure_curr = if(!is.null(W_outcome)){cbind(G[,-i], W_exposure)}else{G[,-i]}
    
    effect_modifier = G[,i]
    
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
    if(progress == T)
    {
      p()
    }
    
    return(tmle_mod)
    
  }, future.seed = T)
  
  plan(sequential)
  
  result_frame = do.call(cbind, tmle_results)
  colnames(result_frame) = if(is.null(colnames(G))){paste0("SNP", 1:G)}else{colnames(G)}
  
  class(result_frame) = "tmle_gxe"
  
  return(result_frame)
}











