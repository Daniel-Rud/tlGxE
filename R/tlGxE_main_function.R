

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

#'
#'Performs tlGxE scan to estimate causal effect modification over a set of candidate SNPs.
#'
#' @param Y Outcome vector, should be numeric. Either \eqn{\{0,1\}} vector for \code{family = "binomial"} or a numeric vector for
#' continuous outcome for \code{family = "gaussian"}
#' @param E Binary exposure vector, generally a binary \eqn{\{0,1\}} vector.
#' @param G Data Frame of size \eqn{n \times g}, where \eqn{g} is the number of candidate SNPs.
#' \code{G} should contain either the genotyped SNPs in an additive allele encoding
#' (ie, SNPs data should be 0, 1, or 2 of minor allele) or in an imputed dosage format (probabilistic minor allele count).
#' SNPs included in the analysis should have adaquete representation across the three levels of minor allele count.  If representation
#' is insufficient for a SNP, results for the SNP are not outputted.
#' @param W Data Frame of \eqn{n \times p}, where p is the number of confounders/adjusting variables.
#' Columns of \code{W} should only include numeric vectors.  If categorical variables are used in the analysis and are not already recoded
#' into dummy variables, consider using the \code{model.matrix()} function and removing the intercept column (usually the first column)
#' @param family One of either \code{"binomial"} or \code{"gaussian"} for binary or continuous outcomes respectively
#' @param case_control_design Boolean, if outcome is binary and data comes from a case control sampling design, should be set to \code{TRUE}.
#' In addition, the user must supply the estimated disease prevalence in the \code{disease_prevalence} argument.  This is because
#' the TMLE procedure is not robust to biased sampling designs and must instead account for this bias through a case control weighted
#' TMLE
#' @param disease_prevalence Single numeric estimated prevalence of disease outcome.  Only necessary if case control design is set to \code{TRUE}.
#' @param obs.weights Vector of observation weights.  Default is to set all weights to 1 unless it is indicated that the sampling design is
#' case control, can be used for custom weighing of observations and overrides any other weights.
#' @param propensity_scores Optional argument to include pre-computed propensity scores, should be a vector of pre-computed propensity
#' scores for the exposure of \code{E}.
#' @param TMLE_args_list A list of arguments controlling the TMLE estimation process. See details below for full description
#' @param propensity_SL.library  List of learners to include in the propensity model SuperLearner (note that the propensity is fitted using a SuperLearner by default unless \code{propensity_formula} is supplied)
#' List of available learners can be viewed using \code{SuperLearner::listWrappers()}.  May require downloading other R packages.
#' @param propensity_SL.cvControl List of options to provide to propensity model SuperLearner.  The options include
#' \code{V} (the number of cross validation folds), \code{stratifyCV} (should the cross validation be stratified on the outcome?  Default is \code{TRUE} for binary outcome),
#' \code{shuffle} (should the data be shuffled?  default is \code{TRUE}), and \code{validRows} (do we want to supply validation data observations?).
#' @param include_G_propensity Option to include all SNPs in \code{G} in the propensity model.  Default is \code{FALSE}.  If set to \code{TRUE}, one should use
#' learners in the SuperLearner that can accommodate the dimension of the data.
#' @param include_W_outcome Option to include all confounders/adjusting variables defined in \code{W} in the outcome model.  Default is \code{TRUE}.
#' @param propensity_formula Option to include a formula for the propensity model to fit a generalized linear model.  \strong{NOTE:} when creating a
#' formula, refer to the exposure variable as \code{A} irregardless of its true name.  Confounders can be references through their columnnames
#' in the \code{W} dataframe, and elements from \code{G} can be included similarly if \code{include_G_propensity = TRUE}.
#' @param outcome_formula Option to include a formula for the outcome model to fit a generalized linear model.  \strong{NOTE:} when creating a
#' formula, refer to the outcome variable as \code{Y} and the exposure variable as \code{A} irregardless of their true names.  Confounders can be references through their columnnames
#' in the \code{G} dataframe, and elements from \code{W} can be included similarly if \code{include_W_propensity = TRUE}.
#' @param SNP_results Option to specify a subset of SNPs to iterate over.  If desired, a vector of column indices corresponding to the SNPs of interest
#' should be provided.  Default is to iterate over all SNPs in \code{G}
#' @param parallel Can \code{tlGxE} utilize parallel computing?  Default is \code{TRUE}
#' @param ncores If the \code{parallel} option is set to \code{TRUE}, how many cores can be used?  Default is the number of available cores on device.
#' @param progress Do we want to display progress bars for the results of the \code{tlGxE} fitting?  Progress will be shown as the function iterates over
#' the SNPs in \code{G} (as determined by the \code{SNP_results} argument).  Default is \code{TRUE}.
#' @param verbose Do we want to see messages regarding the fitting outputted?  Default is \code{TRUE}
#' @description
#' Performs \emph{tlGxE} scan over a set of SNPs.  \emph{tlGxE} estimates causal effect modification of the average treatment effect (ATE) and marginal odds ratio (MOR) of
#' exposure using targeted learning across disjoint subgroups of observations that have 0, 1, or 2 of a minor allele for a particular SNP. \emph{tlGxE} iterates over
#' a set of candidate SNPs that are included in the input \code{G}, where the argument \code{SNP_results} can be used to confine the analysis to only on a subset of
#' \code{G} (while adjusting for confounding by all SNPs in G).
#'
#' @details
#' Please see the guided examples of \code{tlGxE} usage in the vignette at \url{https://github.com/Daniel-Rud/tlgxe}.
#'
#' Below, please find the specification for the \code{TMLE_args_list} argument.
#'
#' \code{TMLE_args_list}: a list containing the following elements:
#' \itemize{
#' \item \code{outcome_method}: Specification of the outcome model for tlGxE, should be one of
#' either "glmnet", "glmnet_int", "gesso", or "SL". "glmnet" corresponds to the elastic net model implemented in the
#' \code{glmnet} package, "glmnet_int" corresponds to a non-hierarchical interaction model that includes all main effects for \emph{G}, \emph{E}, and
#' confounders, along with interactions between \emph{E} and all covariates in \eqn{\{G, W\}}.  "gesso" fits the hierarchical lasso GxE model described
#' in Zemlianskaia et al 2022.  "SL" corresponds to using SuperLearning for the outcome model.  If the "SL" option is specified, the
#' SuperLearner can be configured through the \code{outcome_SL.library} and \code{outcome_SL.cvControl} options.
#' \item \code{npv_thresh}: numeric value between \eqn{[0, 0.2]} which specifies the near positivity violator threshold, essentially setting
#' a bounds on defining observations with extreme propensity for exposure or no-exposure.
#' \item \code{near_positivity_method}: Either one of "trim" or "rebound".  Defines how to deal with near-positivity violators, as
#' determined by observations that have extreme propensity (with respect to the \code{npv_thresh}).  Near positivity violators are
#' either removed from the analysis if \code{near_positivity_method = "trim"} or rebounded to the maximum allowable propensity
#' (either npv_thresh or 1-npv_thresh) if \code{near_positivity_method = "rebound"}.
#' \item \code{nfolds_cv_Q_init}: Number of cross validation folds for the CV-TMLE procedure.  Controls number of cv folds for the outcome model,
#' where the outcome model is refitted \code{nfolds_cv_Q_init} times on \code{nfolds_cv_Q_init -1} data folds and the expected outcome is predicted on the
#' out-of-fold observations.
#' \item \code{nfolds_cv_glmnet_outcome}: Number of cv folds for the outcome model if a
#' \code{glmnet}-based model is used.  If \code{outcome_method = "SL"}, argument is ignored.
#' \item \code{alpha_outcome}: Numeric alpha parameter if \code{glmnet} model is used.  Alpha should be between \eqn{[0,1]} where \eqn{0} corresponds
#' to the Ridge Regression,  \eqn{.5} corresponds to the standard Elastic Net, and \eqn{1} corresponds to the LASSO.
#' \item \code{clever_cov_propensity_wt}: Option to either include clever covariate in TMLE procedure as a weighted regression instead of
#' including the clever covariates in the logistic fluctuation model.  Default is \code{TRUE}.
#' \item \code{outcome_SL.library}: List of learners to include in the outcome SuperLearner \strong{IF} \code{outcome_method = "SL"}, otherwise ignored.
#' List of available learners can be viewed using   \code{SuperLearner::listWrappers()}.  May require downloading other R packages.
#' \item \code{outcome_SL.cvControl}: List of options to provide to outcome SuperLearner \strong{IF} \code{outcome_method = "SL"}.  The options include
#' \code{V} (the number of cross validation folds), \code{stratifyCV} (should the cross validation be stratified on the outcome?  Default is \code{TRUE} for binary outcome),
#' \code{shuffle} (should the data be shuffled?  default is \code{TRUE}), and \code{validRows} (do we want to supply validation data observations?).
#' }
#' @export
tlGxE = function(Y, E, G, W = NULL, family = "binomial",
                    case_control_design = F, disease_prevalence = NULL,
                    obs.weights = NULL,
                    propensity_scores = NULL,
                    TMLE_args_list = list(
                      outcome_method = c("glmnet", "glmnet_int", "gesso", "SL"),
                      npv_thresh = (5/sqrt(length(Y)))/log(length(Y)),
                      near_positivity_method = c("trim", "rebound"),
                      nfolds_cv_Q_init = 10,
                      nfolds_cv_glmnet_outcome = 5,
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
                    include_G_propensity = F,
                    include_W_outcome = T,
                    propensity_formula = NULL,
                    outcome_formula = NULL,
                    SNP_results = 1:(ifelse(is.null(ncol(G)), 1, ncol(G))),
                    parallel = T,
                    ncores = ifelse(parallel, future::availableCores(), 1),
                    progress = T,
                    verbose = T)
{

  # save for later
  args_list = list(family = family,
                    case_control_design = case_control_design,
                    disease_prevalence = disease_prevalence,
                    obs.weights = obs.weights,
                    propensity_scores = propensity_scores,
                    TMLE_args_list = TMLE_args_list,
                    propensity_SL.library = propensity_SL.library,
                    propensity_SL.cvControl = propensity_SL.cvControl,
                    include_G_propensity = include_G_propensity,
                    include_W_outcome = include_W_outcome,
                    propensity_formula = propensity_formula,
                    outcome_formula = outcome_formula,
                    SNP_results = SNP_results,
                    parallel = parallel,
                    ncores = ncores,
                    progress = progress,
                    verbose = verbose
  )

  ###########################################
  # Perform Arg Checks
  ###########################################
  Y = as.numeric(Y); E = as.numeric(E)
  set_list_args = set_arg_lists(TMLE_args_list, propensity_SL.cvControl, family, Y)
  TMLE_args_list = set_list_args[[1]]; propensity_SL.cvControl = set_list_args[[2]]

  perform_arg_checks(Y = Y,
                     E = E,
                     G = G,
                     W = W,
                     family = family,
                     case_control_design = case_control_design,
                     disease_prevalence = disease_prevalence,
                     obs.weights = obs.weights,
                     propensity_scores = propensity_scores,
                     TMLE_args_list = TMLE_args_list,
                     propensity_SL.library = propensity_SL.library,
                     propensity_SL.cvControl = propensity_SL.cvControl,
                     include_G_propensity = include_G_propensity,
                     include_W_outcome = include_W_outcome,
                     propensity_formula = propensity_formula,
                     outcome_formula = outcome_formula,
                     SNP_results = SNP_results,
                     parallel = parallel,
                     ncores = ncores,
                     progress = progress,
                     verbose = verbose)
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

      full_exposure_data = cbind(E = E, G)
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
        message("Fitting global propensity model using SuperLearner...")
      }
      suppress_output(propensity_scores <- generate_propensity_SL(exposure_data = data.frame(if(is.null(W_propensity)){E}else{cbind(E, W_propensity)}),
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
  tlGxE_results = NULL

  if(verbose == T)
  {
    message("Iterating tlGxE over SNPs in G...")
  }

  if(progress == T)
  {
    progressr::with_progress(tlGxE_results <- tlGxE_EM_iterator(Y = Y,
                                                   E = E,
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
    tlGxE_results <- tlGxE_EM_iterator(Y = Y,
                                     E = E,
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
                                     ncores = ncores)
  }


  result_frame = do.call(cbind, tlGxE_results)
  if(!is.null(result_frame))
  {
    if(is.null(colnames(G)))
    {
      colnames(result_frame) = paste0("SNP_", SNP_results)
    }else
    {
      colnames(result_frame) = colnames(G)[SNP_results]
    }
  }

  result_list = list(
    tlGxE_scan_results = result_frame,
    global_propensity_scores = propensity_scores,
    passed_arguments = args_list
  )

  if(!is.null(result_list))
  {
    class(result_list) = "tlGxE"
  }

  return(result_list)
}

tlGxE_EM_iterator = function(Y, E, G, W = NULL,propensity_scores = NULL, family = "binomial",
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
                            ncores = future::availableCores())
{

  p <- NULL

  num_G = length(SNP_results)
  if(progress)
  {
    p <- progressr::progressor(steps = num_G)
  }

  tmle_results = NULL

  if(parallel)
  {
    future::plan(future::multisession, workers = ncores)
    on.exit(future::plan(future::sequential), add = TRUE)
  }

    tmle_results = future.apply::future_lapply(SNP_results, FUN = function(i)
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
                                   E = E,
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
      future::plan(future::sequential)
    }

  return(tmle_results)
}











