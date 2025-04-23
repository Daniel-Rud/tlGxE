
perform_arg_checks = function(Y, E, G, W, family, case_control_design,
                              disease_prevalence,  obs.weights,  propensity_scores,
                              TMLE_args_list,  propensity_SL.library,  propensity_SL.cvControl, include_G_propensity, include_G_outcome,
                              include_W_outcome,  propensity_formula,  outcome_formula,  SNP_results,  parallel,  ncores,  progress,  verbose)
{
  ####################################
  # family
  if(!(family %in% c("binomial", "gaussian")))
  {
    stop("`family` must be one of either 'binomial' or 'gaussian'.")
  }
  ####################################
  # Y
  if(family == "binomial")
  {
    if(sum(!(Y %in% c(0,1))) > 0)
    {
      stop("If 'family' = 'binomial', all elements of Y should be either 0 or 1.")
    }
  }
  n = length(Y)

  ####################################
  # E
  # check if all elemnts of exposure E is 0 1
  if(sum(!(E %in% c(0,1))) > 0)
  {
    stop("Exposure variable 'E' MUST be a binary 0 1 vector.")
  }
  # check if E has same length as Y
  if(length(E) != n)
  {
    stop("Size mismatch between 'Y' and 'E' vectors.")
  }
  ####################################
  # G
  # check if size mismath between Y and G
  if(ifelse(is.null(nrow(G)), length(G), nrow(G)) != n)
  {
    stop("Size mismatch between 'Y' and 'G'.")
  }
  ####################################
  # W
  if(!is.null(W))
  {
    # check for size mismatch
    if(ifelse(is.null(nrow(W)), length(W), nrow(W)) != n)
    {
      stop("Size mismatch between 'Y' and 'W'.")
    }
    # check that all columns are numeric
    if(!all(sapply(W, is.numeric)))
    {
      stop("All columns in 'W' MUST be numeric.  Are any factor variables included?  Consider using 'model.matrix' with 'W'.")
    }
  }
  ####################################

  ####################################
  # case_control_design

  if(!is.logical(case_control_design))
  {
    stop("'case_control_design' must be a logical.")
  }
  ####################################
  # disease_prevalence

  if(is.null(disease_prevalence))
  {
    if(case_control_design)
    {
      stop("'disease_prevalence' must be supplied if 'case_control_design = T'.")
    }
  }else
  {
    if(disease_prevalence < 0 || disease_prevalence > 1)
    {
      stop("'disease_prevalence' must be between 0 and 1.")
    }
  }

  ####################################
  # obs.weights
  if (!is.null(obs.weights))
  {
    if (length(obs.weights) != n) {
      stop("Size mismatch between 'Y' and 'obs.weights'.")
    }
    if (any(obs.weights < 0)) {
      stop("'obs.weights' must be non-negative.")
    }
  }

  ####################################
  # propensity_scores
  if (!is.null(propensity_scores)) {
    if (length(propensity_scores) != n) {
      stop("Size mismatch between 'Y' and 'propensity_scores'.")
    }
    if (any(propensity_scores <= 0 | propensity_scores >= 1)) {
      stop("'propensity_scores' must be strictly between 0 and 1.")
    }
  }

  ####################################
  # TMLE_args_list
  # required_tmle_args <- c(
  #   "outcome_method", "npv_thresh", "near_positivity_method",
  #   "nfolds_cv_Q_init", "nfolds_cv_glmnet_outcome", "alpha_outcome",
  #   "clever_cov_propensity_wt", "outcome_SL.library", "outcome_SL.cvControl"
  # )

  if (!is.list(TMLE_args_list)) {
    stop("'TMLE_args_list' must be a list.")
  }

  # missing_args = setdiff(required_tmle_args, names(TMLE_args_list))
  # if (length(missing_args) > 0) {
  #   stop(paste("The following required elements are missing from 'TMLE_args_list':",
  #              paste(missing_args, collapse = ", ")))
  # }

  # Ensure outcome_SL.cvControl is a list
  if (!is.list(TMLE_args_list$outcome_SL.cvControl)) {
    stop("'outcome_SL.cvControl' must be a list within 'TMLE_args_list'.")
  }
  # # Check required elements in outcome_SL.cvControl
  # required_cv_control_args <- c("V", "stratifyCV", "shuffle", "validRows")
  # missing_cv_args <- setdiff(required_cv_control_args, names(TMLE_args_list$outcome_SL.cvControl))
  #
  # if (length(missing_cv_args) > 0) {
  #   stop(paste("The following required elements are missing from 'outcome_SL.cvControl':",
  #              paste(missing_cv_args, collapse = ", ")))
  # }

  ####################################
  # propensity_SL.library
  if (!is.character(propensity_SL.library) || length(propensity_SL.library) == 0) {
    stop("'propensity_SL.library' must be a non-empty character vector.")
  }

  ####################################
  # propensity_SL.cvControl
  if (!is.list(propensity_SL.cvControl)) {
    stop("'propensity_SL.cvControl' must be a list.")
  }

  # Check required elements in propensity_SL.cvControl
  # required_propensity_cv_control_args <- c("V", "stratifyCV", "shuffle", "validRows")
  # missing_propensity_cv_args <- setdiff(required_propensity_cv_control_args, names(propensity_SL.cvControl))
  #
  # if (length(missing_propensity_cv_args) > 0) {
  #   stop(paste("The following required elements are missing from 'propensity_SL.cvControl':",
  #              paste(missing_propensity_cv_args, collapse = ", ")))
  # }


  ####################################
  # include_G_propensity, include_W_outcome
  if (!is.logical(include_G_propensity)) {
    stop("'include_G_propensity' must be logical.")
  }
  if (!is.logical(include_G_outcome)) {
    stop("'include_G_outcome' must be logical.")
  }
  if (!is.logical(include_W_outcome)) {
    stop("'include_W_outcome' must be logical.")
  }
  ####################################
  # SNP_results
  if (!is.numeric(SNP_results)) {
    stop("'SNP_results' must be a numeric vector.")
  }

  if (is.matrix(G) || is.data.frame(G))
  {
    n_cols_G <- ncol(G)  # Number of columns in G (if it's a matrix or data frame)
  } else if (is.numeric(G))
  {
    n_cols_G <- 1  # If G is a numeric vector, we treat it as having one column
  } else
  {
    stop("'G' must be a numeric matrix/data frame or a numeric vector.")
  }

  # Check that all elements in SNP_results are valid indices for columns in G
  if (any(SNP_results < 1 | SNP_results > n_cols_G)) {
    stop("'SNP_results' indices must be between 1 and the number of columns in 'G'.")
  }
  if (anyDuplicated(SNP_results)) {
    stop("'SNP_results' should not contain duplicate indices.")
  }

  ####################################
  # parallel
  if (!is.logical(parallel)) {
    stop("'parallel' must be logical.")
  }

  ####################################
  # ncores
  if (!is.numeric(ncores) || ncores < 1 || floor(ncores) != ncores) {
    stop("'ncores' must be a positive integer.")
  }

  ####################################
  # progress, verbose
  if (!is.logical(progress)) {
    stop("'progress' must be logical.")
  }
  if (!is.logical(verbose)) {
    stop("'verbose' must be logical.")
  }
}



set_arg_lists <- function(TMLE_args_list, propensity_SL.cvControl, family, Y) {

  # Default TMLE arguments
  default_TMLE_args_list <- list(
    outcome_method = c("glmnet", "glmnet_int", "gesso", "SL"),
    npv_thresh = (5 / sqrt(length(Y))) / log(length(Y)),
    near_positivity_method = c("trim", "rebound"),
    nfolds_cv_Q_init = 10,
    nfolds_cv_glmnet_outcome = 5,
    alpha_outcome = 0.5,
    clever_cov_propensity_wt = TRUE,
    outcome_SL.library = c("SL.glmnet", "SL.rpart", "SL.bartMachine"),
    outcome_SL.cvControl = list(
      V = 10L,
      stratifyCV = ifelse(family == "gaussian", FALSE, TRUE),
      shuffle = TRUE,
      validRows = NULL
    )
  )

  # Update TMLE arguments
  arg_names <- names(TMLE_args_list)
  for (i in seq_along(arg_names)) {

    if (hasName(default_TMLE_args_list, arg_names[i]))
    {

      if (arg_names[i] == "outcome_SL.cvControl" && is.list(TMLE_args_list[[arg_names[i]]]))
      {

        cv_names <- names(TMLE_args_list[[arg_names[i]]])

        for (j in seq_along(cv_names))
        {
          if (hasName(default_TMLE_args_list[[arg_names[i]]], cv_names[j]))
          {
            default_TMLE_args_list[[arg_names[i]]][[cv_names[j]]] <- TMLE_args_list[[arg_names[i]]][[cv_names[j]]]
          }
        }
      } else
      {

        default_TMLE_args_list[[arg_names[i]]] <- TMLE_args_list[[arg_names[i]]]
      }
    }
  }

  # Default propensity control list
  default_propensity_SL.cvControl <- list(
    V = 10L,
    stratifyCV = TRUE,
    shuffle = TRUE,
    validRows = NULL
  )

  # Update propensity control arguments
  cv_names <- names(propensity_SL.cvControl)
  for (i in seq_along(cv_names))
  {

    if (hasName(default_propensity_SL.cvControl, cv_names[i]))
    {

      default_propensity_SL.cvControl[[cv_names[i]]] <- propensity_SL.cvControl[[cv_names[i]]]

    }
  }

  return(list(
    changed_TMLE_args_list = default_TMLE_args_list,
    changed_propensity_SL.cvControl = default_propensity_SL.cvControl
  ))
}
