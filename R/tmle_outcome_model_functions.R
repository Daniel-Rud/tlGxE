compute_gesso = function(outcome_data, train_ids, test_ids, obs.weights, args_list)
{
  # gesso does not have implentation with weights!!! So we do not use...

  Y =  outcome_data[train_ids, 1]
  A = outcome_data[train_ids, 2]
  X = outcome_data[train_ids, -c(1:2)] %>% as.matrix

  new_A = outcome_data[test_ids, 2]
  new_X = outcome_data[test_ids, -c(1:2)] %>% as.matrix

  nfolds = args_list$nfolds_cv_glmnet
  family = "binomial" # note that we call function binomial (not 'binomial') for cv.glmnet


  Q0_mod = gesso::gesso.cv(G = X, E = A, Y = Y,
                    family = family, parallel = F,
                    verbose = F,
                    nfolds = nfolds)

  coefficients = gesso::gesso.coef(Q0_mod$fit, Q0_mod$lambda_min)

  beta_0 = coefficients$beta_0; beta_e = coefficients$beta_e
  beta_g = coefficients$beta_g; beta_gxe = coefficients$beta_gxe

  Q0AW = gesso::gesso.predict(beta_0, beta_e, beta_g, beta_gxe,
                       new_G = new_X, new_E = new_A, family = family)
  Q1W = gesso::gesso.predict(beta_0, beta_e, beta_g, beta_gxe,
                      new_G = new_X, new_E = rep(1, nrow(new_X)), family = family)
  Q0W = gesso::gesso.predict(beta_0, beta_e, beta_g, beta_gxe,
                      new_G = new_X, new_E = rep(0, nrow(new_X)), family = family)

  return_list = list(Q0AW = Q0AW, Q1W = Q1W, Q0W = Q0W)
  return(return_list)
}

compute_glmnet_interaction = function(outcome_data, train_ids, test_ids, obs.weights,  args_list)
{
  family = args_list$family
  alpha = args_list$alpha
  nfolds = args_list$nfolds_cv_glmnet
  type.measure = args_list$type.measure
  lambda =  args_list$lambda

  glm_model_matrix_train = cbind(outcome_data[train_ids, ],
                                 outcome_data[train_ids, -c(1:2)] * outcome_data[train_ids,2])
  colnames(glm_model_matrix_train) = c("Y", "A",  c(colnames(outcome_data)[-c(1:2)],
                                                    paste0("A*", colnames(outcome_data)[-c(1:2)])))

  # create fold ids
  stratification_var = interaction(glm_model_matrix_train[,"Y"], glm_model_matrix_train[,"A"])
  fold_ids = caret::createFolds(stratification_var, k = nfolds, list = FALSE)


  glm_model_matrix_test = cbind(outcome_data[test_ids, ],
                                outcome_data[test_ids, -c(1:2)] * outcome_data[test_ids,2])
  colnames(glm_model_matrix_test) = c("Y", "A",  c(colnames(outcome_data)[-c(1:2)],
                                                   paste0("A*", colnames(outcome_data)[-c(1:2)])))

  glm_model_matrix1 = cbind(outcome_data[test_ids, ], outcome_data[test_ids, -c(1:2)] )
  glm_model_matrix1[,2] = 1
  colnames(glm_model_matrix1) = c("Y", "A",  c(colnames(outcome_data)[-c(1:2)], paste0("A*", colnames(outcome_data)[-c(1:2)])))

  glm_model_matrix0 = cbind(outcome_data[test_ids, ],
                            matrix(0, nrow = length(test_ids), ncol = ncol(outcome_data) - 2))
  glm_model_matrix0[,2] = 0
  colnames(glm_model_matrix0) = c("Y", "A",  c(colnames(outcome_data)[-c(1:2)], paste0("A*", colnames(outcome_data)[-c(1:2)])))

  Q0_mod = glmnet::cv.glmnet(y = glm_model_matrix_train[,1], x = as.matrix(glm_model_matrix_train[,-1]),
                     family = family,
                     weights = obs.weights[train_ids],
                     alpha = alpha,
                     nfolds = nfolds,
                     type.measure = type.measure,
                     foldid = fold_ids)

  s = Q0_mod$lambda.1se
  if(lambda == "lambda.min")
  {
    s = Q0_mod$lambda.min
  }

  Q0AW = stats::predict(Q0_mod, newx = glm_model_matrix_test[,-1] %>% as.matrix, type = "response", s = s)
  Q1W = stats::predict(Q0_mod, newx = glm_model_matrix1[,-1]%>% as.matrix, type = "response", s = s)
  Q0W = stats::predict(Q0_mod, newx = glm_model_matrix0[,-1]%>% as.matrix, type = "response", s = s)
  return_list = list(Q0AW = Q0AW, Q1W = Q1W, Q0W = Q0W)
  return(return_list)
}

compute_glmnet = function(outcome_data, train_ids, test_ids, obs.weights, args_list)
{
  family = args_list$family
  alpha = args_list$alpha
  nfolds = args_list$nfolds_cv_glmnet
  type.measure = args_list$type.measure
  lambda =  args_list$lambda

  X = as.matrix(outcome_data[train_ids,-1])
  X_test = as.matrix(outcome_data[test_ids,-1])

  stratification_var = interaction(outcome_data$Y[train_ids], X[,"A"])
  fold_ids = caret::createFolds(stratification_var, k = nfolds, list = FALSE)

  Q0_mod = glmnet::cv.glmnet(y = outcome_data$Y[train_ids], x = X,
                     family = family,
                     weights = obs.weights[train_ids],
                     alpha = alpha,
                     nfolds = nfolds,
                     type.measure = type.measure,
                     foldid = fold_ids)

  s = Q0_mod$lambda.1se
  if(lambda == "lambda.min")
  {
    s = Q0_mod$lambda.min
  }

  Q0AW = stats::predict(Q0_mod, newx = X_test, type = "response", s = s )
  Q1W = stats::predict(Q0_mod, newx = cbind(rep(1, nrow(X_test)), X_test[,-1]), type = "response", s = s)
  Q0W = stats::predict(Q0_mod, newx = cbind(rep(0, nrow(X_test)), X_test[,-1]), type = "response", s = s)
  return_list = list(Q0AW = Q0AW, Q1W = Q1W, Q0W = Q0W)
  return(return_list)
}

compute_Superlearner = function(outcome_data, train_ids, test_ids, obs.weights, args_list)
{
  family = args_list$family
  outcome_SL.library = args_list$outcome_SL.library
  outcome_SL.cvControl = args_list$outcome_SL.cvControl

  X = outcome_data[train_ids,-1] %>% data.frame
  X_test = outcome_data[test_ids,-1] %>% data.frame
  Y = outcome_data[train_ids,1]
  Y_test = outcome_data[test_ids,1]

  suppressWarnings(Q0_mod <- SuperLearner::SuperLearner(Y = Y, X = X,
                                                        newX = NULL,
                                                        family = family,
                                                        SL.library = outcome_SL.library,
                                                        method = "method.NNLS",
                                                        obsWeights = obs.weights[train_ids],
                                                        cvControl = outcome_SL.cvControl))


  Q0AW = stats::predict(Q0_mod, newx = X_test, onlySL = TRUE)$pred %>% bound_limits(0,1)
  Q1W = stats::predict(Q0_mod, newx = cbind(rep(1, nrow(X_test)), X_test[,-1]), onlySL = TRUE)$pred %>% bound_limits(0,1)
  Q0W = stats::predict(Q0_mod, newx = cbind(rep(0, nrow(X_test)), X_test[,-1]), onlySL = TRUE)$pred %>% bound_limits(0,1)
  return_list = list(Q0AW = Q0AW, Q1W = Q1W, Q0W = Q0W)
  return(return_list)
}



# if no CV, set nfolds_cv_Q_init = 1
generate_Q0_cv = function(outcome_data, family, obs.weights, alpha = 0.5, lambda = "lambda.1se",
                          nfolds_cv_Q_init = 10, nfolds_cv_glmnet = 3,
                          outcome_method = c("glmnet_int", "glmnet", "gesso", "SL"),
                          outcome_SL.library = NULL,
                          outcome_SL.cvControl = NULL,
                          type.measure = "deviance")
{
  return_list = NULL
  # if no confounders specified in outcome model -- can only happen if include_W_outcome and include_G_outcome are both FALSE.
  if(ncol(outcome_data) == 2)
  {
    glm_mod = stats::glm(Y ~ A, data = outcome_data, family = "binomial", weights = obs.weights)
    Q0AW = stats::predict(glm_mod, type = "response")

    # Counterfactual predictions
    Q1W = stats::predict(glm_mod, newdata = transform(outcome_data, A = 1), type = "response")
    Q0W = stats::predict(glm_mod, newdata = transform(outcome_data, A = 0), type = "response")

    return_list = list(Q0AW = Q0AW, Q1W = Q1W, Q0W = Q0W)
  }
  else
  {
    # initialize CV splits

    # use factor call for stratified CV
    folds = NULL
    if(family == "gaussian") # do not create splits based on factor
    {
      folds = caret::createFolds(outcome_data$Y, k = nfolds_cv_Q_init)
    }else
    {
      # note that the cross validation folds are stratified by both outcome and and exposure status
      stratification_var = interaction(outcome_data$Y, outcome_data$A)
      folds = caret::createFolds(stratification_var, k = nfolds_cv_Q_init)
    }

    if(outcome_method != "SL")
    {
      # since we use logistic fluctuation, family is always binomial
      # note that we use `binomial` and not `"binomial"`, the former works and the
      # latter does not with a bounded continuous outcome in (0,1), I think cause it
      # calls binomial function in R and the other might be an implementation in glmnet

      # we make the condition whether the SL is being called because
      # super learner did not work when using probability outcome with binomial family

      family = stats::binomial

    }

    args_list = list(family = family,
                     alpha = alpha,
                     lambda = lambda,
                     nfolds_cv_glmnet = nfolds_cv_glmnet,
                     type.measure = type.measure,
                     outcome_SL.library = outcome_SL.library,
                     outcome_SL.cvControl = outcome_SL.cvControl)

    outcome_function = switch(outcome_method,
                              glmnet = compute_glmnet,
                              glmnet_int = compute_glmnet_interaction,
                              gesso = compute_gesso,
                              SL = compute_Superlearner
    )

    results = lapply(1:length(folds), FUN = function(i)
    {

      train_ids = (1:nrow(outcome_data))[-folds[[i]]]
      test_ids = folds[[i]]
      if(nfolds_cv_Q_init == 1)
      {
        train_ids = test_ids
      }

      return_list = outcome_function(outcome_data = outcome_data, train_ids = train_ids,
                                     test_ids = test_ids, obs.weights = obs.weights,
                                     args_list = args_list)

      return(return_list)

    })

    # put predictions in order

    fold_ids = do.call(c, folds)
    order_obs = order(fold_ids)

    Q0AW = do.call(c, lapply(results, "[[", 1) )
    Q0AW = Q0AW[order_obs]

    Q1W = do.call(c, lapply(results, "[[", 2))
    Q1W = Q1W[order_obs]

    Q0W = do.call(c, lapply(results, "[[", 3))
    Q0W = Q0W[order_obs]
    return_list = list(Q0AW = Q0AW, Q1W = Q1W, Q0W = Q0W)
  }
  return(return_list)
}
