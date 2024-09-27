


# propensity model with superlearner 
generate_propensity_SL = function(exposure_data, 
                                  obs.weights,
                                  SL.library = c("SL.gam", 
                                                 "SL.glm.interaction", 
                                                 "SL.glmnet", 
                                                 "SL.ranger", 
                                                 "SL.xgboost"), 
                                  SL.cvControl = list(V = 10L, 
                                                      stratifyCV = T, 
                                                      shuffle = TRUE, 
                                                      validRows = NULL), 
                                  parallel = FALSE, 
                                  ncores = NULL)
{
  # length of PS vector, if exposure data is vector or data frame
  l = nrow(exposure_data)
  PS = numeric(length = l)
  
  if(ncol(exposure_data) == 1) # if no confounders of exposure
  {
    
    PS_mod = glm(formula = A ~ 1, data = exposure_data, family = "binomial", 
                 weights = obs.weights)
    
    PS = rep(plogis(coef(PS_mod)), l)
    
  }else if(ncol(exposure_data) == 2) # if there is just a single confounder 
  {
    # we just run a glm because some SL methods do not accomodate a single predictor
    # glmnet fails with only 1 confounder -- x should be a matrix with 2 or more columns
    PS_mod = glm(formula = A ~ ., data = exposure_data, family = "binomial", 
                 weights = obs.weights)
    
    PS = rep(plogis(coef(PS_mod)), l)
  }
  
  else # if there are confounders for exposure 
  {
    Y = exposure_data[,1]
    X = exposure_data[,-1, drop = F]
    
    if(!parallel)
    {
      suppressWarnings(PS_mod <- SuperLearner(Y = Y, X = X,
                                              newX = X, 
                                              family = binomial, 
                                              SL.library = SL.library, 
                                              method = "method.NNLS", 
                                              obsWeights = obs.weights, 
                                              cvControl = SL.cvControl))
      
      PS = PS_mod$SL.predict
    }else
    {
      # set cluster with number of cores 
      cluster = parallel::makeCluster(ncores)
      # call superlearner library in each core 
      invisible(parallel::clusterEvalQ(cluster, {library(SuperLearner)}))
      # set random seed 
      parallel::clusterSetRNGStream(cluster, 1)
      # call snow parallelization
      suppressWarnings(PS_mod <- snowSuperLearner(Y = Y, X = X,
                                                  newX = X, 
                                                  family = binomial, 
                                                  SL.library = SL.library, 
                                                  cluster = cluster,
                                                  method = "method.NNLS", 
                                                  obsWeights = obs.weights, 
                                                  cvControl = SL.cvControl))
      PS = PS_mod$SL.predict
      
      parallel::stopCluster(cluster)
    }
  }
  return(PS)
}
