################################################################################
## environment setup
################################################################################
rm(list=ls())
setwd("C:/Users/ywang/OneDrive - Corcept Therapeutics, Inc/Documents/JSM2025/")

library(readxl)
library(clipr)
library(dplyr)
library(haven)
library(Hmisc)
library(tidyr)
library(survival)
library(olsrr)
library(glmnet)
library(ggplot2)
library(reshape2)
library(randomForest)
library(Boruta)
library(leaps)
library(foreach)
library(doParallel)
library(RColorBrewer)
library(ggradar)
library(scales)
library(patchwork)
library(caret)
library(ncvreg)
library(doSNOW)

sim.xy <- function(n, p, nval, rho=0, s=5, beta.type=1, snr=1) {
  # Generate predictors
  x = matrix(rnorm(n*p),n,p)
  xval = matrix(rnorm(nval*p),nval,p)
  
  # Introduce autocorrelation, if needed
  if (rho != 0) {
    inds = 1:p
    Sigma = rho^abs(outer(inds, inds, "-"))
    obj = svd(Sigma)
    Sigma.half = obj$u %*% (sqrt(diag(obj$d))) %*% t(obj$v)
    x = x %*% Sigma.half
    xval = xval %*% Sigma.half
  }
  else Sigma = diag(1,p)
  
  # Generate underlying coefficients
  s = min(s,p)
  beta = rep(0,p)
  if (beta.type==1) {
    beta[round(seq(1,p,length=s))] = 1
  } else if (beta.type==2) {
    beta[1:s] = 1
  } else if (beta.type==3) {
    beta[1:s] = seq(10,0.5,length=s)
  } else if (beta.type==4) {
    beta[1:6] = c(-10,-6,-2,2,6,10)
  } else {
    beta[1:s] = 1
    beta[(s+1):p] = 0.5^(1:(p-s))
  }
  
  # Set snr based on sample variance on infinitely large test set
  vmu = as.numeric(t(beta) %*% Sigma %*% beta)
  sigma = sqrt(vmu/snr)
  
  # Generate responses
  y = as.numeric(x %*% beta + rnorm(n)*sigma)
  yval = as.numeric(xval %*% beta + rnorm(nval)*sigma)
  
  list(x,y,xval,yval,Sigma,beta,sigma)
}


calc_mcc <- function(train_data, model_coeff, beta) {
  all_var <- colnames(train_data)[-length(colnames(train_data))]
  true_signals <- colnames(train_data)[which(beta>coeff_cutoff)]
  true_noise <- setdiff(all_var, true_signals)
  
  tp <- sum(model_coeff %in% true_signals)
  fp <- sum(model_coeff %in% true_noise)
  fn <- sum(!(true_signals %in% model_coeff))
  tn <- sum(!(true_noise %in% model_coeff))

  # print(tp)
  # print(fp)
  # print(fn)
  # print(tn)
  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)
  mcc <- ((tp * tn) - (fp * fn))/( sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) )
  
  # Precision and Recall
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  
  # F1 score
  F1 <- 2 * (precision * recall) / (precision + recall)
  F1
  
  if (length(model_coeff)>0) accuracy <- tp/length(model_coeff)
  else accuracy <- 0
  list(tp, fp, fn, tn, accuracy, sensitivity, specificity, mcc, F1)
}

rr <- function(beta, Sigma, sigma, betahat, betahat0) {
  risk.null = diag(t(beta) %*% Sigma %*% beta)
  delta = betahat - beta
  risk = diag(t(delta) %*% Sigma %*% delta) + as.numeric(betahat0)^2 ## relative risk
  err.test =  risk + sigma^2
  
  list(risk/risk.null, err.test/sigma^2)
}

calc_rmse <- function(yval, predictions) {
  
  # sqrt(mean((yval - predictions)^2))/(max(yval) - min(yval))
  sqrt(mean((yval - predictions)^2))/sd(yval)
  
}

calc_r2 <- function(yval, predictions) {
  
  1 - sum((yval - predictions)^2) / sum((yval - mean(yval))^2)
}

## generate output data frame
coeff_cutoff <- 1e-2 ## minimum coeff weights to be considered as selected
nfold <- 10

n_values <- c(100) #c(100, 500, 1000)
p_values <- c(30) #c(15, 50, 100, 1000)
rho_values <- c(0.2, 0.5, 0.8)
s_values <- c(5)
snr_values <- exp(seq(-3,1.5,0.5))
beta_type <- 5
nsim <- 50 # Number of replicates in each scenarios
method_list <- c("stepwise", "lasso.1se", "lasso.min", "ridge",
                 "elasticnet", "RF", "leaps", "scad")
nmethod <- length(method_list)


grid <- expand.grid(n = n_values, p = p_values, rho=rho_values, 
                    s=s_values, snr = snr_values)
# grid <- grid[1, ]
dim(grid)

# detectCores()
cl = makeCluster(10)
# registerDoParallel(cl)
registerDoSNOW(cl)

pkgs <- c("dplyr", "tidyr", "glmnet", "olsrr", "ncvreg",
          "randomForest", "Boruta", "leaps", "caret")

combinations <- expand.grid(i = 1:nsim, # i is number of simulations
                            j = 1:nrow(grid)) # j is each scenario
set.seed(123)

# Create a text progress bar
pb <- txtProgressBar(max = nrow(combinations), style = 3)

# Progress function
progress <- function(n) setTxtProgressBar(pb, n)

start_time <- Sys.time()
sim_res <- foreach(k = 1:nrow(combinations), .packages = pkgs,
                   .options.snow = list(progress = progress),
                   .combine = rbind) %dopar% {
    i <- combinations$i[k]
    j <- combinations$j[k]
    row <- grid[j, ]
    
    ## set up data
    out <- sim.xy(n=row$n, p=row$p, nval=row$n, rho=row$rho, s=row$s, 
                  beta.type=beta_type, snr=row$snr)
    
    x <- out[[1]]
    y <- out[[2]]
    xval <- out[[3]]
    yval <- out[[4]]
    Sigma <- out[[5]]
    beta <- out[[6]]
    sigma <- out[[7]]

    colnames(x) <- paste0("X", seq(ncol(x)))
    colnames(xval) <- paste0("X", seq(ncol(xval)))
    train_data <- cbind(x, y)
    train_data <- as.data.frame(train_data)
    test_data <- cbind(xval, yval)
    test_data <- as.data.frame(test_data)
    
    each_sim_res <- data.frame(row = rep(j, nmethod),
                               sim = rep(i, nmethod),
                               n = rep(row$n, nmethod),
                               p = rep(row$p, nmethod),
                               rho = rep(row$rho, nmethod),
                               s= rep(row$s, nmethod),
                               snr= rep(row$snr, nmethod),
                               method = method_list,
                               nfeature = rep(NA, nmethod),
                               rmse = rep(NA, nmethod),
                               r2 = rep(NA, nmethod),
                               tp = rep(NA, nmethod),
                               fp = rep(NA, nmethod),
                               fn = rep(NA, nmethod),
                               tn = rep(NA, nmethod),
                               accuracy = rep(NA, nmethod),
                               sensitivity = rep(NA, nmethod),
                               specificity = rep(NA, nmethod),
                               mcc = rep(NA, nmethod),
                               f1 = rep(NA, nmethod),
                               relative_risk = rep(NA, nmethod),
                               relative_test_error = rep(NA, nmethod))
    

    
    ## stepwise based on p-value

    if (row$p <= 100) {

      full_model <- lm(y ~ ., data = train_data)
      stepwise_model <- tryCatch({
       ols_step_both_p(
          full_model,
          pent = 0.10,  # p-value threshold for variable entry
          prem = 0.05,  # p-value threshold for variable removal
          details = F)$model
      }, error = function(e) {
         return(lm(y~1, data = train_data))
      })

    } else { # add a prefilter step
      ## choose by pvalue
      # top_vars <- names(sort(sapply(train_data[, -ncol(train_data)], function(x) summary(glm(train_data$y ~ x))$coefficients[2, 4]),
      #                        decreasing = FALSE))[1:100]
      
      ## choose by correlation
      marg_corr <- abs(cor(x, y, use = "pairwise.complete.obs"))
      ranked_vars <- order(marg_corr, decreasing = TRUE)
      top_vars <- rownames(marg_corr)[ranked_vars[1:100]]
      
      train_data1 <- train_data[,c(top_vars, "y")]

      full_model <- lm(y ~ ., data = train_data1)
      stepwise_model <- tryCatch({
        ols_step_both_p(
          full_model,
          pent = 0.10,  # p-value threshold for variable entry
          prem = 0.05,  # p-value threshold for variable removal
          details = F)$model
      }, error = function(e) {
        return(lm(y~1, data = train_data1))
      })

    }

    predictions <- predict(stepwise_model, newdata = test_data)
    each_sim_res[1,"rmse"] <- calc_rmse(yval, predictions)
    each_sim_res[1,"r2"] <- calc_r2(yval, predictions)
    each_sim_res[1,"nfeature"] <- length(names(coef(stepwise_model)[-1]))
    each_sim_res[1,c("tp","fp", "fn", "tn", "accuracy", "sensitivity", "specificity", "mcc",
                     "f1")] <- calc_mcc(train_data, names(coef(stepwise_model)[-1]),
                                                                       beta)
    all_vars <- colnames(x)
    selected_coefs <- coef(stepwise_model)
    full_coefs <- setNames(rep(0, length(all_vars)), all_vars)
    matching_vars <- intersect(names(selected_coefs), all_vars)
    full_coefs[matching_vars] <- selected_coefs[matching_vars]
    each_sim_res[1,c("relative_risk", "relative_test_error")] <- rr(beta, Sigma, sigma, full_coefs,
                                          selected_coefs[1])

    ## LASSO 1se
    cv <- cv.glmnet(x, y, alpha = 1, nfold=nfold, standardize = T)
    predictions <- predict(cv, newx = xval, s = "lambda.1se")
    each_sim_res[2,"rmse"] <- calc_rmse(yval, predictions)
    each_sim_res[2,"r2"] <- calc_r2(yval, predictions) 
    lasso_coeff_all <- coef(cv, s = "lambda.1se")[,"s0"][-1]
    lasso_coeff_sel <- names(lasso_coeff_all[which(abs(lasso_coeff_all)>0)])
    each_sim_res[2,"nfeature"] <- length(lasso_coeff_sel)
    each_sim_res[2,c("tp","fp", "fn", "tn", "accuracy", "sensitivity", "specificity", "mcc", 
                     "f1")] <- calc_mcc(train_data, lasso_coeff_sel, beta)
    each_sim_res[2,c("relative_risk", "relative_test_error")] <- rr(beta, Sigma, sigma, lasso_coeff_all,
                                          coef(cv, s = "lambda.1se")[,"s0"][1])
    
    ## LASSO min
    predictions <- predict(cv, newx = xval, s = "lambda.min")
    each_sim_res[3,"rmse"] <- calc_rmse(yval, predictions)
    each_sim_res[3,"r2"] <- calc_r2(yval, predictions) 
    lasso_coeff_all <- coef(cv, s = "lambda.min")[,"s0"][-1]
    lasso_coeff_sel <- names(lasso_coeff_all[which(abs(lasso_coeff_all)>0)])
    each_sim_res[3,"nfeature"] <- length(lasso_coeff_sel)
    each_sim_res[3,c("tp","fp", "fn", "tn", "accuracy", "sensitivity", "specificity", "mcc",
                     "f1")] <- calc_mcc(train_data, lasso_coeff_sel, beta)
    each_sim_res[3,c("relative_risk", "relative_test_error")] <- rr(beta, Sigma, sigma, lasso_coeff_all,
                                          coef(cv, s = "lambda.min")[,"s0"][1])
    
    ## Ridge
    cv <- cv.glmnet(x, y, alpha = 0, nfold=nfold, standardize = T)
    predictions <- predict(cv, newx = xval, s = "lambda.1se")
    each_sim_res[4,"rmse"] <- calc_rmse(yval, predictions)
    each_sim_res[4,"r2"] <- calc_r2(yval, predictions) 
    ridge_coeff_all <- coef(cv, s = "lambda.1se")[,"s0"][-1]
    ridge_coeff_sel <- names(ridge_coeff_all[which(abs(ridge_coeff_all)>0)])
    each_sim_res[4,"nfeature"] <- length(ridge_coeff_sel)
    each_sim_res[4,c("tp","fp", "fn", "tn", "accuracy", "sensitivity", "specificity", "mcc", 
                     "f1")] <- calc_mcc(train_data, ridge_coeff_sel, beta)
    each_sim_res[4,c("relative_risk", "relative_test_error")] <- rr(beta, Sigma, sigma, ridge_coeff_all,
                                          coef(cv, s = "lambda.1se")[,"s0"][1])
    
    ## Elastic-net
    cv <- cv.glmnet(x, y, alpha = 0.5, nfold=nfold, standardize = T)
    predictions <- predict(cv, newx = xval, s = "lambda.1se")
    each_sim_res[5,"rmse"] <- calc_rmse(yval, predictions)
    each_sim_res[5,"r2"] <- calc_r2(yval, predictions) 
    elastic_coeff_all <- coef(cv, s = "lambda.1se")[,"s0"][-1]
    elastic_coeff_sel <- names(elastic_coeff_all[which(abs(elastic_coeff_all)>0)])
    each_sim_res[5,"nfeature"] <- length(elastic_coeff_sel)
    each_sim_res[5,c("tp","fp", "fn", "tn", "accuracy", "sensitivity", "specificity", "mcc",
                     "f1")] <- calc_mcc(train_data, elastic_coeff_sel, beta)
    each_sim_res[5,c("relative_risk", "relative_test_error")] <- rr(beta, Sigma, sigma, elastic_coeff_all,
                                          coef(cv, s = "lambda.1se")[,"s0"][1])
    
    ## Boruta and RF
    boruta_result <- Boruta(y ~ ., data = train_data, doTrace = 0)
    selected_vars <- getSelectedAttributes(boruta_result, withTentative = F)
    selected_tent_vars <- getSelectedAttributes(boruta_result, withTentative = T)
    each_sim_res[6,"nfeature"] <- length(selected_vars)
    train_selected <- train_data[, c(selected_tent_vars, "y"), drop=F]
    if (length(selected_tent_vars) > 0) {
      # model <- randomForest(y ~ ., data = train_selected)
      
      ## add parameter tuning
      ctrl <- trainControl(method = "cv", number = nfold)
      model <- train(
        y ~ .,
        data = train_selected,
        method = "rf",
        trControl = ctrl,
        tuneLength=5,
        ntree=500
      )
      predictions <- predict(model, newdata = test_data)
      
    } else { ## no feature selected, use an intercept model
      model <- lm(y~1, data = train_data)
      predictions <- predict(model, newdata = test_data)
    }
    each_sim_res[6,"rmse"] <- calc_rmse(yval, predictions)
    each_sim_res[6,"r2"] <- calc_r2(yval, predictions) 
    each_sim_res[6,c("tp","fp", "fn", "tn", "accuracy", "sensitivity", "specificity", "mcc",
                     "f1")] <- calc_mcc(train_data, selected_vars, beta)
    
    ## best subset leaps
    if (row$p <= 50) {
    regsubsets.out <- regsubsets(y ~ ., data = train_data,
                                 nvmax = NULL,
                                 really.big=T,
                                 method = "seqrep")
    } else { # add a prefilter step
      
      ## choose by pvalue
      # top_vars <- names(sort(sapply(train_data[, -ncol(train_data)], function(x) summary(glm(train_data$y ~ x))$coefficients[2, 4]),
      #                        decreasing = FALSE))[1:50]
      
      ## choose by correlation
      marg_corr <- abs(cor(x, y, use = "pairwise.complete.obs"))
      ranked_vars <- order(marg_corr, decreasing = TRUE)
      top_vars <- rownames(marg_corr)[ranked_vars[1:50]]
      
      regsubsets.out <- regsubsets(y ~ ., data = train_data[,c(top_vars, "y"),drop=F],
                                   nvmax = NULL,
                                   really.big=T,
                                   method = "seqrep")
    }

    summary_best_subset <- summary(regsubsets.out)
    best_size <- which.min(summary_best_subset$bic)
    # best_size <- which.max(summary_best_subset$adjr2)
    selected_vars <- names(which(summary_best_subset$which[best_size, -1]))
    each_sim_res[7,"nfeature"] <- length(selected_vars)
    train_selected <- train_data[, c(selected_vars, "y"), drop=F]
    final_model <- lm(y~., data = train_selected)
    predictions <- predict(final_model, newdata = test_data)
    each_sim_res[7,"rmse"] <- calc_rmse(yval, predictions)
    each_sim_res[7,"r2"] <- calc_r2(yval, predictions)
    each_sim_res[7,c("tp","fp", "fn", "tn", "accuracy", "sensitivity", "specificity", "mcc",
                     "f1")] <- calc_mcc(train_data, selected_vars, beta)
    all_vars <- colnames(x)
    selected_coefs <- coef(final_model)
    full_coefs <- setNames(rep(0, length(all_vars)), all_vars)
    matching_vars <- intersect(names(selected_coefs), all_vars)
    full_coefs[matching_vars] <- selected_coefs[matching_vars]
    each_sim_res[7,c("relative_risk", "relative_test_error")] <- rr(beta, Sigma, sigma, full_coefs,
                                          selected_coefs[1])
    
    ## SCAD 
    cv_scad <- cv.ncvreg(x, y, penalty = "SCAD", nfolds=nfold)
    predictions <- predict(cv_scad, X = xval)
    each_sim_res[8,"rmse"] <- calc_rmse(yval, predictions)
    each_sim_res[8,"r2"] <- calc_r2(yval, predictions)
    scad_coeff_all <- coef(cv_scad)[-1]
    scad_coeff_sel <- names(scad_coeff_all[which(abs(scad_coeff_all)>0)])
    each_sim_res[8,"nfeature"] <- length(scad_coeff_sel)
    each_sim_res[8,c("tp","fp", "fn", "tn", "accuracy", "sensitivity", "specificity", "mcc",
                     "f1")] <- calc_mcc(train_data, scad_coeff_sel, beta)
    each_sim_res[8,c("relative_risk", "relative_test_error")] <- rr(beta, Sigma, sigma, scad_coeff_all,
                                                                    coef(cv_scad)[1])
    
    
    each_sim_res

} ## j
stopCluster(cl)
stop_time <- Sys.time()
stop_time - start_time


## checking run numbers
dim(sim_res)
sim_res %>% group_by(n, p, rho, s, snr, method) %>% summarise(sim=n()) %>%
  filter(sim < nsim)
sim_res %>% distinct(n, p, rho, s, snr) %>% dim
saveRDS(sim_res, file = paste0("n", n_values, "p",
                               p_values, "beta", beta_type, ".rds"))