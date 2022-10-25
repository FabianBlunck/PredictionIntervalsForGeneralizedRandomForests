########################################################################################################################
##                                                                                                                    ##
## Code for Chapter 10 - "Emprirical Study"                                                                           ##
##                                                                                                                    ##
########################################################################################################################


# Install and load packages
packages <- c("iterators", "foreach", "ggplot2", "copula", "doParallel", "doSNOW", "parallel", "snow", "grf", "ggpubr")

# install.packages(packages)

lapply(packages, require, character.only=T)


# Random numbers from R version 4.2.1
RNGversion("4.2.1") 





########################################################################################################################
# Change manual settings according to machine

manual.opts.emp <- list("path"=getwd() 
                    #string, directory with all files from the github repository
)


########################################################################################################################
# Other options 

.global.opts.emp <- list(
  manual.opts = manual.opts.emp,
  general.opts = list("alpha"=0.05), # miscoverage level
  tune.opts = list("tune_reps"=10, # number of additional training sets used for tuning (see Chapter 12.2)
                   "excess_error_threshold"=0.001, # fraction of OOB error that the Monte Carlo error is allowed to take
                   "num.trees.start"=formals(grf::regression_forest)$num.trees, # starting number of trees when selecting
                                                                                # optimal number of trees
                   "num.trees.increment"=200), # tree increment
  interval.opts = list("intervals"=c("grf", "scp", "sdcp", "wscp", "boot", "oob", "woob"), # interval types
                       "cp_subsampling_rate"=0.5, # split conformal true training size share of total training data
                       "bootstrap_reps"=499, # number of bootstrap replications 
                       "dcp_taus"=seq(0,1,length=202)[2:201]) # estimated quantiles for SDCP
)



########################################################################################################################
########################################################################################################################

set.seed(2022)


# Load functions
source(file.path(.global.opts.emp$manual.opts$path, "functions.R"))


# create directory and download data if they dont already exist
if(!dir.exists(file.path(.global.opts.emp$manual.opts$path, "00_data")))
  dir.create(file.path(.global.opts.emp$manual.opts$path, "00_data"))

if(!file.exists(file.path(.global.opts.emp$manual.opts$path, "00_data", "CommViolPredUnnormalizedData.txt"))){
  download.file("https://archive.ics.uci.edu/ml/machine-learning-databases/00211/CommViolPredUnnormalizedData.txt", 
                file.path(.global.opts.emp$manual.opts$path, "00_data", "CommViolPredUnnormalizedData.txt"))
}


# read data
data <- read.csv(file.path(.global.opts.emp$manual.opts$path, "00_data", "CommViolPredUnnormalizedData.txt"), header=F)

#reformat missing values
data[which(data=="?", arr.ind=T)] <- NA 


# drop "non-predictive" columns
data[,c(1:5)] <- NA


# drop not used "goal" columns
data[,c(130:145,147)] <- NA


# remove columns with many NAs (non-predictive, unused goals, data about )
col_rm <- unname(is.na(data) |>
  apply(2, sum) > 1500) 

data <- data[,!col_rm]
Y <- data[,ncol(data)]
X <- data[,-ncol(data)]


# remove NAs in rows
na_ind <- which(is.na(Y))
X <- X[-na_ind,]
Y <- Y[-na_ind]


# formatting
X <- apply(X, 2, as.numeric) |>
  unname()
X <- scale(X)
attr(X,"scaled:scale") <- NULL
attr(X,"scaled:center") <- NULL
Y <- as.numeric(Y)
Y <- scale(Y)
attr(Y,"scaled:scale") <- NULL
attr(Y,"scaled:center") <- NULL


# Cross validataion indexes
cv_ind <- sample(c(1:5), length(Y), r=T)


results <- foreach(i=iterators::icount(5)) %do% {
  
  # train/test data
  X_train <- X[which(cv_ind!=i),]
  Y_train <- Y[which(cv_ind!=i)]
  X_new <- X[which(cv_ind==i),]
  Y_new <- Y[which(cv_ind==i)]
  
  # indices for conformal prediction intervals
  ind <- seq_along(Y_train)
  ind_train <- sample(ind, ceiling(.global.opts.emp$interval.opts$cp_subsampling_rate * length(Y_train)))
  ind_calib <- which(!ind %in% ind_train)
  
  s <- length(which(cv_ind!=i))^-0.1
  
  # tune and train models (this is kind of inefficient, beacuse some forests are trained twice, but I just reused most
  # of the simulation code)
  
  # Boot, OOB, WOOB
  model <- grf::regression_forest(X_train, Y_train,
                         tune.parameters = c("mtry", "min.node.size", "honesty.fraction"),
                         sample.fraction = s,
                         num.trees=2000)
  
  res <- cbind.data.frame("iteration"=i,
                          "mtry"=model$tuning.output$params["mtry"],
                          "min.node.size"=model$tuning.output$params["min.node.size"],
                          "honesty.fraction"=model$tuning.output$params["honesty.fraction"],
                          "num.trees"=2000,
                          "error"=Inf)
  
  tree_res <- tune_trees(grf::regression_forest, grf:::predict.regression_forest,
                         model=model, X=X_train, Y=Y_train, params = res,
                         threshold = .global.opts.emp$tune.opts$excess_error_threshold,
                         increment = .global.opts.emp$tune.opts$num.trees.increment)
  
  res$num.trees <- tree_res$num.trees
  res$error <- tree_res$mce
  
  if(res$num.trees>2000){
    model_all <- grf::regression_forest(X_train, Y_train, sample.fraction = s, num.trees=res$num.trees,
                                        mtry=res$mtry, min.node.size = res$min.node.size,
                                        honesty.fraction = res$honesty.fraction,
                                        ci.group.size = 1)
  } else model_all <- model
  boot_param <- res
  
  # SCP, WSCP
  model <- grf::regression_forest(X_train[ind_train,], Y_train[ind_train],
                                     tune.parameters = c("mtry", "min.node.size", "honesty.fraction"),
                                     sample.fraction = s,
                                     num.trees=2000)
  
  res <- cbind.data.frame("iteration"=i,
                          "mtry"=model$tuning.output$params["mtry"],
                          "min.node.size"=model$tuning.output$params["min.node.size"],
                          "honesty.fraction"=model$tuning.output$params["honesty.fraction"],
                          "num.trees"=2000,
                          "error"=Inf)
  
  tree_res <- tune_trees(grf::regression_forest, grf:::predict.regression_forest,
                         model=model, X=X_train[ind_train,], Y=Y_train[ind_train], params = res,
                         threshold = .global.opts.emp$tune.opts$excess_error_threshold,
                         increment = .global.opts.emp$tune.opts$num.trees.increment)
  
  res$num.trees <- tree_res$num.trees
  res$error <- tree_res$mce
  
  if(res$num.trees>2000){
    model_cp <- grf::regression_forest(X_train, Y_train, sample.fraction = s, num.trees=res$num.trees,
                                        mtry=res$mtry, min.node.size = res$min.node.size,
                                        honesty.fraction = res$honesty.fraction,
                                        ci.group.size = 1)
  } else model_cp <- model
  
  
  # GRF
  qf_tuning_results <- tune_forest(X_train, Y_train, grf::quantile_forest, interval.types = "grf", opts=.global.opts.emp,
                                   sample.fraction=s, cv.folds = 5,
                                   tune.parameters=c("mtry", "min.node.size", "honesty.fraction"),
                                   tune.num.reps = 100, length.penalty = 10, num.trees=500,
                                   list("quantiles"=c(.global.opts.emp$general.opts$alpha, 1-.global.opts.emp$general.opts$alpha)))
  
  model_qf <- grf::quantile_forest(X_train, Y_train, num.trees = 2500, mtry=qf_tuning_results$grf["mtry"],
                                   min.node.size = qf_tuning_results$grf["min.node.size"], 
                                   honesty.fraction = qf_tuning_results$grf["honesty.fraction"],
                                   quantiles=c(.global.opts.emp$general.opts$alpha/2, 1-.global.opts.emp$general.opts$alpha/2)) 
  
  
  # SDCP
  qf_dcp_tuning_results <- tune_forest(X_train, Y_train, grf::quantile_forest, interval.types = "sdcp", opts=.global.opts.emp,
                                   sample.fraction=s, cv.folds = 5,
                                   tune.parameters=c("mtry", "min.node.size", "honesty.fraction"),
                                   tune.num.reps = 100, length.penalty = 10, num.trees=500,
                                   list("quantiles"=.global.opts.emp$interval.opts$dcp_taus))
  
  model_qf_dcp <- grf::quantile_forest(X_train, Y_train, num.trees = 2500, mtry=qf_dcp_tuning_results$sdcp["mtry"],
                                   min.node.size = qf_dcp_tuning_results$sdcp["min.node.size"], 
                                   honesty.fraction = qf_dcp_tuning_results$sdcp["honesty.fraction"],
                                   quantiles=.global.opts.emp$interval.opts$dcp_taus)
  
  
  
  # get interval results
  GRF_results <- GRF_interval(model_qf, X_new, Y_new)
  
  SCP_results <- SCP_interval(model_cp, X_train, Y_train, ind_calib, X_new, Y_new, .global.opts.emp$general.opts$alpha)
  
  SDCP_results <- SDCP_interval(model_qf_dcp, X_train, Y_train, ind_calib, X_new, Y_new, 
                                .global.opts.emp$general.opts$alpha, .global.opts.emp$interval.opts$dcp_taus)
  
  WSCP_results <- WSCP_interval(model_cp, X_train, Y_train, ind_calib, X_new, Y_new, .global.opts.emp$general.opts$alpha)
  
  bootstrap_results <- bootstrap_interval(model_all, X_train, Y_train, X_new, Y_new, alpha=.global.opts.emp$general.opts$alpha,
                                          bootstrap_reps=.global.opts.emp$interval.opts$bootstrap_reps, tune_params=boot_param)
  
  OOB_results <- OOB_interval(model_all, Y_train, X_new, Y_new, .global.opts.emp$general.opts$alpha)
  
  WOOB_results <- WOOB_interval(model_all, Y_train, X_new, Y_new, .global.opts.emp$general.opts$alpha)
  
  
  out <- cbind.data.frame("oob_hit"=OOB_results$hit,
                          "oob_weighted_hit"=WOOB_results$hit,
                          "cp_hit"=SCP_results$hit,
                          "cp_weighted_hit"=WSCP_results$hit,
                          "grf_hit"=GRF_results$hit,
                          "dcp_hit"=SDCP_results$hit,
                          "boot_hit"=bootstrap_results$hit,
                          "oob_length"=OOB_results$length,
                          "oob_weighted_length"=WOOB_results$length,
                          "cp_length"=SCP_results$length,
                          "cp_weighted_length"=WSCP_results$length,
                          "grf_length"=GRF_results$length,
                          "dcp_length"=SDCP_results$length,
                          "boot_length"=bootstrap_results$length,
                          "n"=length(Y_train))
  
  rm("boot_param", "bootstrap_results", "GRF_results", "model", "model_all", "model_cp", "model_qf", "model_qf_dcp",
     "OOB_results", "qf_dcp_tuning_results", "qf_tuning_results", "res", "SCP_results", "SDCP_results", "tree_res",
     "WOOB_results", "WSCP_results", "X_new", "X_train", "s", "Y_new", "Y_train")
  
  out
}



# Results to table:
res <- lapply(results, function(x) apply(x, 2, mean))

# Table 1 and 2 contents in empirical study section
res2table(res, "hit")
res2table(res, "length")






