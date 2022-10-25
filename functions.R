########################################################################################################################
#                                                                                                                      #
# Utility functions                                                                                                    #
#                                                                                                                      #
########################################################################################################################


# Tune trees (corresponding to steps (a6)--(a8) described in section 8.2)
tune_trees <- function(forest.fun, pred.fun, model, X, Y, params, threshold, increment){
  
  mce <- pred.fun(model, estimate.variance=T) |>
    mc.error()
  
  while(mce > threshold){
    model_add <- forest.fun(X, Y, num.trees = increment,
                            sample.fraction = length(Y)^-0.2,
                            mtry=params$mtry,
                            min.node.size = params$min.node.size,
                            honesty.fraction = params$honesty.fraction)
    model <- grf::merge_forests(list(model, model_add))
    mce <- pred.fun(model, estimate.variance=T) |>
      mc.error()
  }
  list("num.trees"=model$`_num_trees`,
       "mce"=mce)
}


# Tune forest based on interval performance (corresponding to steps (b1)--(b7) described in section 8.2), but written
# more generally such that other intervals COULD be tuned with it as well 
tune_forest <- function(X, Y, forest.fun, cv.folds=5, tune.parameters=c("mtry", "min.node.size", "honesty.fraction"), 
                        interval.types=c("grf", "scp", "dscp", "boot", "oob", "woob"), tune.num.reps = 100, opts, 
                        length.penalty=20, num.trees=500, ...){
  
  fit.parameters <- list("num.trees" = num.trees,
                         "compute.oob.predictions" = T,
                         "num.threads" = 2)
  
  num.params <- length(tune.parameters)
  data <- grf:::create_train_matrices(X, outcome = Y, sample.weights = NULL)
  
  draw.parameters <- runif(tune.num.reps * num.params) |>
    matrix(tune.num.reps, num.params, dimnames=list(NULL, tune.parameters)) |>
    apply(1, function(d) grf:::get_params_from_draw(nrow(X), ncol(X), d))
  
  cv_ind <- sample(1:cv.folds, length(Y), r=T)
  
  
  cv_results <- foreach(i = icount(cv.folds)) %do% {
    
    X_train <- X[which(cv_ind != i),]
    X_new <- X[which(cv_ind == i),]
    
    Y_train <- Y[which(cv_ind != i)]
    Y_new <- Y[which(cv_ind == i)]
    
    if(any(c("sdcp", "scp", "wscp") %in% interval.types)){
      
      ind <- seq_along(Y_train)
      ind_train <- sample(ind, ceiling(opts$interval.opts$cp_subsampling_rate * length(Y_train)))
      ind_calib <- which(!ind %in% ind_train)
      
      models <- apply(draw.parameters, 2, function(d) do.call(forest.fun, c(fit.parameters, d, "X"=list(X_train[ind_train,]), "Y"=list(Y_train[ind_train]), ...)))
      
    } else {
      models <- apply(draw.parameters, 2, function(d) do.call(forest.fun, c(fit.parameters, d, "X"=list(X_train), "Y"=list(Y_train), ...)))
    }
    
    interval_results <- foreach(m = iterators::icount(length(models)), .combine="rbind") %do% { 
      
      results <- array(,dim=c(1, 2*length(interval.types))) |>
        as.data.frame()
      name_grid <- expand.grid(interval.types, c("coverage", "length"))
      colnames(results) <- sprintf("%s_%s", name_grid[,1], name_grid[,2])
      
      if("grf" %in% interval.types) {
        res <- GRF_interval(models[[m]], X_new, Y_new)
        results$grf_coverage <- mean(res$hit, na.rm=T)-1+opts$general.opts$alpha
        results$grf_length <- mean(res$length, na.rm=T)
      }
      
      if("scp" %in% interval.types){
        res <- SCP_interval(models[[m]], X_train, Y_train, ind_calib, X_new, Y_new, opts$general.opts$alpha)
        results$scp_coverage <- mean(res$hit, na.rm=T)-1+opts$general.opts$alpha
        results$scp_length <- mean(res$length, na.rm=T)
      }
      
      if("sdcp" %in% interval.types) {
        res <- SDCP_interval(models[[m]], X_train, Y_train, ind_calib, X_new, Y_new,
                             opts$general.opts$alpha, opts$interval.opts$dcp_taus)
        results$sdcp_coverage <- mean(res$hit, na.rm=T)-1+opts$general.opts$alpha
        results$sdcp_length <- mean(res$length, na.rm=T)
      }
      
      if("wscp" %in% interval.types){
        res <- WSCP_interval(models[[m]], X_train, Y_train, ind_calib, X_new, Y_new, opts$general.opts$alpha)
        results$wscp_coverage <- mean(res$hit, na.rm=T)-1+opts$general.opts$alpha
        results$wscp_length <- mean(res$length, na.rm=T)
      }
      
      if("boot" %in% interval.types){
        stop("This takes too much time to be worth it unless you lower the settings so much that its worthless.")
      }
      
      if("oob" %in% interval.types) {
        res <- OOB_interval(models[[m]], Y_train, X_new, Y_new, opts$general.opts$alpha)
        results$oob_coverage <- mean(res$hit, na.rm=T)-1+opts$general.opts$alpha
        results$oob_length <- mean(res$length, na.rm=T)
      }
      
      if("woob" %in% interval.types){
        res <- WOOB_interval(models[[m]], Y_train, X_new, Y_new, opts$general.opts$alpha)
        results$woob_coverage <- mean(res$hit, na.rm=T)-1+opts$general.opts$alpha
        results$woob_length <- mean(res$length, na.rm=T)
      }
      results
    }
    interval_results
  }
  
  c_names <- colnames(cv_results[[1]])
  results <- array(unlist(cv_results), dim =c(length(models), 2*length(interval.types), cv.folds))
  results <- apply(results, c(1,2), mean)
  colnames(results) <- c_names
  
  best <- vector("list", length(interval.types))
  names(best) <- interval.types
  
  for(i in 1:length(interval.types)){
    res <- results[,c(i,length(interval.types)+i)]
    res[,2] <- res[,2]/(length.penalty*var(Y))
    res[which(res[,1]>0),1] <- res[which(res[,1]>0),1]/1.5
    
    best[[i]] <- c(draw.parameters[,which.min(apply(res, 1, dist))], min(apply(res, 1, dist)))
    names(best[[i]]) <- c(tune.parameters, "dist")
    
  }
  best
}



# Weighted quantiles
w.quantile <- function(x, w, probs = c(0.025, 0.975)){
  ind <- 1:length(x)
  ord <- order(x)
  x <- x[ord]
  
  reps <- max(1, nrow(w))
  out <- foreach(i=iterators::icount(reps), .combine = "rbind") %do% {
    w_ord <- w[i,ord]
    w_cdf <- cumsum(w_ord)/sum(w_ord)
    
    lower_bound_ind <- unlist(lapply(probs, function(p) max(max(ind[w_cdf<=p],1)))) 
    upper_bound_ind <- unlist(lapply(probs, function(p) min(ind[w_cdf>p]))) 
    (w_ord[lower_bound_ind]*x[lower_bound_ind]+w_ord[upper_bound_ind]*x[upper_bound_ind])/(w_ord[lower_bound_ind]+w_ord[upper_bound_ind])
  }
  rownames(out) <- NULL
  out
}



# generate multivariate distributions with copula
dist_generator <- function(p, .dispstr, .par=0.6, .margins, .paramMargins){
  if(length(.margins)==1 && 
     class(unlist(.paramMargins, recursive=F))!="list") {
    .margins <- rep(.margins, p)
    .paramMargins <- rep(list(.paramMargins), p)
  }
  if(.dispstr=="uncor") cop <- indepCopula(p) else{
    cop <- normalCopula(param=.par, dim = p, dispstr = .dispstr)
  }
  
  mvdc(copula=cop, margins=.margins,
       paramMargins=.paramMargins)
}



# generate a data generating process using an object returned by dist_generator()
DGP_generator <- function(...){
  args <- list(...)
  function(n){
    foreach::foreach(i=iterators::icount(length(args)), .combine="cbind") %do% {rMvdc(n, args[[i]])}
  }
}



# generate a mean function
mean_fun_generator <- function(.body, .args="X"){
  eval(parse(text = paste('function(', .args, ') { if(!"matrix" %in% class(X)) {X <- matrix(X, nrow=1)} 
                          return(' , .body , ')}', sep='')))
}



# generate an error distribution
err_fun_generator <- function(preset, .args, .dist, adjustment=1){
  
  function(X, type=c("r", "q"), dist=.dist, ...){
    
    .fun <- eval(parse(text=paste0(type, dist)))
    
    if(substr(type, 1, 1)=="q"){
      return(foreach(i=icount(nrow(X)), .combine="c") %do% {
        args <- .args
        for(j in 1:length(args)) {if(class(args[[j]])=="character") args[[j]] <- eval(parse(text=args[[j]]))}
        do.call(.fun, c(list(...), args))*adjustment
      })
      
    } else if(substr(type, 1, 1)=="r"){
      return(foreach(i=icount(nrow(X)), .combine="c") %do% {
        args <- .args
        for(j in 1:length(args)) {if(class(args[[j]])=="character") args[[j]] <- eval(parse(text=args[[j]]))}
        do.call(.fun, c(list(n=1), args))*adjustment
      })
    }
  }
}



# progress bar because its fancy
progress_bar <- function(.max){
  progress_bar <- txtProgressBar(max=.max, style=3)
  progress <- function(iter) setTxtProgressBar(progress_bar, iter)
  
  list(progress=progress)
}



# makeshift function to produce latex table contents. easier than copying everything by hand
res2table <- function(res, type=c("hit", "length"), 
                      order=c("grf", "cp", "dcp", "cp_weighted", "boot", "oob", "oob_weighted"),
                      print.mean=T, .digits=4, .format="f", ...){
  out_num <- vector("list", length(res))
  out_pretty <- vector("list", length(res))
  for(i in 1:length(res)){
    dat <- res[[i]][grepl(paste0("_", type, "$"), names(res[[1]]))]
    
    if(length(order)==length(dat)){
      ord <- lapply(order, function(j) which(grepl(paste0("^", j, "_", type), names(dat)))) |>
        unlist()
      
      if(length(unique(ord))!=length(order)) stop("non-unique order")
    } else {
      stop("order and data have unequal lengths")
    }
    
    dat_pretty <- formatC(dat, digits=.digits, format=.format, ...)
    
    out_num[[i]] <- dat[ord]
    out_pretty[[i]] <- paste(paste(i, paste(dat_pretty[ord], collapse=" & "), sep=" & "), "\\\\")
  }
  
  n <- unlist(lapply(res, function(r) r["n"]))
  
  means <- simplify2array(out_num) |>
    t() |>
    apply(2, function(x) weighted.mean(x, n))
  
  out_pretty <- unlist(out_pretty) |>
    paste(collapse="\n ") 
  
  if(print.mean){
    means_pretty <- formatC(means, digits=.digits, format=.format, ...)
    means_pretty <- paste(paste("Mean", paste(means_pretty, collapse=" & "), sep=" & "), "\\\\")
    
    out <- paste(out_pretty, "\\midrule", means_pretty, sep="\n ")
  } else {
    out <- out_pretty
  }
  
  cat(out)
}



# compute monte carlo error threshold
mc.error <- function(preds){
  mpe <- mean(preds$debiased.error+preds$excess.error)
  mean(preds$excess.error)/mpe
}



# make plots pretty =)
plot.aesthetics <- function(type=c("coverage", "length"), data, opts, max.y.override=NULL){
  
  if(type=="coverage"){
    out <- list(geom_hline(yintercept=1-.global.opts$general.opts$alpha, size=1.3, col="gray"),
                scale_y_continuous(limits=c(min(data$hit)*0.99,1),
                                   expand=c(0,1/exp(30)),
                                   breaks=pretty(c(0,1),20),
                                   name="Coverage rate"),
                theme(plot.margin = unit(c(0,0.5,0,0), "cm")))
  } else if(type=="length"){
    out <- list(geom_hline(yintercept=0, size=1.3, col="gray"),
                scale_y_continuous(limits=c(min(-0.1, min(data$ln_rel_length)), ifelse(is.null(max.y.override), max(data$ln_rel_length), max.y.override)),
                                   expand=c(0,1/exp(30)),
                                   breaks=pretty(c(min(data$ln_rel_length),max(data$ln_rel_length)),9),
                                   labels=formatC(pretty(c(min(data$ln_rel_length),max(data$ln_rel_length)),9), digits=2, format="f"),
                                   name="Natural logarithm of relative interval length"),
                theme(plot.margin = unit(c(0,0,0,0.5), "cm")))
  } else stop("undefined type")
  
  list(out,
       scale_x_discrete(expand=c(0,0.5),
                        name="",
                        guide = guide_axis(angle = -45)),
       scale_fill_discrete(name="Number of observations:"),
       theme_bw(),
       theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
             legend.box.background = element_rect(colour="black"),
             legend.background = element_blank(),
             legend.spacing.y = unit(0, "mm"),
             legend.position = "bottom", 
             legend.direction = "horizontal",
             text=element_text(size=26),
             legend.key.size = unit(1, "cm"),
             legend.spacing.x = unit(0, "mm")),
       geom_boxplot())
}





########################################################################################################################
#                                                                                                                      #
# functions to construct all intervals                                                                                 #
#                                                                                                                      #
# Note that models arent trained inside the function but given as an argument to avoid training the same model         #
# multiple times if it can be shared aross multiple models (e.g. OOB and WOOB interval require an identical model and  #
# by giving it as an argument we only need to train it once, not twice)                                                #
#                                                                                                                      #
########################################################################################################################


GRF_interval <- function(model_qf, X_new, Y_new){
  
  qf_preds <- grf:::predict.quantile_forest(model_qf, newdata=X_new)[[1]]
  PI_grf_hit <- ifelse(qf_preds[,1] <= Y_new & qf_preds[,2] >= Y_new, T, F)
  PI_grf_length <- abs(qf_preds[,2] - qf_preds[,1])
  
  list("hit"=PI_grf_hit,
       "length"=PI_grf_length)
}



SCP_interval <- function(model_cp, X_train, Y_train, ind_calib, X_new, Y_new, alpha){
  
  pred_calib <- grf:::predict.regression_forest(model_cp, newdata=X_train[ind_calib,])[[1]]
  resid_calib <- Y_train[ind_calib] - pred_calib
  
  pred_new_cp <- grf:::predict.regression_forest(model_cp, newdata=X_new)[[1]]
  
  boundary_cp <- unname(quantile(resid_calib, probs=c(alpha/2, 1-alpha/2)))
  PI_cp_hit <- ifelse(pred_new_cp+boundary_cp[1] <= Y_new & pred_new_cp+boundary_cp[2] >= Y_new, T, F)
  PI_cp_length <- abs(boundary_cp[2] - boundary_cp[1]) 
  
  list("hit"=PI_cp_hit,
       "length"=PI_cp_length)
}



SDCP_interval <- function(model_qf_dcp, X_train, Y_train, ind_calib, X_new, Y_new, alpha, dcp_taus){
  
  pred_qf_calib <- grf:::predict.quantile_forest(model_qf_dcp, newdata=X_train[ind_calib,])[[1]]
  pred_qf_test <- grf:::predict.quantile_forest(model_qf_dcp, newdata=X_new)[[1]]
  
  u.hat <- rowMeans((pred_qf_calib <= matrix(Y_train[ind_calib],length(Y_train[ind_calib]),length(dcp_taus))))
  
  bhat <- rep(NA,length(Y_train[ind_calib]))
  b.grid <- dcp_taus[dcp_taus<=alpha]
  
  for (t in 1:length(Y_train[ind_calib])){
    leng <- rep(NA,length(b.grid))
    leng.test <- rep(NA,length(b.grid))
    for (b in 1:length(b.grid)){
      Q.yx.u <- approx(x=dcp_taus,y=pred_qf_calib[t,],xout=(b.grid[b]+1-alpha),rule=2)$y
      leng[b] <- Q.yx.u -pred_qf_calib[t,b]
    }
    bhat[t] <- b.grid[which.min(leng)]
  }
  
  cs.opt <- abs(u.hat-bhat-(1-alpha)/2)
  k           <- ceiling((1-alpha)*(1+length(Y_train[ind_calib])))
  threshold   <- sort(cs.opt)[k]
  cov.opt   <- cs.opt <= threshold
  
  PI_dcp_length <- vector("numeric", length(Y_new))
  PI_dcp_hit <- vector("logical", length(Y_new))
  
  for (t in 1:length(Y_new)){
    ci.grid <- abs(dcp_taus - bhat[t]-(1-alpha)/2)
    ci <- pred_qf_calib[t,(ci.grid <= threshold)]
    ub <- max(ci)
    lb <- min(ci)
    PI_dcp_length[t] <- abs(ub-lb)
    PI_dcp_hit[t] <- ifelse(lb <= Y_new[t] & ub >= Y_new[t], T, F)
  }
  PI_dcp_length[which(PI_dcp_length==-Inf)] <- NA

  list("hit"=PI_dcp_hit,
       "length"=PI_dcp_length)
}



WSCP_interval <- function(model_cp, X_train, Y_train, ind_calib, X_new, Y_new, alpha){
  
  pred_calib <- grf:::predict.regression_forest(model_cp, newdata=X_train[ind_calib,])[[1]]
  resid_calib <- Y_train[ind_calib] - pred_calib
  
  pred_new_cp <- grf:::predict.regression_forest(model_cp, newdata=X_new)[[1]]
  weights_cp <- grf::get_forest_weights(model_cp, newdata=X_new)
  
  boundary_cp_weighted <- w.quantile(resid_calib, weights_cp, probs=c(alpha/2, 1-alpha/2))
  PI_cp_weighted_hit <- ifelse(pred_new_cp+boundary_cp_weighted[,1] <= Y_new & pred_new_cp+boundary_cp_weighted[,2] >= Y_new, T, F)
  PI_cp_weighted_length <- abs(boundary_cp_weighted[,2] - boundary_cp_weighted[,1])
  
  list("hit"=PI_cp_weighted_hit,
       "length"=PI_cp_weighted_length)
}



bootstrap_interval <- function(model_all, X_train, Y_train, X_new, Y_new, alpha, bootstrap_reps, tune_params){
  
  pred_oob <- grf:::predict.regression_forest(model_all)
  resid_oob <- Y_train-pred_oob[,1]
  
  pred_new <- grf:::predict.regression_forest(model_all, newdata=X_new)
  
  bootstrap_ind <- lapply(1:bootstrap_reps, function(i) sample(1:length(Y_train), replace=T))
  bootstrap_ind_leftout <- lapply(bootstrap_ind, function(i) setdiff(1:length(Y_train), i))
  
  boot_x <- lapply(bootstrap_ind, function(i) X_train[i,])
  boot_y <- lapply(bootstrap_ind, function(i) Y_train[i])
  
  boot_preds <- vector("list", bootstrap_reps)
  for(i in 1:bootstrap_reps){
    boot_preds[[i]] <- grf:::predict.regression_forest(grf::regression_forest(boot_x[[i]], boot_y[[i]],
                                                                              mtry = tune_params$mtry,
                                                                              min.node.size = tune_params$min.node.size,
                                                                              honesty.fraction = tune_params$honesty.fraction,
                                                                              sample.fraction = length(Y_train)^-0.2,
                                                                              num.threads=2,
                                                                              num.trees=500,
                                                                              ci.group.size = 1), newdata=X_new)[[1]]
  }
  
  boot_preds_array <- do.call(cbind, boot_preds)
  
  boot_pred_avg <- apply(boot_preds_array, 1, mean)
  boot_pred_centered <- sweep(boot_preds_array, 1, boot_pred_avg)
  boot_pred_centered_list <- asplit(boot_pred_centered, 1)
  
  boot_quantiles <- lapply(boot_pred_centered_list, function(m) unname(quantile(apply(expand.grid(m, resid_oob),1,sum), 
                                                                                probs = c(alpha/2, 1-alpha/2))))
  
  bootstrap_intervals <- do.call(rbind, boot_quantiles)
  PI_boot_hit <- ifelse(pred_new[[1]] + bootstrap_intervals[,1] <= Y_new & pred_new[[1]] + bootstrap_intervals[,2] >= Y_new, T, F)
  PI_boot_length <- abs(bootstrap_intervals[,2] - bootstrap_intervals[,1])
  
  
  list("hit"=PI_boot_hit,
       "length"=PI_boot_length)
  
}



OOB_interval <- function(model_all, Y_train, X_new, Y_new, alpha){
  
  pred_oob <- grf:::predict.regression_forest(model_all)
  resid_oob <- Y_train-pred_oob[,1]
  
  pred_new <- grf:::predict.regression_forest(model_all, newdata=X_new)[[1]]
  
  boundary_oob <- unname(quantile(resid_oob, probs = c(alpha/2, 1-alpha/2)))
  PI_oob_hit <- ifelse(pred_new+boundary_oob[1] <= Y_new & pred_new+boundary_oob[2] >= Y_new, T, F)
  PI_oob_length <- abs(boundary_oob[2] - boundary_oob[1]) 
  
  list("hit"=PI_oob_hit,
       "length"=PI_oob_length)
}



WOOB_interval <- function(model_all, Y_train, X_new, Y_new, alpha){
  
  pred_oob <- grf:::predict.regression_forest(model_all)
  resid_oob <- Y_train-pred_oob[,1]
  
  pred_new <- grf:::predict.regression_forest(model_all, newdata=X_new)[[1]]
  weights <- grf::get_forest_weights(model_all, newdata=X_new)
  
  boundary_weighted_oob <- w.quantile(resid_oob, weights, probs=c(alpha/2, 1-alpha/2))
  PI_oob_weighted_hit <- ifelse(pred_new+boundary_weighted_oob[,1] <= Y_new & pred_new+boundary_weighted_oob[,2] >= Y_new, T, F)
  PI_oob_weighted_length <- abs(boundary_weighted_oob[,2] - boundary_weighted_oob[,1])
  
  list("hit"=PI_oob_weighted_hit,
       "length"=PI_oob_weighted_length)
}





########################################################################################################################
#                                                                                                                      #
# Function for the simulations                                                                                         #
#                                                                                                                      #
########################################################################################################################


# Tune forests for the simulation
simulation_tuning <- function(DataGeneratingProcess, opts) {
  tuning_results <- foreach(i=iterators::icount(opts$tune.opts$tune_reps), 
                            .options.snow=progress_bar(opts$tune.opts$tune_reps), 
                            .packages = opts$parallel.opts$packages,
                            .combine="rbind") %dopar% {
                              
                              set.seed(-i)
                              
                              if(!opts$manual.opts$parallel){
                                if(i==1) progress <- progress_bar(opts$tune$tune_reps) else 
                                  progress$progress(i)
                              }
                              
                              foreach(obs=iterators::icount(length(opts$general.opts$train_obs)), 
                                      .combine="rbind.data.frame") %do% {
                                        
                                        # Generate data
                                        N <- opts$general.opts$train_obs[obs]
                                        X_train <- DataGeneratingProcess(N)
                                        e_train <- err_fun(X_train, type="r")
                                        m_train <- mean_fun(X_train)
                                        Y_train <- m_train + e_train
                                        
                                        rm("e_train", "m_train")
                                        
                                        # samlpe split for CP
                                        ind <- seq_along(Y_train)
                                        ind_train <- sample(ind, ceiling(opts$interval.opts$cp_subsampling_rate * N))
                                        ind_calib <- which(!ind %in% ind_train)
                                        
                                        s <- length(Y_train)^-0.2
                                        
                                        ##############################################################################################
                                        # Actual tuning 
                                        
                                        # bootstrap, OOB, WOOB intervals
                                        model_all <- grf::regression_forest(X_train, Y_train,
                                                                            tune.parameters = c("mtry", "min.node.size", "honesty.fraction"),
                                                                            sample.fraction = s,
                                                                            num.trees=2000)
                                        
                                        res <- cbind.data.frame("iteration"=i,
                                                                "type"="",
                                                                "N"=opts$general.opts$train_obs[obs],
                                                                "mtry"=model_all$tuning.output$params["mtry"],
                                                                "min.node.size"=model_all$tuning.output$params["min.node.size"],
                                                                "honesty.fraction"=model_all$tuning.output$params["honesty.fraction"],
                                                                "num.trees"=2000,
                                                                "error"=Inf)
                                        
                                        tree_res <- tune_trees(grf::regression_forest, grf:::predict.regression_forest,
                                                               model=model_all, X=X_train, Y=Y_train, params = res,
                                                               threshold = opts$tune.opts$excess_error_threshold,
                                                               increment = opts$tune.opts$num.trees.increment)
                                        
                                        res$num.trees <- tree_res$num.trees
                                        res$error <- tree_res$mce
                                        
                                        res_boot <- res
                                        res_boot$type <- "boot"
                                        res_oob <- res
                                        res_oob$type <- "oob"
                                        res_woob <- res
                                        res_woob$type <- "woob"
                                        
                                        # SCP, WSCP intervals
                                        model_cp <- grf::regression_forest(X_train[ind_train,], Y_train[ind_train],
                                                                           tune.parameters = c("mtry", "min.node.size", "honesty.fraction"),
                                                                           sample.fraction = s,
                                                                           num.trees=2000)
                                        
                                        res <- cbind.data.frame("iteration"=i,
                                                                "type"="",
                                                                "N"=opts$general.opts$train_obs[obs],
                                                                "mtry"=model_cp$tuning.output$params["mtry"],
                                                                "min.node.size"=model_cp$tuning.output$params["min.node.size"],
                                                                "honesty.fraction"=model_cp$tuning.output$params["honesty.fraction"],
                                                                "num.trees"=2000,
                                                                "error"=Inf)
                                        
                                        tree_res <- tune_trees(grf::regression_forest, grf:::predict.regression_forest,
                                                               model=model_cp, X=X_train[ind_train,], Y=Y_train[ind_train], params = res,
                                                               threshold = opts$tune.opts$excess_error_threshold,
                                                               increment = opts$tune.opts$num.trees.increment)
                                        
                                        res$num.trees <- tree_res$num.trees
                                        res$error <- tree_res$mce
                                        
                                        res_scp <- res
                                        res_scp$type <- "scp"
                                        res_wscp <- res
                                        res_wscp$type <- "wscp"
                                        
                                        # GRF interval
                                        qf_tuning_results <- tune_forest(X_train, Y_train, grf::quantile_forest, interval.types = "grf", opts=opts,
                                                                         sample.fraction=s, cv.folds = 5,
                                                                         tune.parameters=c("mtry", "min.node.size", "honesty.fraction"),
                                                                         tune.num.reps = 100, length.penalty = 10, num.trees=500,
                                                                         list("quantiles"=c(opts$general.opts$alpha/2, 1-opts$general.opts$alpha/2)))
                                        
                                        res_grf <- cbind.data.frame("iteration"=i,
                                                                    "type"="grf",
                                                                    "N"=opts$general.opts$train_obs[obs],
                                                                    "mtry"=qf_tuning_results$grf["mtry"],
                                                                    "min.node.size"=qf_tuning_results$grf["min.node.size"],
                                                                    "honesty.fraction"=qf_tuning_results$grf["honesty.fraction"],
                                                                    "num.trees"=2500,
                                                                    "error"=qf_tuning_results$grf["dist"])
                                        
                                        # SDCP interval
                                        qf_dcp_tuning_results <- tune_forest(X_train, Y_train, grf::quantile_forest, interval.types = "sdcp", opts=opts,
                                                                             sample.fraction=s, cv.folds = 5,
                                                                             tune.parameters=c("mtry", "min.node.size", "honesty.fraction"),
                                                                             tune.num.reps = 100, length.penalty = 10, num.trees=500,
                                                                             list("quantiles"=opts$interval.opts$dcp_taus))
                                        
                                        res_sdcp <- cbind.data.frame("iteration"=i,
                                                                     "type"="sdcp",
                                                                     "N"=opts$general.opts$train_obs[obs],
                                                                     "mtry"=qf_dcp_tuning_results$sdcp["mtry"],
                                                                     "min.node.size"=qf_dcp_tuning_results$sdcp["min.node.size"],
                                                                     "honesty.fraction"=qf_dcp_tuning_results$sdcp["honesty.fraction"],
                                                                     "num.trees"=2500,
                                                                     "error"=qf_dcp_tuning_results$sdcp["dist"])
                                        
                                        
                                        # output
                                        out <- rbind(res_boot, res_oob, res_woob, res_scp, res_wscp, res_grf, res_sdcp)
                                        rownames(out) <- NULL
                                        out
                                      }
                            }
  
  # Formatting
  results <- split(tuning_results, tuning_results$type) |>
    lapply(function(x) split(x, x$N)) |>
    lapply(function(r) lapply(r, function(rr) rr[which.min(rr$error),]))
  
  results
}


# Actual simulation
simulation <- function(DataGeneratingProcess, opts, tune_param=NULL) {
  sim_results <- foreach(rep=iterators::icount(opts$general.opts$reps), 
                         .packages = opts$parallel.opts$packages,
                         .options.snow=progress_bar(opts$general.opts$reps)) %dopar% {
                           
                           set.seed(rep)
                           
                           if(!opts$manual.opts$parallel){
                             if(rep==1) progress <- progress_bar(opts$general.opts$reps) else progress$progress(rep)
                           }
                           
                           foreach(obs=iterators::icount(length(opts$general.opts$train_obs))) %do% { 
                             
                             # Generate data
                             N <- opts$general.opts$train_obs[obs]
                             X_train <- DataGeneratingProcess(N)
                             e_train <- err_fun(X_train, type="r")
                             m_train <- mean_fun(X_train)
                             Y_train <- m_train + e_train
                             
                             rm("e_train", "m_train")
                             
                             # samlpe split for CP
                             ind <- seq_along(Y_train)
                             ind_train <- sample(ind, ceiling(opts$interval.opts$cp_subsampling_rate * N))
                             ind_calib <- which(!ind %in% ind_train)
                             
                             # Training parameters
                             train_params <- vector("list", length(opts$interval.opts$intervals))
                             names(train_params) <- opts$interval.opts$intervals
                             for(int in opts$interval.opts$intervals){
                               train_params[[int]] <- list("num.trees"=ifelse(is.null(tune_param), 
                                                                              2000, 
                                                                              tune_param[[int]][[paste(opts$general.opts$train_obs[obs])]]["num.trees"][[1]]),
                                                           "mtry"=ifelse(is.null(tune_param), 
                                                                         min(ceiling(sqrt(ncol(X_train)) + 20), ncol(X_train)), 
                                                                         tune_param[[int]][[paste(opts$general.opts$train_obs[obs])]]["mtry"][[1]]),
                                                           "min.node.size"=ifelse(is.null(tune_param), 
                                                                                  5, 
                                                                                  tune_param[[int]][[paste(opts$general.opts$train_obs[obs])]]["min.node.size"][[1]]),
                                                           "honesty.fraction"=ifelse(is.null(tune_param), 
                                                                                     0.5, 
                                                                                     tune_param[[int]][[paste(opts$general.opts$train_obs[obs])]]["honesty.fraction"][[1]]))
                             }
                             
                             s <- N^-0.2
                             
                             # bootstrap, OOB, WOOB intervals
                             model_all <- grf::regression_forest(X_train, Y_train,
                                                                 num.trees = train_params$boot$num.trees,
                                                                 mtry = train_params$boot$mtry,
                                                                 min.node.size = train_params$boot$min.node.size,
                                                                 honesty.fraction = train_params$boot$honesty.fraction,
                                                                 sample.fraction = s,
                                                                 ci.group.size = 1)
                             
                             
                             # SCP, WSCP intervals
                             model_cp <- grf::regression_forest(X_train[ind_train,], Y_train[ind_train],
                                                                num.trees = train_params$scp$num.trees,
                                                                mtry = train_params$scp$mtry,
                                                                min.node.size = train_params$scp$min.node.size,
                                                                honesty.fraction = train_params$scp$honesty.fraction,
                                                                sample.fraction = s,
                                                                ci.group.size = 1)
                             
                             # GRF interval
                             model_qf <- grf::quantile_forest(X_train, Y_train,
                                                              num.trees = train_params$grf$num.trees,
                                                              mtry = train_params$grf$mtry,
                                                              min.node.size = train_params$grf$min.node.size,
                                                              honesty.fraction = train_params$grf$honesty.fraction,
                                                              sample.fraction = s,
                                                              quantiles = c(opts$general.opts$alpha/2, 1-opts$general.opts$alpha/2))
                             
                             # SDCP interval
                             model_qf_dcp <- grf::quantile_forest(X_train[ind_train,], Y_train[ind_train],
                                                                  num.trees = train_params$sdcp$num.trees,
                                                                  mtry = train_params$sdcp$mtry,
                                                                  min.node.size = train_params$sdcp$min.node.size,
                                                                  honesty.fraction = train_params$sdcp$honesty.fraction,
                                                                  sample.fraction = s,
                                                                  quantiles = opts$interval.opts$dcp_taus)
                             
                             
                             #####################################
                             # generating test data
                             X_new <- DataGeneratingProcess(opts$general.opts$test_obs)
                             e_new <- err_fun(X_new, type="r")
                             m_new <- mean_fun(X_new)
                             Y_new <- m_new + e_new
                             
                             
                             # True interval
                             true_lower <- err_fun(X_new, type="q", p=opts$general.opts$alpha/2)
                             true_upper <- err_fun(X_new, type="q", p=1-opts$general.opts$alpha/2)
                             PI_t_hit <- ifelse(m_new+true_lower <= Y_new & m_new+true_upper >= Y_new, T, F)
                             PI_t_length <- abs(true_upper - true_lower)
                             
                             rm("m_new", "e_new")
                             
                             
                             
                             GRF_results <- GRF_interval(model_qf, X_new, Y_new)
                             
                             SCP_results <- SCP_interval(model_cp, X_train, Y_train, ind_calib, X_new, Y_new, opts$general.opts$alpha)
                             
                             SDCP_results <- SDCP_interval(model_qf_dcp, X_train, Y_train, ind_calib, X_new, Y_new, 
                                                           opts$general.opts$alpha, opts$interval.opts$dcp_taus)
                             
                             WSCP_results <- WSCP_interval(model_cp, X_train, Y_train, ind_calib, X_new, Y_new, opts$general.opts$alpha)
                             
                             if(rep>(opts$general.opts$reps*opts$interval.opts$boot_skip)) {
                               bootstrap_results <- bootstrap_interval(model_all, X_train, Y_train, X_new, Y_new, alpha=opts$general.opts$alpha,
                                                                       bootstrap_reps=opts$interval.opts$bootstrap_reps, 
                                                                       tune_params=train_params$boot)
                             } else {
                               bootstrap_results <- list("hit"=NA,
                                                         "length"=NA)
                             }
                             
                             OOB_results <- OOB_interval(model_all, Y_train, X_new, Y_new, opts$general.opts$alpha)
                             
                             WOOB_results <- WOOB_interval(model_all, Y_train, X_new, Y_new, opts$general.opts$alpha)
                             
                             
                             
                             ###########################################################
                             out <- cbind.data.frame("true_hit"=PI_t_hit,
                                                     "oob_hit"=OOB_results$hit,
                                                     "oob_weighted_hit"=WOOB_results$hit,
                                                     "cp_hit"=SCP_results$hit,
                                                     "cp_weighted_hit"=WSCP_results$hit,
                                                     "grf_hit"=GRF_results$hit,
                                                     "dcp_hit"=SDCP_results$hit,
                                                     "boot_hit"=bootstrap_results$hit,
                                                     "true_length"=PI_t_length,
                                                     "oob_length"=OOB_results$length,
                                                     "oob_weighted_length"=WOOB_results$length,
                                                     "cp_length"=SCP_results$length,
                                                     "cp_weighted_length"=WSCP_results$length,
                                                     "grf_length"=GRF_results$length,
                                                     "dcp_length"=SDCP_results$length,
                                                     "boot_length"=bootstrap_results$length,
                                                     "n"=N)
                             
                             out
                           }
                         }
  
  
  # Formatting
  out1 <- foreach(i=iterators::icount(length(sim_results)), .combine="rbind") %do% {
    foreach(j=iterators::icount(length(sim_results[[i]])), .combine="rbind")  %do% {
      rbind(cbind.data.frame("hit"=mean(sim_results[[i]][[j]]$true_hit), 
                             "length"=mean(sim_results[[i]][[j]]$true_length),
                             "nobs"=as.factor(sim_results[[i]][[j]]$n[1]), 
                             "type"="true"),
            cbind.data.frame("hit"=mean(sim_results[[i]][[j]]$oob_hit), 
                             "length"=mean(sim_results[[i]][[j]]$oob_length),
                             "nobs"=as.factor(sim_results[[i]][[j]]$n[1]), 
                             "type"="oob"),
            cbind.data.frame("hit"=mean(sim_results[[i]][[j]]$oob_weighted_hit), 
                             "length"=mean(sim_results[[i]][[j]]$oob_weighted_length),
                             "nobs"=as.factor(sim_results[[i]][[j]]$n[1]), 
                             "type"="oob_weighted"),
            cbind.data.frame("hit"=mean(sim_results[[i]][[j]]$cp_hit), 
                             "length"=mean(sim_results[[i]][[j]]$cp_length),
                             "nobs"=as.factor(sim_results[[i]][[j]]$n[1]), 
                             "type"="cp"),
            cbind.data.frame("hit"=mean(sim_results[[i]][[j]]$cp_weighted_hit), 
                             "length"=mean(sim_results[[i]][[j]]$cp_weighted_length),
                             "nobs"=as.factor(sim_results[[i]][[j]]$n[1]), 
                             "type"="cp_weighted"),
            cbind.data.frame("hit"=mean(sim_results[[i]][[j]]$grf_hit), 
                             "length"=mean(sim_results[[i]][[j]]$grf_length),
                             "nobs"=as.factor(sim_results[[i]][[j]]$n[1]), 
                             "type"="grf"),
            cbind.data.frame("hit"=mean(sim_results[[i]][[j]]$dcp_hit), 
                             "length"=mean(sim_results[[i]][[j]]$dcp_length),
                             "nobs"=as.factor(sim_results[[i]][[j]]$n[1]), 
                             "type"="dcp"),
            cbind.data.frame("hit"=mean(sim_results[[i]][[j]]$boot_hit),
                             "length"=mean(sim_results[[i]][[j]]$boot_length),
                             "nobs"=as.factor(sim_results[[i]][[j]]$n[1]),
                             "type"="boot")
      )
      
    }
  }
  
  # log interval length ratio
  out <- foreach(i=iterators::icount(length(unique(out1$type))), .combine="rbind") %do% {
    out_par <- out1[which(out1$type==out1$type[i]),]
    cbind(out_par, "ln_rel_length"=log(out_par$length/out1[which(out1$type=="true"),"length"]))
  }
  out
}



