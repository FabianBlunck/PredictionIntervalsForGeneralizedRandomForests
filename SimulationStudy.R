########################################################################################################################
##                                                                                                                    ##
## Code for Chapter 9 - "Simulation Study"                                                                            ##
##                                                                                                                    ##
########################################################################################################################


# Install and load packages
packages <- c("iterators", "foreach", "ggplot2", "copula", "doParallel", "doSNOW", "parallel", "snow", "grf", "ggpubr")

# install.packages(packages)

lapply(packages, require, character.only=T)


# Random numbers from R version 4.2.1
RNGversion("4.2.1") 





########################################################################################################################
# Change manual settings according to machine. Remainder of the script can then be run at once

manual.opts <- list("path"=getwd(), 
                    #string, directory with all files from the github repository
                    
                    "parallel"=T, 
                    #boolean, should the simulation be done on multiple workers?
                    # (Training the models will be - if possible - parallelised regardless)
                    
                    "ncores"=ceiling(parallel::detectCores()/2) 
                    #integer, number of cores used for the simulation. 
                    # ignored if "parallel"=F
)


########################################################################################################################
# Other options for simulations

.global.opts <- list(
  manual.opts = manual.opts,
  general.opts = list("reps"=200, # repetitions of simulation
                      "train_obs"=c(100,500,1000,5000,10000), # number of training observations
                      "test_obs"=500, # number of test observations
                      "alpha"=0.05), # miscoverage level
  tune.opts = list("tune_reps"=10, # number of additional training sets used for tuning (see Chapter 12.2)
                   "excess_error_threshold"=0.001, # fraction of OOB error that the Monte Carlo error is allowed to take
                   "num.trees.start"=formals(grf::regression_forest)$num.trees, # starting number of trees when selecting
                   # optimal number of trees
                   "num.trees.increment"=200), # tree increment
  interval.opts = list("intervals"=c("grf", "scp", "sdcp", "wscp", "boot", "oob", "woob"), # interval types
                       "cp_subsampling_rate"=0.5, # split conformal true training size share of total training data
                       "bootstrap_reps"=499, # number of bootstrap replications 
                       "dcp_taus"=seq(0,1,length=202)[2:201], # estimated quantiles for SDCP
                       "boot_skip"=0.5), # fraction of how many bootstrap intervals should be skipped for computational
  # reasons
  parallel.opts = list("packages"=c("foreach", "copula", "iterators")) # packages exported to each worker
)



########################################################################################################################
########################################################################################################################


# Create folders for output
if(!dir.exists(file.path(.global.opts$manual.opts$path, "01_simulation_results")))
  dir.create(file.path(.global.opts$manual.opts$path, "01_simulation_results"))

if(!dir.exists(file.path(.global.opts$manual.opts$path, "02_figures")))
  dir.create(file.path(.global.opts$manual.opts$path, "02_figures"))


# Register parallel backend
if(.global.opts$manual.opts$parallel){
  cl <- parallel::makeCluster(.global.opts$manual.opts$ncores)
  doSNOW::registerDoSNOW(cl)
  on.exit(parallel::stopCluster(cl))
} else {
  if(exists("cl")) parallel::stopCluster(cl); rm("cl")
  foreach::registerDoSEQ()
}


# Load functions
source(file.path(.global.opts$manual.opts$path, "functions.R"))


########################################################################################################################
# Actual simulation


# Simulation Scenarios (See Chapter 8.1) 
.simulation.settings <- list("mean_function"=c("linear", "nonlinear", "nonlinear_jump"),
                             "error_distribution"=c("homoscedastic", "heavy_tailed", "heteroscedastic"),
                             "correlated"=c(F, T)) 


# Empty list for results
results <- vector("list", length(.simulation.settings$mean_function)*
                    length(.simulation.settings$error_distribution)*
                    length(.simulation.settings$correlated))


# Simulation (this takes forever to run)
for(i in 1:length(.simulation.settings$mean_function)){
  if(.simulation.settings$mean_function[i]=="linear"){
    mean_fun <- mean_fun_generator(.body="X[,1] + X[,2] + X[,3]")
  } else if (.simulation.settings$mean_function[i]=="nonlinear"){
    mean_fun <- mean_fun_generator(.body="X[,1] * exp(sin(sqrt(5/4)*X[,2]) + abs(X[,3])/2)")
  } else if (.simulation.settings$mean_function[i]=="nonlinear_jump"){
    mean_fun <- mean_fun_generator(.body="0.9*sign(tan(X[,1]))*(floor(3*sin(X[,2])) + 0.1/cosh(X[,3]))")
  } else stop("undefined mean function")
  
  for(j in 1:length(.simulation.settings$error_distribution)){
    if(.simulation.settings$error_distribution[j]=="homoscedastic"){
      err_fun <- err_fun_generator(.args=list("sd"="1", "mean"=0), .dist="norm")
    } else if (.simulation.settings$error_distribution[j]=="heavy_tailed"){
      err_fun <- err_fun_generator(.args=list("df"=4), .dist="t", adjustment=(1/sqrt(2)))
    } else if (.simulation.settings$error_distribution[j]=="heteroscedastic"){
      err_fun <- err_fun_generator(.args=list("sd"="0.9+sin(mean_fun(X[i,]))/sqrt(exp(1))", "mean"=0), .dist="norm")
    } else stop("undefined error distribution")
    
    
    for(k in 1:length(.simulation.settings$correlated)){
      if(.simulation.settings$correlated[k]){
        dist <- dist_generator(10, "ar1", 0.6, "norm", list(mean=0, sd=1))
      } else {
        dist <- dist_generator(10, "uncor", 0, "norm", list(mean=0, sd=1))
      }
      
      if(i==1 & j==1 & k==1) {
        print(paste(length(.simulation.settings$mean_function)*
                      length(.simulation.settings$error_distribution)*
                      length(.simulation.settings$correlated), "simulations started. Note that this will take a while with the default settings. Starting simulation 1:"))
        print("Note that the bootstrap interval makes up roughly 99.9% of the total training time, and thus the loading bar will not necessarily progress equally fast at all times.")
        print("With the default parameters this is very noticable, as the first 50% of the simulations are done way faster as the second 50% as only then the bootstrap interval is calculated.")
      }
      
      DGP <- DGP_generator(dist)
      
      snow::clusterExport(cl, list("mean_fun", "err_fun", "DGP", "tune_trees", "tune_forest", "mc.error", "w.quantile",
                                   "GRF_interval", "SCP_interval", "SDCP_interval", "WSCP_interval",
                                   "bootstrap_interval", "OOB_interval", "WOOB_interval"))
      
      tuning_results <- simulation_tuning(DataGeneratingProcess = DGP, 
                                          opts = .global.opts)
      
      print("Tuning of this scenario done.")
      
      sim_results <- simulation(DataGeneratingProcess = DGP, 
                                opts=.global.opts, 
                                tune_param = tuning_results
      )
      
      results[[(i-1)*length(.simulation.settings$error_distribution)*
                 length(.simulation.settings$correlated)+(j-1)*
                 length(.simulation.settings$correlated)+k]] <-
        cbind(sim_results, 
              "mean_function"=.simulation.settings$mean_function[i],
              "error_distribution"=.simulation.settings$error_distribution[j],
              "correlated"=.simulation.settings$correlated[k])
      
      print(paste((i-1)*length(.simulation.settings$error_distribution)*
                    length(.simulation.settings$correlated)+(j-1)*
                    length(.simulation.settings$correlated)+k, 
                  "out of of", 
                  length(.simulation.settings$mean_function)*
                    length(.simulation.settings$error_distribution)*
                    length(.simulation.settings$correlated),
                  "done -",
                  Sys.time()))
      
    }
    save.image(file=file.path(.global.opts$manual.opts$path, 
                              "01_simulation_results",
                              paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),".RData")))
  }
}





########################################################################################################################
# FIGURES ##############################################################################################################
########################################################################################################################

# Remove parallel backend
if(.global.opts$manual.opts$parallel){
  parallel::stopCluster(cl); rm("cl"); foreach::registerDoSEQ()
} 


# Formatting
results <- do.call("rbind", results)

results$type <- factor(results$type,
                       levels=c("grf", "cp", "dcp", "cp_weighted", "boot", "oob", "oob_weighted", "true"),
                       labels=c("GRF", "SCP", "SDCP", "WSCP", "Bootstrap", "OOB", "WOOB", "True"))

results$mean_function <- factor(results$mean_function, 
                                levels=.simulation.settings$mean_function,
                                labels=c("Linear", "Non-linear", "Non-linear with jump points")) 

results$error_distribution <- factor(results$error_distribution, 
                                     levels=.simulation.settings$error_distribution,
                                     labels=c("Homoscedastic", "Heavy-tailed", "Heteroscedastic"))

levels(results$nobs) <- paste0(levels(results$nobs), "   ")

results <- results[!is.na(results$hit),]

ar_data <- results[which(results$correlated),]
uncorr_data <- results[which(!results$correlated),]


########################################################################################################################
# Plots (figure numbers are only referring to figures within this chapter)


# Data for each plot
figure1data <- uncorr_data[which(uncorr_data$mean_function=="Linear"),]
figure2data <- uncorr_data[which(uncorr_data$mean_function=="Non-linear"),]
figure3data <- uncorr_data[which(uncorr_data$mean_function=="Non-linear with jump points"),]


# plots
f1a <- ggplot(figure1data, aes(x=type, y=hit, fill=nobs))+
  facet_grid(rows=vars(error_distribution), cols=vars(mean_function))+
  plot.aesthetics("coverage", uncorr_data, .global.opts)+
  theme(plot.margin = unit(c(3,20,3,3), "points"))


f1b <- ggplot(figure1data[which(figure1data$type!="True"),], aes(x=type, y=ln_rel_length, fill=nobs))+
  facet_grid(rows=vars(error_distribution), cols=vars(mean_function))+
  plot.aesthetics("length", uncorr_data, .global.opts, max.y.override = 1.05)+
  theme(plot.margin = unit(c(3,3,3,20), "points"))


f2a <- ggplot(figure2data, aes(x=type, y=hit, fill=nobs))+
  facet_grid(rows=vars(error_distribution), cols=vars(mean_function))+
  plot.aesthetics("coverage", uncorr_data, .global.opts)+
  theme(plot.margin = unit(c(3,20,3,3), "points"))


f2b <- ggplot(figure2data[which(figure2data$type!="True"),], aes(x=type, y=ln_rel_length, fill=nobs))+
  facet_grid(rows=vars(error_distribution), cols=vars(mean_function))+
  plot.aesthetics("length", uncorr_data, .global.opts)+
  theme(plot.margin = unit(c(3,3,3,20), "points"))


f3a <- ggplot(figure3data, aes(x=type, y=hit, fill=nobs))+
  facet_grid(rows=vars(error_distribution), cols=vars(mean_function))+
  plot.aesthetics("coverage", uncorr_data, .global.opts)+
  theme(plot.margin = unit(c(3,20,3,3), "points"))


f3b <- ggplot(figure3data[which(figure3data$type!="True"),], aes(x=type, y=ln_rel_length, fill=nobs))+
  facet_grid(rows=vars(error_distribution), cols=vars(mean_function))+
  plot.aesthetics("length", uncorr_data, .global.opts, max.y.override = 1.05)+
  theme(plot.margin = unit(c(3,3,3,20), "points"))

# save plots
pdf.options(reset = TRUE, onefile = FALSE)

pdf(file = file.path(.global.opts$manual.opts$path, "02_figures/f1.pdf"), width = 16, height = 14) 
ggpubr::ggarrange(f1a, f1b, ncol=2, common.legend = TRUE, legend="bottom", widths=c(8,7))
dev.off()

pdf(file = file.path(.global.opts$manual.opts$path, "02_figures/f2.pdf"), width = 16, height = 14) 
ggpubr::ggarrange(f2a, f2b, ncol=2, common.legend = TRUE, legend="bottom", widths=c(8,7))
dev.off()

pdf(file = file.path(.global.opts$manual.opts$path, "02_figures/f3.pdf"), width = 16, height = 14) 
ggpubr::ggarrange(f3a, f3b, ncol=2, common.legend = TRUE, legend="bottom", widths=c(8,7))
dev.off()

########################################################################################################################
# AR plots (in appendix) 

# Data for each plot
figure4data <- ar_data[which(ar_data$mean_function=="Linear"),]
figure5data <- ar_data[which(ar_data$mean_function=="Non-linear"),]
figure6data <- ar_data[which(ar_data$mean_function=="Non-linear with jump points"),]


# plots
f4a <- ggplot(figure4data, aes(x=type, y=hit, fill=nobs))+
  facet_grid(rows=vars(error_distribution), cols=vars(mean_function))+
  plot.aesthetics("coverage", ar_data, .global.opts)+
  theme(plot.margin = unit(c(3,20,3,3), "points"))


f4b <- ggplot(figure4data[which(figure4data$type!="True"),], aes(x=type, y=ln_rel_length, fill=nobs))+
  facet_grid(rows=vars(error_distribution), cols=vars(mean_function))+
  plot.aesthetics("length", ar_data, .global.opts, max.y.override = 1.05)+
  theme(plot.margin = unit(c(3,3,3,20), "points"))


f5a <- ggplot(figure5data, aes(x=type, y=hit, fill=nobs))+
  facet_grid(rows=vars(error_distribution), cols=vars(mean_function))+
  plot.aesthetics("coverage", ar_data, .global.opts)+
  theme(plot.margin = unit(c(3,20,3,3), "points"))


f5b <- ggplot(figure5data[which(figure5data$type!="True"),], aes(x=type, y=ln_rel_length, fill=nobs))+
  facet_grid(rows=vars(error_distribution), cols=vars(mean_function))+
  plot.aesthetics("length", ar_data, .global.opts)+
  theme(plot.margin = unit(c(3,3,3,20), "points"))


f6a <- ggplot(figure6data, aes(x=type, y=hit, fill=nobs))+
  facet_grid(rows=vars(error_distribution), cols=vars(mean_function))+
  plot.aesthetics("coverage", ar_data, .global.opts)+
  theme(plot.margin = unit(c(3,20,3,3), "points"))


f6b <- ggplot(figure6data[which(figure6data$type!="True"),], aes(x=type, y=ln_rel_length, fill=nobs))+
  facet_grid(rows=vars(error_distribution), cols=vars(mean_function))+
  plot.aesthetics("length", ar_data, .global.opts, max.y.override = 1.05)+
  theme(plot.margin = unit(c(3,3,3,20), "points"))


# save plots
pdf.options(reset = TRUE, onefile = FALSE)

pdf(file = file.path(.global.opts$manual.opts$path, "02_figures/f4.pdf"), width = 16, height = 14) 
ggpubr::ggarrange(f4a, f4b, ncol=2, common.legend = TRUE, legend="bottom", widths=c(8,7))
dev.off()

pdf(file = file.path(.global.opts$manual.opts$path, "02_figures/f5.pdf"), width = 16, height = 14) 
ggpubr::ggarrange(f5a, f5b, ncol=2, common.legend = TRUE, legend="bottom", widths=c(8,7))
dev.off()

pdf(file = file.path(.global.opts$manual.opts$path, "02_figures/f6.pdf"), width = 16, height = 14) 
ggpubr::ggarrange(f6a, f6b, ncol=2, common.legend = TRUE, legend="bottom", widths=c(8,7))
dev.off()



