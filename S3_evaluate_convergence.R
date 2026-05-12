# Set the base directory using your favorite method
setwd("XX")

##################################################################################################
# INPUT AND OUTPUT OF THIS SCRIPT (BEGINNING)
##################################################################################################
#	INPUT. Fitted models

#	OUTPUT. MCMC convergence statistics for selected model parameters,
# illustrated (for all RUNs performed thus far in S3) in the file "results/MCMC_convergence.pdf",
# and the text file "results/MCMC_convergence.txt".
##################################################################################################
# INPUT AND OUTPUT OF THIS SCRIPT (END)
##################################################################################################


##################################################################################################
# MAKE THE SCRIPT REPRODUCIBLE (BEGINNING)
##################################################################################################
set.seed(1)
##################################################################################################
## MAKE THE SCRIPT REPRODUCIBLE (END)
##################################################################################################


##################################################################################################
# SETTING COMMONLY ADJUSTED PARAMETERS TO NULL WHICH CORRESPONDS TO DEFAULT CHOICE (BEGINNING)
##################################################################################################
showBeta = NULL #Default: showBeta = TRUE, convergence shown for beta-parameters
showGamma = NULL #Default: showGamma = FALSE, convergence not shown for gamma-parameters
showOmega = NULL #Default: showOmega = FALSE, convergence not shown for Omega-parameters
maxOmega = NULL #Default: convergence of Omega shown for 50 randomly selected species pairs
showRho = NULL #Default: showRho = FALSE, convergence not shown for rho-parameters
showAlpha = NULL #Default: showAlpha = FALSE, convergence not shown for alpha-parameters
##################################################################################################
# SETTING COMMONLY ADJUSTED PARAMETERS TO NULL WHICH CORRESPONDS TO DEFAULT CHOICE (END)
##################################################################################################

##################################################################################################
# CHANGE DEFAULT OPTIONS BY REMOVING COMMENT AND SETTING VALUE (BEGINNING)
# NOTE THAT THIS IS THE ONLY SECTION OF THE SCRIPT THAT YOU TYPICALLY NEED TO MODIFY
##################################################################################################
showBeta = TRUE
showGamma = TRUE
showOmega = TRUE
maxOmega = 100
showRho = TRUE
showAlpha = TRUE
##################################################################################################
# CHANGE DEFAULT OPTIONS BY REMOVING COMMENT AND SETTING VALUE (END)
# NOTE THAT THIS IS THE ONLY SECTION OF THE SCRIPT THAT YOU TYPICALLY NEED TO MODIFY
##################################################################################################

##################################################################################################
# SET DIRECTORIES (BEGINNING)
##################################################################################################
localDir = "."
modelDir = file.path(localDir, "models")
resultDir = file.path(localDir, "results")
if (!dir.exists(resultDir)) dir.create(resultDir)
##################################################################################################
# SET DIRECTORIES (END)
##################################################################################################

if(is.null(showBeta)) showBeta = TRUE
if(is.null(showGamma)) showGamma = FALSE
if(is.null(showOmega)) showOmega = FALSE
if(is.null(maxOmega)) maxOmega = 50
if(is.null(showRho)) showRho = FALSE
if(is.null(showAlpha)) showAlpha = FALSE

library(Hmsc)
library(colorspace)
library(vioplot)

samples_list = c(5,250,250,250,250) 
thin_list = c(1,1,10,100,1000)  
nst = length(thin_list)
nChains = 4

text.file = file.path(resultDir,"/MCMC_convergence.txt")
cat("MCMC Convergennce statistics\n\n",file=text.file,sep="")

ma.beta = NULL
na.beta = NULL
ma.gamma = NULL
na.gamma = NULL
ma.omega= NULL
na.omega = NULL
ma.alpha = NULL
na.alpha = NULL  
ma.rho = NULL
na.rho = NULL
Lst = 1
while(Lst <= nst){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  
  
  filename = file.path(modelDir,paste("models_thin_", as.character(thin),
                                      "_samples_", as.character(samples),
                                      "_chains_",as.character(nChains),
                                      ".Rdata",sep = ""))
  if(file.exists(filename)){
    load(filename)
    cat(c("\n",filename,"\n\n"),file=text.file,sep="",append=TRUE)
    nm = length(models)
    for(j in 1:nm){
      mpost = convertToCodaObject(models[[j]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
      nr = models[[j]]$nr
      cat(c("\n",names(models)[j],"\n\n"),file=text.file,sep="",append=TRUE)
      if(showBeta){
        psrf = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
        tmp = summary(psrf)
        cat("\nbeta\n\n",file=text.file,sep="",append=TRUE)
        cat(tmp[,1],file=text.file,sep="\n",append=TRUE)
        if(is.null(ma.beta)){
          ma.beta = psrf[,1]
          na.beta = paste0(as.character(thin),",",as.character(samples))
        } else {
          ma.beta = cbind(ma.beta,psrf[,1])
          if(j==1){
            na.beta = c(na.beta,paste0(as.character(thin),",",as.character(samples)))
          } else {
            na.beta = c(na.beta,"")
          }
        }
      }
      if(showGamma){
        psrf = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf
        tmp = summary(psrf)
        cat("\ngamma\n\n",file=text.file,sep="",append=TRUE)
        cat(tmp[,1],file=text.file,sep="\n",append=TRUE)
        if(is.null(ma.gamma)){
          ma.gamma = psrf[,1]
          na.gamma = paste0(as.character(thin),",",as.character(samples))
        } else {
          ma.gamma = cbind(ma.gamma,psrf[,1])
          if(j==1){
            na.gamma = c(na.gamma,paste0(as.character(thin),",",as.character(samples)))
          } else {
            na.gamma = c(na.gamma,"")
          }
        }
      }
      if(showRho & !is.null(mpost$Rho)){
        psrf = gelman.diag(mpost$Rho,multivariate=FALSE)$psrf
        cat("\nrho\n\n",file=text.file,sep="",append=TRUE)
        cat(psrf[1],file=text.file,sep="\n",append=TRUE)
      }
      if(showOmega & nr>0){
        cat("\nomega\n\n",file=text.file,sep="",append=TRUE)
        for(k in 1:nr){
          cat(c("\n",names(models[[j]]$ranLevels)[k],"\n\n"),file=text.file,sep="",append=TRUE)
          tmp = mpost$Omega[[k]]
          z = dim(tmp[[1]])[2]
          if(z > maxOmega){
            sel = sample(1:z, size = maxOmega)
            for(i in 1:length(tmp)){
              tmp[[i]] = tmp[[i]][,sel]
            }
          }
          psrf = gelman.diag(tmp, multivariate = FALSE)$psrf
          tmp = summary(psrf)
          cat(tmp[,1],file=text.file,sep="\n",append=TRUE)
          if(is.null(ma.omega)){
            ma.omega = psrf[,1]
            na.omega = paste0(as.character(thin),",",as.character(samples))
          } else {
            ma.omega = cbind(ma.omega,psrf[,1])
            if(j==1){
              na.omega = c(na.omega,paste0(as.character(thin),",",as.character(samples)))
            } else {
              na.omega = c(na.omega,"")
            }
          }
        }
      }
      if(showAlpha & nr>0){
        for(k in 1:nr){
          if(models[[j]]$ranLevels[[k]]$sDim>0){
            cat("\nalpha\n\n",file=text.file,sep="\n",append=TRUE)
            cat(c("\n",names(models[[j]]$ranLevels)[k],"\n\n"),file=text.file,sep="",append=TRUE)
            psrf = gelman.diag(mpost$Alpha[[k]],multivariate = FALSE)$psrf
            cat(psrf[,1],file=text.file,sep="\n",append=TRUE)            
          }
        }
      }
    }
  }
  Lst = Lst + 1
}

pdf(file= file.path(resultDir,"/MCMC_convergence.pdf"))
if(showBeta){
  par(mfrow=c(2,1))
  vioplot(ma.beta,col=rainbow_hcl(nm),names=na.beta,ylim=c(0,max(ma.beta)),main="psrf(beta)")
  legend("topright",legend = names(models), fill=rainbow_hcl(nm))
  vioplot(ma.beta,col=rainbow_hcl(nm),names=na.beta,ylim=c(0.9,1.1),main="psrf(beta)")
}
if(showGamma){
  par(mfrow=c(2,1))
  vioplot(ma.gamma,col=rainbow_hcl(nm),names=na.gamma,ylim=c(0,max(ma.gamma)),main="psrf(gamma)")
  legend("topright",legend = names(models), fill=rainbow_hcl(nm))
  vioplot(ma.gamma,col=rainbow_hcl(nm),names=na.gamma,ylim=c(0.9,1.1),main="psrf(gamma)")
}
if(showOmega & !is.null(ma.omega)){
  par(mfrow=c(2,1))
  vioplot(ma.omega,col=rainbow_hcl(nm),names=na.omega,ylim=c(0,max(ma.omega)),main="psrf(omega)")
  legend("topright",legend = names(models), fill=rainbow_hcl(nm))
  vioplot(ma.omega,col=rainbow_hcl(nm),names=na.omega,ylim=c(0.9,1.1),main="psrf(omega)")
}
dev.off()


##########################################################################################################################
###############################################################################################Drop bad chain for non-naitve.NE model
set.seed(1)

showBeta = TRUE
showGamma = TRUE
showOmega = TRUE
maxOmega = 100
showRho = TRUE
showAlpha = TRUE

localDir = "."
modelDir = file.path(localDir, "models")
resultDir = file.path(localDir, "results")
if (!dir.exists(resultDir)) dir.create(resultDir)

library(Hmsc)
library(colorspace)
library(vioplot)
library(coda)

##################################################################################################
# Load fitted model
##################################################################################################

filename = file.path(modelDir,"models_thin_1000_samples_250_chains_4.Rdata")

load(filename)

##################################################################################################
# Select model to test
##################################################################################################

model_full = models$non_native.NE

##################################################################################################
# Create model list: full model + drop each chain
##################################################################################################

drop_chain <- function(model, chain_to_drop) {
  model_new <- model
  model_new$postList <- model$postList[-chain_to_drop]
  return(model_new)
}

model_list = list()
model_list[["full_4_chains"]] = model_full
model_list[["drop_chain_1"]] = drop_chain(model_full, 1)
model_list[["drop_chain_2"]] = drop_chain(model_full, 2)
model_list[["drop_chain_3"]] = drop_chain(model_full, 3)
model_list[["drop_chain_4"]] = drop_chain(model_full, 4)

nm = length(model_list)

##################################################################################################
# Output text file
##################################################################################################

text.file = file.path(resultDir, "MCMC_convergence_drop_chains_non_native_NE.txt")
cat("MCMC Convergence statistics after dropping each chain\n\n", file = text.file, sep = "")

##################################################################################################
# Containers
##################################################################################################

ma.beta = NULL
na.beta = NULL

ma.gamma = NULL
na.gamma = NULL

ma.omega = NULL
na.omega = NULL

ma.alpha = NULL
na.alpha = NULL

ma.rho = NULL
na.rho = NULL

##################################################################################################
# Calculate PSRF for full model and models with one chain dropped
##################################################################################################

for (j in 1:nm) {
  
  current_model_name = names(model_list)[j]
  current_model = model_list[[j]]
  
  cat(c("\n", current_model_name, "\n\n"), file = text.file, sep = "", append = TRUE)
  
  mpost = convertToCodaObject(
    current_model,
    spNamesNumbers = c(TRUE, FALSE),
    covNamesNumbers = c(TRUE, FALSE)
  )
  
  nr = current_model$nr
  
  ################################################################################################
  # Beta
  ################################################################################################
  
  if (showBeta) {
    psrf = gelman.diag(mpost$Beta, multivariate = FALSE)$psrf
    tmp = summary(psrf)
    
    cat("\nbeta\n\n", file = text.file, sep = "", append = TRUE)
    cat(tmp[, 1], file = text.file, sep = "\n", append = TRUE)
    cat("\n\n", file = text.file, append = TRUE)
    
    if (is.null(ma.beta)) {
      ma.beta = matrix(psrf[, 1], ncol = 1)
      na.beta = current_model_name
    } else {
      ma.beta = cbind(ma.beta, psrf[, 1])
      na.beta = c(na.beta, current_model_name)
    }
  }
  
  ################################################################################################
  # Gamma
  ################################################################################################
  
  if (showGamma) {
    psrf = gelman.diag(mpost$Gamma, multivariate = FALSE)$psrf
    tmp = summary(psrf)
    
    cat("\ngamma\n\n", file = text.file, sep = "", append = TRUE)
    cat(tmp[, 1], file = text.file, sep = "\n", append = TRUE)
    cat("\n\n", file = text.file, append = TRUE)
    
    if (is.null(ma.gamma)) {
      ma.gamma = matrix(psrf[, 1], ncol = 1)
      na.gamma = current_model_name
    } else {
      ma.gamma = cbind(ma.gamma, psrf[, 1])
      na.gamma = c(na.gamma, current_model_name)
    }
  }
  
  ################################################################################################
  # Rho
  ################################################################################################
  
  if (showRho & !is.null(mpost$Rho)) {
    psrf = gelman.diag(mpost$Rho, multivariate = FALSE)$psrf
    
    cat("\nrho\n\n", file = text.file, sep = "", append = TRUE)
    cat(psrf[, 1], file = text.file, sep = "\n", append = TRUE)
    cat("\n\n", file = text.file, append = TRUE)
    
    if (is.null(ma.rho)) {
      ma.rho = matrix(psrf[, 1], ncol = 1)
      na.rho = current_model_name
    } else {
      ma.rho = cbind(ma.rho, psrf[, 1])
      na.rho = c(na.rho, current_model_name)
    }
  }
  
  ################################################################################################
  # Omega
  ################################################################################################
  
  if (showOmega & nr > 0) {
    cat("\nomega\n\n", file = text.file, sep = "", append = TRUE)
    
    for (k in 1:nr) {
      
      ran_name = names(current_model$ranLevels)[k]
      cat(c("\n", ran_name, "\n\n"), file = text.file, sep = "", append = TRUE)
      
      tmp = mpost$Omega[[k]]
      z = dim(tmp[[1]])[2]
      
      if (z > maxOmega) {
        sel = sample(1:z, size = maxOmega)
        for (i in 1:length(tmp)) {
          tmp[[i]] = tmp[[i]][, sel]
        }
      }
      
      psrf = gelman.diag(tmp, multivariate = FALSE)$psrf
      tmp_summary = summary(psrf)
      
      cat(tmp_summary[, 1], file = text.file, sep = "\n", append = TRUE)
      cat("\n\n", file = text.file, append = TRUE)
      
      omega_name = paste0(current_model_name, "_", ran_name)
      
      if (is.null(ma.omega)) {
        ma.omega = matrix(psrf[, 1], ncol = 1)
        na.omega = omega_name
      } else {
        ma.omega = cbind(ma.omega, psrf[, 1])
        na.omega = c(na.omega, omega_name)
      }
    }
  }
  
  ################################################################################################
  # Alpha
  ################################################################################################
  
  if (showAlpha & nr > 0) {
    for (k in 1:nr) {
      if (current_model$ranLevels[[k]]$sDim > 0) {
        
        ran_name = names(current_model$ranLevels)[k]
        
        cat("\nalpha\n\n", file = text.file, sep = "\n", append = TRUE)
        cat(c("\n", ran_name, "\n\n"), file = text.file, sep = "", append = TRUE)
        
        psrf = gelman.diag(mpost$Alpha[[k]], multivariate = FALSE)$psrf
        cat(psrf[, 1], file = text.file, sep = "\n", append = TRUE)
        cat("\n\n", file = text.file, append = TRUE)
        
        alpha_name = paste0(current_model_name, "_", ran_name)
        
        if (is.null(ma.alpha)) {
          ma.alpha = matrix(psrf[, 1], ncol = 1)
          na.alpha = alpha_name
        } else {
          ma.alpha = cbind(ma.alpha, psrf[, 1])
          na.alpha = c(na.alpha, alpha_name)
        }
      }
    }
  }
}

##################################################################################################
# Save summary table for Beta
##################################################################################################

beta_summary = data.frame(
  model = na.beta,
  max_PSRF = apply(ma.beta, 2, max, na.rm = TRUE),
  mean_PSRF = apply(ma.beta, 2, mean, na.rm = TRUE),
  median_PSRF = apply(ma.beta, 2, median, na.rm = TRUE),
  n_over_1.10 = apply(ma.beta, 2, function(x) sum(x > 1.10, na.rm = TRUE)),
  n_over_1.05 = apply(ma.beta, 2, function(x) sum(x > 1.05, na.rm = TRUE))
)

write.csv(
  beta_summary,
  file.path(resultDir, "PSRF_beta_drop_chains_non_native_NE.csv"),
  row.names = FALSE
)

print(beta_summary)

##################################################################################################
# Plot PDF
##################################################################################################

safe_vioplot <- function(mat, names_vec, main_title) {
  
  # Check whether there are valid finite values
  if (is.null(mat)) return(NULL)
  if (all(is.na(mat))) return(NULL)
  if (!any(is.finite(mat))) return(NULL)
  
  par(mfrow = c(2, 1), mar = c(7, 4, 4, 2))
  
  cols = rainbow_hcl(ncol(mat))
  
  ymax = max(mat, na.rm = TRUE)
  
  vioplot(
    mat,
    col = cols,
    names = names_vec,
    las = 2,
    ylim = c(0, ymax),
    main = main_title
  )
  abline(h = 1.1, lty = 2)
  
  vioplot(
    mat,
    col = cols,
    names = names_vec,
    las = 2,
    ylim = c(0.9, 1.1),
    main = paste0(main_title, " zoomed")
  )
  abline(h = 1.1, lty = 2)
}

##################################################################################################
# Plot PDF
##################################################################################################

pdf(
  file = file.path(resultDir, "MCMC_convergence_drop_chains_non_native_NE.pdf"),
  width = 10,
  height = 8
)

if (showBeta) {
  safe_vioplot(ma.beta, na.beta, "psrf(beta)")
}

if (showGamma) {
  safe_vioplot(ma.gamma, na.gamma, "psrf(gamma)")
}

if (showOmega) {
  safe_vioplot(ma.omega, na.omega, "psrf(omega)")
}

if (showRho) {
  safe_vioplot(ma.rho, na.rho, "psrf(rho)")
}

if (showAlpha) {
  safe_vioplot(ma.alpha, na.alpha, "psrf(alpha)")
}

dev.off()

##################################################################################################
# Backup original file
##################################################################################################

file.copy(
  file.path(modelDir, "models_thin_1000_samples_250_chains_4.Rdata"),
  file.path(modelDir, "models_thin_1000_samples_250_chains_4_original.Rdata"),
  overwrite = TRUE
)

##################################################################################################
# Remove chain 3 from non_native.NE
##################################################################################################

models$non_native.NE$postList <- models$non_native.NE$postList[-3]

##################################################################################################
# Check remaining number of chains
##################################################################################################

length(models$non_native.NE$postList)

##################################################################################################
# Save modified model
##################################################################################################

save(models,file = file.path(modelDir, "models_thin_1000_samples_250_chains_4.Rdata"))
