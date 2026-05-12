library(Hmsc)
load("models/models_thin_1000_samples_250_chains_4.Rdata")

pdf("EOP_predictions.pdf")
for(mtype in names(models)){
  print(mtype)
  m = models[[mtype]]
  covariates = all.vars(m$XFormula)
  for(var in covariates){
    print(var)
    Gradient = constructGradient(hM = m,focalVariable = var)
    predY = predict(m, Gradient=Gradient, expected = TRUE)
    n = length(predY)
    for(i in 1:n){
      p = exp(predY[[i]])
      rs = rowSums(p)
      for(j in 1:3) p[,j] = p[,j]/rs
      if(i==1){
        mp = p
      } else {
        mp = mp + p
      }
    }
    mp = mp/n
    plot(Gradient$XDataNew[[var]],mp[,1],type="l",ylim=c(0,1),
         xlab = var,ylab = c("E-O-P"), main = mtype)
    lines(Gradient$XDataNew[[var]],mp[,1]+mp[,2])
  }
}
dev.off()


###################################################################
load("models/models_thin_1000_samples_250_chains_4.Rdata")

pdf("EOP_predictions_colored.pdf")

# Plot
par(
    cex.axis = 1.5,         
    cex.lab  = 1.8,        
    cex.main = 2.0)   

cols_base <- c("#ED620F", "#1C5A99", "#146C36")  # E / O / P

# transparency
cols <- adjustcolor(cols_base, alpha.f = 0.8)

for(mtype in names(models)){
  print(mtype)
  ## change the mtype to the title text
  mtype_label <- if (mtype == "native") {
    "Native"
  } else if (mtype == "non_native") {
    "Non-native"
  } else {
    mtype 
  }
  m = models[[mtype]]
  covariates = all.vars(m$XFormula)
  for(var in covariates){
    print(var)
    Gradient = constructGradient(hM = m, focalVariable = var)
    predY = predict(m, Gradient = Gradient, expected = TRUE)
    n = length(predY)
    
    for(i in 1:n){
      p = exp(predY[[i]])
      rs = rowSums(p)
      for(j in 1:3) p[, j] = p[, j] / rs  # rescale to 0–1
      if(i == 1){
        mp = p
      } else {
        mp = mp + p
      }
    }
    mp = mp / n   #mean proportion（E, O, P）
    
    # x variate values
    x <- Gradient$XDataNew[[var]]
    
    # polygon
    plot(x, mp[,1], type = "n", ylim = c(0, 1),
         xlab = var, ylab = "Mean proportion", main = mtype)
    
    # （E）
    polygon(
      x = c(x, rev(x)),
      y = c(rep(0, length(x)), rev(mp[,1])),
      col = cols[1], border = NA
    )
    
    # （O）
    polygon(
      x = c(x, rev(x)),
      y = c(mp[,1], rev(mp[,1] + mp[,2])),
      col = cols[2], border = NA
    )
    
    # （P）
    polygon(
      x = c(x, rev(x)),
      y = c(mp[,1] + mp[,2], rev(rep(1, length(x)))),
      col = cols[3], border = NA
    )
    
    ## Legend: Horizontal, placed in the upper right corner outside the panel
    legend("topright",
           inset = c(-0.25, -0.02),   
           legend = c("Equilibrium", "Opportunistic", "Periodic"),
           fill   = cols,
           bty    = "n",
           horiz  = TRUE,             
           cex    = 1.6,            
           xpd    = TRUE)
  }
}

dev.off()
