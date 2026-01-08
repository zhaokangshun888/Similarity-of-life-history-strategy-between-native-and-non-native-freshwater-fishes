
################################################################################# Loading data and preprocessing
rm(list=ls())  
setwd("XX")
# Load data 读取数据
data <- read.csv("Basin_Strategy_with_predictors_model.csv")

summary(data)
head(data)

# Remove missing values in explanatory variables
data <- na.omit(data)

# Split Native and Non-native data 
data_native <- data[data$Native.Exotic.Status == "Native", ]
data_non_native <- data[data$Native.Exotic.Status == "Non-native", ]

# Standardize fixed effect variables (z-score)
native_predictors <- c("MAT", "TS", "DIS", "RAD", "GDP", "CSI", "NSR", "RBA")
non_native_predictors <- c("MAT", "TS", "DIS", "RAD", "GDP", "CSI", "NSR", "RBA")

data_native[native_predictors] <- scale(data_native[native_predictors])
data_non_native[non_native_predictors] <- scale(data_non_native[non_native_predictors])



################################################################################# CLR MAKE THE p->x TRANSFORMATION
library(zCompositions)
library(compositions)

# Select the columns for CLR transformation
clr_vars <- c("Equilibrium", "Periodic", "Opportunistic")

# ---- Native ----
# Zero replacement (CZM method) + CLR transformation
native_ready <- cmultRepl(data_native[clr_vars], method="CZM")
native_clr <- clr(acomp(native_ready))
colnames(native_clr) <- paste0(clr_vars, "_CLR")
data_native_clr <- cbind(data_native, native_clr)

# ---- Non-native ----
non_native_ready <- cmultRepl(data_non_native[clr_vars], method="CZM")
non_native_clr <- clr(acomp(non_native_ready))
colnames(non_native_clr) <- paste0(clr_vars, "_CLR")
data_non_native_clr <- cbind(data_non_native, non_native_clr)



################################################################################# Data distribution after CLR
library(ggplot2)
library(gridExtra)

clr_vars <- c("Equilibrium_CLR", "Periodic_CLR", "Opportunistic_CLR")

# ---- plot ----
plot_clr_distribution <- function(df, vars, group_name) {
  plots <- lapply(vars, function(v) {
    ggplot(df, aes_string(x = v)) +
      geom_histogram(aes(y=..density..), bins=30, fill="skyblue", color="black", alpha=0.6) +
      geom_density(color="red", size=1) +
      labs(title = paste(group_name, "-", gsub("_CLR","",v)),
           x = "CLR", y = "Density") +
      theme_minimal()
  })
  return(plots)
}

# Native
plots_native_dist <- plot_clr_distribution(data_native_clr, clr_vars, "Native")
# Non-native
plots_non_native_dist <- plot_clr_distribution(data_non_native_clr, clr_vars, "Non-native")

grid.arrange(grobs = plots_native_dist, ncol=3)
grid.arrange(grobs = plots_non_native_dist, ncol=3)



################################################################################################### model-Native
library(Hmsc)
# Plot point distribution
plot(data_native_clr$Median.Longitude, data_native_clr$Median.Latitude,
     xlab="Longitude", ylab="Latitude", main="Native Basins")

# Study design
studyDesign <- data.frame(
  basin = as.factor(data_native_clr$Basin.Name),
  ecoregion = as.factor(data_native_clr$Ecoregion))

# Random effects
rL.ecoregion <- HmscRandomLevel(units = levels(studyDesign$ecoregion))
rL.ecoregion <- setPriors(rL.ecoregion,nfMin = 1, nfMax=1)

lonlat <- cbind(data_native_clr$Median.Longitude,
                data_native_clr$Median.Latitude)

colnames(lonlat) <- c("longitude","latitude")
rownames(lonlat) <- data_native_clr$Basin.Name

rL.basin <- HmscRandomLevel(sData = lonlat, longlat = TRUE, sMethod = "NNGP")
rL.basin <- setPriors(rL.basin, nfMin = 1,nfMax = 1)

rL.basin.NS <- HmscRandomLevel(units = levels(studyDesign$basin))
rL.basin.NS <- setPriors(rL.basin.NS, nfMin = 1,nfMax = 1)

# Response variables (strategies after CLR transformation)
Y <- cbind(
  data_native_clr$Equilibrium_CLR,
  data_native_clr$Opportunistic_CLR,
  data_native_clr$Periodic_CLR)

colnames(Y) = c("E","O","P")

# Explanatory variables (including interactions)
XData <- data.frame(
  MAT = data_native_clr$MAT,
  TS = data_native_clr$TS,
  DIS = data_native_clr$DIS,
  RAD = data_native_clr$RAD,
  NSR = data_native_clr$NSR,
  RBA = data_native_clr$RBA,
  GDP = data_native_clr$GDP,
  CSI = data_native_clr$CSI)

XFormula <- ~ MAT + TS + DIS + RAD + NSR + RBA + GDP + CSI

# Build HMSC model
m_native <- Hmsc(
  Y = Y,
  XData = XData,
  XFormula = XFormula,
  studyDesign = studyDesign,
  ranLevels = list(ecoregion = rL.ecoregion, basin = rL.basin)
)

m_native.NS <- Hmsc(
  Y = Y,
  XData = XData,
  XFormula = XFormula,
  studyDesign = studyDesign,
  ranLevels = list(ecoregion = rL.ecoregion, basin = rL.basin.NS)
)

models = list()
models$native = m_native
models$native.NS = m_native.NS

################################################################################################### model-Non-native
# Plot point distribution
plot(data_non_native_clr$Median.Longitude, data_non_native_clr$Median.Latitude,
     xlab="Longitude", ylab="Latitude", main="non_native Basins")

# Study design
studyDesign <- data.frame(
  basin = as.factor(data_non_native_clr$Basin.Name),
  ecoregion = as.factor(data_non_native_clr$Ecoregion))

# Random effects
rL.ecoregion <- HmscRandomLevel(units = levels(studyDesign$ecoregion))
rL.ecoregion <- setPriors(rL.ecoregion,nfMin = 1, nfMax=1)

lonlat <- cbind(data_non_native_clr$Median.Longitude,
                data_non_native_clr$Median.Latitude)

colnames(lonlat) <- c("longitude","latitude")
rownames(lonlat) <- data_non_native_clr$Basin.Name

rL.basin <- HmscRandomLevel(sData = lonlat, longlat = TRUE, sMethod = "NNGP")
rL.basin <- setPriors(rL.basin, nfMin = 1,nfMax = 1)

rL.basin.NS <- HmscRandomLevel(units = levels(studyDesign$basin))
rL.basin.NS <- setPriors(rL.basin.NS, nfMin = 1,nfMax = 1)

# Response variables (strategies after CLR transformation)
Y <- cbind(
  data_non_native_clr$Equilibrium_CLR,
  data_non_native_clr$Opportunistic_CLR,
  data_non_native_clr$Periodic_CLR)

colnames(Y) = c("E","O","P")

# Explanatory variables (including interactions)
XData <- data.frame(
  MAT = data_non_native_clr$MAT,
  TS = data_non_native_clr$TS,
  DIS = data_non_native_clr$DIS,
  RAD = data_non_native_clr$RAD,
  NSR = data_non_native_clr$NSR,
  RBA = data_non_native_clr$RBA,
  GDP = data_non_native_clr$GDP,
  CSI = data_non_native_clr$CSI)

XFormula <- ~ MAT + TS + DIS + RAD + NSR + RBA + GDP + CSI

# Build HMSC model
m_non_native <- Hmsc(
  Y = Y,
  XData = XData,
  XFormula = XFormula,
  studyDesign = studyDesign,
  ranLevels = list(ecoregion = rL.ecoregion, basin = rL.basin)
)

m_non_native.NS <- Hmsc(
  Y = Y,
  XData = XData,
  XFormula = XFormula,
  studyDesign = studyDesign,
  ranLevels = list(ecoregion = rL.ecoregion, basin = rL.basin.NS)
)

models$non_native = m_non_native
models$non_native.NS = m_non_native.NS

save(models,file = "models/unfitted_models.RData")
