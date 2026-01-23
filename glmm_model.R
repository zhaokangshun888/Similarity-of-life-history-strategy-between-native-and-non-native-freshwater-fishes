
###############################################################################################       
############################################################################################### Model-Similarity 
############################################################################################### 
setwd("XXXX")

############################################################# glmm
library(glmmTMB)
library(nlme)
library(DHARMa)
library(car)

#######Load data 
data <- read.csv("Similarity_with_predictors_model_1720.csv")
######Remove missing values in explanatory variables 
data <- na.omit(data)
head(data)
hist(data$Similarity_Centroid)

######Standardize fixed effect variables (z-score) 
predictors <- c("MAT", "TS", "DIS", "RAD", "GDP", "NSR", "RBA", "CSI")
data[predictors] <- scale(data[predictors])

###### VIF check
#models examining the effects of environmental and anthropogenic factors on similarity between native and non-native freshwater fishes at basin scale
vif_Model1 <- lm(Similarity_Centroid ~ MAT + TS + DIS + RAD + NSR + RBA + GDP + CSI, data = data)
vif(vif_Model1)

vif_Model2 <- lm(Similarity_Centroid ~ MAT + TS + DIS + RAD + NSR + RBA + GDP + CSI + 
                   MAT*TS + DIS*RAD + MAT*DIS + MAT*RAD + TS*RAD, data = data) #TS*DIS was removed due to its high VIF
vif(vif_Model2)

vif_Model3 <- lm(Similarity_Centroid  ~ MAT + TS + DIS + RAD + NSR + RBA + GDP + CSI + 
                   DIS*RAD + MAT*DIS + MAT*RAD + TS*RAD, data = data)  #TS*DIS was removed due to its high VIF
vif(vif_Model3)

vif_Model4 <- lm(Similarity_Centroid  ~ MAT + TS + DIS + RAD + NSR + RBA + GDP + CSI + 
                   MAT*TS + MAT*DIS + MAT*RAD + TS*RAD, data = data) #TS*DIS was removed due to its high VIF
vif(vif_Model4)

#model examining whether the main anthropogenic driver modifies ecological relationships between main environmental factors and the similarity between native and non-native freshwater fishes at the basin scale
vif_Model <- lm(Similarity_Centroid ~ MAT + TS + DIS + RAD + NSR + RBA + GDP + CSI + 
                  MAT*TS + MAT*DIS+ MAT*RAD + TS*RAD + MAT*GDP + TS*GDP, data = data)#TS*DIS was removed due to its high VIF
vif(vif_Model)



############################################################# glmm + spatial effect
data$pos <- numFactor(scale(data$Med_Longit), scale(data$Med_Latit))
data$group <- factor(rep(1, nrow(data)))  # Dummy group for global spatial structure 
  
#Model-1 only main factors
glmm_s1 <- glmmTMB(
  Similarity_Centroid  ~ MAT+TS + DIS+RAD + NSR + RBA + GDP + CSI  + (1|Ecoregion) + mat(pos + 0 | group), data = data, family = gaussian(link = "logit"))

summary(glmm_s1)
sims1 <- simulateResiduals(glmm_s1)
plot(sims1)

#Model-4 remove DIS*RAD and retain MAT*TS
glmm_s4 <- glmmTMB(
  Similarity_Centroid  ~ MAT+TS + DIS+RAD + NSR + RBA + GDP + CSI + MAT*TS + MAT*DIS+ MAT*RAD + TS*RAD + (1|Ecoregion) + mat(pos + 0 | group), 
  data = data, family = gaussian(link = "logit"))

summary(glmm_s4)
sims4 <- simulateResiduals(glmm_s4)
plot(sims4)

#Model with interaction terms consider the impact of human activities-GDP 
glmm_s <- glmmTMB(
  Similarity_Centroid  ~ MAT+TS + DIS+RAD + NSR + RBA + GDP + CSI + MAT*TS + MAT*DIS+ MAT*RAD + TS*RAD + MAT*GDP + TS*GDP+ (1|Ecoregion) + mat(pos + 0 | group), 
  data = data, family = gaussian(link = "logit"))

summary(glmm_s)
sims <- simulateResiduals(glmm_s)
plot(sims)



#######################################################check the spatial autocorrelation and R2
######Moran's I
library(spdep)
coords <- cbind(data$Med_Longit, data$Med_Latit)
nb <- knn2nb(knearneigh(coords, k=4))
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
resid_vals <- residuals(glmm_s4, type = "pearson")
moran.test(resid_vals, lw, zero.policy = TRUE)

######R2
library(MuMIn)
r2_values <- r.squaredGLMM(glmm_s4)
r2_values

# Plot marginal effect
library(graphics)
library(sjPlot)
library(ggplot2)

p<-plot_model(glmm_s6, type = "pred", terms = c("MAT"),colors ="#FF9900")
p+geom_line(size=1.5,colour="#FF9900")+theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(size =16,colour="black"),
        axis.text.y = element_text(size =16,colour="black"),
        axis.title.x = element_text(size =16,colour="black"),
        axis.title.y = element_text(size =16,colour="black"))+
  theme()+
  labs(x="MAT",y="Predicted values of similarity")+theme(axis.title = element_text(size = 12))+
  theme(legend.position = "top",plot.title = element_blank())

ggsave("Marginal effect MAT and similarity.pdf", width = 5.5, height = 5)

 
# Plot
p<-plot_model(glmm_s6, type = "pred", terms = c("TS"),colors =c("#146C36"))
p+geom_line(size=1.5,colour="#146C36")+theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(size =16,colour="black"),
        axis.text.y = element_text(size =16,colour="black"),
        axis.title.x = element_text(size =16,colour="black"),
        axis.title.y = element_text(size =16,colour="black"))+
  labs(x="TS",y="Predicted values of Similarity")+theme(axis.title = element_text(size = 12))+
  theme(legend.position = "top",plot.title = element_blank())

ggsave("Marginal effect TS and similarity.pdf", width = 5.5, height = 5)



