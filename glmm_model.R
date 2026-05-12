###############################################################################################       
############################################################################################### Model-Similarity 
############################################################################################### 
setwd("XXXX")

############################################################# glmm
library(glmmTMB)
library(nlme)
library(DHARMa)
library(car)

# Load data 
data <- read.csv("Similarity_with_predictors_model_1720-2.csv")
# Remove missing values in explanatory variables 
data <- na.omit(data)
head(data)
hist(data$Similarity_Centroid)
# Standardize fixed effect variables (z-score) 
predictors <- c("MAT", "TS", "PRD","SRD","CSI", "NSR", "RBA", "GDP")
data[predictors] <- scale(data[predictors])

# vif check
lm_check1 <- lm(Similarity_Centroid ~ MAT+TS+PRD+SRD+CSI+NSR+RBA+GDP, data = data)
vif(lm_check1)

lm_check2 <- lm(Similarity_Centroid ~ MAT+TS+PRD+SRD+CSI+NSR+RBA+GDP+ MAT*TS+ PRD*SRD+ MAT*PRD+ MAT*SRD+ TS*PRD+ TS*SRD, data = data) 
vif(lm_check2)

lm_check3 <- lm(Similarity_Centroid ~ MAT+TS+PRD+SRD+CSI+NSR+RBA+GDP +MAT*PRD+ MAT*SRD+ TS*PRD+ TS*SRD, data = data)
vif(lm_check3)

lm_check4 <- lm(Similarity_Centroid ~ MAT+TS+PRD+SRD+CSI+NSR+RBA+GDP+ MAT*TS+ MAT*PRD+ MAT*SRD+ TS*PRD+ TS*SRD, data = data)
vif(lm_check4)

lm_check5 <- lm(Similarity_Centroid ~ MAT+TS+PRD+SRD+CSI+NSR+RBA+GDP+ PRD*SRD+ MAT*PRD+ MAT*SRD+ TS*PRD+ TS*SRD, data = data)
vif(lm_check5)


lm_check6 <- lm(Similarity_Centroid ~ MAT+TS+PRD+SRD+CSI+NSR+RBA+GDP+ TS*GDP + PRD*GDP+ TS*RBA + PRD*RBA, data = data)
vif(lm_check6)

############################################################# glmm + spatial effect
data$pos <- numFactor(scale(data$Med_Longit), scale(data$Med_Latit))
data$group <- factor(rep(1, nrow(data)))

#Model-1 only main factors
glmm_s1 <- glmmTMB(
  Similarity_Centroid  ~ MAT+TS+PRD+SRD+CSI+NSR+RBA+GDP + mat(pos + 0 | group), data = data, family = gaussian(link = "logit"))
summary(glmm_s1)

#Model-3-5 with interaction terms
glmm_s2 <- glmmTMB(Similarity_Centroid  ~ MAT+TS+PRD+SRD+CSI+NSR+RBA+GDP+ MAT*TS+ PRD*SRD+ MAT*PRD+ MAT*SRD+ TS*PRD+ TS*SRD + mat(pos + 0 | group), 
                   data = data, family = gaussian(link = "logit"))
summary(glmm_s2)

glmm_s3 <- glmmTMB(Similarity_Centroid  ~ MAT+TS+PRD+SRD+CSI+NSR+RBA+GDP +MAT*PRD+ MAT*SRD+ TS*PRD+ TS*SRD + mat(pos + 0 | group), 
                   data = data, family = gaussian(link = "logit"))
summary(glmm_s3)

glmm_s4 <- glmmTMB(Similarity_Centroid  ~ MAT+TS+PRD+SRD+CSI+NSR+RBA+GDP+ MAT*TS+ MAT*PRD+ MAT*SRD+ TS*PRD+ TS*SRD + mat(pos + 0 | group), 
                   data = data, family = gaussian(link = "logit"))
summary(glmm_s4)

glmm_s5 <- glmmTMB(Similarity_Centroid  ~ MAT+TS+PRD+SRD+CSI+NSR+RBA+GDP+ PRD*SRD+ MAT*PRD+ MAT*SRD+ TS*PRD+ TS*SRD + mat(pos + 0 | group), 
                   data = data, family = gaussian(link = "logit"))
summary(glmm_s5)


#Model-6 with interaction terms consider the impact of GDP and RBA
glmm_s6 <- glmmTMB(Similarity_Centroid  ~ MAT+TS+PRD+SRD+CSI+NSR+RBA+GDP+ TS*GDP + PRD*GDP+ TS*RBA + PRD*RBA + mat(pos + 0 | group), 
                   data = data, family = gaussian(link = "logit"))
summary(glmm_s6)


#Check model residuals
sims5 <- simulateResiduals(glmm_s5)
plot(sims5)
res_raw <- residuals(glmm_s5, type = "response")
hist(res_raw,
     breaks = 20,
     col = "grey80",
     border = "black",
     main = "Histogram of model residuals",
     xlab = "Residuals")

# Calculate Moran's I
library(spdep)
coords <- cbind(data$Med_Longit, data$Med_Latit)
nb <- knn2nb(knearneigh(coords, k=4))
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
resid_vals <- residuals(glmm_s5, type = "pearson")
moran.test(resid_vals, lw, zero.policy = TRUE)

# Calculate R²
library(MuMIn)
r2_values <- r.squaredGLMM(glmm_s5)
r2_values


# Plot
library(graphics)
library(sjPlot)
library(ggplot2)

p<-plot_model(glmm_s5, type = "pred", terms = c("TS"),colors ="#FF9900")
p+geom_line(size=1.5,colour="#FF9900")+theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(size =16,colour="black"),
        axis.text.y = element_text(size =16,colour="black"),
        axis.title.x = element_text(size =16,colour="black"),
        axis.title.y = element_text(size =16,colour="black"))+
  theme()+
  labs(x="TS",y="Predicted values of similarity")+theme(axis.title = element_text(size = 12))+
  theme(legend.position = "top",plot.title = element_blank())

ggsave("Marginal effect TS and similarity.pdf", width = 5.5, height = 5)


# Plot
p<-plot_model(glmm_s5, type = "pred", terms = c("PRD"),colors =c("#146C36"))
p+geom_line(size=1.5,colour="#146C36")+theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(size =16,colour="black"),
        axis.text.y = element_text(size =16,colour="black"),
        axis.title.x = element_text(size =16,colour="black"),
        axis.title.y = element_text(size =16,colour="black"))+
  labs(x="PRD",y="Predicted values of Similarity")+theme(axis.title = element_text(size = 12))+
  theme(legend.position = "top",plot.title = element_blank())

ggsave("Marginal effect PRD and similarity.pdf", width = 5.5, height = 5)
