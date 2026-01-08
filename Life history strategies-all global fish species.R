
########################################################################################## Archetypal analysis (AA)
setwd("xx")

library(ggplot2)
library(factoextra)
library(psych)
library(GPArotation)
library(archetypes)  # Archetypal Analysis (AA)
library(dplyr)

# Load data
trait_data <- read.csv("imputed_means_new_tree-BM.csv", header = TRUE)

# Reset the row name
rownames(trait_data) <- trait_data$species
head(trait_data)
traits <- names(trait_data)[2:ncol(trait_data)] 
par(mfrow=c(3, 3))  

for (trait in traits) {
  hist(trait_data[[trait]], main = trait, xlab = trait)
}

# log transformation
log_columns <- c("Life.span", "Fecundity", "Age.at.maturation", "Offspring.size", "Max.total.length","K.value")
trait_data[log_columns] <- log(trait_data[log_columns])

# Select all columns except the first one to draw a histogram and check the distribution
mapply(function(x, name) hist(x, main = name, xlab = name), 
       trait_data[, -1], 
       names(trait_data)[-1])

# Standardized data
trait_data_scaled <- scale(trait_data[, -1])  

# Convert it to data.frame for mapply() 
trait_data_scaled_df <- as.data.frame(trait_data_scaled)

# Select the standardized data and draw a histogram
mapply(function(x, name) hist(x, main = paste("Scaled", name), xlab = name), 
       trait_data_scaled_df, 
       names(trait_data_scaled_df))


# Archetypal Analysis （the number of Archetypes k = 1: 10)
set.seed(123)
aa_result <- stepArchetypes(trait_data_scaled, k = 1:10, nrep = 10,family = archetypesFamily("original"))
rss(aa_result)
screeplot(aa_result)

# It can be adjusted according to the elbow method. Here, k=3 is selected
aa_3 <- bestModel(aa_result[[3]])

# The Archetype parameters are retained to three decimal places
archetype_traits <- as.data.frame(round(t(parameters(aa_3)), 3))
print(archetype_traits)

# Obtain the contribution of each species on each Archetype
aa_3_profile <- coef(aa_3)

# Set row name
rownames(aa_3_profile) <- rownames(trait_data_scaled)

# Identify the main Archetype (the one with the greatest contribution) of each species
aa_3_cluster <- max.col(aa_3_profile)
table(aa_3_cluster)

# Match Archetype type with species name
archetype_results <- data.frame(
  species = rownames(aa_3_profile),  # species name
  archetype = aa_3_cluster,        # Archetype of each species
  aa_3_profile                     # The composition of each species of different Archetypes
)
# Rename the contribution column
colnames(archetype_results)[-(1:2)] <- paste0("Archetype_", 1:3)

# Add the original trait data to the result
final_results <- left_join(archetype_results, trait_data, by = "species")

# Result
print(head(final_results))
final_results$strategies <- factor(final_results$archetype,
                                levels = c(1, 2, 3),
                                labels = c("Equilibrium","Opportunistic", "Periodic"))
# Save result
write.csv(final_results, "archetype_results_robust.csv", row.names = FALSE)
dev.off()

png("screeplot_aa_result.png", width = 600, height = 600, res = 150)
screeplot(aa_result)
dev.off()



##################################Principal component analysis (PCA) of the eight life-history traits and life-history strategy endpoints identified based on archetypal analysis (AA)
library(ggrepel)
# PCA
pca_result <- prcomp(trait_data_scaled)
summary(pca_result)

# Calculate the variance of each PC axis
explained_var <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100  
pc1_var <- round(explained_var[1], 1)  
pc2_var <- round(explained_var[2], 1)
pc3_var <- round(explained_var[3], 1)

# Extract the first three principal components of the PCA results and add the species names
pca_scores_df <- as.data.frame(pca_result$x[, 1:3]) 
pca_scores_df$species <- rownames(aa_3_profile) 
pca_scores_df$cluster <- factor(aa_3_cluster) 

# Calculate the coordinates of archetype in the PCA space
archetype_pca_scores <- predict(pca_result, newdata = parameters(aa_3))
archetype_pca_df <- as.data.frame(archetype_pca_scores[, 1:3])
archetype_pca_df$Archetype <- paste0(1:nrow(archetype_pca_df))  # Tag archetype

# PCA for extracting feature traits
loadings_df <- as.data.frame(pca_result$rotation[, 1:3])  # the first three PC axis
loadings_df$trait <- recode(rownames(loadings_df),
                            "Life.span" = "Life span",
                            "Fecundity" = "Fecundity",
                            "Parental.care" = "Parental care",
                            "Age.at.maturation" = "Age at maturation",
                            "Offspring.size" = "Offspring size",
                            "Max.total.length" = "Max total length",
                            "Trophic.level" = "Trophic level",
                            "K.value" = "growth coefficient")
# Plot PCA result
# Set cluster name 
pca_scores_df$cluster <- factor(pca_scores_df$cluster,
                                levels = c(1, 2, 3),
                                labels = c("Equilibrium","Opportunistic", "Periodic"))

# Tag Archetype
archetype_pca_df$Archetype <- c("Equilibrium","Opportunistic", "Periodic")

ggplot() +
  # Species data points
  geom_point(data = pca_scores_df, aes(x = PC1, y = PC2, color = "grey10"), alpha = 0.3) +
  # Archetype 
  geom_point(data = archetype_pca_df, aes(x = PC1, y = PC2, color = Archetype), size = 5) +
  # Tag Archetype
  geom_text_repel(data = archetype_pca_df, aes(x = PC1, y = PC2, label = Archetype), 
                  color = "black", size = 6,nudge_x = -0.1,nudge_y = -0.6, segment.color = "NA")+
  # Add trait variables
  geom_segment(data = loadings_df, 
               aes(x = 0, y = 0, xend = PC1 * 7, yend = PC2 * 7), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "grey10", size = 0.5) +
  geom_text_repel(data = loadings_df, 
                  aes(x = PC1 * 5.8, y = PC2 * 7.1, label = trait), 
                  color = "grey10", size = 4, segment.color = "NA") +
  # Set colour
  scale_color_manual(name = "Archetype Group",
    values = c("Opportunistic" = "#1f78b4", "Periodic" = "#33a02c", "Equilibrium" = "#ff7f00"))+
  
  # Set the X-axis and Y-axis labels
  labs(title = "",
       x = paste0("PC1 (", pc1_var, "%)"),
       y = paste0("PC2 (", pc2_var, "%)"),
       color = "Archetype Group") +
  
  # Theme
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black", size = 1.2),
    panel.grid.major = element_line(color = "NA"),
    panel.grid.minor = element_blank(), 
    legend.position = "",
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title = element_text(size = 20))

ggsave("PCA_Biplot_PC12.png", width = 7, height = 7, dpi = 300)
