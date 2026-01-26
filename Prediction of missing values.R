
###########################################################################################################Prediction of missing values-Rphylopars
setwd("xx")
library(ape)
library(phytools)
library(dplyr)
library(Rphylopars)

# load data
trait_data <- read.csv("Final trait data of all global freshwater fish-for prediction-Rphylopars.csv", header = TRUE)
na_proportion <- colMeans(is.na(trait_data))  # Calculate the proportion of NA in each column
print(na_proportion)
mean(na_proportion[-1])

new_tree <- read.newick("fish_phylogenetic_tree.tre")

####Ensure that the data matches the species of the tree
#Obtain the names of the species in the data and the tree
species_in_trait_data <- trait_data$species
species_in_tree <- new_tree$tip.label

# Identify the species in the data but not in the phylogenetic tree
missing_species <- setdiff(species_in_trait_data, species_in_tree)
print(missing_species) 

# pruning tree
tree_trimmed <- drop.tip(new_tree, setdiff(new_tree$tip.label, trait_data$species))
print(tree_trimmed)

####################################################Filter the species and check if all of them are in the pruned tree
#Eliminate species that are not in the tree from the data
trait_data_filtered <- trait_data[trait_data$species %in% tree_trimmed$tip.label, ]
head(trait_data_filtered)

#Check if all the species are in the pruned tree
missing_species <- setdiff(trait_data_filtered$species, tree_trimmed$tip.label)
print(missing_species) 
write.csv(missing_species, "missing species in phylogenetic tree.csv", row.names = FALSE) 

# ********** Take the natural log transformation for all traits **********
log_trait_data <- trait_data_filtered
log_trait_data[, -1] <- log(log_trait_data[, -1])

# Imputation
imputed_results <- phylopars(trait_data=log_trait_data, tree = tree_trimmed, model = "BM") #change to different model
imputed_results
AIC(imputed_results)
BIC(imputed_results)

imputed_means <- imputed_results$anc_recon[1:14679,]
head(imputed_means)
imputed_var <- imputed_results$anc_var[1:14679,]

# ********** Perform the inverse exp transform on the imputation result **********
imputed_means <- exp(imputed_means) 
imputed_var <- (imputed_means^2) * imputed_var

# ***  Handle the rounding of the parental.care and limit it to between 1 and 3***
if ("Parental.care" %in% colnames(imputed_means)) {
  imputed_means[, "Parental.care"] <- round(imputed_means[, "Parental.care"])
  imputed_means[, "Parental.care"] <- pmax(pmin(imputed_means[, "Parental.care"], 3), 1)
} else {
  warning("The column Parental.care was not found in imputed_means. Please confirm if the column name is correct")
}
head(imputed_means)

# Saving result
write.csv(imputed_means, "imputed_means_new_tree-BM.csv", row.names = TRUE)
write.csv(imputed_var, "imputed_var_new_tree-BM.csv", row.names = TRUE)



##############################################################################Evaluate the performance of the model - use all non-NA traits for evaluation
library(ape)
library(phytools)
library(dplyr)
library(Rphylopars)
library(tibble)

# load data
trait_data <- read.csv("Final trait data of all global freshwater fish-for prediction-Rphylopars.csv", header = TRUE)
new_tree <- read.newick("fish_phylogenetic_tree.tre")
# pruning tree
tree_trimmed <- drop.tip(new_tree, setdiff(new_tree$tip.label, trait_data$species))

# Filter the data and only retain the species on the tree
trait_data_filtered <- trait_data[trait_data$species %in% tree_trimmed$tip.label, ]
rownames(trait_data_filtered) <- trait_data_filtered$species
# natural log transformation
log_trait_data <- trait_data_filtered
log_trait_data[, -1] <- log(log_trait_data[, -1])

# ---- Evaluate the imputation performance ----
# Only select the complete subset data without NA
complete_cases <- log_trait_data %>% filter(if_all(-species, ~ !is.na(.)))
rownames(complete_cases) <- complete_cases$species
# Only retain trees within a complete subset of the data
tree_complete <- drop.tip(tree_trimmed, setdiff(tree_trimmed$tip.label, complete_cases$species))

#Set the number of repeats
n_iter <- 100
spearman_results <- matrix(NA, nrow = n_iter, ncol = ncol(complete_cases) - 1)
colnames(spearman_results) <- colnames(complete_cases)[-1]
# Record the AIC and BIC of each model fit
aic_bic_results <- matrix(NA, nrow = n_iter, ncol = 2)
colnames(aic_bic_results) <- c("AIC", "BIC")

# Exclude traits that should not be created as missing
excluded_traits <- c("Max.total.length", "Trophic.level")
trait_names <- colnames(complete_cases)[-1]
target_trait_cols <- which(!trait_names %in% excluded_traits)

set.seed(123)
for (i in 1:n_iter) {
  # Create a 25% missing (only in the target traits)
  simulated_data <- complete_cases
  trait_matrix <- as.matrix(simulated_data[, -1])
  n_total <- nrow(trait_matrix) * length(target_trait_cols)
  n_missing <- round(0.25 * n_total)
  missing_rows <- sample(1:nrow(trait_matrix), n_missing, replace = TRUE)
  missing_cols <- sample(target_trait_cols, n_missing, replace = TRUE)
  missing_indices <- cbind(missing_rows, missing_cols)
  
  for (k in 1:nrow(missing_indices)) {
    trait_matrix[missing_indices[k,1], missing_indices[k,2]] <- NA
  }
  simulated_data[, -1] <- trait_matrix
  # imputation
  imputed <- phylopars(trait_data = simulated_data, tree = tree_complete, model = "BM") #change to different model
  imputed_means <- imputed$anc_recon
  # record AIC and BIC
  aic_bic_results[i, "AIC"] <- AIC(imputed)
  aic_bic_results[i, "BIC"] <- BIC(imputed)

  # exp transformation
  imputed_means <- exp(imputed_means)
  true_values <- exp(as.matrix(complete_cases[, -1]))
  imputed_means <- imputed_means[rownames(true_values), colnames(true_values)]
  
  for (j in 1:ncol(true_values)) {
    trait_name <- colnames(true_values)[j]
    if (trait_name %in% excluded_traits) next
    
    missing_species <- rownames(true_values)[missing_indices[missing_indices[,2] == j, 1]]
    
    if (length(missing_species) > 1) {
      predicted <- imputed_means[missing_species, trait_name]
      observed  <- true_values[missing_species, trait_name]

      if (trait_name == "Parental.care") {
        predicted <- round(predicted)
        predicted[predicted < 1] <- 1
        predicted[predicted > 3] <- 3
      }
      
      spearman_cor <- cor(predicted, observed, method = "spearman")
      spearman_results[i, j] <- spearman_cor
    }
  }
}
# summary
spearman_summary <- data.frame(
  Trait = colnames(spearman_results),
  Mean_Spearman = apply(spearman_results, 2, mean, na.rm = TRUE),
  SD_Spearman = apply(spearman_results, 2, sd, na.rm = TRUE)
)
# results
print(spearman_summary)
print(aic_bic_results)
mean_spearman_across_traits <- mean(spearman_summary$Mean_Spearman, na.rm = TRUE)
cat("Mean Spearman correlation across traits:", mean_spearman_across_traits, "\n")

# saving results
write.csv(spearman_summary, "spearman_summary_model.csv", row.names = FALSE)
write.csv(aic_bic_results, "aic_bic_results-BM.csv", row.names = FALSE)


##################################################################################Evaluate the  performance of the model - use the mean of the same genus for evaluation
setwd("xx")

library(dplyr)
# Load data
trait_data <- read.csv("Final trait data of all global freshwater fish-for prediction-mean.csv", header = TRUE)
# Load missing species in the phylogenetic tree
missing_species <- read.csv("missing species in phylogenetic tree.csv", header = TRUE)
# Delete the missing species in trait_data
missing_species$x <- gsub("_", " ", missing_species$x)
trait_data_filtered <- trait_data[!trait_data$species %in% missing_species$x, ]

# Screen for species with all six traits without NA
complete_cases <- trait_data_filtered %>%
  filter(if_all(6:ncol(.), ~ !is.na(.)))
rownames(complete_cases) <- complete_cases$species

trait_cols <- 6:ncol(complete_cases)
trait_names <- names(complete_cases)[trait_cols]

# Exclude the traits that do not participate in the missing simulation
excluded_traits <- c("Max.total.length", "Trophic.level")
included_trait_idx <- which(!(trait_names %in% excluded_traits))
included_cols <- trait_cols[included_trait_idx]

# Genus information
genus_vector <- complete_cases$Genus

# Initialize the result matrix
n_iter <- 100
n_traits <- length(trait_cols)
spearman_genusMean <- matrix(NA, nrow = n_iter, ncol = n_traits)
colnames(spearman_genusMean) <- trait_names

set.seed(123)

for (i in 1:n_iter) {
  cat("Simulation", i, "\n")
  
  # Create a simulated missing matrix
  trait_simulated <- as.matrix(complete_cases[, trait_cols])
  n_total <- length(included_cols) * nrow(trait_simulated)
  n_missing <- round(0.25 * n_total)
  
  # Avoid sampling at duplicate positions
  possible_positions <- expand.grid(row = 1:nrow(trait_simulated), col = included_cols)
  sampled_rowcol <- possible_positions[sample(nrow(possible_positions), n_missing), ]
  
  for (k in 1:n_missing) {
    col_idx <- sampled_rowcol$col[k] - 5
    trait_simulated[sampled_rowcol$row[k], col_idx] <- NA
  }
  
  # -------- imputation with genus mean value --------
  imputed_gm <- trait_simulated
  
  for (j in included_cols) {
    j_trait_idx <- j - 5
    na_rows <- which(is.na(imputed_gm[, j_trait_idx]))
    
    for (r in na_rows) {
      genus_r <- genus_vector[r]
      genus_rows <- which(genus_vector == genus_r & !is.na(imputed_gm[, j_trait_idx]))
      if (length(genus_rows) > 0) {
        imputed_gm[r, j_trait_idx] <- mean(imputed_gm[genus_rows, j_trait_idx], na.rm = TRUE)
      }
    }
  }
  
  true_matrix <- as.matrix(complete_cases[, trait_cols])
  
  for (j in 1:n_traits) {
    missing_rows <- sampled_rowcol$row[sampled_rowcol$col == trait_cols[j]]
    if (length(missing_rows) > 1) {
      spearman_genusMean[i, j] <- suppressWarnings(cor(
        imputed_gm[missing_rows, j],
        true_matrix[missing_rows, j],
        method = "spearman",
        use = "complete.obs"
      ))
    }
  }
}

# Summary
summary_gm <- data.frame(
  Trait = colnames(spearman_genusMean),
  Method = "GenusMean",
  Mean_Spearman = apply(spearman_genusMean, 2, mean, na.rm = TRUE),
  SD_Spearman = apply(spearman_genusMean, 2, sd, na.rm = TRUE)
)

print(summary_gm)
mean_spearman_across_traits <- mean(summary_gm$Mean_Spearman, na.rm = TRUE)
cat("Mean Spearman correlation across traits:", mean_spearman_across_traits, "\n")

# Save result
write.csv(summary_gm, "spearman_summary_genusMean.csv", row.names = FALSE)


