###########################################################################################################Prediction of missing values-missForest
###########################################################################################################
############################################### Prediction of missing values - missForest (with phylogeny)
rm(list=ls())
setwd("xxx")

library(dplyr)
library(missForest)
library(ape)

# Load data
trait_data <- read.csv("Final trait data of all global freshwater fish-for prediction-2.csv", header = TRUE)

# Load phylogenetic tree
tree <- read.tree("fish_phylogenetic_tree.tre")

# Only retain the species shared by the tree and the data
trait_data2 <- trait_data[trait_data$species %in% tree$tip.label, ]
tree_trimmed <- drop.tip(tree, setdiff(tree$tip.label, trait_data2$species))

# Align in tree order
trait_data2 <- trait_data2[match(tree_trimmed$tip.label, trait_data2$species), ]

# Calculate the phylogenetic distance and perform PCoA
phylo_dist <- cophenetic.phylo(tree_trimmed)
phylo_pcoa <- cmdscale(phylo_dist, k = 10, eig = TRUE)

phylo_axes <- as.data.frame(phylo_pcoa$points)
colnames(phylo_axes) <- paste0("PhyloPC", 1:ncol(phylo_axes))

# Split data
taxa_info  <- trait_data2[, 1:5]
trait_only <- trait_data2[, 6:ncol(trait_data2)]

trait_names <- names(trait_only)

# log transform: exclude Parental care 
log_trait <- trait_only
log_transform_cols <- which(trait_names != "Parental.care")
log_trait[, log_transform_cols] <- log(log_trait[, log_transform_cols])

# trait + phylo axes
mf_input <- cbind(log_trait, phylo_axes)

# Data imputation 
set.seed(123)
imputed_mf <- missForest(mf_input, verbose = FALSE)

# Only extract the trait part
imputed_traits <- imputed_mf$ximp[, trait_names]

# Restore the original scale
imputed_traits[, log_transform_cols] <- exp(imputed_traits[, log_transform_cols])

# Parental.care treatment
pc_idx <- which(names(imputed_traits) == "Parental.care")
if (length(pc_idx) == 1) {
  imputed_traits[, pc_idx] <- round(imputed_traits[, pc_idx])
  imputed_traits[, pc_idx][imputed_traits[, pc_idx] < 1] <- 1
  imputed_traits[, pc_idx][imputed_traits[, pc_idx] > 3] <- 3
}

# Output result
imputed_data_full <- cbind(taxa_info, imputed_traits)
write.csv(imputed_data_full, "imputed_trait_data_missForest_with_phylogeny-2.csv", row.names = FALSE)
print(imputed_mf$OOBerror)


#################################################################
# Evaluate the predictive performance of missForest

library(dplyr)
library(missForest)
library(ape)

# Load data
trait_data <- read.csv("Final trait data of all global freshwater fish-for prediction-2.csv", header = TRUE)

# Load phylogenetic tree
tree <- read.tree("fish_phylogenetic_tree.tre")

# Only retain the species shared by the tree and the data
trait_data2 <- trait_data[trait_data$species %in% tree$tip.label, ]
tree_trimmed <- drop.tip(tree, setdiff(tree$tip.label, trait_data2$species))

# Align in tree order
trait_data2 <- trait_data2[match(tree_trimmed$tip.label, trait_data2$species), ]

# Calculate the phylogenetic distance and perform PCoA
phylo_dist <- cophenetic.phylo(tree_trimmed)
phylo_pcoa <- cmdscale(phylo_dist, k = 10, eig = TRUE)

phylo_axes <- as.data.frame(phylo_pcoa$points)
colnames(phylo_axes) <- paste0("PhyloPC", 1:ncol(phylo_axes))

# Combine data
trait_data2 <- cbind(trait_data2, phylo_axes)

# Filter species with trait completeness
trait_cols  <- 6:(ncol(trait_data2) - ncol(phylo_axes))
trait_names <- names(trait_data2)[trait_cols]

complete_cases <- trait_data2 %>%
  filter(if_all(all_of(trait_names), ~ !is.na(.)))

rownames(complete_cases) <- complete_cases$species

# trait matrix
trait_matrix_full <- complete_cases[, trait_cols]

# log transform: exclude Parental care 
log_trait_matrix <- trait_matrix_full
log_transform_cols <- which(trait_names != "Parental.care")
log_trait_matrix[, log_transform_cols] <- log(log_trait_matrix[, log_transform_cols])

# phylo matrix
phylo_cols <- grep("^PhyloPC", names(complete_cases))
phylo_matrix <- complete_cases[, phylo_cols]

# Initialize the result matrix
n_iter <- 100
n_traits <- length(trait_names)
spearman_missForest_phylo <- matrix(NA, nrow = n_iter, ncol = n_traits)
colnames(spearman_missForest_phylo) <- trait_names


# No trait is excluded
excluded_traits <- NULL 
target_trait_cols <- which(!trait_names %in% excluded_traits)

set.seed(123)
for (i in 1:n_iter) {
  cat("Simulation", i, "\n")
  
  trait_simulated <- as.data.frame(log_trait_matrix)
  
  possible_positions <- expand.grid(
    row = 1:nrow(trait_simulated),
    col = target_trait_cols
  )
  
  n_total <- nrow(possible_positions)
  n_missing <- round(0.25 * n_total)
  sampled_idx <- sample(1:n_total, n_missing, replace = FALSE)
  sampled_rowcol <- possible_positions[sampled_idx, ]
  
  for (k in 1:nrow(sampled_rowcol)) {
    trait_simulated[sampled_rowcol$row[k], sampled_rowcol$col[k]] <- NA
  }
  
  # missForest input：trait + phylo axes
  mf_input <- cbind(trait_simulated, phylo_matrix)
  imputed_all <- missForest(mf_input, verbose = FALSE)$ximp
  
  # Only take the trait part
  imputed_mf <- imputed_all[, trait_names]
  imputed_mf[, log_transform_cols] <- exp(imputed_mf[, log_transform_cols])
  
  true_matrix <- as.matrix(trait_matrix_full)
  
  for (j in 1:n_traits) {
    missing_rows <- sampled_rowcol$row[sampled_rowcol$col == j]
    
    if (length(missing_rows) > 1) {
      predicted <- imputed_mf[missing_rows, j]
      observed  <- true_matrix[missing_rows, j]
      
      if (trait_names[j] == "Parental.care") {
        predicted <- round(as.numeric(predicted))
        predicted[predicted < 1] <- 1
        predicted[predicted > 3] <- 3
      }
      
      spearman_missForest_phylo[i, j] <- suppressWarnings(cor(
        as.numeric(predicted),
        as.numeric(observed),
        method = "spearman",
        use = "complete.obs"
      ))
    }
  }
}

summary_mf_phylo <- data.frame(
  Trait = colnames(spearman_missForest_phylo),
  Method = "missForest_phylogeny",
  Mean_Spearman = apply(spearman_missForest_phylo, 2, mean, na.rm = TRUE),
  SD_Spearman = apply(spearman_missForest_phylo, 2, sd, na.rm = TRUE)
)

print(summary_mf_phylo)

mean_spearman_across_traits <- mean(summary_mf_phylo$Mean_Spearman, na.rm = TRUE)
cat("Mean Spearman correlation across traits:", mean_spearman_across_traits, "\n")

write.csv(summary_mf_phylo, "spearman_summary_missForest_with_phylogeny-2.csv", row.names = FALSE)




##################################################################################Evaluate the predictive performance of Rphylopars
library(ape)
library(phytools)
library(dplyr)
library(Rphylopars)
library(tibble)

# Load data
trait_data <- read.csv("Final trait data of all global freshwater fish-for prediction-Rphylopars-2.csv", header = TRUE)
new_tree <- read.newick("fish_phylogenetic_tree.tre")

# Trim tree
tree_trimmed <- drop.tip(new_tree, setdiff(new_tree$tip.label, trait_data$species))

# Only retain the species on the tree
trait_data_filtered <- trait_data[trait_data$species %in% tree_trimmed$tip.label, ]
rownames(trait_data_filtered) <- trait_data_filtered$species

# natural log
log_trait_data <- trait_data_filtered

### log transform: exclude Parental care 
trait_names_all <- colnames(log_trait_data)[-1]
log_cols <- setdiff(trait_names_all, "Parental.care")
log_trait_data[, log_cols] <- log(log_trait_data[, log_cols])

# Only select the complete subset without NA
complete_cases <- log_trait_data %>% filter(if_all(-species, ~ !is.na(.)))
rownames(complete_cases) <- complete_cases$species

# Trim tree
tree_complete <- drop.tip(tree_trimmed, setdiff(tree_trimmed$tip.label, complete_cases$species))

# Set the number of cycles
n_iter <- 100
spearman_results <- matrix(NA, nrow = n_iter, ncol = ncol(complete_cases) - 1)
colnames(spearman_results) <- colnames(complete_cases)[-1]

# AIC and BIC
aic_bic_results <- matrix(NA, nrow = n_iter, ncol = 2)
colnames(aic_bic_results) <- c("AIC", "BIC")

# No trait is excluded
excluded_traits <- NULL
trait_names <- colnames(complete_cases)[-1]
target_trait_cols <- seq_along(trait_names)

set.seed(123)
for (i in 1:n_iter) {
  # Create a 25% deficiency
  simulated_data <- complete_cases
  trait_matrix <- as.matrix(simulated_data[, -1])
  n_total <- nrow(trait_matrix) * length(target_trait_cols)
  n_missing <- round(0.25 * n_total)
  
  possible_positions <- expand.grid(
    row = 1:nrow(trait_matrix),
    col = target_trait_cols
  )
  sampled_pos <- possible_positions[sample(1:nrow(possible_positions), n_missing, replace = FALSE), ]
  missing_indices <- as.matrix(sampled_pos)
  
  for (k in 1:nrow(missing_indices)) {
    trait_matrix[missing_indices[k,1], missing_indices[k,2]] <- NA
  }
  simulated_data[, -1] <- trait_matrix
  # Imputation
  imputed <- phylopars(trait_data = simulated_data, tree = tree_complete, model = "BM")
  imputed_means <- imputed$anc_recon
  
  # AIC and BIC
  aic_bic_results[i, "AIC"] <- AIC(imputed)
  aic_bic_results[i, "BIC"] <- BIC(imputed)
  
  true_values <- as.matrix(complete_cases[, -1])
  imputed_means <- imputed_means[rownames(true_values), colnames(true_values)]
  
  imputed_means[, log_cols] <- exp(imputed_means[, log_cols])
  true_values[, log_cols]   <- exp(true_values[, log_cols])
  
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
# Summary
spearman_summary <- data.frame(
  Trait = colnames(spearman_results),
  Mean_Spearman = apply(spearman_results, 2, mean, na.rm = TRUE),
  SD_Spearman = apply(spearman_results, 2, sd, na.rm = TRUE)
)

# Output result
print(spearman_summary)
print(aic_bic_results)
mean_spearman_across_traits <- mean(spearman_summary$Mean_Spearman, na.rm = TRUE)
cat("Mean Spearman correlation across traits:", mean_spearman_across_traits, "\n")

write.csv(spearman_summary, "spearman_summary-BM-2.csv", row.names = FALSE)


########################################################################################################Evaluate the predictive performance of the genus-mean imputation approach
library(dplyr)

# Load data
trait_data <- read.csv("Final trait data of all global freshwater fish-for prediction-2.csv", header = TRUE)

# Read the list of missing species in the phylogenetic tree
missing_species <- read.csv("missing species in phylogenetic tree-2.csv", header = TRUE)

# Delete the missing species in trait_data
missing_species$x <- gsub("_", " ", missing_species$x)
trait_data_filtered <- trait_data[!trait_data$species %in% missing_species$x, ]

# Screen complete species
complete_cases <- trait_data_filtered %>%
  filter(if_all(6:ncol(.), ~ !is.na(.)))
rownames(complete_cases) <- complete_cases$species

# trait matrix
trait_cols <- 6:ncol(complete_cases)
trait_names <- names(complete_cases)[trait_cols]

### No trait is excluded
excluded_traits <- NULL
included_trait_idx <- seq_along(trait_names)
included_cols <- trait_cols[included_trait_idx]

# Genus, family and order
genus_vector <- complete_cases$Genus
family_vector <- complete_cases$Family
order_vector  <- complete_cases$Order

# Initialize the result matrix
n_iter <- 100
n_traits <- length(trait_cols)
spearman_genusMean <- matrix(NA, nrow = n_iter, ncol = n_traits)
colnames(spearman_genusMean) <- trait_names

set.seed(123)

for (i in 1:n_iter) {
  cat("Simulation", i, "\n")
 
   #Create a 25% deficiency
  trait_simulated <- as.matrix(complete_cases[, trait_cols])
  n_total <- length(included_cols) * nrow(trait_simulated)
  n_missing <- round(0.25 * n_total)
  
  possible_positions <- expand.grid(row = 1:nrow(trait_simulated), col = included_cols)
  sampled_rowcol <- possible_positions[sample(nrow(possible_positions), n_missing, replace = FALSE), ]
  
  for (k in 1:n_missing) {
    col_idx <- sampled_rowcol$col[k] - 5
    trait_simulated[sampled_rowcol$row[k], col_idx] <- NA
  }
  
  # Imputation
  imputed_gm <- trait_simulated
  
  for (j in included_cols) {
    j_trait_idx <- j - 5
    na_rows <- which(is.na(imputed_gm[, j_trait_idx]))
    
    for (r in na_rows) {
      genus_r  <- genus_vector[r]
      family_r <- family_vector[r]
      order_r  <- order_vector[r]
      
      # 1) genus
      genus_rows <- which(genus_vector == genus_r & !is.na(imputed_gm[, j_trait_idx]))
      if (length(genus_rows) > 0) {
        imputed_gm[r, j_trait_idx] <- mean(imputed_gm[genus_rows, j_trait_idx], na.rm = TRUE)
      } else {
        
        ### family
        family_rows <- which(family_vector == family_r & !is.na(imputed_gm[, j_trait_idx]))
        if (length(family_rows) > 0) {
          imputed_gm[r, j_trait_idx] <- mean(imputed_gm[family_rows, j_trait_idx], na.rm = TRUE)
        } else {
          
          ### order
          order_rows <- which(order_vector == order_r & !is.na(imputed_gm[, j_trait_idx]))
          if (length(order_rows) > 0) {
            imputed_gm[r, j_trait_idx] <- mean(imputed_gm[order_rows, j_trait_idx], na.rm = TRUE)
          }
          
        }
        
      }
    }
  }
  
  true_matrix <- as.matrix(complete_cases[, trait_cols])
  
  for (j in 1:n_traits) {
    missing_rows <- sampled_rowcol$row[sampled_rowcol$col == trait_cols[j]]
    if (length(missing_rows) > 1) {
      
      predicted <- imputed_gm[missing_rows, j]
      observed  <- true_matrix[missing_rows, j]
      
      trait_name <- trait_names[j]
      if (trait_name == "Parental.care") {
        predicted <- round(predicted)
        predicted[predicted < 1] <- 1
        predicted[predicted > 3] <- 3
      }
      
      spearman_genusMean[i, j] <- suppressWarnings(cor(
        predicted,
        observed,
        method = "spearman",
        use = "complete.obs"
      ))
    }
  }
}

# Summary
summary_gm <- data.frame(
  Trait = colnames(spearman_genusMean),
  Method = "Genus_Family_Order_Mean",
  Mean_Spearman = apply(spearman_genusMean, 2, mean, na.rm = TRUE),
  SD_Spearman = apply(spearman_genusMean, 2, sd, na.rm = TRUE)
)

#Output result
print(summary_gm)
mean_spearman_across_traits <- mean(summary_gm$Mean_Spearman, na.rm = TRUE)
cat("Mean Spearman correlation across traits:", mean_spearman_across_traits, "\n")

write.csv(summary_gm, "spearman_summary_genusMean_only_nolog_filtered-2.csv", row.names = FALSE)


