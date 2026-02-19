library(sommer)
library(dplyr)

# 1. LOAD DATA

data <- read.csv("microbe_abundance.csv", stringsAsFactors = FALSE)
data$sex <- as.factor(data$sex)
data$cage <- as.factor(data$cage)
data$cohort <- as.factor(data$cohort)
data$genotype <- as.factor(data$genotype)
data$rfid_char <- as.character(data$rfid)

# Set reference genotype
genotype_levels <- levels(data$genotype)
ref_patterns <- c("Homo_REF", "Homo_Ref", "homo_ref", "REF", "Ref", "ref")
ref_level <- genotype_levels[1]
for (pattern in ref_patterns) {
  if (pattern %in% genotype_levels) {
    ref_level <- pattern
    break
  }
}
data$genotype <- relevel(data$genotype, ref = ref_level)
cat("Reference genotype:", ref_level, "\n")

# 2. LOAD GRM 

use_grm <- FALSE
tryCatch({
  read_gcta_grm <- function(grm_prefix) {
    grm_id <- read.table(paste0(grm_prefix, ".grm.id"), header = FALSE)
    colnames(grm_id) <- c("FID", "IID")
    n_ids <- nrow(grm_id)
    
    grm_bin <- file(paste0(grm_prefix, ".grm.bin"), "rb")
    grm_values <- readBin(grm_bin, what = "numeric", 
                          n = n_ids * (n_ids + 1) / 2, size = 4)
    close(grm_bin)
    
    grm_matrix <- matrix(0, n_ids, n_ids)
    k <- 1
    for (i in 1:n_ids) {
      for (j in 1:i) {
        grm_matrix[i, j] <- grm_values[k]
        grm_matrix[j, i] <- grm_values[k]
        k <- k + 1
      }
    }
    
    rownames(grm_matrix) <- as.character(grm_id$IID)
    colnames(grm_matrix) <- as.character(grm_id$IID)
    return(list(grm = grm_matrix, ids = grm_id))
  }
  
  grm_data <- read_gcta_grm("my_subset_no_chr10")
  grm_matrix <- grm_data$grm
  common_ids <- intersect(data$rfid_char, as.character(grm_data$ids$IID))
  data <- data[data$rfid_char %in% common_ids, ]
  grm_matched <- grm_matrix[common_ids, common_ids]
  use_grm <- TRUE
  cat("✓ GRM loaded successfully\n")
}, error = function(e) {
  cat("! Using simple random effects (no GRM)\n")
})

# 3. FIT MODELS

if (use_grm) {
  model_full <- mmer(microbe_abundance ~ genotype + sex + cohort,
                     random = ~ cage + vs(rfid_char, Gu = grm_matched),
                     data = data,
                     verbose = FALSE)
} else {
  model_full <- mmer(microbe_abundance ~ genotype + sex + cohort,
                     random = ~ cage + rfid_char,
                     data = data,
                     verbose = FALSE)
}
cat("Full model AIC:", model_full$AIC, "\n")

cat("\n--- Fitting NULL model (WITHOUT genotype) ---\n")
if (use_grm) {
  model_null <- mmer(microbe_abundance ~ sex + cohort,
                     random = ~ cage + vs(rfid_char, Gu = grm_matched),
                     data = data,
                     verbose = FALSE)
} else {
  model_null <- mmer(microbe_abundance ~ sex + cohort,
                     random = ~ cage + rfid_char,
                     data = data,
                     verbose = FALSE)
}
cat("Null model AIC:", model_null$AIC, "\n")

# 4. LIKELIHOOD RATIO TEST (from AIC)

# AIC = -2*logLik + 2*k, where k is number of parameters
# So: logLik = (AIC - 2*k) / -2

# Count parameters
n_params_full <- length(model_full$Beta) + length(model_full$sigma)
n_params_null <- length(model_null$Beta) + length(model_null$sigma)

# Calculate log-likelihoods
loglik_full <- (model_full$AIC - 2*n_params_full) / -2
loglik_null <- (model_null$AIC - 2*n_params_null) / -2

# Calculate LR statistic
lr_stat <- 2 * (loglik_full - loglik_null)

# Degrees of freedom
n_genotype_levels <- length(unique(data$genotype))
df <- n_genotype_levels - 1

# Calculate p-value
p_value <- 1 - pchisq(lr_stat, df)

# 5. RESULTS

cat("\nQUESTION: Is there a significant difference in microbe abundance\n")
cat("          between the three genotypes (controlling for sex, cohort,\n")
cat("          cage, and", ifelse(use_grm, "GRM", "animal ID"), ")?\n")
cat("\n")
cat("LIKELIHOOD RATIO TEST:\n")
cat("Chi-square statistic:", round(lr_stat, 4), "\n")
cat("Degrees of freedom:", df, "\n")
cat("P-VALUE:", formatC(p_value, format = "e", digits = 4), "\n")
cat("\n")


# 6. PAIRWISE COMPARISONS

cat("PAIRWISE COMPARISONS\n")
cat("(Which specific genotypes differ from the reference?)\n\n")

fixed_effects <- model_full$Beta
genotype_effects <- fixed_effects[grep("genotype", rownames(fixed_effects)), , drop = FALSE]

if (nrow(genotype_effects) > 0) {
  pairwise_results <- data.frame(
    Comparison = rownames(genotype_effects),
    Estimate = genotype_effects[, "Estimate"],
    Std_Error = genotype_effects[, "Std.Error"],
    t_value = genotype_effects[, "Estimate"] / genotype_effects[, "Std.Error"]
  )
  
  # Calculate approximate p-values from t-statistics
  pairwise_results$P_value_approx <- 2 * (1 - pnorm(abs(pairwise_results$t_value)))
  pairwise_results$Significant <- ifelse(pairwise_results$P_value_approx < 0.05, "YES", "NO")
  
  print(pairwise_results, row.names = FALSE)
  
  cat("\nINTERPRETATION:\n")
  for (i in 1:nrow(pairwise_results)) {
    cat("\n", pairwise_results$Comparison[i], ":\n", sep="")
    cat("  Difference from reference:", round(pairwise_results$Estimate[i], 2), "units\n")
    if (pairwise_results$Significant[i] == "YES") {
      cat("  → SIGNIFICANT (p < 0.05)\n")
    } else {
      cat("  → NOT significant (p >= 0.05)\n")
    }
  }
}

# 7. ADDITIONAL INFO: AIC Comparison

cat("MODEL COMPARISON (AIC)\n")

delta_aic <- model_full$AIC - model_null$AIC
cat("\nFull model AIC:", round(model_full$AIC, 2), "\n")
cat("Null model AIC:", round(model_null$AIC, 2), "\n")
cat("Difference (Δ AIC):", round(delta_aic, 2), "\n\n")

if (delta_aic < -2) {
  cat("→ Full model is BETTER (lower AIC)\n")
  cat("  Including genotype improves model fit\n")
} else if (delta_aic > 2) {
  cat("→ Null model is BETTER (lower AIC)\n")
  cat("  Including genotype does NOT improve model fit\n")
} else {
  cat("→ Models are EQUIVALENT (similar AIC)\n")
  cat("  Including genotype provides minimal improvement\n")
}

cat("\nNote: ΔAIC < -2 suggests the full model is substantially better\n")

# 8. SAVE RESULTS

# Overall test
lrt_results <- data.frame(
  Test = "Overall Genotype Effect",
  Chi_square = lr_stat,
  DF = df,
  P_value = p_value,
  AIC_full = model_full$AIC,
  AIC_null = model_null$AIC,
  Delta_AIC = delta_aic,
  Significant = ifelse(p_value < 0.05, "YES", "NO")
)
write.csv(lrt_results, "genotype_test_results.csv", row.names = FALSE)
cat("✓ Saved: genotype_test_results.csv\n")

# Pairwise comparisons
if (nrow(genotype_effects) > 0) {
  write.csv(pairwise_results, "pairwise_genotype_comparisons.csv", row.names = FALSE)
  cat("✓ Saved: pairwise_genotype_comparisons.csv\n")
}

# Genotype means
genotype_means <- data %>%
  group_by(genotype) %>%
  summarise(
    n = n(),
    mean = mean(microbe_abundance, na.rm = TRUE),
    sd = sd(microbe_abundance, na.rm = TRUE),
    median = median(microbe_abundance, na.rm = TRUE),
    min = min(microbe_abundance, na.rm = TRUE),
    max = max(microbe_abundance, na.rm = TRUE)
  )
write.csv(genotype_means, "genotype_summary_statistics.csv", row.names = FALSE)
