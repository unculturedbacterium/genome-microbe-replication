# genome-microbe-replication
Code for analyzing host genetic effects on the gut microbiome in outbred rats. Includes LIMIX and GCTA mixed-model GWAS, GRM construction, alpha and beta diversity analyses, differential abundance testing, and publication-ready plotting scripts for microbiome and genetic association results.


"limix_gwas".py - This script performs a chromosome 10 specific GWAS of microbial abundance using a linear mixed model framework implemented with LIMIX/glimix. It aligns phenotype metadata, genotypes, and a GRM (excluding chr10), while accounting for fixed effects (sex, cohort) and random effects (genetic relatedness and cage). A null mixed model is first fit to estimate variance components, after which each SNP on chromosome 10 is tested using a likelihood ratio test comparing models with and without the SNP effect. The pipeline outputs SNP effect sizes and p-values for downstream interpretation and visualization.

"gwas_plots.py" - This script visualizes chromosome 10 GWAS results for microbial abundance. It loads SNP-level association statistics, identifies the top-associated variant, and computes Bonferroni and suggestive significance thresholds.

"sommer.R" - This R script tests whether host genotype is associated with microbial abundance using linear mixed models. It fits a full model including genotype, sex, and cohort, and a null model excluding genotype, while accounting for cage effects and either genetic relatedness via a GRM or individual ID as random effects. The models are compared using a likelihood ratio test to assess the overall genotype effect. The script also performs pairwise genotype comparisons, evaluates model fit using AIC.
