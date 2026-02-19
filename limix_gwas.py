
# LIMIX GWAS - Chromosome 10 only
# Phenotype: microbe_abundance
# Fixed effects: sex, cohort
# Random effects: GRM (chr10 excluded), cage


import numpy as np
import pandas as pd
from pandas_plink import read_plink
from scipy.linalg import cho_factor, cho_solve
from numpy_sugar.linalg import economic_qs
from glimix_core.lmm import LMM
from scipy.stats import chi2 as chi2_dist

# 1.  LOAD METADATA / PHENOTYPE
meta = pd.read_csv("microbe_abundance.csv")
meta = meta.dropna(subset=["microbe_abundance"])  # drop samples with missing phenotype
meta = meta.reset_index(drop=True)

print(f"  Samples after dropping NaN phenotype: {len(meta)}")

# 2.  LOAD GRM (LOCO – chr10 excluded)
grm_id = pd.read_csv(
    "my_subset_no_chr10.grm.id",
    sep="\t", header=None, names=["fid", "iid"]
)

n_grm = len(grm_id)
# Read binary GRM matrix
grm_bin = np.fromfile("my_subset_no_chr10.grm.bin", dtype=np.float32)
# Lower-triangle packed format → full symmetric matrix
K_grm = np.zeros((n_grm, n_grm), dtype=np.float64)
idx = 0
for i in range(n_grm):
    for j in range(i + 1):
        K_grm[i, j] = grm_bin[idx]
        K_grm[j, i] = grm_bin[idx]
        idx += 1

# 3.  ALIGN SAMPLES (meta ↔ GRM)
print("Aligning samples between metadata and GRM...")
grm_id["rfid"] = grm_id["iid"].astype(str)
meta["rfid"] = meta["rfid"].astype(str)

# Keep only individuals present in BOTH meta and GRM
common_rfids = list(set(meta["rfid"]) & set(grm_id["rfid"]))
print(f"  Individuals in common: {len(common_rfids)}")

# Re-index
meta = meta[meta["rfid"].isin(common_rfids)].copy()
meta = meta.set_index("rfid")

grm_order = grm_id[grm_id["rfid"].isin(common_rfids)]["rfid"].tolist()
meta = meta.loc[grm_order].reset_index()  # align meta to GRM order

# Subset GRM to common individuals
grm_mask = grm_id["rfid"].isin(common_rfids).values
K_grm = K_grm[np.ix_(grm_mask, grm_mask)]

n = len(meta)
print(f"  Final sample size: {n}")

# 4.  BUILD COVARIATE MATRICES

sex_dummy = pd.get_dummies(meta["sex"], drop_first=True, prefix="sex").astype(float)
cohort_dummies = pd.get_dummies(meta["cohort"], drop_first=True, prefix="cohort").astype(float)

intercept = np.ones((n, 1))
W = np.hstack([intercept, sex_dummy.values, cohort_dummies.values])  # (n, p)
print(f"  Fixed-effect covariate matrix shape: {W.shape}")

# Random effect 2: Cage relatedness matrix 
# Build cage kinship: K_cage[i,j] = 1 if same cage, 0 otherwise
cage_arr = meta["cage"].values
K_cage = (cage_arr[:, None] == cage_arr[None, :]).astype(float)
# Normalise so diagonal is mean 1 (standard practice)
K_cage = K_cage / K_cage.diagonal().mean()
print(f"  Cage kinship matrix shape: {K_cage.shape}")


# 5.  PHENOTYPE VECTOR
y = meta["microbe_abundance"].values.astype(float).reshape(-1, 1)


# 6.  LOAD GENOTYPES PLINK FILE
#     Extracting only MY rats and only CHR 10

# read_plink returns (bim, fam, bed) — note different order from npplink
# bed is a dask array of shape (n_variants, n_samples)
(f, b, bed) = read_plink('/tscc/projects/ps-palmer/gwas/databases/rounds/r11.2.1')

# f = bim dataframe: columns chrom, snp, cm, pos, a0, a1, i
# b = fam dataframe: columns fid, iid, father, mother, gender, trait, i
# Filter SNPs to chromosome 10
chr10_mask = f["chrom"].astype(str) == "10"
f_chr10 = f[chr10_mask].copy()
print(f"  SNPs on chromosome 10: {len(f_chr10)}")

# Identify my rats in the shared file by rfid (matched to iid)
b["rfid"] = b["iid"].astype(str)
my_rfids = set(meta["rfid"].astype(str))
sample_mask = b["rfid"].isin(my_rfids)

# Get per-sample index mapping: shared file order → my sample order (aligned to GRM)
b_sub = b[sample_mask].copy()
b_sub = b_sub.set_index("rfid").loc[meta["rfid"].values].reset_index()
sample_indices = b_sub["i"].values  # use the 'i' column (variant/sample integer index)

print(f"  Matched samples in plink file: {len(sample_indices)}")


# 7.  FIT NULL MODEL WITH TWO RANDOM EFFECTS
#     limix.qtl.scan only accepts a single K matrix.
#     Solution: fit a null LMM with two variance components (GRM + cage)
#     using GLMM, extract the fitted variance parameters, then form a single
#     combined K = sigma2_grm * K_grm + sigma2_cage * K_cage + sigma2_e * I
#     and pass that combined covariance to qtl.scan.

print("Fitting null model with two random effects (GRM + Cage)...")

# QS decomposition of the combined covariance (sum of the two Ks)
# We start with equal weighting; the null LMM will re-estimate.
K_init = K_grm + K_cage   # combined kinship as starting point for QS decomp

QS = economic_qs(K_init)

# Fit null LMM (no SNP term) to estimate variance components
null_lmm = LMM(y.ravel(), W, QS)
null_lmm.fit(verbose=True)

# Extract fitted variance components
v_g = null_lmm.v0   # variance explained by K_init (GRM + cage combined)
v_e = null_lmm.v1   # residual variance

print(f"  Null model fitted: v_genetic={v_g:.4f}, v_residual={v_e:.4f}")

# Form the effective single K for the scan:
# K_eff = v_g  K_init + v_e  I  (already implicitly in the QS structure)
# We pass K_init to limix.qtl.scan which will re-estimate internally,
# but this is consistent with how the null model was parameterised.
K_combined = K_grm + K_cage   # single effective K

# 8.  RUN LIMIX LMM GWAS – CHR 10

snp_indices = f_chr10["i"].values  # integer indices for bed accessor

# bed in pandas_plink is (n_variants, n_samples) as a dask array
# Subset to chr10 SNPs and my samples, then compute to numpy
G_chr10 = bed[snp_indices, :][:, sample_indices].compute().T  # → (n_samples, n_snps)
G_chr10 = G_chr10.astype(float)

# Mean-impute missing genotypes (NaN → column mean)
col_means = np.nanmean(G_chr10, axis=0)
nan_locs = np.isnan(G_chr10)
G_chr10[nan_locs] = np.take(col_means, np.where(nan_locs)[1])

print(f"  Genotype matrix shape: {G_chr10.shape}")

# 8.  GWAS SCAN – CHR 10
#     limix.qtl.scan is broken with pandas>=2.0.
#     Instead we use glimix_core.LMM directly:
#       - Null model (already fit above): W as covariates, K_combined as K
#       - Alt model per SNP: [W | g] as covariates, same K
#       - LRT p-value = chi2(1) on 2*(logL_alt - logL_null)



n_snps = G_chr10.shape[1]
pvals  = np.zeros(n_snps)
betas  = np.zeros(n_snps)

null_logL = null_lmm.lml()   # log marginal likelihood of null model
print(f"  Null log-likelihood: {null_logL:.4f}")
print(f"  Scanning {n_snps} SNPs on chr10...")

for j in range(n_snps):
    if j % 5000 == 0:
        print(f"    SNP {j}/{n_snps}")

    g = G_chr10[:, j]

    # Skip monomorphic SNPs
    if np.std(g) == 0:
        pvals[j] = 1.0
        betas[j] = 0.0
        continue

    # Standardise genotype (mean 0, std 1) — common practice
    g = (g - g.mean()) / g.std()

    # Build alt covariate matrix: [W | g]
    M_alt = np.hstack([W, g.reshape(-1, 1)])

    # Fit alt LMM with same QS decomposition
    alt_lmm = LMM(y.ravel(), M_alt, QS)
    alt_lmm.fit(verbose=False)

    # LRT statistic
    lrt = 2.0 * (alt_lmm.lml() - null_logL)
    lrt = max(lrt, 0.0)   # numerical safety

    pvals[j] = chi2_dist.sf(lrt, df=1)
    # beta = coefficient of the last covariate (the SNP)
    betas[j] = alt_lmm.beta[-1]


# 9.  EXTRACT AND SAVE RESULTS


results_df = pd.DataFrame({
    "snp":     f_chr10["snp"].values,
    "chrom":   f_chr10["chrom"].values,
    "pos":     f_chr10["pos"].values,
    "a0":      f_chr10["a0"].values,
    "a1":      f_chr10["a1"].values,
    "beta":    betas,
    "pval":    pvals,
    "-log10p": -np.log10(np.clip(pvals, 1e-300, 1.0))
})

results_df = results_df.sort_values("pval").reset_index(drop=True)

out_file = "gwas_chr10_results.csv"
results_df.to_csv(out_file, index=False)
print(f"\nDone! Results saved to: {out_file}")
print("Top hits:")
print(results_df.head(10).to_string(index=False))
