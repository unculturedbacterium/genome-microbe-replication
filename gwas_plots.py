import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import Normalize


RESULTS_FILE = "gwas_chr10_results.csv"
LOCUS_WINDOW = 1_000_000   # ±1 Mb around top SNP
SUGGESTIVE_P = 1e-5
GENOME_WIDE_P_BONF = None   # None → auto Bonferroni (0.05 / n_snps)

# 1.  LOAD RESULTS
res = pd.read_csv(RESULTS_FILE)
res = res.dropna(subset=["pval"])
res["logp"] = -np.log10(np.clip(res["pval"], 1e-300, 1.0))
res = res.sort_values("pos").reset_index(drop=True)

top_snp   = res.loc[res["pval"].idxmin()]
top_pos   = int(top_snp["pos"])
bonf_p    = GENOME_WIDE_P_BONF or (0.05 / len(res))
bonf_logp = -np.log10(bonf_p)
sug_logp  = -np.log10(SUGGESTIVE_P)

# 2.  MANHATTAN PLOT
fig, ax = plt.subplots(figsize=(14, 5))
fig.patch.set_facecolor("#0d0d1a")
ax.set_facecolor("#0d0d1a")

# Background SNPs
ax.scatter(res["pos"], res["logp"],
           s=6, c="#2ec4b6", alpha=0.65, linewidths=0, rasterized=True)

# Suggestive SNPs (above suggestive but below Bonferroni)
sug_mask = (res["logp"] >= sug_logp) & (res["logp"] < bonf_logp)
ax.scatter(res.loc[sug_mask, "pos"], res.loc[sug_mask, "logp"],
           s=16, c="#ffbe0b", alpha=1.0, linewidths=0, zorder=5)

# Bonferroni-significant SNPs
sig_mask = res["logp"] >= bonf_logp
ax.scatter(res.loc[sig_mask, "pos"], res.loc[sig_mask, "logp"],
           s=20, c="#ff006e", alpha=1.0, linewidths=0, zorder=6)

# Top SNP gold star
ax.scatter(top_snp["pos"], top_snp["logp"],
           s=120, c="#ffbe0b", marker="*",
           edgecolors="white", linewidths=0.5, zorder=10,
           label=f"Lead SNP: {top_snp['snp']}")

# Threshold lines
ax.axhline(bonf_logp, color="#ff006e", linewidth=0.9, linestyle="--", alpha=0.85,
           label=f"Bonferroni (p={bonf_p:.1e})")
ax.axhline(sug_logp,  color="#ffbe0b", linewidth=0.8, linestyle=":",  alpha=0.65,
           label=f"Suggestive (p={SUGGESTIVE_P:.0e})")

# Styling
ax.set_xlabel("Chromosome 10 position", color="#c9d1d9", fontsize=11, labelpad=8)
ax.set_ylabel("−log₁₀(p)", color="#c9d1d9", fontsize=11, labelpad=8)
ax.set_title("GWAS Manhattan Plot · Chromosome 10 · Microbe Abundance",
             color="white", fontsize=13, fontweight="bold", pad=12)
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x/1e6:.0f} Mb"))
ax.set_xlim(res["pos"].min() - 0.5e6, res["pos"].max() + 0.5e6)
ax.set_ylim(0, max(res["logp"].max() * 1.18, bonf_logp * 1.3))
ax.tick_params(colors="#8b949e", labelsize=9)
for spine in ax.spines.values():
    spine.set_edgecolor("#30363d")
ax.legend(facecolor="#161b22", edgecolor="#30363d",
          labelcolor="white", fontsize=8.5, loc="upper left", markerscale=1.3)

plt.tight_layout()
plt.savefig("manhattan_plot.png", dpi=180, bbox_inches="tight",
            facecolor=fig.get_facecolor())
plt.close()

# 3.  SUBSET LOCUS  (top SNP ± 1 Mb)

win_lo = top_pos - LOCUS_WINDOW
win_hi = top_pos + LOCUS_WINDOW
locus_res = res[(res["pos"] >= win_lo) & (res["pos"] <= win_hi)].copy()
print(f"\n  SNPs in locus window (±1 Mb): {len(locus_res)}")

# 4.  LOCUS ZOOM PLOT
#     Points colored by −log10(p) via plasma colormap.
#     Top SNP = gold diamond with arrow annotation.


logp_max = locus_res["logp"].max()
cmap     = plt.cm.get_cmap("plasma")
norm     = Normalize(vmin=0, vmax=logp_max)

fig, (ax_main, ax_chr) = plt.subplots(
    2, 1, figsize=(13, 7),
    gridspec_kw={"height_ratios": [6, 1]},
    sharex=True
)
fig.patch.set_facecolor("#0d0d1a")
for ax in (ax_main, ax_chr):
    ax.set_facecolor("#0d0d1a")

# Scatter colored by -log10(p)
sc = ax_main.scatter(
    locus_res["pos"] / 1e6,
    locus_res["logp"],
    c=locus_res["logp"],
    cmap=cmap,
    norm=norm,
    s=30,
    alpha=0.90,
    linewidths=0.3,
    edgecolors="#ffffff20",
    zorder=4
)

# Top SNP diamond
ax_main.scatter(
    top_pos / 1e6, top_snp["logp"],
    s=180, marker="D", c="#ffbe0b",
    edgecolors="white", linewidths=1.2,
    zorder=10, label=f"Lead SNP: {top_snp['snp']}"
)

# Threshold lines
ax_main.axhline(bonf_logp, color="#ff006e", linewidth=0.9,
                linestyle="--", alpha=0.85, label=f"Bonferroni (p={bonf_p:.1e})")
ax_main.axhline(sug_logp,  color="#aaaaaa", linewidth=0.8,
                linestyle=":",  alpha=0.60, label=f"Suggestive (p={SUGGESTIVE_P:.0e})")

# Colorbar
cbar = fig.colorbar(sc, ax=ax_main, pad=0.01, fraction=0.018)
cbar.set_label("−log₁₀(p)", color="#c9d1d9", fontsize=9)
cbar.ax.yaxis.set_tick_params(color="#8b949e", labelcolor="#8b949e", labelsize=8)
cbar.outline.set_edgecolor("#30363d")

# Annotate top SNP
offset_x = LOCUS_WINDOW * 0.28 / 1e6
ax_main.annotate(
    f"{top_snp['snp']}\np = {top_snp['pval']:.2e}",
    xy=(top_pos / 1e6, top_snp["logp"]),
    xytext=(top_pos / 1e6 + offset_x, top_snp["logp"] * 0.88),
    color="#ffbe0b", fontsize=8.5, fontweight="bold",
    ha="left", va="top",
    arrowprops=dict(arrowstyle="-|>", color="#ffbe0b", lw=0.9)
)

# Axis styling
ax_main.set_ylabel("−log₁₀(p)", color="#c9d1d9", fontsize=11, labelpad=8)
ax_main.set_title(
    f"Locus Zoom · Chr10 ±1 Mb around {top_snp['snp']}  "
    f"({top_pos/1e6:.2f} Mb)\nMicrobe Abundance GWAS",
    color="white", fontsize=12, fontweight="bold", pad=10
)
ax_main.set_ylim(0, logp_max * 1.22)
ax_main.tick_params(colors="#8b949e", labelsize=9)
for spine in ax_main.spines.values():
    spine.set_edgecolor("#30363d")
ax_main.legend(facecolor="#161b22", edgecolor="#30363d",
               labelcolor="white", fontsize=8.5, loc="upper left")

# Chromosome position track
ax_chr.barh(0.5, (win_hi - win_lo) / 1e6, left=win_lo / 1e6,
            height=0.35, color="#3a86ff", alpha=0.35)
ax_chr.axvline(top_pos / 1e6, color="#ffbe0b", linewidth=1.3,
               linestyle="--", alpha=0.9)
ax_chr.text(top_pos / 1e6, 0.05,
            f"{top_pos/1e6:.3f} Mb",
            color="#ffbe0b", fontsize=7.5, ha="center", va="bottom")
ax_chr.set_ylim(0, 1)
ax_chr.set_yticks([])
ax_chr.set_ylabel("Chr 10", color="#8b949e", fontsize=8,
                  rotation=0, labelpad=32, va="center")
ax_chr.set_xlabel("Position (Mb)", color="#c9d1d9", fontsize=11, labelpad=8)
ax_chr.tick_params(colors="#8b949e", labelsize=9)
for spine in ax_chr.spines.values():
    spine.set_edgecolor("#30363d")

plt.tight_layout(h_pad=0.2)
plt.savefig("locus_zoom_chr10.png", dpi=180, bbox_inches="tight",
            facecolor=fig.get_facecolor())
plt.close()
