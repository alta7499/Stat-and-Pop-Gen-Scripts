import allel
import numpy as np
import pandas as pd
from scipy.stats import binom_test
import matplotlib
matplotlib.use("Agg")  # ✅ headless-safe backend
import matplotlib.pyplot as plt
import os

# ---------- Load population information ----------
sample_info = pd.read_csv("sample_info_v3", sep="\t", header=None, names=["sample", "population"])
sample_to_pop = dict(zip(sample_info['sample'], sample_info['population']))

# ---------- Load VCF ----------
callset = allel.read_vcf(
    "Malay_Negrito_200s.vcf",
    fields=['samples', 'variants/CHROM', 'variants/POS', 'calldata/GT'],
    alt_number=1
)
samples = callset['samples']
positions = np.array(callset['variants/POS'])
chroms = np.array(callset['variants/CHROM'])
genotypes = allel.GenotypeArray(callset['calldata/GT'])

# ---------- Define populations ----------
pop1 = [i for i, s in enumerate(samples) if sample_to_pop.get(s) == 'Negrito']
pop2 = [i for i, s in enumerate(samples) if sample_to_pop.get(s) == 'Malay']

print(f"Negrito samples: {len(pop1)}")
print(f"Malay samples: {len(pop2)}")

# ---------- Optional: check unmatched samples ----------
missing = [s for s in samples if s not in sample_to_pop]
if missing:
    print("Warning: some VCF samples not found in Sample_v3.txt:", missing)

# ---------- Calculate Fst per SNP ----------
fst_all = []
for i in range(genotypes.shape[0]):
    gt = genotypes.take([i], axis=0)
    ac1 = gt[:, pop1].count_alleles()
    ac2 = gt[:, pop2].count_alleles()
    num, den = allel.hudson_fst(ac1, ac2)
    fst_val = (num / den)[0] if den > 0 else np.nan
    fst_all.append(fst_val)

fst_all = np.array(fst_all)

# ---------- Clean arrays ----------
assert len(fst_all) == len(positions) == len(chroms), "Array length mismatch"

valid_mask = ~np.isnan(fst_all)
fst_clean = fst_all[valid_mask]
positions_clean = positions[valid_mask]
chroms_clean = chroms[valid_mask]

# ---------- Compute empirical percentiles ----------
percentiles = 1 - (np.argsort(np.argsort(fst_clean)) / len(fst_clean))

# ---------- Perform binomial test per 100kb window ----------
window_size = 100_000
results = []

for chrom in np.unique(chroms_clean):
    mask_chrom = chroms_clean == chrom
    pos_chrom = positions_clean[mask_chrom]
    perc_chrom = percentiles[mask_chrom]

    for start in range(pos_chrom.min(), pos_chrom.max(), window_size):
        end = start + window_size
        mask = (pos_chrom >= start) & (pos_chrom < end)
        n_total = np.sum(mask)
        n_high = np.sum(perc_chrom[mask] < 0.01)  # top 1%
        if n_total == 0:
            pval = np.nan
        else:
            pval = binom_test(n_high, n_total, p=0.01, alternative='greater')
        results.append((chrom, start, end, n_total, n_high, pval))

# ---------- Save results ----------
df = pd.DataFrame(results, columns=['chrom', 'start', 'end', 'n_snps', 'n_top1pct', 'binomial_pval'])
df.to_csv("fst_enrichment_binomial_windows.csv", index=False)

# ---------- Create output folder for plots ----------
os.makedirs("fst_binomial_plots", exist_ok=True)

# ---------- Plot each chromosome ----------
for chrom in df['chrom'].unique():
    df_chr = df[df['chrom'] == chrom]
    plt.figure()
    plt.plot(df_chr['start'], -np.log10(df_chr['binomial_pval']), marker='o')
    plt.axhline(-np.log10(0.05), color='gray', linestyle='--')
    plt.title(f"Binomial test for high Fst SNP enrichment ({chrom})")
    plt.xlabel("Genomic Position")
    plt.ylabel("-log10(P-value)")
    plt.tight_layout()
    plt.savefig(f"fst_binomial_plots/fst_binomial_{chrom}.png", dpi=300)
    plt.close()

print("✅ Analysis complete.")
print("📄 Results saved to: fst_enrichment_binomial_windows.csv")
print("📊 Plots saved in:   fst_binomial_plots/")

