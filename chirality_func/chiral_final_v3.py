#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Final analysis script for comparing semi-automatic vs automatic chirality scores.
Creates a timestamped output folder and saves all results inside it.
Ready for publication/PhD thesis.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from datetime import datetime

# ---------------------------------------------------------
# CREATE TIMESTAMPED OUTPUT DIRECTORY
# ---------------------------------------------------------

timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
output_dir = f"chiral_results_{timestamp}"

os.makedirs(output_dir, exist_ok=True)

print(f"📁 Output directory created: {output_dir}")

# ---------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------

input_file  = "salida.xlsx"
col_id    = "CMPD_CHEMBLID"
col_semia = "f(Q1/Q0)inchi"
col_auto  = "f_chiral"

THRESH_LARGE = 0.10  # absolute difference cutoff

# ---------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------

print(f"Reading file: {input_file}")
df = pd.read_excel(input_file)

# Keep only relevant columns
keep_cols = [col_id, col_semia, col_auto]
if "observations" in df.columns:
    keep_cols.append("observations")

df = df[keep_cols].copy()
df = df.drop_duplicates(subset=[col_id])

# Rename for clarity
df = df.rename(columns={
    col_semia: "f_chiral_semiautomatica",
    col_auto:  "f_chiral_automatica"
})

# Ensure numeric
df["f_chiral_semiautomatica"] = pd.to_numeric(df["f_chiral_semiautomatica"], errors="coerce")
df["f_chiral_automatica"]     = pd.to_numeric(df["f_chiral_automatica"], errors="coerce")

df = df.dropna(subset=["f_chiral_semiautomatica", "f_chiral_automatica"], how="all")

# ---------------------------------------------------------
# DIFFERENCES
# ---------------------------------------------------------

df["diff_abs"] = df["f_chiral_semiautomatica"] - df["f_chiral_automatica"]
df["diff_rel"] = df["diff_abs"] / df["f_chiral_semiautomatica"].replace(0, np.nan)

df["flag_gran_diferencia"] = df["diff_abs"].abs() > THRESH_LARGE
df["flag_quiralidad_latente"] = (
    (df["f_chiral_semiautomatica"] > 0) &
    (df["f_chiral_automatica"] == 0)
)

# ---------------------------------------------------------
# STATISTICS
# ---------------------------------------------------------

corr_global = df["f_chiral_semiautomatica"].corr(df["f_chiral_automatica"])

n_total = len(df)
n_equal = (df["diff_abs"].abs() < 1e-9).sum()
n_diff_001 = (df["diff_abs"].abs() < 0.01).sum()
n_diff_005 = (df["diff_abs"].abs() < 0.05).sum()
n_large = df["flag_gran_diferencia"].sum()
n_latent = df["flag_quiralidad_latente"].sum()

pct_equal  = 100 * n_equal / n_total
pct_001    = 100 * n_diff_001 / n_total
pct_005    = 100 * n_diff_005 / n_total
pct_large  = 100 * n_large  / n_total
pct_latent = 100 * n_latent / n_total

desc_global = df[["f_chiral_semiautomatica", "f_chiral_automatica", "diff_abs"]].describe()

# ---------------------------------------------------------
# FILTERED DATASET (f > 0 in both)
# ---------------------------------------------------------

df_filt = df[(df["f_chiral_semiautomatica"] > 0) & (df["f_chiral_automatica"] > 0)].copy()

corr_filt = df_filt["f_chiral_semiautomatica"].corr(df_filt["f_chiral_automatica"])
desc_filt = df_filt[["f_chiral_semiautomatica", "f_chiral_automatica", "diff_abs"]].describe()

# ---------------------------------------------------------
# SAVE COMPARISON TABLE
# ---------------------------------------------------------

df.to_excel(os.path.join(output_dir, "comparacion_chiral.xlsx"), index=False)
print("✔ comparacion_chiral.xlsx saved.")

# ---------------------------------------------------------
# SAVE STATISTICS TO EXCEL
# ---------------------------------------------------------

summary_df = pd.DataFrame({
    "metric": [
        "n_total",
        "n_equal",
        "n_diff_abs<0.01",
        "n_diff_abs<0.05",
        "n_diff_abs>0.10",
        "n_latent_chirality",
        "corr_global",
        "corr_filtered",
        "n_filtered"
    ],
    "value": [
        n_total,
        n_equal,
        n_diff_001,
        n_diff_005,
        n_large,
        n_latent,
        corr_global,
        corr_filt,
        len(df_filt)
    ]
})

with pd.ExcelWriter(os.path.join(output_dir, "estadisticas_chiral.xlsx")) as writer:
    desc_global.to_excel(writer, sheet_name="global_stats")
    desc_filt.to_excel(writer, sheet_name="filtered_stats")
    summary_df.to_excel(writer, sheet_name="summary", index=False)

print("✔ estadisticas_chiral.xlsx saved.")

# ---------------------------------------------------------
# PLOTS
# ---------------------------------------------------------

# Utility: save inside output directory
def savefig(name):
    plt.savefig(os.path.join(output_dir, name), dpi=200)
    plt.close()

# 1. Scatter global
plt.figure(figsize=(8,6))
plt.scatter(df["f_chiral_semiautomatica"], df["f_chiral_automatica"], alpha=0.4)
plt.xlabel("f_chiral (semi-automatic)")
plt.ylabel("f_chiral (automatic)")
plt.title(f"Global correlation (r = {corr_global:.4f})")
plt.grid(True)
plt.tight_layout()
savefig("scatter_completo.png")

# 2. Scatter filtered
plt.figure(figsize=(8,6))
plt.scatter(df_filt["f_chiral_semiautomatica"], df_filt["f_chiral_automatica"], alpha=0.4)
plt.xlabel("f_chiral (semi-automatic)")
plt.ylabel("f_chiral (automatic)")
plt.title(f"Correlation f>0 in both (r = {corr_filt:.4f})")
plt.grid(True)
plt.tight_layout()
savefig("scatter_filtrado.png")

# 3. Identity plot
plt.figure(figsize=(8,6))
plt.scatter(df_filt["f_chiral_semiautomatica"], df_filt["f_chiral_automatica"], alpha=0.4)
minv = min(df_filt["f_chiral_semiautomatica"].min(), df_filt["f_chiral_automatica"].min())
maxv = max(df_filt["f_chiral_semiautomatica"].max(), df_filt["f_chiral_automatica"].max())
plt.plot([minv, maxv], [minv, maxv], 'r--', linewidth=2)
plt.xlabel("f_chiral (semi-automatic)")
plt.ylabel("f_chiral (automatic)")
plt.title("Identity Plot")
plt.grid(True)
plt.tight_layout()
savefig("identity_plot_filtrado.png")

# 4. Bland–Altman plot
semi = df_filt["f_chiral_semiautomatica"].values
auto = df_filt["f_chiral_automatica"].values
mean_vals = (semi + auto) / 2
diffs = semi - auto
mean_diff = np.mean(diffs)
sd_diff = np.std(diffs)

plt.figure(figsize=(8,6))
plt.scatter(mean_vals, diffs, alpha=0.4)
plt.axhline(mean_diff, color='red', linestyle='--', label=f"Mean diff = {mean_diff:.4f}")
plt.axhline(mean_diff + 1.96*sd_diff, color='gray', linestyle='--', label='Upper 95% limit')
plt.axhline(mean_diff - 1.96*sd_diff, color='gray', linestyle='--', label='Lower 95% limit')
plt.xlabel("Mean of both methods")
plt.ylabel("Difference (semi - auto)")
plt.title("Bland–Altman Plot")
plt.legend()
plt.grid(True)
plt.tight_layout()
savefig("bland_altman_filtrado.png")

# 5. Histogram of all differences
plt.figure(figsize=(8,6))
plt.hist(df["diff_abs"].dropna(), bins=40)
plt.xlabel("Difference (semi - auto)")
plt.ylabel("Frequency")
plt.title("Absolute Differences (All data)")
plt.grid(True)
plt.tight_layout()
savefig("hist_diferencias_total.png")

# 6. Histogram of filtered differences
plt.figure(figsize=(8,6))
plt.hist(df_filt["diff_abs"].dropna(), bins=40)
plt.xlabel("Difference (semi - auto)")
plt.ylabel("Frequency")
plt.title("Absolute Differences (Filtered)")
plt.grid(True)
plt.tight_layout()
savefig("hist_diferencias_filtrado.png")

# 7. Histogram of absolute differences (filtered)
plt.figure(figsize=(8,6))
plt.hist(np.abs(df_filt["diff_abs"].dropna()), bins=40)
plt.xlabel("|Difference|")
plt.ylabel("Frequency")
plt.title("Absolute Differences (Filtered, |semi - auto|)")
plt.grid(True)
plt.tight_layout()
savefig("hist_abs_diff_filtrado.png")

print("\n🎉 ALL RESULTS SAVED INSIDE:")
print(f"   {output_dir}")
print("Done!")
