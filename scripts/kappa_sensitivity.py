import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# CONFIG
# ============================================================
results_dir = Path("results")

# These must match the profiles used in your scan
kappa_values = {
    "central": 5e4,
    "fk_minus30pct": 5e4 * 0.7,
    "fk_plus30pct": 5e4 * 1.3,
}

csv_files = {
    profile: results_dir / f"scan_Deltakappa_40x40_{profile}.csv"
    for profile in kappa_values
}

# ============================================================
# Extract optimal points for each profile
# ============================================================
m_opt = []
dk_opt = []
kappa_list = []

for profile, kappa in kappa_values.items():
    df = pd.read_csv(csv_files[profile])

    # Select viable region |Δκ| ≤ 2%
    mask = np.abs(df["Deltakappa_pct"]) <= 2.0
    df_viable = df[mask]

    # Find point with maximum σ×BR
    idx = df_viable["sigmaBR_WW_fb"].idxmax()
    row = df_viable.loc[idx]

    m_opt.append(row["mrho_GeV"] / 1000.0)  # TeV
    dk_opt.append(row["Deltakappa_pct"])
    kappa_list.append(kappa)

# Convert to arrays
kappa_arr = np.array(kappa_list)
m_opt = np.array(m_opt)
dk_opt = np.array(dk_opt)

# ============================================================
# Plot
# ============================================================
plt.figure(figsize=(8,6))

# Left axis: m_rho,opt
plt.plot(kappa_arr, m_opt, marker="o", color="navy", label=r"$m_{\rho,\mathrm{opt}}$")

# Right axis: |Δκ|_opt
plt.twinx()
plt.plot(kappa_arr, np.abs(dk_opt), marker="s", color="darkred",
         linestyle="--", label=r"$|\Delta\kappa|_{\mathrm{opt}}$")

plt.xscale("log")
plt.xlabel(r"$\kappa_{\mathrm{spec}}$", fontsize=14)
plt.title(r"Sensitivity of $m_{\rho,\mathrm{opt}}$ and $|\Delta\kappa|_{\mathrm{opt}}$ to $\kappa_{\mathrm{spec}}$", fontsize=14)

plt.grid(True, linestyle=":", alpha=0.5)

# Save
out = results_dir / "kappa_sensitivity.pdf"
plt.savefig(out, dpi=300, bbox_inches="tight")
print(f"✓ Figure saved to {out}")
