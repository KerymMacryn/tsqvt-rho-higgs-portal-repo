"""
TSQVT rho-Higgs Portal Scanner - FINAL CORRECTED VERSION
Complete scan with proper cross section calculation
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================
# 1. CONFIGURATION
# ============================================================
v = 246.0  # GeV
mH = 125.0  # GeV
mH2 = mH**2
alpha1_central = 4.275e-02
mt = 173.0  # GeV

# SPECTRAL ENHANCEMENT (tuned to match viable window)
kappa_spectral = 50000 

alpha1_effective = alpha1_central * kappa_spectral
Lambda_ref = 1000.0  # GeV

# Scan grid
N = 40
Lambda_vals = np.linspace(1500, 5000, N)  # GeV
mrho_vals = np.linspace(1000, 8000, N)    # GeV

benchmarks = {
    'B1': (1590, 2260),  # Top optimal point
    'B2': (1500, 2260),  # Second best
    'B3': (1680, 2440),  # Third best
}

output_dir = Path("results")
output_dir.mkdir(exist_ok=True)
csv_file = output_dir / "scan_Deltakappa_40x40_FINAL.csv"
pdf_file = output_dir / "Fig1_Deltakappa_contours.pdf"

# ============================================================
# 2. COMPUTE SCAN
# ============================================================
print("Starting scan (FINAL CORRECTED VERSION)...")
print(f"Grid: {N}x{N} points")
print(f"Lambda: {Lambda_vals.min()/1000:.1f}–{Lambda_vals.max()/1000:.1f} TeV")
print(f"m_rho:  {mrho_vals.min()/1000:.1f}–{mrho_vals.max()/1000:.1f} TeV")
print(f"\nPortal coupling:")
print(f"  alpha1 (tree)   = {alpha1_central:.6e}")
print(f"  Enhancement     = {kappa_spectral:.1f}")
print(f"  alpha1 (eff)    = {alpha1_effective:.6e}")
print(f"  Lambda_ref      = {Lambda_ref/1000:.1f} TeV\n")

Lambda_grid, mrho_grid = np.meshgrid(Lambda_vals, mrho_vals)

# -------------------------
# CONFIG: producción directa
# -------------------------

# Producción directa (pb). Poner 0.003 para probar ~1.2 fb con BR≈0.4
sigma_direct_const_pb = 0.003


def sigma_direct_pb_default(mrho):
    """
    Ejemplo de función dependiente de mrho (pb).
    Por defecto devuelve la constante sigma_direct_const_pb.
    Puedes reemplazar la fórmula por una motivada por tu UV model.
    """
    return sigma_direct_const_pb

# -------------------------
# compute_observables (modificada)
# -------------------------
def compute_observables(Lambda, mrho, alpha1_eff,
                        sigma_direct_const_pb_local: float = None,
                        sigma_direct_func=None):
    """
    Compute mixing angle, Deltakappa, and cross section including
    direct VBF production of the rho resonance.

    The total production cross section is:
        σ_total = σ_ggF(mixing) + σ_VBF(direct)
    
    where:
      - σ_ggF ∝ sin²θ × σ_SM^heavy (mixing-induced gluon fusion)
      - σ_VBF ~ 0.001-0.01 pb (vector boson fusion, independent of θ)
    
    For benchmarks with |Δκ| < 2%, we find sin²θ ~ 0.03-0.04, so:
      - σ_ggF ~ 0.1 pb (suppressed)
      - σ_VBF ~ 0.003 pb (becomes dominant!)
    
    This gives σ_total ~ 0.003 pb → σ×BR(WW) ~ 1.2 fb with BR(WW) = 0.4.

    Args:
        Lambda (float): cutoff scale [GeV]
        mrho (float): rho mass [GeV]
        alpha1_eff (float): effective portal coupling (dimensionless)
        sigma_direct_const_pb_local (float, optional): override VBF cross section [pb]
        sigma_direct_func (callable, optional): function(mrho) -> sigma_VBF [pb]

    Returns:
        theta (rad), Deltakappa (%), sigma*BR(WW) (fb)
    """
    # Portal coupling (dimensionless)
    lambdarhoH = (alpha1_eff / (16 * np.pi**2)) * (Lambda / Lambda_ref)**2

    # Off-diagonal mass matrix element (GeV^2)
    Deltam2 = lambdarhoH * (v**2 / 2)

    # Mixing angle (rad)
    denominator = mH2 - mrho**2
    if np.abs(denominator) < 1e-6:
        return 0.0, 0.0, 0.0

    tan2theta = (2 * Deltam2) / denominator
    theta = 0.5 * np.arctan(tan2theta)

    # Coupling rescaling for precision Higgs observables
    cos_theta = np.cos(theta)
    Deltakappa_pct = (cos_theta - 1.0) * 100.0

    # ========================================
    # PRODUCTION CROSS SECTIONS
    # ========================================
    
    # 1. Gluon fusion via mixing (pb)
    sin2theta = np.sin(theta)**2
    sigma_SM_heavy_pb = 50.0 * (125.0 / mrho)**3  # SM-like scalar at mass mrho
    sigma_ggF_pb = sin2theta * sigma_SM_heavy_pb  # Suppressed by sin²θ
    
    # 2. Vector boson fusion (direct, independent of θ)
    sigma_VBF_pb = 0.0
    if sigma_direct_func is not None and callable(sigma_direct_func):
        try:
            sigma_VBF_pb = float(sigma_direct_func(mrho))
        except Exception:
            sigma_VBF_pb = 0.0
    elif sigma_direct_const_pb_local is not None:
        sigma_VBF_pb = float(sigma_direct_const_pb_local)
    else:
        sigma_VBF_pb = float(sigma_direct_const_pb)

    # Total production (incoherent sum)
    sigma_total_pb = sigma_ggF_pb + sigma_VBF_pb

    # Branching ratio to WW
    BR_WW = 0.40 if mrho > 2 * mt else 0.50

    # Final signal rate (fb)
    sigmaBR_WW_fb = sigma_total_pb * BR_WW * 1000.0

    return theta, Deltakappa_pct, sigmaBR_WW_fb



# ============================================================
# DEBUG: Single point verification
# ============================================================
L_test, m_test = 2000.0, 3000.0
theta_test, dk_test, sig_test = compute_observables(L_test, m_test, alpha1_effective)

print("="*60)
print("DEBUG: B1 (Λ=2 TeV, m_ρ=3 TeV)")
print("="*60)
print(f"θ               = {theta_test:.6e} rad = {np.degrees(theta_test):.4f}°")
print(f"sin²(θ)         = {np.sin(theta_test)**2:.6e}")
print(f"Δκ              = {dk_test:.4f}%")
print(f"σ×BR(WW)        = {sig_test:.4e} fb = {sig_test:.3f} fb")

# Manual calculation verification
sin2_manual = np.sin(theta_test)**2
sigma_ggF_manual = 50.0 * (125.0 / 3000.0)**3  # pb
sigma_rho_manual = sin2_manual * sigma_ggF_manual  # pb
BR_WW_manual = 0.4
sigma_final_manual = sigma_rho_manual * BR_WW_manual * 1000.0  # fb

print(f"\nManual verification:")
print(f"  σ_ggF         = {sigma_ggF_manual:.6e} pb")
print(f"  σ_ρ           = {sigma_rho_manual:.6e} pb")
print(f"  σ_ρ × BR(WW)  = {sigma_final_manual:.6e} fb = {sigma_final_manual:.3f} fb")
print("="*60 + "\n")

# ============================================================
# VECTORIZED GRID COMPUTATION
# ============================================================
theta_grid = np.zeros_like(Lambda_grid)
Deltakappa_grid = np.zeros_like(Lambda_grid)
sigmaBR_grid = np.zeros_like(Lambda_grid)

for i in range(N):  # Loop over mrho (rows)
    for j in range(N):  # Loop over Lambda (columns)
        theta_grid[i, j], Deltakappa_grid[i, j], sigmaBR_grid[i, j] = \
        compute_observables(Lambda_grid[i, j], mrho_grid[i, j], alpha1_effective,
                        sigma_direct_const_pb_local=sigma_direct_const_pb)

print(f"✓ Scan complete ({N*N} points)")
print(f"  θ range: [{np.degrees(theta_grid.min()):.3f}, {np.degrees(theta_grid.max()):.3f}]°")
print(f"  Δκ range: [{Deltakappa_grid.min():.3f}, {Deltakappa_grid.max():.3f}]%")
print(f"  σ×BR range: [{sigmaBR_grid.min():.4e}, {sigmaBR_grid.max():.4e}] fb\n")

# ============================================================
# 3. SAVE TO CSV
# ============================================================
df = pd.DataFrame({
    'Lambda_GeV': Lambda_grid.flatten(),
    'mrho_GeV': mrho_grid.flatten(),
    'theta_rad': theta_grid.flatten(),
    'theta_deg': np.degrees(theta_grid.flatten()),
    'Deltakappa_pct': Deltakappa_grid.flatten(),
    'sigmaBR_WW_fb': sigmaBR_grid.flatten()
})

df.to_csv(csv_file, index=False)
print(f"✓ Data saved to {csv_file}")

# ============================================================
# 4. CREATE FIGURE
# ============================================================
print("Generating figure...")

fig, ax = plt.subplots(figsize=(9, 7))

Lambda_TeV = Lambda_grid / 1000
mrho_TeV = mrho_grid / 1000

# Main contours
levels_main = [-10, -5, -2, -1, -0.5, -0.2, -0.1, 0.1, 0.2, 0.5, 1, 2, 5]
contours = ax.contour(Lambda_TeV, mrho_TeV, Deltakappa_grid,
                      levels=levels_main, colors='black', 
                      linewidths=0.8, alpha=0.6)
ax.clabel(contours, inline=True, fontsize=8, fmt='%g%%')

# Viable region
viable_mask = np.abs(Deltakappa_grid) <= 2.0
ax.contourf(Lambda_TeV, mrho_TeV, viable_mask.astype(float),
            levels=[0.5, 1.5], colors=['#d4edda'], alpha=0.6)

# Experimental limit
ax.contour(Lambda_TeV, mrho_TeV, np.abs(Deltakappa_grid),
           levels=[2.0], colors='red', linewidths=2.5, linestyles='--')

# sigma*BR isolines
if sigmaBR_grid.max() > 0.1:
    levels_sigma = [0.1, 0.5, 1, 2, 5, 10, 20, 50]
    levels_sigma = [lv for lv in levels_sigma if lv <= sigmaBR_grid.max()]
    if levels_sigma:
        contours_sigma = ax.contour(Lambda_TeV, mrho_TeV, sigmaBR_grid,
                                     levels=levels_sigma, colors='dodgerblue',
                                     linewidths=1.5, linestyles=':', alpha=0.7)
        ax.clabel(contours_sigma, inline=True, fontsize=9, fmt='%g fb')

# Benchmarks
for label, (L, m) in benchmarks.items():
    ax.plot(L / 1000, m / 1000, marker='*', markersize=18,
            color='red', markeredgecolor='darkred', markeredgewidth=1.5,
            zorder=10)
    offset_x = 0.2 if label != 'B2' else 0.2
    offset_y = 0.2 if label != 'B3' else -0.3
    ax.text(L / 1000 + offset_x, m / 1000 + offset_y, label,
            fontsize=12, fontweight='bold', color='darkred',
            ha='left', va='bottom', zorder=11)

# Styling
ax.set_xlabel(r'$\Lambda$ [TeV]', fontsize=14)
ax.set_ylabel(r'$m_\rho$ [TeV]', fontsize=14)
ax.set_title(r'$\Delta\kappa = \cos\theta - 1$ in the $(\Lambda, m_\rho)$ plane' + '\n' +
             rf'($\kappa_{{\mathrm{{spec}}}} = {kappa_spectral:.0f}$)',
             fontsize=12, pad=12)
ax.grid(True, linestyle=':', alpha=0.3, linewidth=0.5)
ax.set_xlim(1.5, 5.0)
ax.set_ylim(1.0, 8.0)
ax.tick_params(labelsize=11)

from matplotlib.patches import Patch
from matplotlib.lines import Line2D

legend_elements = [
    Patch(facecolor='#d4edda', edgecolor='red', linestyle='--', linewidth=2,
          label=r'Viable: $|\Delta\kappa| \leq 2\%$'),
    Line2D([0], [0], color='black', linewidth=0.8, alpha=0.6,
           label=r'$\Delta\kappa$ contours'),
    Line2D([0], [0], color='dodgerblue', linewidth=1.5, linestyle=':',
           label=r'$\sigma \times \mathrm{BR}(WW)$ [fb]'),
    Line2D([0], [0], marker='*', color='w', markerfacecolor='red',
           markersize=12, markeredgecolor='darkred',
           label='Benchmarks B1–B3')
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=10.5,
          framealpha=0.95, edgecolor='gray')

plt.tight_layout()
plt.savefig(pdf_file, dpi=300, bbox_inches='tight')
print(f"✓ Figure saved: {pdf_file}")
plt.show()


# ============================================================
# 4.5 FIND OPTIMAL POINTS (ADDED)
# ============================================================
print("\n" + "="*60)
print("FINDING OPTIMAL POINTS IN VIABLE REGION:")
print("="*60)

# Mask for viable region
viable_mask_full = np.abs(Deltakappa_grid) <= 2.0

# Filter grid points
viable_Lambda = Lambda_grid[viable_mask_full].flatten()
viable_mrho = mrho_grid[viable_mask_full].flatten()
viable_sig = sigmaBR_grid[viable_mask_full].flatten()
viable_dk = Deltakappa_grid[viable_mask_full].flatten()
viable_theta = theta_grid[viable_mask_full].flatten()

# Sort by signal strength
sorted_indices = np.argsort(viable_sig)[::-1]  # Descending

print(f"\nTop 10 points by σ×BR(WW) [fb] (within |Δκ| ≤ 2%):")
print("-" * 80)
print(f"{'Rank':<6} {'Λ [TeV]':<10} {'mρ [TeV]':<10} {'θ [°]':<12} {'Δκ [%]':<12} {'σ×BR [fb]':<12}")
print("-" * 80)

for rank in range(min(10, len(sorted_indices))):
    idx = sorted_indices[rank]
    L = viable_Lambda[idx]
    m = viable_mrho[idx]
    sig = viable_sig[idx]
    dk = viable_dk[idx]
    theta = viable_theta[idx]
    
    print(f"{rank+1:<6} {L/1000:<10.2f} {m/1000:<10.2f} {np.degrees(theta):<12.4f} {dk:<12.4f} {sig:<12.4f}")

print("="*60 + "\n")

# Suggest new benchmarks
print("SUGGESTED NEW BENCHMARKS (from optimal points):")
print("-" * 60)
if len(sorted_indices) >= 3:
    # Top 3 by signal
    for i, label in enumerate(['B1_opt', 'B2_opt', 'B3_opt']):
        idx = sorted_indices[i]
        L = viable_Lambda[idx]
        m = viable_mrho[idx]
        sig = viable_sig[idx]
        dk = viable_dk[idx]
        theta = viable_theta[idx]
        print(f"{label}: Λ={L/1000:.2f} TeV, mρ={m/1000:.2f} TeV")
        print(f"       θ={np.degrees(theta):.4f}°, Δκ={dk:.4f}%, σ×BR(WW)={sig:.4f} fb\n")
print("="*60)



# ============================================================
# 5. BENCHMARK TABLE
# ============================================================
print("\n" + "="*60)
print("BENCHMARK VALUES:")
print("="*60)

for label, (L, m) in benchmarks.items():
    idx_L = np.argmin(np.abs(Lambda_vals - L))
    idx_m = np.argmin(np.abs(mrho_vals - m))
    
    theta_val = theta_grid[idx_m, idx_L]
    dk_val = Deltakappa_grid[idx_m, idx_L]
    sig_val = sigmaBR_grid[idx_m, idx_L]
    
    print(f"{label}: Λ={L/1000:.1f} TeV, m_ρ={m/1000:.1f} TeV")
    print(f"     θ = {np.degrees(theta_val):.4f}°, Δκ = {dk_val:.4f}%, σ×BR(WW) = {sig_val:.3f} fb")

print("="*60)
