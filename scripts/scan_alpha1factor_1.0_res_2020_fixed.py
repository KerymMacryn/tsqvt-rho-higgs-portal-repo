"""
TSQVT rho-Higgs Portal Scanner - FINAL CORRECTED VERSION
Complete scan with proper cross section calculation, K-factor and fk profile scan
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
kappa_spectral = 5e4

# --- NEW: theory/systematic scan over spectral moments fk (profiles)
# We create a small ensemble: central, -30%, +30% variations
fk_profiles = {
    'central': 1.0,
    'fk_minus30pct': 0.7,
    'fk_plus30pct': 1.3,
}

# Effective portal scaling (will be recomputed per fk profile below)
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
pdf_file = output_dir / "Fig1_Deltakappa_contours.pdf"

# --- NEW: phenomenology/systematics knobs
# K-factor to approximate higher-order QCD corrections to ggF (multiplicative)
K_factor = 1.6  # typical NLO/NNLO K-factor estimate; tune as needed

# Direct production (pb). Poner 0.003 para probar ~1.2 fb con BR≈0.4
sigma_direct_const_pb = 0.003

# Output filename template (one file per fk profile)
csv_template = output_dir / "scan_Deltakappa_40x40_{profile}.csv"

# ============================================================
# 2. COMPUTE OBSERVABLES (with direct production and K-factor)
# ============================================================
def compute_observables(Lambda, mrho, alpha1_eff,
                        sigma_direct_const_pb_local: float = None,
                        sigma_direct_func=None):
    """
    Compute mixing angle, Deltakappa, and cross section including optional
    direct production of the rho resonance.

    Args:
        Lambda (float): cutoff scale [GeV]
        mrho (float): rho mass [GeV]
        alpha1_eff (float): effective portal coupling (dimensionless)
        sigma_direct_const_pb_local (float, optional): override constant direct production [pb]
        sigma_direct_func (callable, optional): function(mrho) -> sigma_direct_pb

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

    # Coupling rescaling
    cos_theta = np.cos(theta)
    Deltakappa_pct = (cos_theta - 1.0) * 100.0

    # Production cross section from mixing
    sin2theta = np.sin(theta)**2

    # SM-like heavy scalar production (gluon fusion, pb)
    sigma_ggF_pb = 50.0 * (125.0 / mrho)**3

    # Apply K-factor for higher-order QCD corrections
    sigma_ggF_pb *= K_factor

    # Mixing-suppressed production (pb)
    sigma_rho_from_mix_pb = sin2theta * sigma_ggF_pb

    # Direct production (pb) independent of mixing
    # Priority: sigma_direct_func > sigma_direct_const_pb_local > global sigma_direct_const_pb
    sigma_rho_direct_pb = 0.0
    if sigma_direct_func is not None and callable(sigma_direct_func):
        try:
            sigma_rho_direct_pb = float(sigma_direct_func(mrho))
        except Exception:
            sigma_rho_direct_pb = 0.0
    elif sigma_direct_const_pb_local is not None:
        sigma_rho_direct_pb = float(sigma_direct_const_pb_local)
    else:
        sigma_rho_direct_pb = float(sigma_direct_const_pb)

    # Total production (pb)
    sigma_rho_total_pb = sigma_rho_from_mix_pb + sigma_rho_direct_pb

    # Branching ratio to WW
    BR_WW = 0.40 if mrho > 2 * mt else 0.50

    # Final signal rate (fb)
    sigmaBR_WW_fb = sigma_rho_total_pb * BR_WW * 1000.0

    return theta, Deltakappa_pct, sigmaBR_WW_fb

# ============================================================
# DEBUG: Single point verification (will be run for central fk later)
# ============================================================
# We'll print a debug for a representative point after computing central profile

# ============================================================
# VECTORIZED GRID COMPUTATION (with outer loop over fk profiles)
# ============================================================
Lambda_grid, mrho_grid = np.meshgrid(Lambda_vals, mrho_vals)

for profile_name, fk_scale in fk_profiles.items():
    print(f"\n--- Running scan for fk profile: {profile_name} (scale={fk_scale:.2f}) ---")

    # Recompute effective portal coupling if it depends on fk (simple scaling here)
    # If you have a more detailed mapping from fk -> alpha1_eff, replace this line.
    alpha1_eff_profile = alpha1_central * kappa_spectral * fk_scale

    print("Starting scan (FINAL CORRECTED VERSION)...")
    print(f"Grid: {N}x{N} points")
    print(f"Lambda: {Lambda_vals.min()/1000:.1f}–{Lambda_vals.max()/1000:.1f} TeV")
    print(f"m_rho:  {mrho_vals.min()/1000:.1f}–{mrho_vals.max()/1000:.1f} TeV")
    print(f"\nPortal coupling (profile = {profile_name}):")
    print(f"  alpha1 (tree)   = {alpha1_central:.6e}")
    print(f"  Enhancement     = {kappa_spectral:.1f}")
    print(f"  fk scale        = {fk_scale:.3f}")
    print(f"  alpha1 (eff)    = {alpha1_eff_profile:.6e}")
    print(f"  Lambda_ref      = {Lambda_ref/1000:.1f} TeV\n")
    print(f"  K_factor        = {K_factor:.2f}")
    print(f"  sigma_direct_pb = {sigma_direct_const_pb:.6e} pb\n")

    theta_grid = np.zeros_like(Lambda_grid)
    Deltakappa_grid = np.zeros_like(Lambda_grid)
    sigmaBR_grid = np.zeros_like(Lambda_grid)

    for i in range(N):  # Loop over mrho (rows)
        for j in range(N):  # Loop over Lambda (columns)
            theta_grid[i, j], Deltakappa_grid[i, j], sigmaBR_grid[i, j] = \
                compute_observables(Lambda_grid[i, j], mrho_grid[i, j], alpha1_eff_profile,
                                    sigma_direct_const_pb_local=sigma_direct_const_pb)

    print(f"✓ Scan complete ({N*N} points) for profile {profile_name}")
    print(f"  θ range: [{np.degrees(theta_grid.min()):.3f}, {np.degrees(theta_grid.max()):.3f}]°")
    print(f"  Δκ range: [{Deltakappa_grid.min():.3f}, {Deltakappa_grid.max():.3f}]%")
    print(f"  σ×BR range: [{sigmaBR_grid.min():.4e}, {sigmaBR_grid.max():.4e}] fb\n")

    # ============================================================
    # SAVE TO CSV (one file per fk profile)
    # ============================================================
    df = pd.DataFrame({
        'Lambda_GeV': Lambda_grid.flatten(),
        'mrho_GeV': mrho_grid.flatten(),
        'theta_rad': theta_grid.flatten(),
        'theta_deg': np.degrees(theta_grid.flatten()),
        'Deltakappa_pct': Deltakappa_grid.flatten(),
        'sigmaBR_WW_fb': sigmaBR_grid.flatten(),
        'fk_profile': profile_name
    })

    csv_file_profile = csv_template.with_name(csv_template.name.format(profile=profile_name))
    df.to_csv(csv_file_profile, index=False)
    print(f"✓ Data saved to {csv_file_profile}")

    # ============================================================
    # CREATE FIGURE (per profile)
    # ============================================================
    print("Generating figure for profile:", profile_name)

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
                 rf'Profile: {profile_name}; $ \kappa_{{\mathrm{{spec}}}} = {kappa_spectral:.0f}$',
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
    pdf_file_profile = output_dir / f"Fig1_Deltakappa_contours_{profile_name}.pdf"
    plt.savefig(pdf_file_profile, dpi=300, bbox_inches='tight')
    print(f"✓ Figure saved: {pdf_file_profile}")
    plt.close(fig)

    # ============================================================
    # 4.5 FIND OPTIMAL POINTS (ADDED)
    # ============================================================
    print("\n" + "="*60)
    print(f"FINDING OPTIMAL POINTS IN VIABLE REGION (profile={profile_name}):")
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
    # 5. BENCHMARK TABLE (per profile)
    # ============================================================
    print("\n" + "="*60)
    print(f"BENCHMARK VALUES (profile={profile_name}):")
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

print("\nAll profiles completed.")
