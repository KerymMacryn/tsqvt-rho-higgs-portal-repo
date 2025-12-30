#!/usr/bin/env python3
"""
Fig2_production_modes.py
===========================
Production cross sections: VBF vs ggF with uncertainty bands.

Output: Fig2_production_modes.pdf
"""

import numpy as np
import matplotlib.pyplot as plt

# ============================================================================
# PLOT STYLE
# ============================================================================
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 11,
    'axes.labelsize': 13,
    'axes.titlesize': 13,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'legend.fontsize': 10,
    'figure.figsize': (7, 5),
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 1.0,
    'text.usetex': False,
})

# ============================================================================
# PHYSICS PARAMETERS
# ============================================================================

# Mass range [TeV]
m_rho = np.linspace(1.5, 4.5, 100)

# Benchmark mixing angle
sin2_theta = 0.035  # |θ| ≈ 0.19 rad ≈ 10.7°
sin2_theta_up = 0.04
sin2_theta_down = 0.03

# ============================================================================
# GGF CROSS SECTION (mixing-suppressed)
# ============================================================================
# σ_ggF = sin²θ × σ_SM^heavy(m)
# Fitted to reproduce benchmark table: σ_ggF(B1) = 0.13 fb at m = 2.26 TeV

def sigma_ggF(m, sin2th):
    """ggF cross section in fb, m in TeV"""
    # σ_SM^heavy fitted to give σ_ggF ~ 0.13 fb at benchmarks
    # σ_SM^heavy ≈ 3.5 fb at m = 2.3 TeV (from LHC HXSWG with K-factors)
    sigma_SM_heavy = 5.3 * (2.0 / m)**3  # fb, m in TeV
    return sin2th * sigma_SM_heavy

sigma_ggF_central = sigma_ggF(m_rho, sin2_theta)
sigma_ggF_up = sigma_ggF(m_rho, sin2_theta_up)
sigma_ggF_down = sigma_ggF(m_rho, sin2_theta_down)

# ============================================================================
# VBF CROSS SECTION (direct ρVV coupling)
# ============================================================================
# σ_VBF = σ_0 × (2 TeV / m)^n with σ_0 = 4.5 fb, n = 2.5

def sigma_VBF(m, sigma0=4.5, n=2.5):
    """VBF cross section in fb"""
    return sigma0 * (2.0 / m)**n

# Central and uncertainty (factor ~2 from spectral coefficients)
sigma_VBF_central = sigma_VBF(m_rho, sigma0=4.5, n=2.5)
sigma_VBF_up = sigma_VBF(m_rho, sigma0=6.5, n=2.3)  # Higher, flatter
sigma_VBF_down = sigma_VBF(m_rho, sigma0=3.0, n=2.7)  # Lower, steeper

# ============================================================================
# BENCHMARK POINTS
# ============================================================================
benchmarks = {
    'B1': {'m': 2.26, 'sigma_VBF': 3.0, 'sigma_ggF': 0.13},
    'B2': {'m': 2.26, 'sigma_VBF': 3.0, 'sigma_ggF': 0.11},
    'B3': {'m': 2.44, 'sigma_VBF': 3.0, 'sigma_ggF': 0.10},
}

# ============================================================================
# CREATE FIGURE
# ============================================================================
fig, ax = plt.subplots(figsize=(7, 5))

# VBF band and line
ax.fill_between(m_rho, sigma_VBF_down, sigma_VBF_up, 
                color='orange', alpha=0.25, label='VBF uncertainty')
ax.plot(m_rho, sigma_VBF_central, 'C1-', lw=2.5, 
        label=r'VBF (direct $\rho VV$)')

# ggF band and line
ax.fill_between(m_rho, sigma_ggF_down, sigma_ggF_up, 
                color='blue', alpha=0.15, label='ggF uncertainty')
ax.plot(m_rho, sigma_ggF_central, 'C0--', lw=2, 
        label=r'ggF (mixing, $\sin^2\theta = 0.035$)')

# Benchmark points
for i, (name, bp) in enumerate(benchmarks.items()):
    marker = 'o' if name != 'B3' else 's'
    # Only label first point to avoid legend clutter
    label = 'Benchmarks' if i == 0 else None
    ax.plot(bp['m'], bp['sigma_VBF'], 'r*', markersize=12, 
            markeredgecolor='darkred', markeredgewidth=0.5,
            label=label, zorder=10)

# Ratio annotation
ax.annotate('', xy=(3.0, 1.8), xytext=(3.0, 0.22),
            arrowprops=dict(arrowstyle='<->', color='gray', lw=1.5))
ax.text(3.12, 0.6, r'$R \sim 10$', fontsize=10, color='gray', va='center')

# ============================================================================
# AXIS CONFIGURATION
# ============================================================================
ax.set_xlabel(r'$m_\rho$ [TeV]', fontsize=13)
ax.set_ylabel(r'$\sigma$ [fb]', fontsize=13)
ax.set_xlim(1.5, 4.5)
ax.set_ylim(0.03, 15)
ax.set_yscale('log')

# Grid
ax.grid(True, which='major', linestyle='-', alpha=0.3)
ax.grid(True, which='minor', linestyle=':', alpha=0.2)

# Legend
ax.legend(loc='upper right', framealpha=0.95)

# Title
ax.set_title(r'Production cross sections at $\sqrt{s} = 14$ TeV: VBF dominance',
             fontsize=12, fontweight='bold')

# Info box
ax.text(0.03, 0.03, r'TSQVT $\rho$-Higgs portal',
        transform=ax.transAxes, fontsize=9, color='gray',
        verticalalignment='bottom')

# ============================================================================
# SAVE
# ============================================================================
plt.tight_layout()
plt.savefig('figures/Fig2_production_modes.pdf', format='pdf', bbox_inches='tight')
plt.savefig('figures/Fig2_production_modes.png', format='png', dpi=300, bbox_inches='tight')
print("Saved: Fig2_production_modes.pdf")
print("Saved: Fig2_production_modes.png")

# Print summary
print("\n" + "="*60)
print("PRODUCTION CROSS SECTIONS SUMMARY")
print("="*60)
print(f"{'m_ρ [TeV]':<12} {'σ_VBF [fb]':<15} {'σ_ggF [fb]':<15} {'R = VBF/ggF':<12}")
print("-"*60)
for m in [2.0, 2.3, 2.5, 3.0, 3.5]:
    vbf = sigma_VBF(m, sigma0=4.5, n=2.5)
    ggf = sigma_ggF(m, sin2_theta)
    print(f"{m:<12.1f} {vbf:<15.2f} {ggf:<15.3f} {vbf/ggf:<12.0f}")
print("="*60)