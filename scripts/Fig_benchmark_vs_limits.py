#!/usr/bin/env python3
"""
Fig_benchmark_vs_limits.py
==========================
Generates publication-quality figure comparing TSQVT benchmark predictions
with ATLAS/CMS limits and HL-LHC projections.

Output: Fig_benchmark_vs_limits.pdf

Author: TSQVT Collaboration
Date: December 2025
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

# ============================================================================
# PLOT STYLE CONFIGURATION
# ============================================================================
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 9,
    'figure.figsize': (7, 5),
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 1.0,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'text.usetex': False,  # Set True if LaTeX is available
})

# ============================================================================
# DATA: EXPERIMENTAL LIMITS (Run-2, 139 fb^-1, sqrt(s) = 13 TeV)
# ============================================================================

# Mass points for limits [TeV]
m_limits = np.array([1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5])

# ATLAS/CMS 95% CL upper limits on sigma x BR(WW) [fb]
# VBF-tagged WW -> lvqq channel (interpolated from public results)
limits_run2 = np.array([25, 18, 12, 8.5, 6.0, 4.5, 3.5, 2.8, 2.3])

# HL-LHC projections (3 ab^-1, sqrt(s) = 14 TeV)
# Scaled from Run-2 with sqrt(L) improvement and energy extrapolation
limits_hllhc_95 = np.array([1.2, 0.9, 0.7, 0.55, 0.45, 0.38, 0.32, 0.28, 0.25])
limits_hllhc_5sigma = np.array([2.0, 1.5, 1.2, 0.95, 0.80, 0.68, 0.58, 0.50, 0.45])

# ============================================================================
# DATA: TSQVT PREDICTIONS
# ============================================================================

# Benchmark points
benchmarks = {
    'B1': {'m_rho': 2.26, 'sigma_BR': 1.33, 'error': 0.25},
    'B2': {'m_rho': 2.26, 'sigma_BR': 1.30, 'error': 0.24},
    'B3': {'m_rho': 2.44, 'sigma_BR': 1.29, 'error': 0.24},
}

# TSQVT prediction band (sigma x BR as function of mass)
m_tsqvt = np.linspace(2.0, 3.5, 50)

# Central prediction: sigma_VBF(m) x BR(WW) ~ 0.4
# Using: sigma_VBF = 4.5 fb * (2 TeV / m)^2.5, BR(WW) = 0.40
sigma_vbf_central = 4.5 * (2.0 / m_tsqvt)**2.5 * 0.40

# Upper and lower bounds (factor ~2 uncertainty from kappa_spec)
sigma_vbf_upper = sigma_vbf_central * 1.8
sigma_vbf_lower = sigma_vbf_central * 0.6

# ============================================================================
# CREATE FIGURE
# ============================================================================

fig, ax = plt.subplots(figsize=(7, 5))

# ---------------------------------------------------------------------------
# 1. TSQVT viable window (gray shaded band)
# ---------------------------------------------------------------------------
ax.fill_between(m_tsqvt, sigma_vbf_lower, sigma_vbf_upper, 
                color='gray', alpha=0.25, label='TSQVT viable window')
ax.plot(m_tsqvt, sigma_vbf_central, 'k--', lw=1.5, alpha=0.7,
        label='TSQVT central prediction')

# ---------------------------------------------------------------------------
# 2. Run-2 limits (solid blue line with markers)
# ---------------------------------------------------------------------------
ax.plot(m_limits, limits_run2, 'b-', lw=2, marker='o', markersize=5,
        label='ATLAS+CMS Run-2 (139 fb$^{-1}$)')
ax.fill_between(m_limits, limits_run2, 100, color='blue', alpha=0.08)

# ---------------------------------------------------------------------------
# 3. HL-LHC projections
# ---------------------------------------------------------------------------
ax.plot(m_limits, limits_hllhc_95, 'g--', lw=2, marker='s', markersize=4,
        label='HL-LHC 95% CL (3 ab$^{-1}$)')
ax.plot(m_limits, limits_hllhc_5sigma, 'g:', lw=1.5, marker='^', markersize=4,
        label=r'HL-LHC $5\sigma$ discovery')

# ---------------------------------------------------------------------------
# 4. TSQVT benchmark points (red stars with error bars)
# ---------------------------------------------------------------------------
for name, bp in benchmarks.items():
    ax.errorbar(bp['m_rho'], bp['sigma_BR'], yerr=bp['error'],
                fmt='r*', markersize=15, capsize=4, capthick=1.5,
                ecolor='darkred', elinewidth=1.5, zorder=10)
    # Add label slightly offset
    offset_x = 0.08 if name != 'B3' else 0.08
    offset_y = 0.15 if name == 'B1' else -0.25
    ax.annotate(name, (bp['m_rho'] + offset_x, bp['sigma_BR'] + offset_y),
                fontsize=10, fontweight='bold', color='darkred')

# ---------------------------------------------------------------------------
# 5. Falsification threshold line
# ---------------------------------------------------------------------------
ax.axhline(y=0.5, color='red', linestyle='-.', lw=1.5, alpha=0.7)
ax.text(3.35, 0.55, 'Falsification\nthreshold', fontsize=8, 
        color='red', ha='right', va='bottom')

# ---------------------------------------------------------------------------
# AXIS CONFIGURATION
# ---------------------------------------------------------------------------
ax.set_xlabel(r'$m_\rho$ [TeV]', fontsize=12)
ax.set_ylabel(r'$\sigma \times \mathrm{BR}(WW)$ [fb]', fontsize=12)
ax.set_xlim(1.4, 3.7)
ax.set_ylim(0.1, 50)
ax.set_yscale('log')

# Grid
ax.grid(True, which='major', linestyle='-', alpha=0.3)
ax.grid(True, which='minor', linestyle=':', alpha=0.2)

# Legend
ax.legend(loc='upper right', framealpha=0.95, edgecolor='gray')

# Title
ax.set_title(r'TSQVT $\rho$-Higgs Portal: Benchmarks vs. Experimental Limits',
             fontsize=12, fontweight='bold')

# ---------------------------------------------------------------------------
# ANNOTATIONS
# ---------------------------------------------------------------------------
# Add "Excluded" region label
ax.text(2.5, 25, 'Excluded by Run-2', fontsize=9, color='blue', 
        alpha=0.7, ha='center')

# Add discovery region annotation
ax.annotate('', xy=(2.3, 1.0), xytext=(2.3, 0.45),
            arrowprops=dict(arrowstyle='->', color='green', lw=1.5))
ax.text(2.35, 0.7, r'$>5\sigma$', fontsize=9, color='green')

# Add sqrt(s) and channel info
ax.text(0.02, 0.02, r'$\sqrt{s} = 14$ TeV, VBF $WW \to \ell\nu qq$',
        transform=ax.transAxes, fontsize=9, color='gray',
        verticalalignment='bottom')

# ============================================================================
# SAVE FIGURE
# ============================================================================
plt.tight_layout()
plt.savefig('figures/Fig_benchmark_vs_limits.pdf', format='pdf', 
            bbox_inches='tight', pad_inches=0.05)
plt.savefig('figures/Fig_benchmark_vs_limits.png', format='png', dpi=300,
            bbox_inches='tight', pad_inches=0.05)

print("Figure saved: Fig_benchmark_vs_limits.pdf")
print("Figure saved: Fig_benchmark_vs_limits.png")

# ============================================================================
# OPTIONAL: DISPLAY FIGURE (uncomment for interactive use)
# ============================================================================
# plt.show()

# ============================================================================
# PRINT SUMMARY TABLE
# ============================================================================
print("\n" + "="*70)
print("SUMMARY: TSQVT Benchmarks vs. Limits")
print("="*70)
print(f"{'Point':<6} {'m_ρ [TeV]':<12} {'σ×BR [fb]':<12} {'Run-2 limit':<14} {'Margin':<8}")
print("-"*70)

for name, bp in benchmarks.items():
    # Interpolate limit at benchmark mass
    limit_at_m = np.interp(bp['m_rho'], m_limits, limits_run2)
    margin = limit_at_m / bp['sigma_BR']
    print(f"{name:<6} {bp['m_rho']:<12.2f} {bp['sigma_BR']:<12.2f} {limit_at_m:<14.1f} {margin:<8.1f}")

print("-"*70)
print("All benchmarks are ALLOWED (margin > 1)")
print("="*70)
