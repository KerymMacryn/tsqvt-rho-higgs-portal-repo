#!/usr/bin/env python3
"""
Fig3_branching_ratios.py
===========================
Decay branching ratios with realistic threshold effects.


Output: Fig3_branching_ratios.pdf
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

# Mass range [TeV] - start above 2*m_t to avoid threshold complications
m_rho = np.linspace(1.0, 5.0, 100)

# Particle masses [GeV]
m_W = 80.4
m_Z = 91.2
m_t = 173.0
m_h = 125.0

# ============================================================================
# PARTIAL WIDTH CALCULATIONS
# ============================================================================

def compute_BRs(m_TeV):
    """
    Compute branching ratios matching TSQVT benchmarks.
    Pattern: BR(WW) > BR(tt) ~ BR(ZZ) > BR(hh)
    Values from Eq. (4.15): BR(WW)~40%, BR(ZZ)~18%, BR(tt)~25%, BR(hh)~12%
    
    Mild mass dependence from phase space near thresholds.
    """
    m_GeV = m_TeV * 1000
    
    # Base values (asymptotic, high mass)
    BR_WW_base = 0.40
    BR_ZZ_base = 0.18
    BR_tt_base = 0.25
    BR_hh_base = 0.12
    # Remaining ~5% goes to other channels (bb, gg, etc.) - absorbed in normalization
    
    # Phase space suppression near thresholds
    # tt threshold effect
    if m_GeV < 2 * m_t:
        ps_tt = 0.0
    else:
        beta_t = np.sqrt(1 - (2*m_t/m_GeV)**2)
        ps_tt = beta_t**3  # P-wave suppression
    
    # hh threshold effect  
    if m_GeV < 2 * m_h:
        ps_hh = 0.0
    else:
        beta_h = np.sqrt(1 - (2*m_h/m_GeV)**2)
        ps_hh = beta_h  # S-wave
    
    # VV channels: always open, mild phase space
    ps_WW = np.sqrt(max(0, 1 - (2*m_W/m_GeV)**2))
    ps_ZZ = np.sqrt(max(0, 1 - (2*m_Z/m_GeV)**2))
    
    # Combine: partial widths proportional to base × phase space
    Gamma_WW = BR_WW_base * ps_WW
    Gamma_ZZ = BR_ZZ_base * ps_ZZ
    Gamma_tt = BR_tt_base * ps_tt
    Gamma_hh = BR_hh_base * ps_hh
    
    Gamma_tot = Gamma_WW + Gamma_ZZ + Gamma_tt + Gamma_hh
    
    if Gamma_tot < 1e-10:
        return 0, 0, 0, 0
    
    return (Gamma_WW/Gamma_tot, Gamma_ZZ/Gamma_tot, 
            Gamma_tt/Gamma_tot, Gamma_hh/Gamma_tot)

# ============================================================================
# COMPUTE BRs VS MASS
# ============================================================================

# Arrays for results
BR_WW = np.zeros_like(m_rho)
BR_ZZ = np.zeros_like(m_rho)
BR_tt = np.zeros_like(m_rho)
BR_hh = np.zeros_like(m_rho)

# Central values
for i, m in enumerate(m_rho):
    BR_WW[i], BR_ZZ[i], BR_tt[i], BR_hh[i] = compute_BRs(m)

# Uncertainty bands (±15% relative variation in each channel)
uncertainty = 0.15
BR_WW_up = BR_WW * (1 + uncertainty)
BR_WW_down = BR_WW * (1 - uncertainty)
BR_ZZ_up = BR_ZZ * (1 + uncertainty)
BR_ZZ_down = BR_ZZ * (1 - uncertainty)
BR_tt_up = BR_tt * (1 + uncertainty)
BR_tt_down = BR_tt * (1 - uncertainty)
BR_hh_up = BR_hh * (1 + uncertainty)
BR_hh_down = BR_hh * (1 - uncertainty)

# Convert to percentages
BR_WW *= 100; BR_WW_up *= 100; BR_WW_down *= 100
BR_ZZ *= 100; BR_ZZ_up *= 100; BR_ZZ_down *= 100
BR_tt *= 100; BR_tt_up *= 100; BR_tt_down *= 100
BR_hh *= 100; BR_hh_up *= 100; BR_hh_down *= 100

# ============================================================================
# CREATE FIGURE
# ============================================================================
fig, ax = plt.subplots(figsize=(7, 5))

# Mask for m > 2*m_t (where all channels are open)
mask = m_rho > 0.35  # 350 GeV, safely above tt threshold

# Colors
c_WW = 'C0'
c_ZZ = 'C1'
c_tt = 'C2'
c_hh = 'C3'

# Plot with uncertainty bands
ax.fill_between(m_rho[mask], BR_WW_down[mask], BR_WW_up[mask], color=c_WW, alpha=0.2)
ax.plot(m_rho[mask], BR_WW[mask], c_WW+'-', lw=2.5, label=r'$WW$')

ax.fill_between(m_rho[mask], BR_ZZ_down[mask], BR_ZZ_up[mask], color=c_ZZ, alpha=0.2)
ax.plot(m_rho[mask], BR_ZZ[mask], c_ZZ+'-', lw=2.5, label=r'$ZZ$')

ax.fill_between(m_rho[mask], BR_tt_down[mask], BR_tt_up[mask], color=c_tt, alpha=0.2)
ax.plot(m_rho[mask], BR_tt[mask], c_tt+'-', lw=2.5, label=r'$t\bar{t}$')

ax.fill_between(m_rho[mask], BR_hh_down[mask], BR_hh_up[mask], color=c_hh, alpha=0.2)
ax.plot(m_rho[mask], BR_hh[mask], c_hh+'-', lw=2.5, label=r'$hh$')

# Mark viable window
ax.axvspan(2.0, 3.5, color='green', alpha=0.08, label='Viable window')

# ============================================================================
# AXIS CONFIGURATION
# ============================================================================
ax.set_xlabel(r'$m_\rho$ [TeV]', fontsize=13)
ax.set_ylabel(r'Branching Ratio (%)', fontsize=13)
ax.set_xlim(1.0, 5.0)
ax.set_ylim(0, 50)

# Grid
ax.grid(True, which='major', linestyle='-', alpha=0.3)

# Legend
ax.legend(loc='upper right', framealpha=0.95, ncol=2)

# Title
ax.set_title(r'Decay branching ratios of $h_2$ (TSQVT benchmarks)',
             fontsize=12, fontweight='bold')

# Annotation
ax.text(0.03, 0.97, r'$\mathrm{BR}(VV) > \mathrm{BR}(hh)$: TSQVT signature',
        transform=ax.transAxes, fontsize=9, color='gray',
        verticalalignment='top')

# ============================================================================
# SAVE
# ============================================================================
plt.tight_layout()
plt.savefig('figures/Fig3_branching_ratios.pdf', format='pdf', bbox_inches='tight')
plt.savefig('figures/Fig3_branching_ratios.png', format='png', dpi=300, bbox_inches='tight')
print("Saved: Fig3_branching_ratios.pdf")
print("Saved: Fig3_branching_ratios.png")

# Print summary at benchmark masses
print("\n" + "="*70)
print("BRANCHING RATIOS AT BENCHMARK MASSES")
print("="*70)
print(f"{'m_ρ [TeV]':<12} {'BR(WW) %':<12} {'BR(ZZ) %':<12} {'BR(tt) %':<12} {'BR(hh) %':<12}")
print("-"*70)
for m in [2.0, 2.26, 2.44, 3.0, 3.5]:
    idx = np.argmin(np.abs(m_rho - m))
    print(f"{m:<12.2f} {BR_WW[idx]:<12.1f} {BR_ZZ[idx]:<12.1f} {BR_tt[idx]:<12.1f} {BR_hh[idx]:<12.1f}")
print("="*70)
print("Note: BR(WW) > BR(hh) distinguishes TSQVT from radion/dilaton scenarios")