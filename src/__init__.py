"""
TSQVT ρ-Higgs Portal: Core Modules
==================================

This package provides the numerical implementation for the TSQVT ρ-Higgs portal
phenomenology as described in the manuscript.

Modules:
--------
- branching_ratios: Partial widths and branching ratios for h₂
- cross_sections: Production cross sections (ggF, VBF)
- tsqvt_pipeline: Complete benchmark computation pipeline

Usage:
------
>>> from src.tsqvt_pipeline import compute_benchmark, SpectralParameters
>>> params = SpectralParameters()
>>> B1 = compute_benchmark("B1", 1.59, 2.26, params)
>>> print(f"θ = {B1.theta_deg:.1f}°, σ×BR(WW) = {B1.sigma_BR_WW_fb:.2f} fb")

Notes:
------
The branching ratios use MANUSCRIPT-CALIBRATED values to ensure consistency
with the "narrow resonance" regime (Γ/m ~ few %) of the paper.
See README_BRANCHING_RATIOS.md for technical details.
"""

from .tsqvt_pipeline import (
    SpectralParameters,
    TSQVTParameters,  # Alias for compatibility
    compute_benchmark,
    compute_benchmark_dict,
    compute_alpha1,
    compute_delta_m2,
    compute_mixing_angle,
    compute_delta_kappa,
    compute_branching_ratios,
    BenchmarkResult,
)

from .cross_sections import (
    sigma_ggF,
    sigma_VBF,
    sigma_SM_heavy,
    sigma_ggF_SM_heavy,  # Alias for compatibility
    compute_cross_sections,
    compute_total_xsec,
)

from .branching_ratios import (
    BranchingRatios,
    compute_branching_ratios as compute_BR_from_first_principles,
    get_manuscript_BR,
)

__version__ = "2.3"
__author__ = "Kerym Makraini"
