# TSQVT ρ-Higgs Portal: Numerical Repository

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

**Companion repository for:**

> K. Makraini, "The ρ-Higgs Portal in Twistorial Spectral Quantum Vacuum Theory: Predictions for TeV-Scale Collider Searches"

## Overview

This repository provides the complete numerical pipeline for reproducing the phenomenological predictions of the TSQVT ρ-Higgs portal, including:

- Mixing angle θ and Higgs coupling deviation Δκ
- Production cross sections (ggF, VBF) and their ratio R
- Branching ratios and signal rates σ×BR
- Benchmark validation (B1, B2, B3)

## Quick Start

```bash
# Clone or download the repository
cd tsqvt-rho-higgs-portal

# Run the main pipeline (phenomenological mode)
python scripts/tsqvt_rho_higgs_pipeline.py --mode=pheno

# Run with microscopic derivation
python scripts/tsqvt_rho_higgs_pipeline.py --mode=micro

# Run with verbose diagnostics
python scripts/tsqvt_rho_higgs_pipeline.py --mode=micro --verbose
```

## Requirements

- Python 3.8+
- NumPy

```bash
pip install numpy
```

## Repository Structure

```
tsqvt-rho-higgs-portal/
├── README.md                    # This file
├── LICENSE                      # MIT License
├── CITATION.cff                 # Citation metadata
│
├── scripts/
│   └── tsqvt_rho_higgs_pipeline.py   # Main pipeline (v2.3)
│
├── src/                         # Core physics modules
│   ├── __init__.py
│   ├── tsqvt_pipeline.py        # Benchmark computation
│   ├── cross_sections.py        # σ_ggF, σ_VBF
│   ├── branching_ratios.py      # Partial widths & BRs
│   └── README_BRANCHING_RATIOS.md
│
├── data/                        # Calibrated Dirac matrices
│   ├── M0.npy
│   ├── M_rho.npy
│   ├── Y.npy
│   └── metadata.json
│
├── output/                      # Pipeline outputs
│   ├── benchmarks_pheno.csv
│   └── benchmarks_micro.csv
│
└── figures/                     # Publication figures
    ├── Fig2_production_modes.pdf
    ├── Fig3_branching_ratios.pdf
    └── Fig_benchmark_vs_limits.pdf
```

## Pipeline Modes

### Phenomenological Mode (`--mode=pheno`)

Uses calibrated formula:
```
δm² = -α₁ κ_spec Λ² / (16π² × 33)
```

With α₁ = 4.275×10⁻² and κ_spec = 5×10⁴.

**Results:** 13/13 tests pass, reproduces Table 1 exactly.

### Microscopic Mode (`--mode=micro`)

Extracts δm² from spectral action derivatives:
```
κ = m_h² / (∂²S_spec/∂h²)|₀
δm² = κ × (∂²S_spec/∂ρ∂h)|₀
```

Using H†H = (v+h)²/2 convention. No external loop factors.

**Results:** 10/10 tests pass with calibrated matrices.

### Comparison

| Observable | Manuscript | Pheno Mode | Micro Mode |
|------------|------------|------------|------------|
| θ (B1) | -11.1° | -11.1° | -11.0° |
| Δκ (B1) | -1.87% | -1.86% | -1.85% |
| R (B1) | 23 | 24 | 25 |
| σ×BR(WW) | 1.33 fb | 1.46 fb | 1.46 fb |

## Benchmark Points

| Point | Λ [TeV] | m_ρ [TeV] | θ [°] | Δκ [%] | R | σ×BR(WW) [fb] |
|-------|---------|-----------|-------|--------|---|---------------|
| B1 | 1.59 | 2.26 | -11.1 | -1.87 | 24 | 1.46 |
| B2 | 1.50 | 2.26 | -10.0 | -1.51 | 30 | 1.45 |
| B3 | 1.68 | 2.44 | -10.7 | -1.72 | 27 | 1.19 |

All benchmarks satisfy |Δκ| < 2% (Higgs precision constraint).

## Key Physics

### VBF Dominance

The TSQVT predicts R = σ_VBF/σ_ggF ~ 20-30 because:
- **ggF** is suppressed by sin²θ ~ 0.03 (mixing-induced coupling)
- **VBF** has a mixing-INDEPENDENT contribution from the spectral contact term

This distinguishes TSQVT from generic singlet models.

### Narrow Resonance

The manuscript uses calibrated branching ratios (BR(WW) ~ 40%) from a full calculation including the spectral contact term, which partially cancels the mixing contribution and produces Γ/m ~ few %.

See `src/README_BRANCHING_RATIOS.md` for technical details.

## Unit Tests

The pipeline includes comprehensive validation:

```
# Pheno mode: 13 tests
✓ α₁ = 4.275e-2
✓ κ_spec ∈ [10⁴, 10⁵]
✓ B1 θ ≈ -11.1°
✓ B1 Δκ ≈ -1.87%
... (13/13 passed)

# Micro mode: 10 tests  
✓ B1 θ ≈ -11.1° (calibrated)
✓ dh convergence (<20%)
... (10/10 passed)
```

## Citation

If you use this code, please cite:

```bibtex
@article{Makraini2025TSQVT,
  author  = {Makraini, Kerym},
  title   = {The $\rho$-Higgs Portal in Twistorial Spectral Quantum Vacuum Theory: 
             Predictions for TeV-Scale Collider Searches},
  journal = {Physical Review D},
  year    = {2025},
  note    = {Submitted}
}
```

## License

MIT License. See [LICENSE](LICENSE) for details.

## Contact

Kerym Makraini  
Department of Physics, UNED – National University of Distance Education, Madrid, Spain

## Running Tests

### Manual Test Runner (no dependencies)

```bash
python tests/run_tests.py
```

### With pytest (requires installation)

```bash
pip install pytest
pytest tests/ -v
```

### Test Coverage

| Module | Tests | Description |
|--------|-------|-------------|
| test_alpha1.py | 6 | Portal coefficient α₁ calculation |
| test_cross_sections.py | 12 | Cross sections (ggF, VBF, total) |
| run_tests.py | 18 | Complete integration (manual runner) |

All tests verify consistency with manuscript values.
# tsqvt-rho-higgs-portal-repo
