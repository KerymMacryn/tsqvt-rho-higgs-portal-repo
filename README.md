# TSQVT ρ-Higgs Portal: Numerical Repository

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18074977.svg)](https://doi.org/10.5281/zenodo.18074977)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Companion repository for:**

> K. Makraini, "The ρ-Higgs Portal in Twistorial Spectral Quantum Vacuum Theory: Predictions for TeV-Scale Collider Searches," *Physical Review D* (2025, submitted).

---

## Overview

This repository provides the complete numerical pipeline for reproducing the phenomenological predictions of the TSQVT ρ-Higgs portal, including:

- Mixing angle θ and Higgs coupling deviation Δκ
- Production cross sections (ggF, VBF) and their ratio R
- Branching ratios and signal rates σ×BR
- Parameter space scans with spectral moment sensitivity (f_k ±30%)
- Benchmark validation (B1, B2, B3)
- Publication-quality figures

---

## Quick Start

```bash
# Clone the repository
git clone https://github.com/KerymMacryn/tsqvt-rho-higgs-portal-repo.git
cd tsqvt-rho-higgs-portal-repo

# Install dependencies
pip install -r requirements.txt

# Run the main pipeline (phenomenological mode)
python scripts/tsqvt_rho_higgs_pipeline.py --mode=pheno

# Run with microscopic derivation
python scripts/tsqvt_rho_higgs_pipeline.py --mode=micro

# Run parameter space scan
python scripts/scan_alpha1factor_1.0_res_2020_fixed.py
```

---

## Requirements

- Python ≥ 3.9
- NumPy, SciPy, Pandas, Matplotlib

```bash
pip install -r requirements.txt
```

For full development environment (including Jupyter notebooks and testing):
```bash
pip install -e ".[dev]"
```

---

## Repository Structure

```
tsqvt-rho-higgs-portal/
│
├── README.md                         # This file
├── LICENSE                           # MIT License
├── CITATION.cff                      # Citation metadata
├── requirements.txt                  # Python dependencies
├── pyproject.toml                    # Package configuration
│
├── scripts/                          # Executable scripts
│   ├── tsqvt_rho_higgs_pipeline.py       # Main pipeline (v2.3)
│   ├── scan_alpha1factor_1.0_res_2020_fixed.py  # Parameter scan (κ_spec=5×10⁴)
│   ├── kappa_sensitivity.py              # κ_spec sensitivity analysis
│   ├── Fig2_production_modes.py          # Figure 2 generator
│   ├── Fig3_branching_ratios.py          # Figure 3 generator
│   └── Fig_benchmark_vs_limits.py        # Benchmarks vs limits figure
│
├── src/                              # Core physics modules
│   ├── __init__.py
│   ├── tsqvt_pipeline.py                 # Benchmark computation engine
│   ├── cross_sections.py                 # σ_ggF, σ_VBF calculations
│   ├── branching_ratios.py               # Partial widths & BRs
│   ├── rg_running.py                     # RG evolution utilities
│   └── README_BRANCHING_RATIOS.md        # Technical documentation
│
├── data/                             # Calibrated Dirac matrices
│   ├── M0.npy                            # Constant mass matrix
│   ├── M_rho.npy                         # ρ-dependent mass matrix
│   ├── Y.npy                             # Yukawa matrix
│   ├── metadata.json                     # Calibration parameters
│   └── spectral_inputs.json              # Spectral action inputs
│
├── notebooks/                        # Jupyter analysis notebooks
│   ├── 01_parameter_scan.ipynb           # (Λ, m_ρ) plane exploration
│   ├── 02_benchmark_analysis.ipynb       # Benchmark point validation
│   ├── 03_figure_generation.ipynb        # Publication figures
│   └── 04_TSQVT_Master_Analysis.ipynb    # Complete analysis workflow
│
├── output/                           # Pipeline outputs (CSV, matrices)
│   ├── benchmarks_pheno.csv              # Phenomenological mode results
│   ├── benchmarks_micro.csv              # Microscopic mode results
│   └── benchmarks_v6.csv                 # Final manuscript values
│
├── results/                          # Scan results and sensitivity figures
│   ├── scan_Deltakappa_40x40_central.csv     # Central f_k scan
│   ├── scan_Deltakappa_40x40_fk_minus30pct.csv
│   ├── scan_Deltakappa_40x40_fk_plus30pct.csv
│   ├── Fig1_Deltakappa_contours*.pdf     # Contour plots
│   └── kappa_sensitivity.pdf             # κ_spec sensitivity figure
│
├── figures/                          # Publication-ready figures
│   ├── Fig1_Deltakappa_contours.pdf
│   ├── Fig2_production_modes.pdf
│   ├── Fig3_branching_ratios.pdf
│   └── Fig4_benchmark_vs_limits.pdf
│
├── benchmarks/                       # FeynRules/MadGraph model
│   └── TSQVTrhoPortal.fr                 # FeynRules model file
│
└── tests/                            # Unit test suite
    ├── __init__.py
    ├── conftest.py                       # pytest fixtures
    ├── run_tests.py                      # Manual test runner (18 tests)
    ├── test_alpha1.py                    # Portal coefficient tests
    └── test_cross_sections.py            # Cross section tests
```

---

## Pipeline Modes

### Phenomenological Mode (`--mode=pheno`)

Uses calibrated formula:
```
δm² = -α₁ κ_spec Λ² / (16π² × 33)
```

With α₁ = 4.275×10⁻² and κ_spec = 5×10⁴.

### Microscopic Mode (`--mode=micro`)

Extracts δm² from spectral action derivatives using calibrated Dirac matrices:
```
κ = m_h² / (∂²S_spec/∂h²)|₀
δm² = κ × (∂²S_spec/∂ρ∂h)|₀
```

### Comparison

| Observable | Manuscript | Pheno Mode | Micro Mode |
|------------|------------|------------|------------|
| θ (B1) | −11.1° | −11.1° | −11.0° |
| Δκ (B1) | −1.87% | −1.86% | −1.85% |
| R (B1) | 23 | 24 | 25 |
| σ×BR(WW) | 1.33 fb | 1.46 fb | 1.46 fb |

---

## Benchmark Points

| Point | Λ [TeV] | m_ρ [TeV] | θ [°] | Δκ [%] | R | σ×BR(WW) [fb] |
|-------|---------|-----------|-------|--------|---|---------------|
| **B1** | 1.59 | 2.26 | −11.1 | −1.87 | 24 | 1.46 |
| **B2** | 1.50 | 2.26 | −10.0 | −1.51 | 30 | 1.45 |
| **B3** | 1.68 | 2.44 | −10.7 | −1.72 | 27 | 1.19 |

All benchmarks satisfy |Δκ| < 2% (Higgs precision constraint).

---

## Key Physics

### VBF Dominance

The TSQVT predicts R = σ_VBF/σ_ggF ~ 20–30 because:
- **ggF** is suppressed by sin²θ ~ 0.03 (mixing-induced coupling)
- **VBF** has a mixing-**independent** contribution from the spectral contact term

This distinguishes TSQVT from generic singlet models (R < 1).

### Spectral Moment Sensitivity

Variations of f_k by ±30% produce:
| Observable | Variation |
|------------|-----------|
| m_ρ,opt | ±15% |
| θ | ±1% |
| Δκ | ±2% |
| R | <1% |

Core phenomenological signatures are **robust** to spectral uncertainties.

---

## Running Tests

### Manual Test Runner (no pytest dependency)

```bash
python tests/run_tests.py
```

### With pytest

```bash
pip install pytest
pytest tests/ -v
```

### Test Summary

| Module | Tests | Description |
|--------|-------|-------------|
| `test_alpha1.py` | 6 | Portal coefficient α₁ calculation |
| `test_cross_sections.py` | 12 | Cross sections (ggF, VBF, total) |
| `run_tests.py` | 18 | Complete integration suite |

**All 18 tests pass** with 100% consistency to manuscript values.

---

## FeynRules/MadGraph Integration

The `benchmarks/TSQVTrhoPortal.fr` file provides a FeynRules model for generating UFO output compatible with MadGraph5_aMC@NLO:

```mathematica
(* In Mathematica with FeynRules loaded *)
<< FeynRules`
LoadModel["SM.fr", "TSQVTrhoPortal.fr"]
WriteUFO[LagTSQVT, Output -> "TSQVT_rho_UFO"]
```

Copy the generated `TSQVT_rho_UFO/` directory to `MadGraph5/models/` for Monte Carlo event generation.

---

## Citation

If you use this code, please cite:

```bibtex
@article{Makraini2025TSQVT,
  author  = {Makraini, Kerym},
  title   = {The $\rho$-Higgs Portal in Twistorial Spectral Quantum Vacuum Theory: 
             Predictions for TeV-Scale Collider Searches},
  journal = {Physical Review D},
  year    = {2025},
  note    = {Submitted},
  doi     = {10.5281/zenodo.18074977}
}
```

See also the foundation paper:
```bibtex
@article{Makraini2025Foundation,
  author  = {Makraini, Kerym},
  title   = {Emergent Lorentzian Spacetime and Gauge Dynamics 
             from Twistorial Spectral Data},
  journal = {Next Research},
  pages   = {101114},
  year    = {2025},
  doi     = {10.1016/j.nexres.2025.101114}
}
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

## Contact

**Kerym Makraini**  
Department of Physics  
UNED – National University of Distance Education  
Madrid, Spain  
📧 mhamed34@alumno.uned.es

---

## Changelog

### v2.3 (2025-12-30)
- Fixed κ_spec = 50,000 (corrected from 200,000)
- Added f_k sensitivity scans (±30%)
- Added kappa_sensitivity.py and figure
- 18/18 tests passing
- Zenodo archival

### v2.0 (2025-12-27)
- Microscopic mode with calibrated Dirac matrices
- FeynRules model file
- Jupyter notebook suite

### v1.0 (2025-12-21)
- Initial phenomenological pipeline
- Benchmark points B1, B2, B3
