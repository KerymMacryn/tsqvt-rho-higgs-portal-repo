#!/usr/bin/env python3
"""
Manual test runner for TSQVT tests (pytest-free).
Run from repository root: python tests/run_tests.py
"""

import sys
from pathlib import Path

# Add src to path
src_path = Path(__file__).parent.parent / 'src'
sys.path.insert(0, str(src_path))

import numpy as np

# Track results
passed = 0
failed = 0
errors = []

def test(name, condition):
    global passed, failed, errors
    try:
        if condition:
            print(f"  ✓ {name}")
            passed += 1
        else:
            print(f"  ✗ {name}")
            failed += 1
            errors.append(f"{name}: assertion failed")
    except Exception as e:
        print(f"  ✗ {name} (ERROR: {e})")
        failed += 1
        errors.append(f"{name}: {e}")

def approx(a, b, rel=0.1):
    """Check if a ≈ b within rel tolerance."""
    if b == 0:
        return abs(a) < 1e-10
    return abs(a - b) / abs(b) < rel


print("=" * 70)
print("TSQVT Unit Tests (Manual Runner)")
print("=" * 70)

# =============================================================================
# TEST ALPHA1
# =============================================================================
print("\n[test_alpha1.py]")

from tsqvt_pipeline import compute_alpha1, TSQVTParameters, SpectralParameters

# Create benchmark matrices
M2 = {
    'up': np.diag([1e-6, 5e-4, 0.08]),
    'down': np.diag([2e-6, 3e-4, 0.02]),
    'lepton': np.diag([5e-7, 1e-4, 0.01])
}
Y = {
    'up': np.diag([1.27e-5, 7.38e-3, 0.995]),
    'down': np.diag([2.76e-5, 5.51e-4, 0.024]),
    'lepton': np.diag([2.87e-6, 6.07e-4, 0.0102])
}

alpha1 = compute_alpha1(M2, Y)
test("α₁ in reasonable range (0.01-1.0)", 0.01 < alpha1 < 1.0)
test("α₁ positive", alpha1 > 0)

# Test top dominance
M2_top = {'up': np.diag([0, 0, M2['up'][2,2]]), 'down': np.zeros((3,3)), 'lepton': np.zeros((3,3))}
Y_top = {'up': np.diag([0, 0, Y['up'][2,2]]), 'down': np.zeros((3,3)), 'lepton': np.zeros((3,3))}
alpha1_top = compute_alpha1(M2_top, Y_top)
test("Top dominance (>90%)", alpha1_top / alpha1 > 0.90)

# Test TSQVTParameters
params = TSQVTParameters()
test("Default α₁ = 4.275e-2", abs(params.alpha1 - 4.275e-2) < 1e-4)
test("Default κ_spec = 50000", params.kappa_spec == 50000.0)
test("ρ_c = 2/3", abs(params.rho_c - 2/3) < 1e-10)

# =============================================================================
# TEST CROSS SECTIONS
# =============================================================================
print("\n[test_cross_sections.py]")

from cross_sections import sigma_ggF_SM_heavy, sigma_ggF, sigma_VBF, compute_total_xsec
from tsqvt_pipeline import compute_benchmark

# Test ggF decreasing with mass
masses = [500, 1000, 2000, 3000]
sigmas = [sigma_ggF_SM_heavy(m) for m in masses]
test("σ_ggF decreases with mass", all(sigmas[i] > sigmas[i+1] for i in range(len(sigmas)-1)))
test("σ_ggF positive", all(s > 0 for s in sigmas))

# Test mixing suppression
m_rho = 2000
sin2_theta = 0.04
sigma_SM = sigma_ggF_SM_heavy(m_rho)
sigma_mix = sigma_ggF(m_rho, sin2_theta)
test("Mixing suppression σ_ggF = sin²θ × σ_SM", approx(sigma_mix, sin2_theta * sigma_SM, 0.01))

# Test VBF
sigma_vbf = sigma_VBF(2300)
test("σ_VBF at 2.3 TeV in range", 0.001 < sigma_vbf/1000 < 0.01)  # Convert fb to pb

# Test VBF scaling
sigma_2TeV = sigma_VBF(2000)
sigma_4TeV = sigma_VBF(4000)
ratio = sigma_4TeV / sigma_2TeV
test("VBF scales as ~1/m^2.5", 0.15 < ratio < 0.35)

# Test VBF dominance
xsec = compute_total_xsec(2300, 0.035)
test("VBF > ggF", xsec['sigma_VBF_pb'] > xsec['sigma_ggF_pb'])
test("R > 10", xsec['R_VBF_ggF'] > 10)

# Test total is sum
expected_tot = xsec['sigma_ggF_pb'] + xsec['sigma_VBF_pb']
test("σ_tot = σ_ggF + σ_VBF", approx(xsec['sigma_tot_pb'], expected_tot, 0.001))

# =============================================================================
# TEST BENCHMARK INTEGRATION
# =============================================================================
print("\n[Benchmark Integration]")

params = TSQVTParameters()
result = compute_benchmark(1590, 2260, "B1", params)

test("B1 |Δκ| < 2%", abs(result['Delta_kappa_percent']) < 2.0)
test("B1 R > 5 (VBF dominance)", result['R_VBF_ggF'] > 5)
test("B1 σ×BR(WW) detectable (0.1-5 fb)", 0.1 < result['sigma_x_BR_WW'] < 5.0)

# Test all benchmarks viable
benchmarks = [(1590, 2260), (1500, 2260), (1680, 2440)]
all_viable = True
for Lambda, m_rho in benchmarks:
    result = compute_benchmark(Lambda, m_rho, params=params)
    if abs(result['Delta_kappa_percent']) >= 2.0:
        all_viable = False
test("All benchmarks satisfy |Δκ| < 2%", all_viable)

# =============================================================================
# SUMMARY
# =============================================================================
print("\n" + "=" * 70)
total = passed + failed
print(f"Results: {passed}/{total} tests passed")

if failed > 0:
    print("\nFailed tests:")
    for e in errors:
        print(f"  - {e}")
    sys.exit(1)
else:
    print("\n✅ All tests passed!")
    sys.exit(0)
