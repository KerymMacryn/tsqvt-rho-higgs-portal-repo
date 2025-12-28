#!/usr/bin/env python3
"""
Unit tests for cross section calculations.
"""

import pytest
import numpy as np
import sys
from pathlib import Path

# Add src to path (works from any directory)
src_path = Path(__file__).parent.parent / 'src'
sys.path.insert(0, str(src_path))

from cross_sections import sigma_ggF_SM_heavy, sigma_ggF, sigma_VBF, compute_total_xsec
from tsqvt_pipeline import compute_benchmark, TSQVTParameters


class TestGluonFusion:
    """Tests for ggF cross sections."""
    
    def test_sigma_ggF_SM_decreasing(self):
        """Test that σ_ggF decreases with mass."""
        masses = [500, 1000, 2000, 3000]
        sigmas = [sigma_ggF_SM_heavy(m) for m in masses]
        
        for i in range(len(sigmas) - 1):
            assert sigmas[i] > sigmas[i+1], "σ_ggF should decrease with mass"
    
    def test_sigma_ggF_positive(self):
        """Test that σ_ggF is positive."""
        for m in [500, 1000, 2000, 3000]:
            sigma = sigma_ggF_SM_heavy(m)
            assert sigma > 0
    
    def test_sigma_ggF_mixing_suppression(self):
        """Test that mixing suppresses ggF."""
        m_rho = 2000
        sin2_theta = 0.04  # ~11 degrees
        
        sigma_SM = sigma_ggF_SM_heavy(m_rho)
        sigma_mix = sigma_ggF(m_rho, sin2_theta)
        
        assert sigma_mix == pytest.approx(sin2_theta * sigma_SM, rel=0.01)


class TestVBF:
    """Tests for VBF cross sections."""
    
    def test_sigma_VBF_benchmark(self):
        """Test VBF at benchmark mass."""
        sigma = sigma_VBF(2300)
        assert 0.001 < sigma < 0.01  # ~3 fb
    
    def test_sigma_VBF_mass_scaling(self):
        """Test that σ_VBF scales as 1/m²."""
        sigma_2TeV = sigma_VBF(2000)
        sigma_4TeV = sigma_VBF(4000)
        
        # Should scale as (2000/4000)² = 0.25
        ratio = sigma_4TeV / sigma_2TeV
        assert 0.2 < ratio < 0.35


class TestTotalCrossSection:
    """Tests for total cross section."""
    
    def test_VBF_dominance(self):
        """Test VBF > ggF in viable region."""
        xsec = compute_total_xsec(2300, 0.035)  # Benchmark
        
        assert xsec['sigma_VBF_pb'] > xsec['sigma_ggF_pb'], "VBF should dominate"
        assert xsec['R_VBF_ggF'] > 10, "R should be > 10"
    
    def test_total_is_sum(self):
        """Test σ_tot = σ_ggF + σ_VBF."""
        xsec = compute_total_xsec(2000, 0.04)
        
        expected_tot = xsec['sigma_ggF_pb'] + xsec['sigma_VBF_pb']
        assert xsec['sigma_tot_pb'] == pytest.approx(expected_tot, rel=0.001)


class TestBenchmarkIntegration:
    """Integration tests with full benchmark calculation."""
    
    def test_benchmark_B1(self):
        """Test benchmark B1 values."""
        params = TSQVTParameters()
        result = compute_benchmark(1590, 2260, "B1", params)
        
        # Check Δκ constraint
        assert abs(result['Delta_kappa_percent']) < 2.0
        
        # Check VBF dominance (R > 5 is sufficient for VBF dominance)
        assert result['R_VBF_ggF'] > 5
        
        # Check signal rate is detectable
        assert 0.1 < result['sigma_x_BR_WW'] < 5.0
    
    def test_all_benchmarks_viable(self):
        """Test all benchmarks satisfy viability."""
        params = TSQVTParameters()
        
        benchmarks = [(1590, 2260), (1500, 2260), (1680, 2440)]
        
        for Lambda, m_rho in benchmarks:
            result = compute_benchmark(Lambda, m_rho, params=params)
            assert abs(result['Delta_kappa_percent']) < 2.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
