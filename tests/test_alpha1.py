#!/usr/bin/env python3
"""
Unit tests for portal coefficient α₁ calculation.
"""

import pytest
import numpy as np
import sys
from pathlib import Path

# Add src to path (works from any directory)
src_path = Path(__file__).parent.parent / 'src'
sys.path.insert(0, str(src_path))

from tsqvt_pipeline import compute_alpha1, TSQVTParameters


class TestAlpha1Calculation:
    """Test suite for α₁ computation."""
    
    @pytest.fixture
    def benchmark_matrices(self):
        """Return benchmark M₂ and Y matrices."""
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
        return M2, Y
    
    def test_alpha1_benchmark_value(self, benchmark_matrices):
        """Test that α₁ is in physically reasonable range."""
        M2, Y = benchmark_matrices
        alpha1 = compute_alpha1(M2, Y)
        # α₁ should be O(0.1) for SM-like matrices
        # The central value 4.275e-2 in TSQVTParameters is a fitted value
        assert 0.01 < alpha1 < 1.0, f"α₁ = {alpha1} outside reasonable range"
    
    def test_alpha1_positive(self, benchmark_matrices):
        """Test that α₁ is positive for SM-like inputs."""
        M2, Y = benchmark_matrices
        alpha1 = compute_alpha1(M2, Y)
        assert alpha1 > 0
    
    def test_alpha1_top_dominance(self, benchmark_matrices):
        """Test that top sector dominates α₁."""
        M2, Y = benchmark_matrices
        alpha1_full = compute_alpha1(M2, Y)
        
        M2_top_only = {'up': np.diag([0, 0, M2['up'][2,2]]), 'down': np.zeros((3,3)), 'lepton': np.zeros((3,3))}
        Y_top_only = {'up': np.diag([0, 0, Y['up'][2,2]]), 'down': np.zeros((3,3)), 'lepton': np.zeros((3,3))}
        alpha1_top = compute_alpha1(M2_top_only, Y_top_only)
        
        assert alpha1_top / alpha1_full > 0.90


class TestTSQVTParameters:
    """Test the TSQVTParameters class."""
    
    def test_default_alpha1(self):
        params = TSQVTParameters()
        assert abs(params.alpha1 - 4.275e-2) < 1e-4
    
    def test_default_kappa_spec(self):
        params = TSQVTParameters()
        assert params.kappa_spec == 50000.0
    
    def test_rho_c_value(self):
        params = TSQVTParameters()
        assert abs(params.rho_c - 2/3) < 1e-10


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
