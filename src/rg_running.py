#!/usr/bin/env python3
"""
Renormalization Group Running for TSQVT ρ-Higgs Portal
======================================================

This module provides functions for computing the RG evolution of
gauge and scalar couplings in the extended Higgs sector.

Author: Kerym Makraini
Date: December 2025
"""

import numpy as np
from scipy.integrate import solve_ivp
from typing import Tuple, Dict, Optional
import warnings

# Physical constants
v_EW = 246.22    # GeV
m_t = 172.69     # GeV
m_Z = 91.1876    # GeV

# SM gauge coupling values at m_Z
alpha1_mZ = 0.01017  # U(1)_Y
alpha2_mZ = 0.03381  # SU(2)_L
alpha3_mZ = 0.1179   # SU(3)_c

# Convert to g couplings: g² = 4π α
g1_mZ = np.sqrt(4 * np.pi * alpha1_mZ) * np.sqrt(5/3)  # GUT normalized
g2_mZ = np.sqrt(4 * np.pi * alpha2_mZ)
g3_mZ = np.sqrt(4 * np.pi * alpha3_mZ)

# Top Yukawa at m_t
y_t_mt = m_t / v_EW

# SM Higgs quartic at m_Z
lambda_H_mZ = 0.13


def beta_gauge_1loop(g: np.ndarray) -> np.ndarray:
    """One-loop gauge beta functions (SM)."""
    b = np.array([41.0/10.0, -19.0/6.0, -7.0])
    return b * g**3 / (16 * np.pi**2)


def beta_gauge_2loop(g: np.ndarray) -> np.ndarray:
    """Two-loop gauge beta functions (SM)."""
    b = np.array([41.0/10.0, -19.0/6.0, -7.0])
    beta_1loop = b * g**3 / (16 * np.pi**2)
    
    B = np.array([
        [199.0/50.0,  27.0/10.0,  44.0/5.0],
        [9.0/10.0,    35.0/6.0,   12.0],
        [11.0/10.0,   9.0/2.0,    -26.0]
    ])
    
    beta_2loop = np.zeros(3)
    for i in range(3):
        for j in range(3):
            beta_2loop[i] += B[i,j] * g[i]**3 * g[j]**2
    beta_2loop /= (16 * np.pi**2)**2
    
    return beta_1loop + beta_2loop


def beta_lambda_H(lambda_H: float, lambda_rhoH: float,
                   g1: float, g2: float, y_t: float) -> float:
    """One-loop beta function for Higgs quartic λ_H."""
    scalar = 24 * lambda_H**2 + 2 * lambda_rhoH**2
    yukawa = -6 * y_t**4 + 12 * y_t**2 * lambda_H
    gauge_quartic = (3.0/8.0) * (2*g2**4 + (g1**2 + g2**2)**2)
    gauge_linear = lambda_H * (-9*g2**2/2 - 3*g1**2/2)
    
    return (scalar + yukawa + gauge_quartic + gauge_linear) / (16 * np.pi**2)


def beta_lambda_rho(lambda_rho: float, lambda_rhoH: float) -> float:
    """One-loop beta function for ρ quartic λ_ρ."""
    return (20 * lambda_rho**2 + 2 * lambda_rhoH**2) / (16 * np.pi**2)


def beta_lambda_rhoH(lambda_rhoH: float, lambda_H: float,
                      lambda_rho: float, g1: float, g2: float, y_t: float) -> float:
    """One-loop beta function for portal coupling λ_ρH."""
    scalar = 4 * lambda_rhoH**2 + 12 * lambda_H * lambda_rhoH + 8 * lambda_rho * lambda_rhoH
    yukawa = 6 * y_t**2 * lambda_rhoH
    gauge = -(9*g2**2/4 + 3*g1**2/4) * lambda_rhoH
    
    return (scalar + yukawa + gauge) / (16 * np.pi**2)


def beta_top_yukawa(y_t: float, g1: float, g2: float, g3: float) -> float:
    """One-loop beta function for top Yukawa."""
    yukawa = (9.0/2.0) * y_t**3
    gauge = -y_t * (17*g1**2/20 + 9*g2**2/4 + 8*g3**2)
    return (yukawa + gauge) / (16 * np.pi**2)


def run_gauge_couplings(mu_init: float, mu_final: float,
                         g_init: np.ndarray = None,
                         use_2loop: bool = False,
                         n_points: int = 100) -> Dict:
    """Integrate gauge coupling RGEs."""
    if g_init is None:
        g_init = np.array([g1_mZ, g2_mZ, g3_mZ])
    
    t_init = 0.0
    t_final = np.log(mu_final / mu_init)
    
    def rhs(t, g):
        return beta_gauge_2loop(g) if use_2loop else beta_gauge_1loop(g)
    
    t_eval = np.linspace(t_init, t_final, n_points)
    sol = solve_ivp(rhs, [t_init, t_final], g_init, 
                    t_eval=t_eval, method='RK45', rtol=1e-8, atol=1e-10)
    
    mu = mu_init * np.exp(sol.t)
    
    return {
        'mu': mu, 't': sol.t,
        'g1': sol.y[0], 'g2': sol.y[1], 'g3': sol.y[2],
        'alpha1': sol.y[0]**2 / (4 * np.pi) * 3/5,
        'alpha2': sol.y[1]**2 / (4 * np.pi),
        'alpha3': sol.y[2]**2 / (4 * np.pi)
    }


def run_quartic_couplings(mu_init: float, mu_final: float,
                           lambda_init: np.ndarray,
                           gauge_result: Dict = None,
                           n_points: int = 100) -> Dict:
    """Integrate scalar quartic RGEs."""
    from scipy.interpolate import interp1d
    
    if gauge_result is None:
        gauge_result = run_gauge_couplings(mu_init, mu_final, n_points=n_points)
    
    g1_func = interp1d(gauge_result['mu'], gauge_result['g1'], 
                       kind='cubic', fill_value='extrapolate')
    g2_func = interp1d(gauge_result['mu'], gauge_result['g2'],
                       kind='cubic', fill_value='extrapolate')
    g3_func = interp1d(gauge_result['mu'], gauge_result['g3'],
                       kind='cubic', fill_value='extrapolate')
    
    t_init = 0.0
    t_final = np.log(mu_final / mu_init)
    
    def rhs(t, couplings):
        lambda_H, lambda_rho, lambda_rhoH, y_t = couplings
        mu = mu_init * np.exp(t)
        g1, g2, g3 = g1_func(mu), g2_func(mu), g3_func(mu)
        
        return [
            beta_lambda_H(lambda_H, lambda_rhoH, g1, g2, y_t),
            beta_lambda_rho(lambda_rho, lambda_rhoH),
            beta_lambda_rhoH(lambda_rhoH, lambda_H, lambda_rho, g1, g2, y_t),
            beta_top_yukawa(y_t, g1, g2, g3)
        ]
    
    t_eval = np.linspace(t_init, t_final, n_points)
    sol = solve_ivp(rhs, [t_init, t_final], lambda_init,
                    t_eval=t_eval, method='RK45', rtol=1e-8, atol=1e-10)
    
    mu = mu_init * np.exp(sol.t)
    stability = check_vacuum_stability_array(sol.y[0], sol.y[1], sol.y[2])
    
    return {
        'mu': mu, 't': sol.t,
        'lambda_H': sol.y[0], 'lambda_rho': sol.y[1],
        'lambda_rhoH': sol.y[2], 'y_t': sol.y[3],
        'stability': stability
    }


def check_vacuum_stability(lambda_H: float, lambda_rho: float, 
                            lambda_rhoH: float) -> Tuple[bool, str]:
    """Check vacuum stability conditions."""
    conditions = []
    if lambda_H <= 0:
        conditions.append(f"λ_H = {lambda_H:.4f} ≤ 0")
    if lambda_rho <= 0:
        conditions.append(f"λ_ρ = {lambda_rho:.4f} ≤ 0")
    if lambda_H > 0 and lambda_rho > 0:
        if lambda_H * lambda_rho < lambda_rhoH**2 / 4:
            conditions.append(f"λ_H λ_ρ < λ_ρH²/4")
    
    is_stable = len(conditions) == 0
    return is_stable, "Stable" if is_stable else "; ".join(conditions)


def check_vacuum_stability_array(lambda_H: np.ndarray, lambda_rho: np.ndarray,
                                   lambda_rhoH: np.ndarray) -> np.ndarray:
    """Vectorized stability check."""
    return (lambda_H > 0) & (lambda_rho > 0) & (lambda_H * lambda_rho > lambda_rhoH**2 / 4)


def check_perturbativity(couplings: Dict, bound: float = 4*np.pi) -> Tuple[bool, str]:
    """Check perturbativity bounds."""
    violations = []
    for name in ['lambda_H', 'lambda_rho', 'lambda_rhoH', 'y_t']:
        if name in couplings:
            max_val = np.max(np.abs(couplings[name]))
            if max_val > bound:
                violations.append(f"|{name}|_max = {max_val:.2f}")
    
    is_pert = len(violations) == 0
    return is_pert, "Perturbative" if is_pert else "; ".join(violations)


def compute_sin2_theta_W(g1: float, g2: float) -> float:
    """Compute sin²θ_W from gauge couplings."""
    g1_prime_sq = (3.0/5.0) * g1**2
    return g1_prime_sq / (g2**2 + g1_prime_sq)


if __name__ == "__main__":
    print("RG Running Module Test")
    print("=" * 60)
    
    result = run_gauge_couplings(m_Z, 1e16, n_points=50)
    print(f"\nGauge running from m_Z to 10^16 GeV:")
    print(f"g1(m_Z) = {result['g1'][0]:.4f} -> g1(10^16) = {result['g1'][-1]:.4f}")
    print(f"g2(m_Z) = {result['g2'][0]:.4f} -> g2(10^16) = {result['g2'][-1]:.4f}")
    print(f"g3(m_Z) = {result['g3'][0]:.4f} -> g3(10^16) = {result['g3'][-1]:.4f}")
