"""
Scaling Analysis for Yang-Mills Mass Gap

This module implements scaling analysis across multiple beta values to:
1. Verify approach to the continuum limit
2. Extract physical quantities in the a -> 0 limit
3. Confirm universal behavior

Author: Yang-Mills Mass Gap Project
Date: January 2026
"""

import numpy as np
from typing import Tuple, List, Optional, Dict
from dataclasses import dataclass, asdict
import json
import os
import sys
import time

# Add parent directory for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import lattice_gauge as lg
import wilson_loops as wl
import mass_gap as mg
import correlation_functions as cf


# =============================================================================
# Lattice Spacing from Beta
# =============================================================================

def lattice_spacing_two_loop(beta: float, N: int = 3) -> float:
    """
    Compute lattice spacing from beta using two-loop beta function.

    a(beta) = (1/Lambda) * exp(-beta/(12*b0*N)) * (corrections)

    Parameters
    ----------
    beta : float
        Inverse coupling.
    N : int
        Number of colors (default: 3 for SU(3)).

    Returns
    -------
    float
        Lattice spacing in units of 1/Lambda_lat.
    """
    b0 = 11.0 / (16.0 * np.pi**2)
    b1 = 102.0 / (16.0 * np.pi**2)**2

    # Leading order
    log_a = -beta / (12.0 * b0 * N)

    # Next-to-leading order correction
    g2 = 6.0 / beta
    log_a += (b1 / (2.0 * b0**2)) * np.log(b0 * g2)

    return np.exp(log_a)


def lattice_spacing_fm(beta: float, sqrt_sigma_phys: float = 0.44) -> float:
    """
    Estimate lattice spacing in fm using physical string tension.

    Uses: sqrt(sigma) = 440 MeV, hbar*c = 197 MeV*fm

    Parameters
    ----------
    beta : float
        Inverse coupling.
    sqrt_sigma_phys : float
        Physical sqrt(sigma) in GeV.

    Returns
    -------
    float
        Lattice spacing in fm.
    """
    # Empirical fit for SU(3) (Necco & Sommer 2002)
    # a * sqrt(sigma) as function of beta
    r0_a = np.exp(-1.6805 - 1.7139 * (beta - 6.0) + 0.8155 * (beta - 6.0)**2)

    # r0 ~ 0.5 fm, so a ~ r0_a * 0.5 fm
    a_fm = r0_a * 0.5

    return a_fm


# =============================================================================
# Scaling Result Data Structures
# =============================================================================

@dataclass
class SingleBetaResult:
    """Results from simulation at a single beta value."""
    beta: float
    lattice_size: Tuple[int, int, int, int]
    n_configs: int

    # Plaquette
    plaquette: float
    plaquette_err: float

    # String tension
    sigma: float
    sigma_err: float

    # Mass gap from string tension
    mass_gap_string: float
    mass_gap_string_err: float

    # Mass gap from correlator (if available)
    mass_gap_corr: Optional[float] = None
    mass_gap_corr_err: Optional[float] = None

    # Lattice spacing estimate
    a_fm: Optional[float] = None


@dataclass
class ScalingAnalysisResult:
    """Results from scaling analysis across multiple beta values."""
    beta_values: List[float]
    results: List[SingleBetaResult]

    # Physical quantities (in appropriate units)
    sigma_phys: Optional[float] = None
    sigma_phys_err: Optional[float] = None
    mass_gap_phys: Optional[float] = None
    mass_gap_phys_err: Optional[float] = None

    # Scaling quality
    scaling_chi2: Optional[float] = None


# =============================================================================
# Single Beta Simulation
# =============================================================================

def run_single_beta(beta: float,
                    lattice_size: Tuple[int, int, int, int],
                    n_therm: int = 100,
                    n_configs: int = 50,
                    n_sep: int = 10,
                    verbose: bool = True) -> SingleBetaResult:
    """
    Run simulation at a single beta value.

    Parameters
    ----------
    beta : float
        Inverse coupling.
    lattice_size : Tuple
        (Nx, Ny, Nz, Nt) lattice dimensions.
    n_therm : int
        Thermalization sweeps.
    n_configs : int
        Number of configurations.
    n_sep : int
        Sweeps between configurations.
    verbose : bool
        Print progress.

    Returns
    -------
    SingleBetaResult
        Results for this beta value.
    """
    Nx, Ny, Nz, Nt = lattice_size

    if verbose:
        print(f"\n{'='*60}")
        print(f"Beta = {beta}, Lattice = {Nx}x{Ny}x{Nz}x{Nt}")
        print(f"{'='*60}")

    # Setup
    config = lg.LatticeConfig(Nx=Nx, Ny=Ny, Nz=Nz, Nt=Nt, beta=beta)
    gauge = lg.GaugeField(config)
    gauge.hot_start()

    # Thermalize
    if verbose:
        print(f"Thermalizing ({n_therm} sweeps)...")

    for i in range(n_therm):
        lg.sweep(gauge, lg.UpdateMethod.METROPOLIS)
        if verbose and (i + 1) % 20 == 0:
            plaq = lg.compute_average_plaquette(gauge)
            print(f"  Sweep {i+1}: <P> = {plaq:.4f}")

    # Generate configurations and measure
    if verbose:
        print(f"Generating {n_configs} configurations...")

    plaq_values = []
    wilson_data_all = []

    for i in range(n_configs):
        # Update
        for _ in range(n_sep):
            lg.sweep(gauge, lg.UpdateMethod.METROPOLIS)

        # Measure plaquette
        plaq = lg.compute_average_plaquette(gauge)
        plaq_values.append(plaq)

        # Measure Wilson loops
        analyzer = wl.WilsonLoopAnalyzer(gauge)
        R_max = min(Nx // 2, 4)
        T_max = min(Nt // 2, 6)
        analyzer.measure_all(R_max, T_max, sample_fraction=0.5, verbose=False)
        wilson_data_all.append(analyzer.get_wilson_loop_dict())

        if verbose and (i + 1) % 10 == 0:
            print(f"  Config {i+1}: <P> = {plaq:.4f}")

    # Analyze
    if verbose:
        print("Analyzing...")

    # Plaquette statistics
    plaq_mean = np.mean(plaq_values)
    plaq_err = np.std(plaq_values) / np.sqrt(len(plaq_values))

    # Average Wilson loops
    avg_wilson = {}
    for key in wilson_data_all[0].keys():
        values = [w[key] for w in wilson_data_all]
        avg_wilson[key] = np.mean(values)

    # Extract potential and string tension
    T_values = list(range(2, T_max + 1))
    R, V, V_err = mg.extract_potential_multiT(avg_wilson, T_values)

    fitter = mg.StaticPotentialFitter(R, V, V_err)
    try:
        fit = fitter.fit_linear(R_min=1.0)
        sigma = fit.sigma
        sigma_err = fit.sigma_err
    except Exception:
        sigma = np.nan
        sigma_err = np.nan

    # Mass gap from string tension
    if sigma > 0:
        mass_gap = np.sqrt(sigma)
        mass_gap_err = sigma_err / (2 * mass_gap)
    else:
        mass_gap = np.nan
        mass_gap_err = np.nan

    # Lattice spacing
    a_fm = lattice_spacing_fm(beta)

    result = SingleBetaResult(
        beta=beta,
        lattice_size=lattice_size,
        n_configs=n_configs,
        plaquette=plaq_mean,
        plaquette_err=plaq_err,
        sigma=sigma,
        sigma_err=sigma_err,
        mass_gap_string=mass_gap,
        mass_gap_string_err=mass_gap_err,
        a_fm=a_fm
    )

    if verbose:
        print(f"\nResults for beta = {beta}:")
        print(f"  <P> = {plaq_mean:.4f} +/- {plaq_err:.4f}")
        print(f"  sigma = {sigma:.4f} +/- {sigma_err:.4f}")
        print(f"  mass_gap = {mass_gap:.4f} +/- {mass_gap_err:.4f}")
        print(f"  a ~ {a_fm:.4f} fm")

    return result


# =============================================================================
# Multi-Beta Scaling Analysis
# =============================================================================

def run_scaling_analysis(beta_values: List[float],
                          lattice_size: Tuple[int, int, int, int],
                          n_therm: int = 100,
                          n_configs: int = 30,
                          n_sep: int = 10,
                          verbose: bool = True) -> ScalingAnalysisResult:
    """
    Run scaling analysis across multiple beta values.

    Parameters
    ----------
    beta_values : List[float]
        List of beta values to simulate.
    lattice_size : Tuple
        Lattice dimensions (same for all beta).
    n_therm, n_configs, n_sep : int
        Simulation parameters.
    verbose : bool
        Print progress.

    Returns
    -------
    ScalingAnalysisResult
        Combined scaling analysis results.
    """
    if verbose:
        print("\n" + "=" * 70)
        print("SCALING ANALYSIS")
        print("=" * 70)
        print(f"Beta values: {beta_values}")
        print(f"Lattice: {lattice_size}")

    results = []

    for beta in beta_values:
        result = run_single_beta(
            beta, lattice_size,
            n_therm=n_therm, n_configs=n_configs, n_sep=n_sep,
            verbose=verbose
        )
        results.append(result)

    # Scaling analysis
    if verbose:
        print("\n" + "=" * 70)
        print("SCALING SUMMARY")
        print("=" * 70)

    # Check scaling: a^2 * sigma should be constant
    a2_sigma = []
    for r in results:
        if r.sigma > 0 and r.a_fm is not None:
            a2_sigma.append(r.a_fm**2 * r.sigma)

    if len(a2_sigma) > 1:
        sigma_phys = np.mean(a2_sigma)
        sigma_phys_err = np.std(a2_sigma)
        scaling_chi2 = np.sum((np.array(a2_sigma) - sigma_phys)**2) / sigma_phys**2
    else:
        sigma_phys = None
        sigma_phys_err = None
        scaling_chi2 = None

    # Physical mass gap
    if sigma_phys is not None and sigma_phys > 0:
        mass_gap_phys = np.sqrt(sigma_phys) / 0.197  # Convert to GeV (hbar*c = 0.197 GeV*fm)
        mass_gap_phys_err = sigma_phys_err / (2 * np.sqrt(sigma_phys)) / 0.197
    else:
        mass_gap_phys = None
        mass_gap_phys_err = None

    if verbose:
        print("\nScaling table:")
        print(f"{'Beta':>6} {'a (fm)':>8} {'sigma':>8} {'a^2*sigma':>10} {'Delta':>8}")
        print("-" * 50)
        for r in results:
            a2s = r.a_fm**2 * r.sigma if r.a_fm and r.sigma > 0 else np.nan
            print(f"{r.beta:6.2f} {r.a_fm:8.4f} {r.sigma:8.4f} {a2s:10.4f} {r.mass_gap_string:8.4f}")

        if sigma_phys is not None:
            print(f"\nPhysical string tension: sqrt(sigma) = {np.sqrt(sigma_phys)/0.197:.0f} MeV")
            print(f"Physical mass gap: Delta = {mass_gap_phys:.0f} +/- {mass_gap_phys_err:.0f} MeV")
            print(f"Scaling chi^2: {scaling_chi2:.2f}")

    return ScalingAnalysisResult(
        beta_values=beta_values,
        results=results,
        sigma_phys=sigma_phys,
        sigma_phys_err=sigma_phys_err,
        mass_gap_phys=mass_gap_phys,
        mass_gap_phys_err=mass_gap_phys_err,
        scaling_chi2=scaling_chi2
    )


# =============================================================================
# Save and Load Results
# =============================================================================

def save_scaling_results(results: ScalingAnalysisResult, filepath: str) -> None:
    """Save scaling analysis results to JSON."""
    data = {
        'beta_values': results.beta_values,
        'sigma_phys': results.sigma_phys,
        'sigma_phys_err': results.sigma_phys_err,
        'mass_gap_phys': results.mass_gap_phys,
        'mass_gap_phys_err': results.mass_gap_phys_err,
        'scaling_chi2': results.scaling_chi2,
        'results': [asdict(r) for r in results.results]
    }

    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2)

    print(f"Results saved to {filepath}")


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("YANG-MILLS SCALING ANALYSIS")
    print("=" * 70)

    # Configuration
    beta_values = [5.7, 6.0, 6.2]
    lattice_size = (4, 4, 4, 8)
    n_configs = 20  # Reduced for quick test

    # Run analysis
    results = run_scaling_analysis(
        beta_values=beta_values,
        lattice_size=lattice_size,
        n_therm=50,
        n_configs=n_configs,
        n_sep=5,
        verbose=True
    )

    # Save results
    os.makedirs("results", exist_ok=True)
    save_scaling_results(results, "results/scaling_analysis.json")

    print("\n" + "=" * 70)
    print("Analysis complete!")
    print("=" * 70)
