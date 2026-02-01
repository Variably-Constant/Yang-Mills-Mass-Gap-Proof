"""
Correlation Function Measurements for Yang-Mills Mass Gap

This module implements correlation function measurements for extracting
the mass gap from exponential decay:

    G(t) = <O(t) O(0)> ~ exp(-Delta * t)

The mass gap Delta is extracted from the large-t behavior.

Author: Yang-Mills Mass Gap Project
Date: January 2026
"""

import numpy as np
from typing import Tuple, List, Optional, Dict
from dataclasses import dataclass
from scipy.optimize import curve_fit
from scipy.stats import linregress
import sys
import os

# Add parent directory for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import lattice_gauge as lg
import wilson_loops as wl


# =============================================================================
# Glueball Operators
# =============================================================================

def plaquette_operator(gauge: lg.GaugeField, site: Tuple[int, int, int, int],
                       mu: int = 0, nu: int = 1) -> float:
    """
    Compute the plaquette operator at a given site and orientation.

    This is the simplest gauge-invariant local operator, corresponding
    to the 0++ glueball.

    Returns
    -------
    float
        Re Tr(U_p) / N
    """
    return lg.compute_plaquette_trace(gauge, site, mu, nu)


def spatial_plaquette_sum(gauge: lg.GaugeField, t: int) -> float:
    """
    Sum of all spatial plaquettes at time slice t.

    This operator has quantum numbers J^PC = 0++.

    Parameters
    ----------
    gauge : GaugeField
        The gauge configuration.
    t : int
        Time slice.

    Returns
    -------
    float
        Sum of Re Tr(U_p) / N over spatial plaquettes.
    """
    config = gauge.config
    total = 0.0

    # Sum over all spatial sites and spatial plaquette orientations
    for x in range(config.Nx):
        for y in range(config.Ny):
            for z in range(config.Nz):
                site = (x, y, z, t)
                # Spatial plaquettes: (0,1), (0,2), (1,2)
                for mu in range(3):
                    for nu in range(mu + 1, 3):
                        total += plaquette_operator(gauge, site, mu, nu)

    return total


def spacetime_plaquette_sum(gauge: lg.GaugeField, t: int, mu: int) -> float:
    """
    Sum of all space-time plaquettes at time slice t in direction mu-t.

    Parameters
    ----------
    gauge : GaugeField
        The gauge configuration.
    t : int
        Time slice.
    mu : int
        Spatial direction (0, 1, or 2).

    Returns
    -------
    float
        Sum of Re Tr(U_p) / N over space-time plaquettes.
    """
    config = gauge.config
    total = 0.0

    for x in range(config.Nx):
        for y in range(config.Ny):
            for z in range(config.Nz):
                site = (x, y, z, t)
                total += plaquette_operator(gauge, site, mu, 3)  # 3 = time direction

    return total


def polyakov_loop_operator(gauge: lg.GaugeField,
                           spatial_site: Tuple[int, int, int]) -> complex:
    """
    Compute the Polyakov loop at a spatial site.

    The Polyakov loop is an order parameter for confinement.

    Returns
    -------
    complex
        Tr(P) / N where P is the product of temporal links.
    """
    P = wl.compute_polyakov_loop(gauge, spatial_site)
    return np.trace(P) / 3.0


# =============================================================================
# Correlation Functions
# =============================================================================

@dataclass
class CorrelationResult:
    """Results from correlation function measurement."""
    t_values: np.ndarray
    G_values: np.ndarray
    G_errors: np.ndarray
    operator_type: str


def compute_plaquette_correlator(gauge: lg.GaugeField,
                                  t_max: Optional[int] = None) -> CorrelationResult:
    """
    Compute the plaquette-plaquette correlation function.

    G(t) = <O(t) O(0)> - <O>^2

    where O is the spatial plaquette sum.

    Parameters
    ----------
    gauge : GaugeField
        The gauge configuration.
    t_max : int, optional
        Maximum time separation (default: Nt // 2).

    Returns
    -------
    CorrelationResult
        Correlation function data.
    """
    config = gauge.config
    if t_max is None:
        t_max = config.Nt // 2

    # Compute O(t) for all t
    O_t = np.array([spatial_plaquette_sum(gauge, t) for t in range(config.Nt)])

    # Mean value
    O_mean = np.mean(O_t)

    # Correlation function (using periodic boundary conditions)
    G = np.zeros(t_max)
    for dt in range(t_max):
        corr = 0.0
        for t0 in range(config.Nt):
            t1 = (t0 + dt) % config.Nt
            corr += O_t[t0] * O_t[t1]
        corr /= config.Nt
        G[dt] = corr - O_mean**2

    return CorrelationResult(
        t_values=np.arange(t_max),
        G_values=G,
        G_errors=np.zeros(t_max),  # Single config, no error estimate
        operator_type="spatial_plaquette"
    )


def compute_polyakov_correlator(gauge: lg.GaugeField,
                                 r_max: Optional[int] = None) -> CorrelationResult:
    """
    Compute the Polyakov loop correlator as a function of distance.

    C(r) = <P(0) P(r)^*> - |<P>|^2

    Parameters
    ----------
    gauge : GaugeField
        The gauge configuration.
    r_max : int, optional
        Maximum distance (default: min(Nx, Ny, Nz) // 2).

    Returns
    -------
    CorrelationResult
        Correlator data as function of distance.
    """
    config = gauge.config
    if r_max is None:
        r_max = min(config.Nx, config.Ny, config.Nz) // 2

    # Compute Polyakov loops at all spatial sites
    P = np.zeros((config.Nx, config.Ny, config.Nz), dtype=complex)
    for x in range(config.Nx):
        for y in range(config.Ny):
            for z in range(config.Nz):
                P[x, y, z] = polyakov_loop_operator(gauge, (x, y, z))

    P_mean = np.mean(P)
    P_mean_sq = np.abs(P_mean)**2

    # Correlator at each distance
    C = np.zeros(r_max)
    counts = np.zeros(r_max)

    for x in range(config.Nx):
        for y in range(config.Ny):
            for z in range(config.Nz):
                for dx in range(config.Nx):
                    for dy in range(config.Ny):
                        for dz in range(config.Nz):
                            r = int(np.sqrt(dx**2 + dy**2 + dz**2))
                            if 0 < r < r_max:
                                x2 = (x + dx) % config.Nx
                                y2 = (y + dy) % config.Ny
                                z2 = (z + dz) % config.Nz
                                C[r] += np.real(P[x,y,z] * np.conj(P[x2,y2,z2]))
                                counts[r] += 1

    # Normalize
    for r in range(1, r_max):
        if counts[r] > 0:
            C[r] /= counts[r]
            C[r] -= P_mean_sq

    return CorrelationResult(
        t_values=np.arange(r_max),
        G_values=C,
        G_errors=np.zeros(r_max),
        operator_type="polyakov"
    )


# =============================================================================
# Ensemble Averaging
# =============================================================================

def average_correlators(correlators: List[CorrelationResult]) -> CorrelationResult:
    """
    Average correlation functions over an ensemble of configurations.

    Parameters
    ----------
    correlators : List[CorrelationResult]
        List of correlators from different configurations.

    Returns
    -------
    CorrelationResult
        Averaged correlator with statistical errors.
    """
    n_configs = len(correlators)
    t_values = correlators[0].t_values
    n_t = len(t_values)

    # Collect all measurements
    G_all = np.array([c.G_values for c in correlators])

    # Mean and standard error
    G_mean = np.mean(G_all, axis=0)
    G_std = np.std(G_all, axis=0)
    G_err = G_std / np.sqrt(n_configs)

    return CorrelationResult(
        t_values=t_values,
        G_values=G_mean,
        G_errors=G_err,
        operator_type=correlators[0].operator_type
    )


# =============================================================================
# Mass Gap Extraction from Correlators
# =============================================================================

@dataclass
class MassFromCorrelator:
    """Mass gap extracted from correlation function decay."""
    mass: float
    mass_err: float
    chi_squared: float
    fit_range: Tuple[int, int]
    method: str


def extract_mass_from_correlator(corr: CorrelationResult,
                                  t_min: int = 1,
                                  t_max: Optional[int] = None) -> MassFromCorrelator:
    """
    Extract mass gap from exponential decay of correlator.

    Fits G(t) = A * exp(-m * t) to extract m.

    Parameters
    ----------
    corr : CorrelationResult
        Correlation function data.
    t_min : int
        Minimum t for fit (to avoid contact terms).
    t_max : int, optional
        Maximum t for fit.

    Returns
    -------
    MassFromCorrelator
        Extracted mass and fit quality.
    """
    t = corr.t_values
    G = corr.G_values
    G_err = corr.G_errors

    if t_max is None:
        t_max = len(t) - 1

    # Select fit range
    mask = (t >= t_min) & (t <= t_max) & (G > 0)
    t_fit = t[mask]
    G_fit = G[mask]

    if len(t_fit) < 2:
        return MassFromCorrelator(
            mass=np.nan, mass_err=np.nan,
            chi_squared=np.nan, fit_range=(t_min, t_max),
            method="exponential_fit"
        )

    # Fit log(G) = log(A) - m*t
    log_G = np.log(G_fit)

    if G_err is not None and np.any(G_err[mask] > 0):
        # Weighted fit
        weights = G_fit / (G_err[mask] + 1e-10)
    else:
        weights = None

    try:
        slope, intercept, r_value, p_value, std_err = linregress(t_fit, log_G)
        mass = -slope
        mass_err = std_err

        # Chi-squared
        log_G_pred = intercept + slope * t_fit
        chi_sq = np.sum((log_G - log_G_pred)**2)

    except Exception:
        mass = np.nan
        mass_err = np.nan
        chi_sq = np.nan

    return MassFromCorrelator(
        mass=mass,
        mass_err=mass_err,
        chi_squared=chi_sq,
        fit_range=(t_min, t_max),
        method="exponential_fit"
    )


def extract_mass_effective(corr: CorrelationResult) -> np.ndarray:
    """
    Compute effective mass at each time slice.

    m_eff(t) = log(G(t) / G(t+1))

    The effective mass should plateau at large t, giving the mass gap.

    Parameters
    ----------
    corr : CorrelationResult
        Correlation function data.

    Returns
    -------
    np.ndarray
        Effective mass at each t.
    """
    G = corr.G_values
    m_eff = np.zeros(len(G) - 1)

    for t in range(len(G) - 1):
        if G[t] > 0 and G[t+1] > 0:
            m_eff[t] = np.log(G[t] / G[t+1])
        else:
            m_eff[t] = np.nan

    return m_eff


# =============================================================================
# Analysis Pipeline
# =============================================================================

class CorrelatorAnalyzer:
    """
    Complete analysis pipeline for extracting mass gap from correlators.
    """

    def __init__(self, gauge: lg.GaugeField):
        """Initialize with a gauge configuration."""
        self.gauge = gauge
        self.plaquette_corr: Optional[CorrelationResult] = None
        self.polyakov_corr: Optional[CorrelationResult] = None
        self.mass_result: Optional[MassFromCorrelator] = None

    def measure_all(self, verbose: bool = True) -> None:
        """Measure all correlation functions."""
        if verbose:
            print("Measuring correlation functions...")

        # Plaquette correlator
        self.plaquette_corr = compute_plaquette_correlator(self.gauge)
        if verbose:
            print(f"  Plaquette correlator: {len(self.plaquette_corr.t_values)} points")

        # Polyakov correlator
        self.polyakov_corr = compute_polyakov_correlator(self.gauge)
        if verbose:
            print(f"  Polyakov correlator: {len(self.polyakov_corr.t_values)} points")

    def extract_mass(self, t_min: int = 2, verbose: bool = True) -> MassFromCorrelator:
        """Extract mass gap from plaquette correlator."""
        if self.plaquette_corr is None:
            self.measure_all(verbose=False)

        self.mass_result = extract_mass_from_correlator(
            self.plaquette_corr, t_min=t_min
        )

        if verbose:
            print(f"\nMass gap from correlator:")
            print(f"  m = {self.mass_result.mass:.4f} +/- {self.mass_result.mass_err:.4f}")
            print(f"  chi^2 = {self.mass_result.chi_squared:.4f}")

        return self.mass_result


# =============================================================================
# Testing
# =============================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("Correlation Function Analysis - Test")
    print("=" * 60)

    # Create and thermalize a small lattice
    config = lg.LatticeConfig(Nx=4, Ny=4, Nz=4, Nt=8, beta=6.0)
    gauge = lg.GaugeField(config)
    gauge.cold_start()

    print("\n1. Thermalizing (30 sweeps)...")
    for i in range(30):
        lg.sweep(gauge, lg.UpdateMethod.METROPOLIS)
    plaq = lg.compute_average_plaquette(gauge)
    print(f"   Plaquette: {plaq:.4f}")

    print("\n2. Measuring correlators...")
    analyzer = CorrelatorAnalyzer(gauge)
    analyzer.measure_all()

    print("\n3. Plaquette correlator G(t):")
    for t, G in zip(analyzer.plaquette_corr.t_values[:5],
                    analyzer.plaquette_corr.G_values[:5]):
        print(f"   G({t}) = {G:.6f}")

    print("\n4. Effective mass:")
    m_eff = extract_mass_effective(analyzer.plaquette_corr)
    for t in range(min(4, len(m_eff))):
        if not np.isnan(m_eff[t]):
            print(f"   m_eff({t}) = {m_eff[t]:.4f}")

    print("\n5. Extracting mass gap...")
    mass_result = analyzer.extract_mass(t_min=1)

    print("\n" + "=" * 60)
    print("Test completed!")
    print("=" * 60)
