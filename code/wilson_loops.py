"""
Wilson Loop Measurements for Lattice Gauge Theory

This module implements Wilson loop observables:
- Rectangular Wilson loops W(R, T)
- Polyakov loops
- Creutz ratios
- Correlation function measurements

Wilson loops are fundamental observables for extracting the static quark potential
and detecting confinement.

Author: Yang-Mills Mass Gap Project
Date: January 2026
"""

import numpy as np
from typing import Tuple, List, Optional, Dict
from dataclasses import dataclass
import lattice_gauge as lg


# =============================================================================
# Wilson Loop Computation
# =============================================================================

def compute_wilson_loop(gauge: lg.GaugeField,
                        start_site: Tuple[int, int, int, int],
                        R: int, T: int,
                        spatial_dir: int = 0,
                        temporal_dir: int = 3) -> np.ndarray:
    """
    Compute a single Wilson loop W(R, T) starting at given site.

    The Wilson loop is the trace of the product of link variables around
    a rectangular path of spatial extent R and temporal extent T.

            T
      +-----<-----+
      |           |
    R v           ^ R
      |           |
      +----->-----+
          T

    Parameters
    ----------
    gauge : GaugeField
        The gauge field configuration.
    start_site : Tuple
        Starting corner of the Wilson loop.
    R : int
        Spatial extent (in units of lattice spacing).
    T : int
        Temporal extent (in units of lattice spacing).
    spatial_dir : int
        Spatial direction (0, 1, or 2).
    temporal_dir : int
        Temporal direction (typically 3).

    Returns
    -------
    np.ndarray
        The Wilson loop matrix (3x3 for SU(3)).
    """
    site = start_site
    W = np.eye(3, dtype=np.complex128)

    # Bottom edge: go R steps in spatial direction
    for _ in range(R):
        W = W @ gauge.get_link(site, spatial_dir)
        site = gauge.shift_site(site, spatial_dir, 1)

    # Right edge: go T steps in temporal direction
    for _ in range(T):
        W = W @ gauge.get_link(site, temporal_dir)
        site = gauge.shift_site(site, temporal_dir, 1)

    # Top edge: go R steps in negative spatial direction
    for _ in range(R):
        site = gauge.shift_site(site, spatial_dir, -1)
        W = W @ gauge.get_link(site, spatial_dir).conj().T

    # Left edge: go T steps in negative temporal direction
    for _ in range(T):
        site = gauge.shift_site(site, temporal_dir, -1)
        W = W @ gauge.get_link(site, temporal_dir).conj().T

    return W


def compute_wilson_loop_trace(gauge: lg.GaugeField,
                               start_site: Tuple[int, int, int, int],
                               R: int, T: int,
                               spatial_dir: int = 0) -> complex:
    """
    Compute the trace of a Wilson loop.

    Returns
    -------
    complex
        Tr(W(R, T)).
    """
    W = compute_wilson_loop(gauge, start_site, R, T, spatial_dir)
    return np.trace(W)


def compute_average_wilson_loop(gauge: lg.GaugeField, R: int, T: int,
                                 sample_fraction: float = 1.0) -> Tuple[float, float]:
    """
    Compute the average Wilson loop over all starting sites and orientations.

    Parameters
    ----------
    gauge : GaugeField
        The gauge field configuration.
    R : int
        Spatial extent.
    T : int
        Temporal extent.
    sample_fraction : float
        Fraction of sites to sample (for efficiency).

    Returns
    -------
    Tuple[float, float]
        (mean, standard error) of Re Tr(W) / N.
    """
    config = gauge.config
    values = []

    for x in range(config.Nx):
        for y in range(config.Ny):
            for z in range(config.Nz):
                for t in range(config.Nt):
                    if sample_fraction < 1.0 and np.random.random() > sample_fraction:
                        continue

                    site = (x, y, z, t)

                    # Average over 3 spatial orientations
                    for spatial_dir in range(3):
                        W_trace = compute_wilson_loop_trace(
                            gauge, site, R, T, spatial_dir
                        )
                        values.append(np.real(W_trace) / 3.0)

    values = np.array(values)
    mean = np.mean(values)
    std_err = np.std(values) / np.sqrt(len(values))

    return mean, std_err


# =============================================================================
# APE Smearing (for noise reduction)
# =============================================================================

def ape_smear_link(gauge: lg.GaugeField, site: Tuple[int, int, int, int],
                   mu: int, alpha: float = 0.5) -> np.ndarray:
    """
    Compute APE-smeared link.

    APE smearing replaces spatial links with a weighted average including
    staples, which reduces UV fluctuations.

    Parameters
    ----------
    gauge : GaugeField
        The gauge field configuration.
    site : Tuple
        Site coordinates.
    mu : int
        Direction (only spatial directions 0, 1, 2 are smeared).
    alpha : float
        Smearing parameter (0 < alpha < 1).

    Returns
    -------
    np.ndarray
        Smeared link matrix.
    """
    if mu == 3:  # Don't smear temporal links
        return gauge.get_link(site, mu)

    U = gauge.get_link(site, mu)
    staple_sum = np.zeros((3, 3), dtype=np.complex128)

    # Sum over spatial staples only
    for nu in range(3):
        if nu == mu:
            continue

        # Forward staple
        site_plus_mu = gauge.shift_site(site, mu, 1)
        site_plus_nu = gauge.shift_site(site, nu, 1)

        U1 = gauge.get_link(site_plus_mu, nu)
        U2 = gauge.get_link(site_plus_nu, mu).conj().T
        U3 = gauge.get_link(site, nu).conj().T
        staple_sum += U1 @ U2 @ U3

        # Backward staple
        site_plus_mu_minus_nu = gauge.shift_site(site_plus_mu, nu, -1)
        site_minus_nu = gauge.shift_site(site, nu, -1)

        U1 = gauge.get_link(site_plus_mu_minus_nu, nu).conj().T
        U2 = gauge.get_link(site_minus_nu, mu).conj().T
        U3 = gauge.get_link(site_minus_nu, nu)
        staple_sum += U1 @ U2 @ U3

    # Smeared link: project (1 - alpha) * U + alpha/4 * staples onto SU(3)
    smeared = (1 - alpha) * U + (alpha / 4.0) * staple_sum
    return lg.su3_project(smeared)


def ape_smear_configuration(gauge: lg.GaugeField, n_smear: int = 10,
                            alpha: float = 0.5) -> lg.GaugeField:
    """
    Apply APE smearing to an entire configuration.

    Parameters
    ----------
    gauge : GaugeField
        Original gauge field.
    n_smear : int
        Number of smearing iterations.
    alpha : float
        Smearing parameter.

    Returns
    -------
    GaugeField
        Smeared gauge field (new object).
    """
    # Create copy
    smeared = lg.GaugeField(gauge.config)
    smeared.U = gauge.U.copy()

    for _ in range(n_smear):
        # Compute all smeared links first, then update
        new_U = np.zeros_like(smeared.U)

        for x in range(gauge.config.Nx):
            for y in range(gauge.config.Ny):
                for z in range(gauge.config.Nz):
                    for t in range(gauge.config.Nt):
                        site = (x, y, z, t)
                        for mu in range(4):
                            new_U[x, y, z, t, mu] = ape_smear_link(
                                smeared, site, mu, alpha
                            )

        smeared.U = new_U

    return smeared


# =============================================================================
# Smeared Wilson Loops
# =============================================================================

def compute_smeared_wilson_loop(gauge: lg.GaugeField,
                                 start_site: Tuple[int, int, int, int],
                                 R: int, T: int,
                                 spatial_dir: int = 0,
                                 n_smear: int = 10,
                                 alpha: float = 0.5) -> np.ndarray:
    """
    Compute Wilson loop with APE-smeared spatial links.

    Only spatial links are smeared; temporal links remain unsmeared
    to preserve the transfer matrix interpretation.

    Parameters
    ----------
    gauge : GaugeField
        The gauge field configuration.
    start_site : Tuple
        Starting corner.
    R, T : int
        Spatial and temporal extents.
    spatial_dir : int
        Spatial direction.
    n_smear : int
        Number of APE smearing iterations.
    alpha : float
        Smearing parameter.

    Returns
    -------
    np.ndarray
        Smeared Wilson loop matrix.
    """
    # Create smeared configuration
    smeared = ape_smear_configuration(gauge, n_smear, alpha)

    # Compute Wilson loop using smeared spatial links but original temporal links
    site = start_site
    W = np.eye(3, dtype=np.complex128)

    # Bottom edge: smeared spatial links
    for _ in range(R):
        W = W @ smeared.get_link(site, spatial_dir)
        site = gauge.shift_site(site, spatial_dir, 1)

    # Right edge: original temporal links
    for _ in range(T):
        W = W @ gauge.get_link(site, 3)
        site = gauge.shift_site(site, 3, 1)

    # Top edge: smeared spatial links (conjugate)
    for _ in range(R):
        site = gauge.shift_site(site, spatial_dir, -1)
        W = W @ smeared.get_link(site, spatial_dir).conj().T

    # Left edge: original temporal links (conjugate)
    for _ in range(T):
        site = gauge.shift_site(site, 3, -1)
        W = W @ gauge.get_link(site, 3).conj().T

    return W


# =============================================================================
# Polyakov Loop
# =============================================================================

def compute_polyakov_loop(gauge: lg.GaugeField,
                          spatial_site: Tuple[int, int, int]) -> np.ndarray:
    """
    Compute the Polyakov loop at a spatial site.

    The Polyakov loop is the trace of the product of temporal links
    around the periodic time direction:

    P(x) = Tr[ prod_{t=0}^{Nt-1} U_4(x, t) ]

    Parameters
    ----------
    gauge : GaugeField
        The gauge field configuration.
    spatial_site : Tuple
        Spatial coordinates (x, y, z).

    Returns
    -------
    np.ndarray
        Polyakov loop matrix.
    """
    x, y, z = spatial_site
    P = np.eye(3, dtype=np.complex128)

    for t in range(gauge.config.Nt):
        P = P @ gauge.get_link((x, y, z, t), 3)

    return P


def compute_average_polyakov_loop(gauge: lg.GaugeField) -> Tuple[complex, float]:
    """
    Compute spatially averaged Polyakov loop.

    Returns
    -------
    Tuple[complex, float]
        (mean, abs(mean)) of Tr(P) / N.
    """
    config = gauge.config
    P_sum = 0.0 + 0.0j

    for x in range(config.Nx):
        for y in range(config.Ny):
            for z in range(config.Nz):
                P = compute_polyakov_loop(gauge, (x, y, z))
                P_sum += np.trace(P) / 3.0

    mean = P_sum / config.spatial_volume
    return mean, np.abs(mean)


def compute_polyakov_correlator(gauge: lg.GaugeField, r: int) -> float:
    """
    Compute Polyakov loop correlator <P(0) P(r)^dag> at distance r.

    This is related to the free energy of a static quark-antiquark pair.

    Parameters
    ----------
    gauge : GaugeField
        The gauge field configuration.
    r : int
        Separation distance.

    Returns
    -------
    float
        Averaged correlator value.
    """
    config = gauge.config
    correlator_sum = 0.0
    count = 0

    for x in range(config.Nx):
        for y in range(config.Ny):
            for z in range(config.Nz):
                P1 = compute_polyakov_loop(gauge, (x, y, z))

                # Average over directions
                for direction in range(3):
                    if direction == 0:
                        x2 = (x + r) % config.Nx
                        spatial_site2 = (x2, y, z)
                    elif direction == 1:
                        y2 = (y + r) % config.Ny
                        spatial_site2 = (x, y2, z)
                    else:
                        z2 = (z + r) % config.Nz
                        spatial_site2 = (x, y, z2)

                    P2 = compute_polyakov_loop(gauge, spatial_site2)

                    # <P P^dag>
                    correlator_sum += np.real(np.trace(P1 @ P2.conj().T)) / 9.0
                    count += 1

    return correlator_sum / count


# =============================================================================
# Creutz Ratios
# =============================================================================

def compute_creutz_ratio(W: Dict[Tuple[int, int], float],
                         R: int, T: int) -> Optional[float]:
    """
    Compute the Creutz ratio chi(R, T).

    chi(R, T) = -ln( W(R,T) W(R-1,T-1) / (W(R,T-1) W(R-1,T)) )

    The Creutz ratio provides a lattice estimate of the string tension.

    Parameters
    ----------
    W : Dict
        Dictionary of Wilson loop values W[(R, T)].
    R, T : int
        Loop dimensions.

    Returns
    -------
    Optional[float]
        Creutz ratio, or None if not computable.
    """
    try:
        W_RT = W[(R, T)]
        W_Rm1_Tm1 = W[(R - 1, T - 1)]
        W_R_Tm1 = W[(R, T - 1)]
        W_Rm1_T = W[(R - 1, T)]

        if W_RT <= 0 or W_Rm1_Tm1 <= 0 or W_R_Tm1 <= 0 or W_Rm1_T <= 0:
            return None

        ratio = (W_RT * W_Rm1_Tm1) / (W_R_Tm1 * W_Rm1_T)
        if ratio <= 0:
            return None

        return -np.log(ratio)

    except (KeyError, ValueError):
        return None


# =============================================================================
# Wilson Loop Measurements Class
# =============================================================================

@dataclass
class WilsonLoopMeasurement:
    """Container for Wilson loop measurement results."""
    R: int
    T: int
    mean: float
    error: float
    n_samples: int


class WilsonLoopAnalyzer:
    """
    Class for systematic Wilson loop measurements.

    Attributes
    ----------
    gauge : GaugeField
        The gauge field configuration.
    measurements : Dict
        Stored measurement results.
    """

    def __init__(self, gauge: lg.GaugeField):
        """Initialize analyzer with gauge configuration."""
        self.gauge = gauge
        self.measurements: Dict[Tuple[int, int], WilsonLoopMeasurement] = {}

    def measure_all(self, R_max: int, T_max: int,
                    sample_fraction: float = 1.0,
                    verbose: bool = True) -> None:
        """
        Measure all Wilson loops up to given dimensions.

        Parameters
        ----------
        R_max : int
            Maximum spatial extent.
        T_max : int
            Maximum temporal extent.
        sample_fraction : float
            Fraction of sites to sample.
        verbose : bool
            Print progress.
        """
        if verbose:
            print(f"Measuring Wilson loops up to R={R_max}, T={T_max}...")

        for R in range(1, R_max + 1):
            for T in range(1, T_max + 1):
                mean, error = compute_average_wilson_loop(
                    self.gauge, R, T, sample_fraction
                )
                n_samples = int(self.gauge.config.volume * 3 * sample_fraction)

                self.measurements[(R, T)] = WilsonLoopMeasurement(
                    R=R, T=T, mean=mean, error=error, n_samples=n_samples
                )

                if verbose:
                    print(f"  W({R},{T}) = {mean:.6e} +/- {error:.6e}")

    def get_wilson_loop_dict(self) -> Dict[Tuple[int, int], float]:
        """Return dictionary of Wilson loop means."""
        return {key: m.mean for key, m in self.measurements.items()}

    def compute_creutz_ratios(self, R_max: int, T_max: int) -> Dict[Tuple[int, int], float]:
        """Compute all Creutz ratios up to given dimensions."""
        W = self.get_wilson_loop_dict()
        creutz = {}

        for R in range(2, R_max + 1):
            for T in range(2, T_max + 1):
                chi = compute_creutz_ratio(W, R, T)
                if chi is not None:
                    creutz[(R, T)] = chi

        return creutz

    def extract_potential_from_loops(self, T_fit: int) -> Tuple[np.ndarray, np.ndarray]:
        """
        Extract V(R) from Wilson loops at fixed T.

        V(R) = -ln(W(R, T)) / T (for large T)

        Parameters
        ----------
        T_fit : int
            Temporal extent to use for fitting.

        Returns
        -------
        Tuple[np.ndarray, np.ndarray]
            Arrays of (R values, V(R) values).
        """
        R_values = []
        V_values = []

        for (R, T), m in self.measurements.items():
            if T == T_fit and m.mean > 0:
                V = -np.log(m.mean) / T
                R_values.append(R)
                V_values.append(V)

        # Sort by R
        idx = np.argsort(R_values)
        return np.array(R_values)[idx], np.array(V_values)[idx]


# =============================================================================
# Testing
# =============================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("Wilson Loop Measurements - Basic Tests")
    print("=" * 60)

    # Create and thermalize a small lattice
    config = lg.LatticeConfig(Nx=4, Ny=4, Nz=4, Nt=8, beta=6.0)
    gauge = lg.GaugeField(config)
    gauge.cold_start()

    print("\n1. Thermalizing configuration (20 sweeps)...")
    lg.thermalize(gauge, n_sweeps=20, verbose=False)
    plaq = lg.compute_average_plaquette(gauge)
    print(f"   Average plaquette: {plaq:.6f}")

    print("\n2. Computing Wilson loops...")
    analyzer = WilsonLoopAnalyzer(gauge)
    analyzer.measure_all(R_max=3, T_max=4, verbose=True)

    print("\n3. Computing Creutz ratios...")
    creutz = analyzer.compute_creutz_ratios(R_max=3, T_max=4)
    for (R, T), chi in sorted(creutz.items()):
        print(f"   chi({R},{T}) = {chi:.6f}")

    print("\n4. Computing Polyakov loop...")
    P_mean, P_abs = compute_average_polyakov_loop(gauge)
    print(f"   <P> = {P_mean:.6f}")
    print(f"   |<P>| = {P_abs:.6f}")

    print("\n5. Extracting potential from W(R, T=3)...")
    R_vals, V_vals = analyzer.extract_potential_from_loops(T_fit=3)
    for R, V in zip(R_vals, V_vals):
        print(f"   V({R}) = {V:.6f}")

    print("\n" + "=" * 60)
    print("Wilson loop tests completed!")
    print("=" * 60)
