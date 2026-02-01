"""
Mass Gap Extraction from Lattice Gauge Theory

This module implements the extraction of the mass gap from Wilson loop data:
- Static quark potential V(R) extraction
- Cornell potential fitting: V(R) = sigma*R + mu + c/R
- String tension and mass gap determination
- Scaling analysis with lattice spacing

The mass gap Delta is the fundamental quantity needed for the
Millennium Prize Yang-Mills problem.

Author: Yang-Mills Mass Gap Project
Date: January 2026
"""

import numpy as np
from typing import Tuple, List, Optional, Dict
from dataclasses import dataclass
from scipy.optimize import curve_fit, minimize
from scipy.stats import linregress
import warnings

import lattice_gauge as lg
import wilson_loops as wl


# =============================================================================
# Potential Fitting Functions
# =============================================================================

def cornell_potential(R: np.ndarray, sigma: float, mu: float, c: float) -> np.ndarray:
    """
    Cornell potential: V(R) = sigma * R + mu + c / R

    Parameters
    ----------
    R : np.ndarray
        Distance values.
    sigma : float
        String tension (confining term).
    mu : float
        Constant self-energy term.
    c : float
        Coulomb coefficient (perturbative term).

    Returns
    -------
    np.ndarray
        Potential values.
    """
    return sigma * R + mu + c / R


def pure_linear_potential(R: np.ndarray, sigma: float, mu: float) -> np.ndarray:
    """
    Pure linear (confining) potential: V(R) = sigma * R + mu

    Parameters
    ----------
    R : np.ndarray
        Distance values.
    sigma : float
        String tension.
    mu : float
        Constant term.

    Returns
    -------
    np.ndarray
        Potential values.
    """
    return sigma * R + mu


def coulomb_plus_linear(R: np.ndarray, sigma: float, alpha: float, mu: float) -> np.ndarray:
    """
    Physical QCD potential: V(R) = sigma * R - (4/3) * alpha / R + mu

    Parameters
    ----------
    R : np.ndarray
        Distance values.
    sigma : float
        String tension.
    alpha : float
        Strong coupling constant.
    mu : float
        Constant term.

    Returns
    -------
    np.ndarray
        Potential values.
    """
    return sigma * R - (4.0/3.0) * alpha / R + mu


# =============================================================================
# Potential Extraction from Wilson Loops
# =============================================================================

def extract_potential_multiT(wilson_data: Dict[Tuple[int, int], float],
                              T_values: List[int],
                              method: str = "log_ratio") -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Extract V(R) from Wilson loops using multiple T values.

    The potential is extracted using:
    - Log method: V(R) = -ln(W(R,T)) / T  (for large T)
    - Log ratio: V(R) = ln(W(R,T) / W(R,T+1))  (more accurate)

    Parameters
    ----------
    wilson_data : Dict
        Wilson loop values W[(R, T)].
    T_values : List[int]
        List of T values to use.
    method : str
        Extraction method ("log" or "log_ratio").

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        (R values, V(R) values, V(R) errors).
    """
    # Find all R values
    R_set = set()
    for (R, T) in wilson_data.keys():
        R_set.add(R)

    R_values = sorted(R_set)
    V_values = []
    V_errors = []

    for R in R_values:
        V_estimates = []

        for T in T_values:
            if (R, T) not in wilson_data:
                continue

            W_RT = wilson_data[(R, T)]
            if W_RT <= 0:
                continue

            if method == "log_ratio" and (R, T + 1) in wilson_data:
                W_RT1 = wilson_data[(R, T + 1)]
                if W_RT1 > 0:
                    V = np.log(W_RT / W_RT1)
                    V_estimates.append(V)
            else:
                V = -np.log(W_RT) / T
                V_estimates.append(V)

        if len(V_estimates) > 0:
            V_mean = np.mean(V_estimates)
            V_err = np.std(V_estimates) / np.sqrt(len(V_estimates)) if len(V_estimates) > 1 else 0.1 * abs(V_mean)
            V_values.append(V_mean)
            V_errors.append(V_err)
        else:
            V_values.append(np.nan)
            V_errors.append(np.nan)

    return np.array(R_values), np.array(V_values), np.array(V_errors)


# =============================================================================
# Fitting Class
# =============================================================================

@dataclass
class PotentialFitResult:
    """Results from fitting the static quark potential."""
    sigma: float           # String tension
    sigma_err: float       # String tension error
    mu: float              # Constant term
    mu_err: float          # Constant term error
    c: float               # Coulomb coefficient
    c_err: float           # Coulomb coefficient error
    chi_squared: float     # Chi-squared of fit
    dof: int               # Degrees of freedom
    chi_squared_per_dof: float
    fit_type: str          # Type of fit used
    R_range: Tuple[float, float]  # R range used for fit


class StaticPotentialFitter:
    """
    Class for fitting the static quark potential from Wilson loop data.

    Attributes
    ----------
    R : np.ndarray
        Distance values.
    V : np.ndarray
        Potential values.
    V_err : np.ndarray
        Potential errors.
    """

    def __init__(self, R: np.ndarray, V: np.ndarray, V_err: Optional[np.ndarray] = None):
        """
        Initialize fitter with potential data.

        Parameters
        ----------
        R : np.ndarray
            Distance values.
        V : np.ndarray
            Potential values.
        V_err : np.ndarray, optional
            Potential errors.
        """
        # Remove NaN values
        valid = ~np.isnan(V)
        self.R = R[valid]
        self.V = V[valid]

        if V_err is not None:
            self.V_err = V_err[valid]
        else:
            self.V_err = np.ones_like(self.V) * 0.01

    def fit_cornell(self, R_min: float = 1.0, R_max: float = np.inf) -> PotentialFitResult:
        """
        Fit Cornell potential: V(R) = sigma * R + mu + c / R

        Parameters
        ----------
        R_min : float
            Minimum R for fit.
        R_max : float
            Maximum R for fit.

        Returns
        -------
        PotentialFitResult
            Fit results including string tension.
        """
        # Select data in range
        mask = (self.R >= R_min) & (self.R <= R_max)
        R_fit = self.R[mask]
        V_fit = self.V[mask]
        V_err_fit = self.V_err[mask]

        if len(R_fit) < 3:
            raise ValueError("Not enough points for Cornell fit")

        # Initial guesses
        p0 = [0.1, 0.0, -0.3]  # sigma, mu, c

        try:
            popt, pcov = curve_fit(
                cornell_potential, R_fit, V_fit,
                p0=p0, sigma=V_err_fit, absolute_sigma=True,
                maxfev=5000
            )
            perr = np.sqrt(np.diag(pcov))
        except RuntimeError:
            # Fall back to simpler fitting
            return self.fit_linear(R_min, R_max)

        sigma, mu, c = popt
        sigma_err, mu_err, c_err = perr

        # Compute chi-squared
        V_pred = cornell_potential(R_fit, *popt)
        chi_sq = np.sum(((V_fit - V_pred) / V_err_fit)**2)
        dof = len(R_fit) - 3

        return PotentialFitResult(
            sigma=sigma,
            sigma_err=sigma_err,
            mu=mu,
            mu_err=mu_err,
            c=c,
            c_err=c_err,
            chi_squared=chi_sq,
            dof=dof,
            chi_squared_per_dof=chi_sq / dof if dof > 0 else np.inf,
            fit_type="Cornell",
            R_range=(R_min, R_max)
        )

    def fit_linear(self, R_min: float = 1.0, R_max: float = np.inf) -> PotentialFitResult:
        """
        Fit linear potential: V(R) = sigma * R + mu

        Parameters
        ----------
        R_min : float
            Minimum R for fit.
        R_max : float
            Maximum R for fit.

        Returns
        -------
        PotentialFitResult
            Fit results.
        """
        mask = (self.R >= R_min) & (self.R <= R_max)
        R_fit = self.R[mask]
        V_fit = self.V[mask]
        V_err_fit = self.V_err[mask]

        if len(R_fit) < 2:
            raise ValueError("Not enough points for linear fit")

        # Weighted linear regression
        weights = 1.0 / V_err_fit**2
        W = np.sum(weights)
        Wx = np.sum(weights * R_fit)
        Wy = np.sum(weights * V_fit)
        Wxx = np.sum(weights * R_fit**2)
        Wxy = np.sum(weights * R_fit * V_fit)

        delta = W * Wxx - Wx**2
        sigma = (W * Wxy - Wx * Wy) / delta
        mu = (Wxx * Wy - Wx * Wxy) / delta

        sigma_err = np.sqrt(W / delta)
        mu_err = np.sqrt(Wxx / delta)

        # Chi-squared
        V_pred = pure_linear_potential(R_fit, sigma, mu)
        chi_sq = np.sum(((V_fit - V_pred) / V_err_fit)**2)
        dof = len(R_fit) - 2

        return PotentialFitResult(
            sigma=sigma,
            sigma_err=sigma_err,
            mu=mu,
            mu_err=mu_err,
            c=0.0,
            c_err=0.0,
            chi_squared=chi_sq,
            dof=dof,
            chi_squared_per_dof=chi_sq / dof if dof > 0 else np.inf,
            fit_type="Linear",
            R_range=(R_min, R_max)
        )


# =============================================================================
# Mass Gap Extraction
# =============================================================================

@dataclass
class MassGapResult:
    """Results from mass gap extraction."""
    mass_gap: float        # Mass gap in lattice units
    mass_gap_err: float    # Mass gap error
    string_tension: float  # String tension sigma
    string_tension_err: float
    method: str            # Extraction method
    lattice_spacing: Optional[float] = None  # a in fm (if known)
    mass_gap_MeV: Optional[float] = None     # Mass gap in MeV


def extract_mass_gap_from_string_tension(sigma: float, sigma_err: float,
                                         lattice_spacing: Optional[float] = None) -> MassGapResult:
    """
    Extract mass gap from string tension.

    The mass gap is related to the string tension by:
    Delta = sqrt(sigma) (in appropriate units)

    For QCD, the fundamental scale is sqrt(sigma) ~ 440 MeV.

    Parameters
    ----------
    sigma : float
        String tension in lattice units.
    sigma_err : float
        String tension error.
    lattice_spacing : float, optional
        Lattice spacing in fm.

    Returns
    -------
    MassGapResult
        Mass gap results.
    """
    if sigma <= 0:
        warnings.warn("Non-positive string tension; mass gap undefined")
        return MassGapResult(
            mass_gap=np.nan,
            mass_gap_err=np.nan,
            string_tension=sigma,
            string_tension_err=sigma_err,
            method="string_tension"
        )

    # Mass gap in lattice units
    Delta = np.sqrt(sigma)
    Delta_err = sigma_err / (2 * Delta)

    result = MassGapResult(
        mass_gap=Delta,
        mass_gap_err=Delta_err,
        string_tension=sigma,
        string_tension_err=sigma_err,
        method="string_tension"
    )

    # Convert to physical units if lattice spacing known
    if lattice_spacing is not None:
        # sqrt(sigma) ~ 440 MeV in QCD
        hbarc = 197.327  # MeV * fm
        result.lattice_spacing = lattice_spacing
        result.mass_gap_MeV = Delta * hbarc / lattice_spacing

    return result


def extract_mass_gap_from_correlator(gauge: lg.GaugeField,
                                      operator: str = "plaquette",
                                      T_max: int = 10) -> MassGapResult:
    """
    Extract mass gap from correlation function decay.

    C(t) ~ exp(-Delta * t) for large t

    Parameters
    ----------
    gauge : GaugeField
        Gauge configuration.
    operator : str
        Operator to correlate ("plaquette" or "polyakov").
    T_max : int
        Maximum time separation.

    Returns
    -------
    MassGapResult
        Mass gap results.
    """
    config = gauge.config
    correlators = np.zeros(T_max)

    if operator == "plaquette":
        # Compute plaquette-plaquette correlator
        for t in range(T_max):
            # Average over all spatial sites and orientations
            for x in range(config.Nx):
                for y in range(config.Ny):
                    for z in range(config.Nz):
                        site0 = (x, y, z, 0)
                        sitet = (x, y, z, t)

                        # Space-time plaquettes
                        for mu in range(3):
                            P0 = lg.compute_plaquette_trace(gauge, site0, mu, 3)
                            Pt = lg.compute_plaquette_trace(gauge, sitet, mu, 3)
                            correlators[t] += P0 * Pt

            correlators[t] /= (config.spatial_volume * 3)

    elif operator == "polyakov":
        # Polyakov loop correlator
        for t in range(T_max):
            correlators[t] = wl.compute_polyakov_correlator(gauge, t)

    # Fit exponential decay
    # log(C(t)) = const - Delta * t
    t_vals = np.arange(1, T_max)
    log_C = np.log(np.abs(correlators[1:]) + 1e-10)

    # Linear fit
    slope, intercept, r_value, p_value, std_err = linregress(t_vals, log_C)

    Delta = -slope
    Delta_err = std_err

    return MassGapResult(
        mass_gap=Delta,
        mass_gap_err=Delta_err,
        string_tension=Delta**2,  # Approximate
        string_tension_err=2 * Delta * Delta_err,
        method=f"correlator_{operator}"
    )


# =============================================================================
# Scaling Analysis
# =============================================================================

@dataclass
class ScalingResult:
    """Results from scaling analysis."""
    beta_values: np.ndarray
    sigma_values: np.ndarray
    sigma_errors: np.ndarray
    a_values: np.ndarray        # Lattice spacings
    sigma_physical: float       # Physical string tension
    sigma_physical_err: float
    scaling_violations: np.ndarray


def compute_lattice_spacing(beta: float, sigma_physical: float = 0.44**2) -> float:
    """
    Estimate lattice spacing from beta using two-loop perturbation theory.

    The lattice spacing a(beta) follows from the beta function:
    a * Lambda_lat = exp(-beta / (2 * b0)) * (b0 * g^2)^(-b1/(2*b0^2))

    For SU(3): b0 = 11/(16 pi^2), b1 = 102/(16 pi^2)^2

    Parameters
    ----------
    beta : float
        Inverse coupling (beta = 6/g^2).
    sigma_physical : float
        Physical string tension in GeV^2 (default: (440 MeV)^2).

    Returns
    -------
    float
        Lattice spacing in fm.
    """
    # Two-loop beta function coefficients for SU(3)
    b0 = 11.0 / (16.0 * np.pi**2)
    b1 = 102.0 / (16.0 * np.pi**2)**2

    g_squared = 6.0 / beta

    # Approximate lattice spacing
    # a * sqrt(sigma) = f(beta) from scaling relation
    # Empirical fit for SU(3)
    log_a = -beta / (12.0 * b0) + b1 / (2.0 * b0**2) * np.log(beta / (12.0 * b0))

    # Normalize to physical units (fm)
    a = np.exp(log_a)

    # Scale to physical string tension
    sqrt_sigma_physical = np.sqrt(sigma_physical)  # GeV
    hbarc = 0.197327  # GeV * fm

    a_fm = a * hbarc / sqrt_sigma_physical

    return a_fm


def perform_scaling_analysis(results: List[Tuple[float, PotentialFitResult]],
                              sigma_physical: float = 0.44**2) -> ScalingResult:
    """
    Analyze scaling behavior across multiple beta values.

    Parameters
    ----------
    results : List[Tuple[float, PotentialFitResult]]
        List of (beta, fit_result) pairs.
    sigma_physical : float
        Physical string tension in GeV^2.

    Returns
    -------
    ScalingResult
        Scaling analysis results.
    """
    beta_vals = []
    sigma_vals = []
    sigma_errs = []
    a_vals = []

    for beta, fit in results:
        beta_vals.append(beta)
        sigma_vals.append(fit.sigma)
        sigma_errs.append(fit.sigma_err)
        a_vals.append(compute_lattice_spacing(beta, sigma_physical))

    beta_vals = np.array(beta_vals)
    sigma_vals = np.array(sigma_vals)
    sigma_errs = np.array(sigma_errs)
    a_vals = np.array(a_vals)

    # sigma_lat * a^2 should be constant = sigma_phys
    sigma_phys_estimates = sigma_vals * a_vals**2

    sigma_phys_mean = np.average(sigma_phys_estimates, weights=1.0/sigma_errs**2)
    sigma_phys_err = np.sqrt(1.0 / np.sum(1.0/sigma_errs**2))

    # Scaling violations
    violations = (sigma_phys_estimates - sigma_phys_mean) / sigma_phys_mean

    return ScalingResult(
        beta_values=beta_vals,
        sigma_values=sigma_vals,
        sigma_errors=sigma_errs,
        a_values=a_vals,
        sigma_physical=sigma_phys_mean,
        sigma_physical_err=sigma_phys_err,
        scaling_violations=violations
    )


# =============================================================================
# Main Analysis Pipeline
# =============================================================================

class MassGapAnalyzer:
    """
    Complete analysis pipeline for mass gap extraction.

    This class orchestrates:
    1. Wilson loop measurements
    2. Potential extraction
    3. Cornell fit
    4. Mass gap determination
    """

    def __init__(self, gauge: lg.GaugeField, n_smear: int = 0):
        """
        Initialize analyzer.

        Parameters
        ----------
        gauge : GaugeField
            Gauge configuration.
        n_smear : int
            Number of APE smearing iterations.
        """
        self.gauge = gauge
        self.n_smear = n_smear
        self.wilson_analyzer: Optional[wl.WilsonLoopAnalyzer] = None
        self.potential_fitter: Optional[StaticPotentialFitter] = None
        self.fit_result: Optional[PotentialFitResult] = None
        self.mass_gap_result: Optional[MassGapResult] = None

    def run_analysis(self, R_max: int = 5, T_max: int = 8,
                     sample_fraction: float = 1.0,
                     fit_R_min: float = 1.0, fit_R_max: float = np.inf,
                     verbose: bool = True) -> MassGapResult:
        """
        Run complete mass gap analysis.

        Parameters
        ----------
        R_max : int
            Maximum R for Wilson loops.
        T_max : int
            Maximum T for Wilson loops.
        sample_fraction : float
            Fraction of sites to sample.
        fit_R_min, fit_R_max : float
            R range for potential fit.
        verbose : bool
            Print progress.

        Returns
        -------
        MassGapResult
            Final mass gap result.
        """
        if verbose:
            print("\n" + "=" * 60)
            print("MASS GAP ANALYSIS")
            print("=" * 60)

        # Step 1: Measure Wilson loops
        if verbose:
            print("\n1. Measuring Wilson loops...")

        if self.n_smear > 0:
            smeared = wl.ape_smear_configuration(self.gauge, self.n_smear)
            self.wilson_analyzer = wl.WilsonLoopAnalyzer(smeared)
        else:
            self.wilson_analyzer = wl.WilsonLoopAnalyzer(self.gauge)

        self.wilson_analyzer.measure_all(R_max, T_max, sample_fraction, verbose)

        # Step 2: Extract potential
        if verbose:
            print("\n2. Extracting static potential...")

        wilson_data = self.wilson_analyzer.get_wilson_loop_dict()
        T_values = list(range(2, T_max + 1))

        R, V, V_err = extract_potential_multiT(wilson_data, T_values, method="log_ratio")

        if verbose:
            print("   R    V(R)     Error")
            for i in range(len(R)):
                print(f"   {R[i]:.0f}    {V[i]:.4f}   {V_err[i]:.4f}")

        # Step 3: Fit potential
        if verbose:
            print("\n3. Fitting Cornell potential...")

        self.potential_fitter = StaticPotentialFitter(R, V, V_err)

        try:
            self.fit_result = self.potential_fitter.fit_cornell(fit_R_min, fit_R_max)
            if verbose:
                print(f"   sigma = {self.fit_result.sigma:.6f} +/- {self.fit_result.sigma_err:.6f}")
                print(f"   mu    = {self.fit_result.mu:.6f} +/- {self.fit_result.mu_err:.6f}")
                print(f"   c     = {self.fit_result.c:.6f} +/- {self.fit_result.c_err:.6f}")
                print(f"   chi^2/dof = {self.fit_result.chi_squared_per_dof:.3f}")
        except Exception as e:
            if verbose:
                print(f"   Cornell fit failed: {e}")
                print("   Falling back to linear fit...")
            self.fit_result = self.potential_fitter.fit_linear(fit_R_min, fit_R_max)
            if verbose:
                print(f"   sigma = {self.fit_result.sigma:.6f} +/- {self.fit_result.sigma_err:.6f}")
                print(f"   mu    = {self.fit_result.mu:.6f} +/- {self.fit_result.mu_err:.6f}")

        # Step 4: Extract mass gap
        if verbose:
            print("\n4. Extracting mass gap...")

        self.mass_gap_result = extract_mass_gap_from_string_tension(
            self.fit_result.sigma,
            self.fit_result.sigma_err,
            lattice_spacing=compute_lattice_spacing(self.gauge.config.beta)
        )

        if verbose:
            print(f"   Mass gap (lattice units): {self.mass_gap_result.mass_gap:.6f} +/- {self.mass_gap_result.mass_gap_err:.6f}")
            if self.mass_gap_result.mass_gap_MeV is not None:
                print(f"   Mass gap (physical):      {self.mass_gap_result.mass_gap_MeV:.1f} MeV")

        if verbose:
            print("\n" + "=" * 60)
            print("Analysis complete!")
            print("=" * 60)

        return self.mass_gap_result


# =============================================================================
# Testing
# =============================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("Mass Gap Extraction - Basic Tests")
    print("=" * 60)

    # Create and thermalize a lattice
    config = lg.LatticeConfig(Nx=6, Ny=6, Nz=6, Nt=12, beta=6.0)
    gauge = lg.GaugeField(config)
    gauge.hot_start()

    print("\n1. Thermalizing configuration (50 sweeps)...")
    lg.thermalize(gauge, n_sweeps=50, verbose=False)
    plaq = lg.compute_average_plaquette(gauge)
    print(f"   Average plaquette: {plaq:.6f}")

    # Run mass gap analysis
    analyzer = MassGapAnalyzer(gauge)
    result = analyzer.run_analysis(R_max=4, T_max=6, verbose=True)

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"String tension:    sigma = {result.string_tension:.6f}")
    print(f"Mass gap:          Delta = {result.mass_gap:.6f}")
    if result.mass_gap_MeV is not None:
        print(f"Mass gap (phys):   Delta ~ {result.mass_gap_MeV:.0f} MeV")
    print("=" * 60)
