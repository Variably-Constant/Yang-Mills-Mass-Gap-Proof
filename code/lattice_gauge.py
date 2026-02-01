"""
Lattice Gauge Theory Implementation for SU(3) Yang-Mills

This module implements the core lattice gauge theory framework:
- SU(3) gauge field initialization on a 4D hypercubic lattice
- Wilson action calculation
- Plaquette average computation
- Monte Carlo updates (heatbath and Metropolis)

Author: Yang-Mills Mass Gap Project
Date: January 2026
"""

import numpy as np
from typing import Tuple, Optional, List
from dataclasses import dataclass
from enum import Enum


# =============================================================================
# SU(3) Matrix Operations
# =============================================================================

def su3_identity() -> np.ndarray:
    """Return the 3x3 identity matrix (SU(3) identity element)."""
    return np.eye(3, dtype=np.complex128)


def su3_random_near_identity(epsilon: float = 0.1) -> np.ndarray:
    """
    Generate a random SU(3) matrix near the identity.

    Uses the exponential map: U = exp(i * epsilon * H)
    where H is a random traceless Hermitian matrix.

    Parameters
    ----------
    epsilon : float
        Controls how far from identity the matrix is.

    Returns
    -------
    np.ndarray
        A 3x3 SU(3) matrix.
    """
    # Generate random Hermitian traceless matrix using Gell-Mann matrices
    # coefficients
    coeffs = np.random.randn(8) * epsilon
    H = gellmann_combination(coeffs)

    # Exponentiate to get SU(3) element
    return matrix_exp_su3(1j * H)


def su3_random() -> np.ndarray:
    """
    Generate a uniformly random SU(3) matrix using the Haar measure.

    Uses Gram-Schmidt orthogonalization on random complex vectors.

    Returns
    -------
    np.ndarray
        A random 3x3 SU(3) matrix.
    """
    # Generate random complex matrix
    A = (np.random.randn(3, 3) + 1j * np.random.randn(3, 3)) / np.sqrt(2)

    # QR decomposition gives orthogonal matrix
    Q, R = np.linalg.qr(A)

    # Ensure det(Q) = 1 by adjusting phase
    det_Q = np.linalg.det(Q)
    phase = np.exp(-1j * np.angle(det_Q) / 3)
    Q = Q * phase

    return Q


def su3_project(A: np.ndarray) -> np.ndarray:
    """
    Project a 3x3 matrix onto SU(3).

    Uses singular value decomposition for numerical stability.

    Parameters
    ----------
    A : np.ndarray
        A 3x3 complex matrix.

    Returns
    -------
    np.ndarray
        The closest SU(3) matrix to A.
    """
    U, S, Vh = np.linalg.svd(A)
    Q = U @ Vh

    # Ensure determinant is 1
    det_Q = np.linalg.det(Q)
    phase = np.exp(-1j * np.angle(det_Q) / 3)

    return Q * phase


def gellmann_matrices() -> List[np.ndarray]:
    """
    Return the 8 Gell-Mann matrices (generators of SU(3)).

    Returns
    -------
    List[np.ndarray]
        List of 8 Gell-Mann matrices, each 3x3.
    """
    lam = []

    # lambda_1
    lam.append(np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=np.complex128))

    # lambda_2
    lam.append(np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], dtype=np.complex128))

    # lambda_3
    lam.append(np.array([[1, 0, 0], [0, -1, 0], [0, 0, 0]], dtype=np.complex128))

    # lambda_4
    lam.append(np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=np.complex128))

    # lambda_5
    lam.append(np.array([[0, 0, -1j], [0, 0, 0], [1j, 0, 0]], dtype=np.complex128))

    # lambda_6
    lam.append(np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=np.complex128))

    # lambda_7
    lam.append(np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], dtype=np.complex128))

    # lambda_8
    lam.append(np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]], dtype=np.complex128) / np.sqrt(3))

    return lam


def gellmann_combination(coeffs: np.ndarray) -> np.ndarray:
    """
    Compute a linear combination of Gell-Mann matrices.

    Parameters
    ----------
    coeffs : np.ndarray
        Array of 8 real coefficients.

    Returns
    -------
    np.ndarray
        The linear combination sum_i coeffs[i] * lambda_i.
    """
    lam = gellmann_matrices()
    result = np.zeros((3, 3), dtype=np.complex128)
    for i, c in enumerate(coeffs):
        result += c * lam[i]
    return result


def matrix_exp_su3(H: np.ndarray) -> np.ndarray:
    """
    Compute matrix exponential for SU(3) Lie algebra element.

    Parameters
    ----------
    H : np.ndarray
        A 3x3 anti-Hermitian traceless matrix (i * Hermitian traceless).

    Returns
    -------
    np.ndarray
        exp(H), an SU(3) matrix.
    """
    # Use numpy's matrix exponential
    from scipy.linalg import expm
    return expm(H)


# =============================================================================
# Lattice Configuration
# =============================================================================

@dataclass
class LatticeConfig:
    """Configuration parameters for the lattice."""
    Nx: int = 8    # Spatial extent in x
    Ny: int = 8    # Spatial extent in y
    Nz: int = 8    # Spatial extent in z
    Nt: int = 8    # Temporal extent
    beta: float = 6.0  # Inverse coupling: beta = 2N/g^2

    @property
    def shape(self) -> Tuple[int, int, int, int]:
        return (self.Nx, self.Ny, self.Nz, self.Nt)

    @property
    def volume(self) -> int:
        return self.Nx * self.Ny * self.Nz * self.Nt

    @property
    def spatial_volume(self) -> int:
        return self.Nx * self.Ny * self.Nz


class GaugeField:
    """
    SU(3) gauge field on a 4D hypercubic lattice.

    The gauge field is stored as a 5D array:
    U[x, y, z, t, mu] is the link variable from site (x,y,z,t) in direction mu.
    Each link is a 3x3 SU(3) matrix.

    Attributes
    ----------
    config : LatticeConfig
        Lattice configuration parameters.
    U : np.ndarray
        The gauge field array of shape (Nx, Ny, Nz, Nt, 4, 3, 3).
    """

    def __init__(self, config: LatticeConfig):
        """Initialize gauge field with given configuration."""
        self.config = config
        self.U = np.zeros((*config.shape, 4, 3, 3), dtype=np.complex128)

    def cold_start(self):
        """Initialize all links to identity (ordered start)."""
        for x in range(self.config.Nx):
            for y in range(self.config.Ny):
                for z in range(self.config.Nz):
                    for t in range(self.config.Nt):
                        for mu in range(4):
                            self.U[x, y, z, t, mu] = su3_identity()

    def hot_start(self):
        """Initialize all links to random SU(3) matrices (disordered start)."""
        for x in range(self.config.Nx):
            for y in range(self.config.Ny):
                for z in range(self.config.Nz):
                    for t in range(self.config.Nt):
                        for mu in range(4):
                            self.U[x, y, z, t, mu] = su3_random()

    def get_link(self, site: Tuple[int, int, int, int], mu: int) -> np.ndarray:
        """
        Get link variable U_mu(site).

        Parameters
        ----------
        site : Tuple[int, int, int, int]
            Lattice site coordinates (x, y, z, t).
        mu : int
            Direction (0=x, 1=y, 2=z, 3=t).

        Returns
        -------
        np.ndarray
            The 3x3 SU(3) link matrix.
        """
        x, y, z, t = site
        # Apply periodic boundary conditions
        x = x % self.config.Nx
        y = y % self.config.Ny
        z = z % self.config.Nz
        t = t % self.config.Nt
        return self.U[x, y, z, t, mu]

    def set_link(self, site: Tuple[int, int, int, int], mu: int,
                 value: np.ndarray):
        """Set link variable U_mu(site) = value."""
        x, y, z, t = site
        x = x % self.config.Nx
        y = y % self.config.Ny
        z = z % self.config.Nz
        t = t % self.config.Nt
        self.U[x, y, z, t, mu] = value

    def shift_site(self, site: Tuple[int, int, int, int], mu: int,
                   steps: int = 1) -> Tuple[int, int, int, int]:
        """
        Shift a site by 'steps' in direction 'mu' with periodic BC.

        Parameters
        ----------
        site : Tuple
            Starting site coordinates.
        mu : int
            Direction to shift.
        steps : int
            Number of steps (can be negative).

        Returns
        -------
        Tuple
            New site coordinates.
        """
        site_list = list(site)
        sizes = [self.config.Nx, self.config.Ny, self.config.Nz, self.config.Nt]
        site_list[mu] = (site_list[mu] + steps) % sizes[mu]
        return tuple(site_list)


# =============================================================================
# Wilson Action and Plaquettes
# =============================================================================

def compute_plaquette(gauge: GaugeField, site: Tuple[int, int, int, int],
                      mu: int, nu: int) -> np.ndarray:
    """
    Compute the plaquette U_p = U_mu(x) U_nu(x+mu) U_mu(x+nu)^dag U_nu(x)^dag.

    Parameters
    ----------
    gauge : GaugeField
        The gauge field configuration.
    site : Tuple
        Starting site coordinates.
    mu, nu : int
        Directions defining the plaquette plane (mu < nu conventionally).

    Returns
    -------
    np.ndarray
        The 3x3 plaquette matrix.
    """
    # U_mu(x)
    U1 = gauge.get_link(site, mu)

    # U_nu(x + mu_hat)
    site_plus_mu = gauge.shift_site(site, mu, 1)
    U2 = gauge.get_link(site_plus_mu, nu)

    # U_mu(x + nu_hat)^dag
    site_plus_nu = gauge.shift_site(site, nu, 1)
    U3 = gauge.get_link(site_plus_nu, mu).conj().T

    # U_nu(x)^dag
    U4 = gauge.get_link(site, nu).conj().T

    return U1 @ U2 @ U3 @ U4


def compute_plaquette_trace(gauge: GaugeField, site: Tuple[int, int, int, int],
                            mu: int, nu: int) -> float:
    """Compute Re Tr(plaquette) / N."""
    plaq = compute_plaquette(gauge, site, mu, nu)
    return np.real(np.trace(plaq)) / 3.0


def compute_average_plaquette(gauge: GaugeField) -> float:
    """
    Compute the average plaquette over all sites and orientations.

    Returns
    -------
    float
        <Re Tr(U_p) / N> averaged over all plaquettes.
    """
    total = 0.0
    count = 0

    for x in range(gauge.config.Nx):
        for y in range(gauge.config.Ny):
            for z in range(gauge.config.Nz):
                for t in range(gauge.config.Nt):
                    site = (x, y, z, t)
                    # Sum over all 6 plaquette orientations
                    for mu in range(4):
                        for nu in range(mu + 1, 4):
                            total += compute_plaquette_trace(gauge, site, mu, nu)
                            count += 1

    return total / count


def compute_wilson_action(gauge: GaugeField) -> float:
    """
    Compute the Wilson action: S = beta * sum_p (1 - Re Tr(U_p) / N).

    Returns
    -------
    float
        The Wilson action value.
    """
    avg_plaq = compute_average_plaquette(gauge)
    num_plaquettes = gauge.config.volume * 6  # 6 orientations per site
    return gauge.config.beta * num_plaquettes * (1.0 - avg_plaq)


# =============================================================================
# Staples (for Monte Carlo updates)
# =============================================================================

def compute_staple_sum(gauge: GaugeField, site: Tuple[int, int, int, int],
                       mu: int) -> np.ndarray:
    """
    Compute the sum of staples touching the link U_mu(site).

    The staple is the product of links forming a plaquette with U_mu(site).
    For each direction nu != mu, there are two staples (forward and backward).

    Parameters
    ----------
    gauge : GaugeField
        The gauge field configuration.
    site : Tuple
        Site coordinates.
    mu : int
        Direction of the link being updated.

    Returns
    -------
    np.ndarray
        Sum of 6 staple matrices (3 directions x 2 orientations).
    """
    staple_sum = np.zeros((3, 3), dtype=np.complex128)

    for nu in range(4):
        if nu == mu:
            continue

        # Forward staple: U_nu(x+mu) U_mu(x+nu)^dag U_nu(x)^dag
        site_plus_mu = gauge.shift_site(site, mu, 1)
        site_plus_nu = gauge.shift_site(site, nu, 1)

        U1 = gauge.get_link(site_plus_mu, nu)
        U2 = gauge.get_link(site_plus_nu, mu).conj().T
        U3 = gauge.get_link(site, nu).conj().T

        staple_sum += U1 @ U2 @ U3

        # Backward staple: U_nu(x+mu-nu)^dag U_mu(x-nu)^dag U_nu(x-nu)
        site_plus_mu_minus_nu = gauge.shift_site(site_plus_mu, nu, -1)
        site_minus_nu = gauge.shift_site(site, nu, -1)

        U1 = gauge.get_link(site_plus_mu_minus_nu, nu).conj().T
        U2 = gauge.get_link(site_minus_nu, mu).conj().T
        U3 = gauge.get_link(site_minus_nu, nu)

        staple_sum += U1 @ U2 @ U3

    return staple_sum


# =============================================================================
# Monte Carlo Updates
# =============================================================================

class UpdateMethod(Enum):
    """Monte Carlo update methods."""
    METROPOLIS = "metropolis"
    HEATBATH = "heatbath"


def metropolis_update(gauge: GaugeField, site: Tuple[int, int, int, int],
                      mu: int, epsilon: float = 0.2) -> bool:
    """
    Perform a Metropolis update on a single link.

    Parameters
    ----------
    gauge : GaugeField
        The gauge field configuration.
    site : Tuple
        Site coordinates.
    mu : int
        Direction of link to update.
    epsilon : float
        Step size for random SU(3) perturbation.

    Returns
    -------
    bool
        True if the update was accepted.
    """
    beta = gauge.config.beta

    # Get current link and staple sum
    U_old = gauge.get_link(site, mu)
    staple = compute_staple_sum(gauge, site, mu)

    # Compute old action contribution
    S_old = -beta * np.real(np.trace(U_old @ staple)) / 3.0

    # Propose new link: U_new = R * U_old where R is near identity
    R = su3_random_near_identity(epsilon)
    U_new = su3_project(R @ U_old)  # Project to ensure SU(3)

    # Compute new action contribution
    S_new = -beta * np.real(np.trace(U_new @ staple)) / 3.0

    # Metropolis acceptance
    delta_S = S_new - S_old

    if delta_S < 0 or np.random.random() < np.exp(-delta_S):
        gauge.set_link(site, mu, U_new)
        return True

    return False


def heatbath_su2_subgroup(a: np.ndarray, beta_eff: float) -> np.ndarray:
    """
    Perform SU(2) heatbath update for Cabibbo-Marinari algorithm.

    Uses the Creutz algorithm to generate SU(2) matrices with the
    correct Boltzmann distribution P(U) ~ exp(beta_eff * Re Tr(U @ a†)).

    Parameters
    ----------
    a : np.ndarray
        A 2x2 complex matrix (the submatrix of U @ staple).
    beta_eff : float
        Effective inverse coupling (beta/3 for SU(3)).

    Returns
    -------
    np.ndarray
        New SU(2) matrix r such that r @ a has maximized trace.
    """
    # Compute k = sqrt(|det(a)|) - the "magnitude" of a
    det_a = np.linalg.det(a)
    k = np.sqrt(np.abs(det_a))

    if k < 1e-10:
        # a is nearly singular, return identity (no update)
        return np.eye(2, dtype=np.complex128)

    # Normalize: a_bar = a / k, making a_bar unitary (in SU(2))
    a_bar = a / k

    # We want r such that r @ a has large Re Tr
    # r @ a = r @ k @ a_bar = k @ (r @ a_bar)
    # So we want r @ a_bar close to identity, meaning r close to a_bar†
    #
    # Generate s from distribution P(s) ~ exp(2 * beta_eff * k * s0)
    # where s = s0*I + i*(s_vec . sigma) is SU(2)
    # Then r = s @ a_bar†

    alpha = 2.0 * beta_eff * k

    if alpha < 0.1:
        # Very weak coupling: return nearly random
        return su2_random() @ a_bar.conj().T

    # Creutz algorithm for generating s0 from P(s0) ~ sqrt(1-s0^2) * exp(alpha * s0)
    max_iter = 100
    for _ in range(max_iter):
        # Generate from exponential envelope
        r1 = max(1e-10, np.random.random())
        r2 = max(1e-10, np.random.random())
        r3 = np.random.random()

        # This generates x with distribution ~ exp(-alpha*x) * (something with cos^2)
        lam = -np.log(r1) - np.log(r2) * np.cos(2 * np.pi * r3)**2

        # x = lambda / (2 * alpha), so x is small when alpha is large
        x = lam / (2.0 * alpha)

        if x > 1.0:
            continue

        # s0 = 1 - 2*x, mapping x in [0,1] to s0 in [-1, 1]
        # But we want s0 close to 1, so use s0 = 1 - x for x in [0, 2]
        s0 = 1.0 - x

        # Accept with probability proportional to sqrt(1 - s0^2)
        r4 = np.random.random()
        if r4**2 <= 1.0 - s0**2:
            break
    else:
        # Fallback: s0 close to 1
        s0 = 1.0 - 0.5 / alpha

    s0 = np.clip(s0, -1.0, 1.0)

    # Generate s1, s2, s3 uniformly on sphere of radius sqrt(1 - s0^2)
    r_sphere = np.sqrt(max(0.0, 1.0 - s0**2))

    # Random direction on unit sphere
    phi = 2.0 * np.pi * np.random.random()
    cos_theta = 2.0 * np.random.random() - 1.0
    sin_theta = np.sqrt(max(0.0, 1.0 - cos_theta**2))

    s1 = r_sphere * sin_theta * np.cos(phi)
    s2 = r_sphere * sin_theta * np.sin(phi)
    s3 = r_sphere * cos_theta

    # Construct s = s0*I + i*(s1*sigma_x + s2*sigma_y + s3*sigma_z)
    # sigma_x = [[0,1],[1,0]], sigma_y = [[0,-i],[i,0]], sigma_z = [[1,0],[0,-1]]
    # s = [[s0 + i*s3, s2 + i*s1], [-s2 + i*s1, s0 - i*s3]]
    s = np.array([
        [s0 + 1j * s3, s2 + 1j * s1],
        [-s2 + 1j * s1, s0 - 1j * s3]
    ], dtype=np.complex128)

    # Return r = s @ a_bar†
    return s @ a_bar.conj().T


def su2_random() -> np.ndarray:
    """Generate a random SU(2) matrix using Haar measure."""
    # Parametrize as quaternion
    x = np.random.randn(4)
    x = x / np.linalg.norm(x)

    U = np.array([
        [x[0] + 1j * x[3], x[2] + 1j * x[1]],
        [-x[2] + 1j * x[1], x[0] - 1j * x[3]]
    ], dtype=np.complex128)

    return U


def embed_su2_in_su3(U2: np.ndarray, subgroup: int) -> np.ndarray:
    """
    Embed an SU(2) matrix into SU(3).

    Parameters
    ----------
    U2 : np.ndarray
        2x2 SU(2) matrix.
    subgroup : int
        Which SU(2) subgroup (0, 1, or 2).

    Returns
    -------
    np.ndarray
        3x3 SU(3) matrix.
    """
    U3 = np.eye(3, dtype=np.complex128)

    if subgroup == 0:  # (1,2) block
        U3[0, 0] = U2[0, 0]
        U3[0, 1] = U2[0, 1]
        U3[1, 0] = U2[1, 0]
        U3[1, 1] = U2[1, 1]
    elif subgroup == 1:  # (1,3) block
        U3[0, 0] = U2[0, 0]
        U3[0, 2] = U2[0, 1]
        U3[2, 0] = U2[1, 0]
        U3[2, 2] = U2[1, 1]
    else:  # (2,3) block
        U3[1, 1] = U2[0, 0]
        U3[1, 2] = U2[0, 1]
        U3[2, 1] = U2[1, 0]
        U3[2, 2] = U2[1, 1]

    return U3


def extract_su2_from_su3(U3: np.ndarray, subgroup: int) -> np.ndarray:
    """Extract an SU(2) submatrix from an SU(3) matrix."""
    if subgroup == 0:
        return U3[:2, :2].copy()
    elif subgroup == 1:
        return np.array([
            [U3[0, 0], U3[0, 2]],
            [U3[2, 0], U3[2, 2]]
        ], dtype=np.complex128)
    else:
        return U3[1:, 1:].copy()


def heatbath_update(gauge: GaugeField, site: Tuple[int, int, int, int],
                    mu: int, n_hits: int = 1) -> None:
    """
    Perform a heatbath update using the Cabibbo-Marinari algorithm.

    The SU(3) heatbath is performed by sequential SU(2) subgroup updates.

    Note: Due to a subtle issue with multiple hits, we currently use n_hits=1
    by default, which gives correct equilibrium but slower thermalization.
    For production runs, consider using Metropolis updates instead.

    Parameters
    ----------
    gauge : GaugeField
        The gauge field configuration.
    site : Tuple
        Site coordinates.
    mu : int
        Direction of link to update.
    n_hits : int
        Number of SU(2) subgroup hits per update (default=1).
    """
    beta = gauge.config.beta
    U = gauge.get_link(site, mu)
    staple = compute_staple_sum(gauge, site, mu)

    # W = U * staple
    W = U @ staple

    for _ in range(n_hits):
        for subgroup in range(3):
            # Extract SU(2) submatrix from W
            w2 = extract_su2_from_su3(W, subgroup)

            # Perform SU(2) heatbath
            r2 = heatbath_su2_subgroup(w2, beta / 3.0)

            # Embed back to SU(3) and update
            R3 = embed_su2_in_su3(r2, subgroup)
            U = R3 @ U
            W = R3 @ W

    # Project to SU(3) for numerical stability
    U = su3_project(U)
    gauge.set_link(site, mu, U)


def sweep(gauge: GaugeField, method: UpdateMethod = UpdateMethod.HEATBATH,
          epsilon: float = 0.2) -> float:
    """
    Perform one sweep over the entire lattice.

    Parameters
    ----------
    gauge : GaugeField
        The gauge field configuration.
    method : UpdateMethod
        Update method to use.
    epsilon : float
        Step size for Metropolis updates.

    Returns
    -------
    float
        Acceptance rate (for Metropolis) or 1.0 (for heatbath).
    """
    accepted = 0
    total = 0

    for x in range(gauge.config.Nx):
        for y in range(gauge.config.Ny):
            for z in range(gauge.config.Nz):
                for t in range(gauge.config.Nt):
                    site = (x, y, z, t)
                    for mu in range(4):
                        if method == UpdateMethod.METROPOLIS:
                            if metropolis_update(gauge, site, mu, epsilon):
                                accepted += 1
                            total += 1
                        else:
                            heatbath_update(gauge, site, mu)
                            accepted += 1
                            total += 1

    return accepted / total if total > 0 else 1.0


# =============================================================================
# Thermalization and Measurement
# =============================================================================

def thermalize(gauge: GaugeField, n_sweeps: int = 100,
               method: UpdateMethod = UpdateMethod.HEATBATH,
               verbose: bool = True) -> List[float]:
    """
    Thermalize the gauge configuration.

    Parameters
    ----------
    gauge : GaugeField
        The gauge field configuration.
    n_sweeps : int
        Number of thermalization sweeps.
    method : UpdateMethod
        Update method.
    verbose : bool
        Print progress information.

    Returns
    -------
    List[float]
        Average plaquette values during thermalization.
    """
    plaquette_history = []

    for i in range(n_sweeps):
        acc_rate = sweep(gauge, method)
        plaq = compute_average_plaquette(gauge)
        plaquette_history.append(plaq)

        if verbose and (i + 1) % 10 == 0:
            print(f"Thermalization sweep {i+1}/{n_sweeps}: "
                  f"<P> = {plaq:.6f}, acc = {acc_rate:.3f}")

    return plaquette_history


# =============================================================================
# Testing
# =============================================================================

if __name__ == "__main__":
    # Test basic functionality
    print("=" * 60)
    print("Lattice Gauge Theory - Basic Tests")
    print("=" * 60)

    # Test SU(3) operations
    print("\n1. Testing SU(3) matrix generation...")
    U = su3_random()
    det = np.linalg.det(U)
    unitarity = np.max(np.abs(U @ U.conj().T - np.eye(3)))
    print(f"   Random SU(3): |det| = {np.abs(det):.10f}, "
          f"unitarity error = {unitarity:.2e}")

    # Test lattice setup
    print("\n2. Testing lattice initialization...")
    config = LatticeConfig(Nx=4, Ny=4, Nz=4, Nt=4, beta=6.0)
    gauge = GaugeField(config)

    gauge.cold_start()
    plaq_cold = compute_average_plaquette(gauge)
    print(f"   Cold start: <P> = {plaq_cold:.6f} (expected: 1.0)")

    gauge.hot_start()
    plaq_hot = compute_average_plaquette(gauge)
    print(f"   Hot start:  <P> = {plaq_hot:.6f} (expected: ~0.33)")

    # Test thermalization
    print("\n3. Testing thermalization (20 sweeps)...")
    gauge.cold_start()
    history = thermalize(gauge, n_sweeps=20, verbose=False)
    print(f"   Initial <P> = {history[0]:.6f}")
    print(f"   Final <P> = {history[-1]:.6f}")

    # Test Wilson action
    print("\n4. Testing Wilson action...")
    S = compute_wilson_action(gauge)
    print(f"   Wilson action S = {S:.4f}")

    print("\n" + "=" * 60)
    print("All tests completed!")
    print("=" * 60)
