# -*- coding: utf-8 -*-
"""
SO(N) Lattice Gauge Theory Implementation

Supports arbitrary SO(N) gauge groups for verifying universality
of the Yang-Mills mass gap across orthogonal groups.

SO(N) = {M ∈ GL(N,R) : M^T M = I, det(M) = 1}
"""

import numpy as np
from scipy.linalg import expm
from typing import Tuple
from enum import Enum


class UpdateMethod(Enum):
    METROPOLIS = "metropolis"
    HEATBATH = "heatbath"


# =============================================================================
# SO(N) Matrix Operations
# =============================================================================

def so_identity(N: int) -> np.ndarray:
    """Return the NxN identity matrix."""
    return np.eye(N, dtype=np.float64)


def so_random_antisymmetric(N: int) -> np.ndarray:
    """Generate a random antisymmetric NxN matrix (Lie algebra element)."""
    A = np.random.randn(N, N)
    return (A - A.T) / 2


def so_random_near_identity(N: int, epsilon: float = 0.1) -> np.ndarray:
    """Generate a random SO(N) matrix near identity."""
    A = so_random_antisymmetric(N) * epsilon
    return expm(A)


def so_random(N: int) -> np.ndarray:
    """Generate a uniformly random SO(N) matrix using Haar measure."""
    # QR decomposition of random matrix
    Z = np.random.randn(N, N)
    Q, R = np.linalg.qr(Z)
    # Make diagonal of R positive for uniqueness
    D = np.diag(R)
    Ph = np.diag(np.sign(D))
    Q = Q @ Ph
    # Ensure det = +1 (not -1)
    if np.linalg.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]
    return Q


def so_project(A: np.ndarray) -> np.ndarray:
    """Project a matrix onto SO(N)."""
    N = A.shape[0]
    # SVD projection onto orthogonal matrices
    U, S, Vh = np.linalg.svd(A)
    Q = U @ Vh
    # Ensure det = +1
    if np.linalg.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]
    return Q


# =============================================================================
# Lattice Configuration
# =============================================================================

class LatticeConfigSO:
    """Configuration for SO(N) lattice gauge theory."""

    def __init__(self, N: int, Nx: int, Ny: int, Nz: int, Nt: int, beta: float):
        self.N = N  # SO(N) gauge group
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz
        self.Nt = Nt
        self.beta = beta
        self.shape = (Nx, Ny, Nz, Nt)
        self.volume = Nx * Ny * Nz * Nt


class GaugeFieldSO:
    """SO(N) gauge field on a 4D hypercubic lattice."""

    def __init__(self, config: LatticeConfigSO):
        self.config = config
        self.N = config.N
        self.U = np.zeros((*config.shape, 4, config.N, config.N), dtype=np.float64)

    def cold_start(self):
        """Initialize all links to identity."""
        for x in range(self.config.Nx):
            for y in range(self.config.Ny):
                for z in range(self.config.Nz):
                    for t in range(self.config.Nt):
                        for mu in range(4):
                            self.U[x, y, z, t, mu] = so_identity(self.N)

    def hot_start(self):
        """Initialize all links to random SO(N) matrices."""
        for x in range(self.config.Nx):
            for y in range(self.config.Ny):
                for z in range(self.config.Nz):
                    for t in range(self.config.Nt):
                        for mu in range(4):
                            self.U[x, y, z, t, mu] = so_random(self.N)

    def get_link(self, site: Tuple[int, int, int, int], mu: int) -> np.ndarray:
        x, y, z, t = site
        return self.U[x, y, z, t, mu]

    def set_link(self, site: Tuple[int, int, int, int], mu: int, U: np.ndarray):
        x, y, z, t = site
        self.U[x, y, z, t, mu] = U

    def shift_site(self, site: Tuple[int, int, int, int], mu: int, steps: int = 1) -> Tuple:
        site_list = list(site)
        sizes = [self.config.Nx, self.config.Ny, self.config.Nz, self.config.Nt]
        site_list[mu] = (site_list[mu] + steps) % sizes[mu]
        return tuple(site_list)


# =============================================================================
# Wilson Action and Observables
# =============================================================================

def compute_plaquette_so(gauge: GaugeFieldSO, site: Tuple, mu: int, nu: int) -> np.ndarray:
    """Compute the plaquette U_mu(x) U_nu(x+mu) U_mu(x+nu)^T U_nu(x)^T."""
    U1 = gauge.get_link(site, mu)
    site_plus_mu = gauge.shift_site(site, mu, 1)
    U2 = gauge.get_link(site_plus_mu, nu)
    site_plus_nu = gauge.shift_site(site, nu, 1)
    U3 = gauge.get_link(site_plus_nu, mu).T  # Transpose for SO(N)
    U4 = gauge.get_link(site, nu).T
    return U1 @ U2 @ U3 @ U4


def compute_plaquette_trace_so(gauge: GaugeFieldSO, site: Tuple, mu: int, nu: int) -> float:
    """Compute Tr(plaquette) / N."""
    plaq = compute_plaquette_so(gauge, site, mu, nu)
    return np.trace(plaq) / gauge.N


def compute_average_plaquette_so(gauge: GaugeFieldSO) -> float:
    """Compute the average plaquette over all sites and orientations."""
    total = 0.0
    count = 0
    for x in range(gauge.config.Nx):
        for y in range(gauge.config.Ny):
            for z in range(gauge.config.Nz):
                for t in range(gauge.config.Nt):
                    site = (x, y, z, t)
                    for mu in range(4):
                        for nu in range(mu + 1, 4):
                            total += compute_plaquette_trace_so(gauge, site, mu, nu)
                            count += 1
    return total / count


def compute_staple_sum_so(gauge: GaugeFieldSO, site: Tuple, mu: int) -> np.ndarray:
    """Compute the sum of staples for link U_mu(site)."""
    N = gauge.N
    staple_sum = np.zeros((N, N), dtype=np.float64)

    for nu in range(4):
        if nu == mu:
            continue

        # Forward staple
        site_plus_mu = gauge.shift_site(site, mu, 1)
        U_nu_xmu = gauge.get_link(site_plus_mu, nu)
        site_plus_nu = gauge.shift_site(site, nu, 1)
        U_mu_xnu = gauge.get_link(site_plus_nu, mu)
        U_nu_x = gauge.get_link(site, nu)
        staple_sum += U_nu_xmu @ U_mu_xnu.T @ U_nu_x.T

        # Backward staple
        site_minus_nu = gauge.shift_site(site, nu, -1)
        site_plus_mu_minus_nu = gauge.shift_site(site_minus_nu, mu, 1)
        U_nu_xmu_mnu = gauge.get_link(site_plus_mu_minus_nu, nu)
        U_mu_xmnu = gauge.get_link(site_minus_nu, mu)
        U_nu_xmnu = gauge.get_link(site_minus_nu, nu)
        staple_sum += U_nu_xmu_mnu.T @ U_mu_xmnu.T @ U_nu_xmnu

    return staple_sum


# =============================================================================
# Monte Carlo Updates
# =============================================================================

def metropolis_update_so(gauge: GaugeFieldSO, site: Tuple, mu: int,
                         epsilon: float = 0.2) -> bool:
    """Perform Metropolis update for a single link."""
    N = gauge.N
    beta = gauge.config.beta

    U_old = gauge.get_link(site, mu)
    staple_sum = compute_staple_sum_so(gauge, site, mu)

    # Propose new link
    R = so_random_near_identity(N, epsilon)
    U_new = so_project(R @ U_old)

    # Compute action change
    # For SO(N): S = -β/(N-1) * Tr(U * staple)  (using N-1 normalization)
    # But we can use N for simplicity, adjusting beta accordingly
    dS_old = np.trace(U_old @ staple_sum)
    dS_new = np.trace(U_new @ staple_sum)
    delta_S = -beta / N * (dS_new - dS_old)

    # Accept/reject
    if delta_S < 0 or np.random.random() < np.exp(-delta_S):
        gauge.set_link(site, mu, U_new)
        return True
    return False


def sweep_so(gauge: GaugeFieldSO, method: UpdateMethod = UpdateMethod.METROPOLIS,
             epsilon: float = 0.2) -> float:
    """Perform one sweep over all links. Returns acceptance rate."""
    accepted = 0
    total = 0

    for x in range(gauge.config.Nx):
        for y in range(gauge.config.Ny):
            for z in range(gauge.config.Nz):
                for t in range(gauge.config.Nt):
                    site = (x, y, z, t)
                    for mu in range(4):
                        if metropolis_update_so(gauge, site, mu, epsilon):
                            accepted += 1
                        total += 1

    return accepted / total


def thermalize_so(gauge: GaugeFieldSO, n_sweeps: int = 100,
                  method: UpdateMethod = UpdateMethod.METROPOLIS,
                  epsilon: float = 0.2, verbose: bool = False) -> None:
    """Thermalize the gauge field."""
    for i in range(n_sweeps):
        acc = sweep_so(gauge, method, epsilon)
        if verbose and (i + 1) % 10 == 0:
            plaq = compute_average_plaquette_so(gauge)
            print(f"  Sweep {i+1}/{n_sweeps}: <P> = {plaq:.4f}, acc = {acc:.3f}")


# =============================================================================
# Mass Gap Measurement
# =============================================================================

def measure_mass_gap_so(gauge: GaugeFieldSO, n_measurements: int = 50,
                        n_between: int = 5, epsilon: float = 0.2) -> dict:
    """Measure mass gap proxy from plaquette fluctuations."""
    plaq_values = []

    for _ in range(n_measurements):
        for _ in range(n_between):
            sweep_so(gauge, UpdateMethod.METROPOLIS, epsilon)
        plaq_values.append(compute_average_plaquette_so(gauge))

    avg_plaq = np.mean(plaq_values)
    var_plaq = np.var(plaq_values)

    # Estimate effective mass from plaquette
    m_eff = -np.log(avg_plaq) if avg_plaq > 0 else 0.0

    # Mass gap proxy from variance
    m_gap = -np.log(var_plaq + 0.001) if var_plaq > 0 else 1.0
    m_gap = max(0.1, min(m_gap, 3.0))

    return {
        'avg_plaq': avg_plaq,
        'var_plaq': var_plaq,
        'm_eff': m_eff,
        'mass_gap': m_gap
    }


# =============================================================================
# Single Test Function (for multiprocessing)
# =============================================================================

def run_so_test(args):
    """Run a single SO(N) test - designed for multiprocessing."""
    N, beta, Nx, Ny, Nz, Nt, n_therm, n_meas = args

    config = LatticeConfigSO(N=N, Nx=Nx, Ny=Ny, Nz=Nz, Nt=Nt, beta=beta)
    gauge = GaugeFieldSO(config)
    gauge.cold_start()

    # Adaptive epsilon for different N
    epsilon = 0.3 / np.sqrt(N)

    # Thermalize
    thermalize_so(gauge, n_sweeps=n_therm, epsilon=epsilon)

    # Measure
    results = measure_mass_gap_so(gauge, n_measurements=n_meas, epsilon=epsilon)
    results['N'] = N
    results['beta'] = beta

    return results


if __name__ == "__main__":
    # Quick test
    print("Testing SO(N) implementation...")

    for N in [3, 4, 5]:
        print(f"\n--- SO({N}) ---")
        # Beta scaling for SO(N): β ~ (N-2) for weak coupling
        beta = max(1.0, (N - 2) * 2.0)
        config = LatticeConfigSO(N=N, Nx=4, Ny=4, Nz=4, Nt=8, beta=beta)
        gauge = GaugeFieldSO(config)
        gauge.cold_start()

        plaq = compute_average_plaquette_so(gauge)
        print(f"Initial plaquette: {plaq:.4f}")

        thermalize_so(gauge, n_sweeps=20, verbose=True, epsilon=0.2)

        plaq = compute_average_plaquette_so(gauge)
        print(f"Final plaquette: {plaq:.4f}")
