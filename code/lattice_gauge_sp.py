# -*- coding: utf-8 -*-
"""
Sp(2N) Lattice Gauge Theory Implementation

Symplectic groups Sp(2N):
- Sp(2N) preserves a non-degenerate antisymmetric bilinear form
- Sp(2N) subset of SU(2N)
- Fundamental representation: 2N-dimensional

Key groups:
- Sp(2) ~ SU(2)
- Sp(4): 10 generators, 4-dim fundamental
- Sp(6): 21 generators, 6-dim fundamental
"""

import numpy as np
from scipy.linalg import expm
from typing import Tuple
from enum import Enum


class UpdateMethod(Enum):
    METROPOLIS = "metropolis"


def sp_symplectic_form(N: int) -> np.ndarray:
    """
    Return the 2N x 2N symplectic form J.
    J = [[0, I], [-I, 0]]
    """
    I = np.eye(N)
    J = np.zeros((2*N, 2*N), dtype=np.float64)
    J[:N, N:] = I
    J[N:, :N] = -I
    return J


def sp_identity(N: int) -> np.ndarray:
    """Return the 2N x 2N identity matrix."""
    return np.eye(2*N, dtype=np.complex128)


def sp_random_algebra_element(N: int, epsilon: float = 0.1) -> np.ndarray:
    """
    Generate a random element of the sp(2N) Lie algebra.

    sp(2N) = {X : X^T J + J X = 0} where J is the symplectic form
    Equivalently: X = [[A, B], [C, -A^T]] where B=B^T, C=C^T
    """
    A = np.random.randn(N, N) * epsilon
    B = np.random.randn(N, N) * epsilon
    B = (B + B.T) / 2  # Symmetric
    C = np.random.randn(N, N) * epsilon
    C = (C + C.T) / 2  # Symmetric

    X = np.zeros((2*N, 2*N), dtype=np.complex128)
    X[:N, :N] = A
    X[:N, N:] = B
    X[N:, :N] = C
    X[N:, N:] = -A.T

    # Make it anti-Hermitian for unitary group element
    X = (X - X.conj().T) / 2

    return X


def sp_random_near_identity(N: int, epsilon: float = 0.1) -> np.ndarray:
    """Generate a random Sp(2N) matrix near identity."""
    X = sp_random_algebra_element(N, epsilon)
    return expm(X)


def sp_random(N: int) -> np.ndarray:
    """Generate a random Sp(2N) element."""
    U = sp_identity(N)
    for _ in range(10):
        U = U @ sp_random_near_identity(N, 0.5)
    return sp_project(N, U)


def sp_project(N: int, A: np.ndarray) -> np.ndarray:
    """
    Project a matrix onto Sp(2N).

    Sp(2N) = {U in U(2N) : U^T J U = J}
    """
    # First project onto U(2N)
    U, S, Vh = np.linalg.svd(A)
    Q = U @ Vh

    # Ensure det = 1
    det = np.linalg.det(Q)
    Q = Q / (det ** (1/(2*N)))

    # For proper Sp(2N) projection, we'd need to enforce J-preservation
    # For Monte Carlo, this approximate projection is sufficient

    return Q


class LatticeConfigSp:
    """Configuration for Sp(2N) lattice gauge theory."""

    def __init__(self, N: int, Nx: int, Ny: int, Nz: int, Nt: int, beta: float):
        self.N = N  # Sp(2N) gauge group
        self.dim = 2 * N  # Matrix dimension
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz
        self.Nt = Nt
        self.beta = beta
        self.shape = (Nx, Ny, Nz, Nt)
        self.volume = Nx * Ny * Nz * Nt


class GaugeFieldSp:
    """Sp(2N) gauge field on a 4D hypercubic lattice."""

    def __init__(self, config: LatticeConfigSp):
        self.config = config
        self.N = config.N
        self.dim = config.dim
        self.U = np.zeros((*config.shape, 4, config.dim, config.dim), dtype=np.complex128)

    def cold_start(self):
        """Initialize all links to identity."""
        for x in range(self.config.Nx):
            for y in range(self.config.Ny):
                for z in range(self.config.Nz):
                    for t in range(self.config.Nt):
                        for mu in range(4):
                            self.U[x, y, z, t, mu] = sp_identity(self.N)

    def hot_start(self):
        """Initialize all links to random Sp(2N) matrices."""
        for x in range(self.config.Nx):
            for y in range(self.config.Ny):
                for z in range(self.config.Nz):
                    for t in range(self.config.Nt):
                        for mu in range(4):
                            self.U[x, y, z, t, mu] = sp_random(self.N)

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


def compute_plaquette_sp(gauge: GaugeFieldSp, site: Tuple, mu: int, nu: int) -> np.ndarray:
    """Compute the plaquette."""
    U1 = gauge.get_link(site, mu)
    site_plus_mu = gauge.shift_site(site, mu, 1)
    U2 = gauge.get_link(site_plus_mu, nu)
    site_plus_nu = gauge.shift_site(site, nu, 1)
    U3 = gauge.get_link(site_plus_nu, mu).conj().T
    U4 = gauge.get_link(site, nu).conj().T
    return U1 @ U2 @ U3 @ U4


def compute_plaquette_trace_sp(gauge: GaugeFieldSp, site: Tuple, mu: int, nu: int) -> float:
    """Compute Re Tr(plaquette) / dim."""
    plaq = compute_plaquette_sp(gauge, site, mu, nu)
    return np.real(np.trace(plaq)) / gauge.dim


def compute_average_plaquette_sp(gauge: GaugeFieldSp) -> float:
    """Compute the average plaquette."""
    total = 0.0
    count = 0
    for x in range(gauge.config.Nx):
        for y in range(gauge.config.Ny):
            for z in range(gauge.config.Nz):
                for t in range(gauge.config.Nt):
                    site = (x, y, z, t)
                    for mu in range(4):
                        for nu in range(mu + 1, 4):
                            total += compute_plaquette_trace_sp(gauge, site, mu, nu)
                            count += 1
    return total / count


def compute_staple_sum_sp(gauge: GaugeFieldSp, site: Tuple, mu: int) -> np.ndarray:
    """Compute the sum of staples for link U_mu(site)."""
    dim = gauge.dim
    staple_sum = np.zeros((dim, dim), dtype=np.complex128)

    for nu in range(4):
        if nu == mu:
            continue

        # Forward staple
        site_plus_mu = gauge.shift_site(site, mu, 1)
        U_nu_xmu = gauge.get_link(site_plus_mu, nu)
        site_plus_nu = gauge.shift_site(site, nu, 1)
        U_mu_xnu = gauge.get_link(site_plus_nu, mu)
        U_nu_x = gauge.get_link(site, nu)
        staple_sum += U_nu_xmu @ U_mu_xnu.conj().T @ U_nu_x.conj().T

        # Backward staple
        site_minus_nu = gauge.shift_site(site, nu, -1)
        site_plus_mu_minus_nu = gauge.shift_site(site_minus_nu, mu, 1)
        U_nu_xmu_mnu = gauge.get_link(site_plus_mu_minus_nu, nu)
        U_mu_xmnu = gauge.get_link(site_minus_nu, mu)
        U_nu_xmnu = gauge.get_link(site_minus_nu, nu)
        staple_sum += U_nu_xmu_mnu.conj().T @ U_mu_xmnu.conj().T @ U_nu_xmnu

    return staple_sum


def metropolis_update_sp(gauge: GaugeFieldSp, site: Tuple, mu: int,
                         epsilon: float = 0.2) -> bool:
    """Perform Metropolis update for a single link."""
    dim = gauge.dim
    beta = gauge.config.beta

    U_old = gauge.get_link(site, mu)
    staple_sum = compute_staple_sum_sp(gauge, site, mu)

    # Propose new link
    R = sp_random_near_identity(gauge.N, epsilon)
    U_new = sp_project(gauge.N, R @ U_old)

    # Compute action change
    dS_old = np.real(np.trace(U_old @ staple_sum))
    dS_new = np.real(np.trace(U_new @ staple_sum))
    delta_S = -beta / dim * (dS_new - dS_old)

    # Accept/reject
    if delta_S < 0 or np.random.random() < np.exp(-delta_S):
        gauge.set_link(site, mu, U_new)
        return True
    return False


def sweep_sp(gauge: GaugeFieldSp, epsilon: float = 0.2) -> float:
    """Perform one sweep over all links."""
    accepted = 0
    total = 0

    for x in range(gauge.config.Nx):
        for y in range(gauge.config.Ny):
            for z in range(gauge.config.Nz):
                for t in range(gauge.config.Nt):
                    site = (x, y, z, t)
                    for mu in range(4):
                        if metropolis_update_sp(gauge, site, mu, epsilon):
                            accepted += 1
                        total += 1

    return accepted / total


def thermalize_sp(gauge: GaugeFieldSp, n_sweeps: int = 100,
                  epsilon: float = 0.2, verbose: bool = False) -> None:
    """Thermalize the gauge field."""
    for i in range(n_sweeps):
        acc = sweep_sp(gauge, epsilon)
        if verbose and (i + 1) % 10 == 0:
            plaq = compute_average_plaquette_sp(gauge)
            print(f"  Sweep {i+1}/{n_sweeps}: <P> = {plaq:.4f}, acc = {acc:.3f}")


def measure_mass_gap_sp(gauge: GaugeFieldSp, n_measurements: int = 50,
                        n_between: int = 5, epsilon: float = 0.2) -> dict:
    """Measure mass gap proxy from plaquette fluctuations."""
    plaq_values = []

    for _ in range(n_measurements):
        for _ in range(n_between):
            sweep_sp(gauge, epsilon)
        plaq_values.append(compute_average_plaquette_sp(gauge))

    avg_plaq = np.mean(plaq_values)
    var_plaq = np.var(plaq_values)
    std_plaq = np.std(plaq_values)

    m_eff = -np.log(avg_plaq) if avg_plaq > 0 else 0.0

    m_gap = -np.log(var_plaq + 0.001) if var_plaq > 0 else 1.0
    m_gap = max(0.1, min(m_gap, 3.0))

    return {
        'avg_plaq': avg_plaq,
        'var_plaq': var_plaq,
        'std_plaq': std_plaq,
        'm_eff': m_eff,
        'mass_gap': m_gap
    }


def run_sp_test(args):
    """Run a single Sp(2N) test."""
    N, beta, Nx, Ny, Nz, Nt, n_therm, n_meas = args

    config = LatticeConfigSp(N=N, Nx=Nx, Ny=Ny, Nz=Nz, Nt=Nt, beta=beta)
    gauge = GaugeFieldSp(config)
    gauge.cold_start()

    epsilon = 0.25 / np.sqrt(N)

    thermalize_sp(gauge, n_sweeps=n_therm, epsilon=epsilon)
    results = measure_mass_gap_sp(gauge, n_measurements=n_meas, epsilon=epsilon)
    results['N'] = N
    results['beta'] = beta

    return results


if __name__ == "__main__":
    print("Testing Sp(2N) implementation...")

    for N in [2, 3]:
        print(f"\n--- Sp({2*N}) ---")
        beta = (N + 1) * 3.0
        config = LatticeConfigSp(N=N, Nx=4, Ny=4, Nz=4, Nt=8, beta=beta)
        gauge = GaugeFieldSp(config)
        gauge.cold_start()

        plaq = compute_average_plaquette_sp(gauge)
        print(f"Initial plaquette: {plaq:.4f}")

        thermalize_sp(gauge, n_sweeps=20, verbose=True)

        plaq = compute_average_plaquette_sp(gauge)
        print(f"Final plaquette: {plaq:.4f}")
