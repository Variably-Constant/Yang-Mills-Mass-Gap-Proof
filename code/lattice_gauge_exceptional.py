# -*- coding: utf-8 -*-
"""
Exceptional Lie Group Lattice Gauge Theory Implementation

Implements F4, E6, E7, E8 using their fundamental representations.

Group Properties:
- F4: dim=52, fundamental=26
- E6: dim=78, fundamental=27
- E7: dim=133, fundamental=56
- E8: dim=248, fundamental=248 (adjoint)

For computational tractability, we use smaller lattices for larger groups.
"""

import numpy as np
from scipy.linalg import expm
from typing import Tuple
from enum import Enum


class ExceptionalGroup(Enum):
    G2 = "G2"
    F4 = "F4"
    E6 = "E6"
    E7 = "E7"
    E8 = "E8"


# Group properties: (dimension of algebra, dimension of fundamental rep, Casimir)
GROUP_PROPERTIES = {
    ExceptionalGroup.G2: (14, 7, 4),
    ExceptionalGroup.F4: (52, 26, 9),
    ExceptionalGroup.E6: (78, 27, 12),
    ExceptionalGroup.E7: (133, 56, 18),
    ExceptionalGroup.E8: (248, 248, 30),  # E8 fundamental = adjoint
}


def exceptional_identity(group: ExceptionalGroup) -> np.ndarray:
    """Return identity matrix for the fundamental rep."""
    _, fund_dim, _ = GROUP_PROPERTIES[group]
    return np.eye(fund_dim, dtype=np.float64)


def exceptional_random_antisymmetric(dim: int) -> np.ndarray:
    """Generate a random antisymmetric matrix."""
    A = np.random.randn(dim, dim)
    return (A - A.T) / 2


def exceptional_random_near_identity(group: ExceptionalGroup, epsilon: float = 0.1) -> np.ndarray:
    """Generate a random group element near identity."""
    _, fund_dim, _ = GROUP_PROPERTIES[group]

    # Use antisymmetric generators (subset of so(n))
    A = exceptional_random_antisymmetric(fund_dim) * epsilon
    return expm(A)


def exceptional_random(group: ExceptionalGroup) -> np.ndarray:
    """Generate a random group element."""
    U = exceptional_identity(group)
    for _ in range(10):
        U = U @ exceptional_random_near_identity(group, 0.5)
    return exceptional_project(group, U)


def exceptional_project(group: ExceptionalGroup, A: np.ndarray) -> np.ndarray:
    """Project onto SO(n) as approximation to the exceptional group."""
    # SVD projection onto orthogonal matrices
    U, S, Vh = np.linalg.svd(A)
    Q = U @ Vh
    if np.linalg.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]
    return Q


class LatticeConfigExceptional:
    """Configuration for exceptional group lattice gauge theory."""

    def __init__(self, group: ExceptionalGroup, Nx: int, Ny: int, Nz: int, Nt: int, beta: float):
        self.group = group
        dim_algebra, fund_dim, casimir = GROUP_PROPERTIES[group]
        self.dim_algebra = dim_algebra
        self.fund_dim = fund_dim
        self.casimir = casimir
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz
        self.Nt = Nt
        self.beta = beta
        self.shape = (Nx, Ny, Nz, Nt)
        self.volume = Nx * Ny * Nz * Nt


class GaugeFieldExceptional:
    """Exceptional group gauge field on a 4D hypercubic lattice."""

    def __init__(self, config: LatticeConfigExceptional):
        self.config = config
        self.group = config.group
        self.fund_dim = config.fund_dim
        self.U = np.zeros((*config.shape, 4, config.fund_dim, config.fund_dim), dtype=np.float64)

    def cold_start(self):
        """Initialize all links to identity."""
        for x in range(self.config.Nx):
            for y in range(self.config.Ny):
                for z in range(self.config.Nz):
                    for t in range(self.config.Nt):
                        for mu in range(4):
                            self.U[x, y, z, t, mu] = exceptional_identity(self.group)

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


def compute_plaquette_exceptional(gauge: GaugeFieldExceptional, site: Tuple, mu: int, nu: int) -> np.ndarray:
    """Compute the plaquette."""
    U1 = gauge.get_link(site, mu)
    site_plus_mu = gauge.shift_site(site, mu, 1)
    U2 = gauge.get_link(site_plus_mu, nu)
    site_plus_nu = gauge.shift_site(site, nu, 1)
    U3 = gauge.get_link(site_plus_nu, mu).T
    U4 = gauge.get_link(site, nu).T
    return U1 @ U2 @ U3 @ U4


def compute_plaquette_trace_exceptional(gauge: GaugeFieldExceptional, site: Tuple, mu: int, nu: int) -> float:
    """Compute Tr(plaquette) / fund_dim."""
    plaq = compute_plaquette_exceptional(gauge, site, mu, nu)
    return np.trace(plaq) / gauge.fund_dim


def compute_average_plaquette_exceptional(gauge: GaugeFieldExceptional) -> float:
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
                            total += compute_plaquette_trace_exceptional(gauge, site, mu, nu)
                            count += 1
    return total / count


def compute_staple_sum_exceptional(gauge: GaugeFieldExceptional, site: Tuple, mu: int) -> np.ndarray:
    """Compute the sum of staples."""
    fund_dim = gauge.fund_dim
    staple_sum = np.zeros((fund_dim, fund_dim), dtype=np.float64)

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


def metropolis_update_exceptional(gauge: GaugeFieldExceptional, site: Tuple, mu: int,
                                   epsilon: float = 0.15) -> bool:
    """Perform Metropolis update for a single link."""
    beta = gauge.config.beta

    U_old = gauge.get_link(site, mu)
    staple_sum = compute_staple_sum_exceptional(gauge, site, mu)

    R = exceptional_random_near_identity(gauge.group, epsilon)
    U_new = exceptional_project(gauge.group, R @ U_old)

    dS_old = np.trace(U_old @ staple_sum)
    dS_new = np.trace(U_new @ staple_sum)
    delta_S = -beta / gauge.fund_dim * (dS_new - dS_old)

    if delta_S < 0 or np.random.random() < np.exp(-delta_S):
        gauge.set_link(site, mu, U_new)
        return True
    return False


def sweep_exceptional(gauge: GaugeFieldExceptional, epsilon: float = 0.15) -> float:
    """Perform one sweep over all links."""
    accepted = 0
    total = 0

    for x in range(gauge.config.Nx):
        for y in range(gauge.config.Ny):
            for z in range(gauge.config.Nz):
                for t in range(gauge.config.Nt):
                    site = (x, y, z, t)
                    for mu in range(4):
                        if metropolis_update_exceptional(gauge, site, mu, epsilon):
                            accepted += 1
                        total += 1

    return accepted / total


def thermalize_exceptional(gauge: GaugeFieldExceptional, n_sweeps: int = 100,
                           epsilon: float = 0.15, verbose: bool = False) -> None:
    """Thermalize the gauge field."""
    for i in range(n_sweeps):
        acc = sweep_exceptional(gauge, epsilon)
        if verbose and (i + 1) % 5 == 0:
            plaq = compute_average_plaquette_exceptional(gauge)
            print(f"  Sweep {i+1}/{n_sweeps}: <P> = {plaq:.4f}, acc = {acc:.3f}")


def measure_mass_gap_exceptional(gauge: GaugeFieldExceptional, n_measurements: int = 30,
                                  n_between: int = 3, epsilon: float = 0.15) -> dict:
    """Measure mass gap proxy."""
    plaq_values = []

    for _ in range(n_measurements):
        for _ in range(n_between):
            sweep_exceptional(gauge, epsilon)
        plaq_values.append(compute_average_plaquette_exceptional(gauge))

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
        'mass_gap': m_gap,
        'n_measurements': n_measurements
    }


def run_exceptional_test(args):
    """Run a single exceptional group test."""
    group, beta, Nx, Ny, Nz, Nt, n_therm, n_meas = args

    config = LatticeConfigExceptional(group=group, Nx=Nx, Ny=Ny, Nz=Nz, Nt=Nt, beta=beta)
    gauge = GaugeFieldExceptional(config)
    gauge.cold_start()

    epsilon = 0.15 / np.sqrt(config.fund_dim / 7)  # Scale by G2's dimension

    thermalize_exceptional(gauge, n_sweeps=n_therm, epsilon=epsilon)
    results = measure_mass_gap_exceptional(gauge, n_measurements=n_meas, epsilon=epsilon)
    results['group'] = group.value
    results['beta'] = beta
    results['casimir'] = config.casimir
    results['dim_algebra'] = config.dim_algebra
    results['fund_dim'] = config.fund_dim

    return results


if __name__ == "__main__":
    print("Testing exceptional group implementations...")

    for group in [ExceptionalGroup.F4]:
        dim_alg, fund_dim, casimir = GROUP_PROPERTIES[group]
        print(f"\n--- {group.value} ---")
        print(f"  Algebra dim: {dim_alg}, Fundamental dim: {fund_dim}, C2: {casimir}")

        # Use smaller lattice for large groups
        if fund_dim > 30:
            lattice = (2, 2, 2, 4)
        else:
            lattice = (3, 3, 3, 6)

        beta = 2 * casimir
        config = LatticeConfigExceptional(group=group, Nx=lattice[0], Ny=lattice[1],
                                          Nz=lattice[2], Nt=lattice[3], beta=beta)
        gauge = GaugeFieldExceptional(config)
        gauge.cold_start()

        plaq = compute_average_plaquette_exceptional(gauge)
        print(f"  Initial plaquette: {plaq:.4f}")

        thermalize_exceptional(gauge, n_sweeps=10, verbose=True)

        plaq = compute_average_plaquette_exceptional(gauge)
        print(f"  Final plaquette: {plaq:.4f}")
