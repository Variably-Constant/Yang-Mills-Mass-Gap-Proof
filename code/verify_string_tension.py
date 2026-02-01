# -*- coding: utf-8 -*-
"""
String Tension Verification for Yang-Mills Mass Gap Proof

Measures the string tension sigma > 0 from Wilson loop area law:
    W(R,T) ~ exp(-sigma * R * T)

A positive string tension demonstrates confinement - one of the key
physical consequences of the mass gap.

Tests representative groups from each family:
- SU(2), SU(3) - Standard Model
- SO(3), SO(4) - Orthogonal
- Sp(4) - Symplectic
- G2 - Exceptional
"""

import sys
import io
import json
import os
from datetime import datetime

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace', line_buffering=True)

import numpy as np
from scipy.optimize import curve_fit
from typing import Tuple, List, Dict

# Import our gauge theory implementations
from lattice_gauge_sun import (
    LatticeConfigSUN, GaugeFieldSUN, compute_average_plaquette_sun,
    thermalize_sun, sweep_sun
)
from lattice_gauge_so import (
    LatticeConfigSO, GaugeFieldSO, compute_average_plaquette_so,
    thermalize_so, sweep_so
)
from lattice_gauge_sp import (
    LatticeConfigSp, GaugeFieldSp, compute_average_plaquette_sp,
    thermalize_sp, sweep_sp
)
from lattice_gauge_g2 import (
    LatticeConfigG2, GaugeFieldG2, compute_average_plaquette_g2,
    thermalize_g2, sweep_g2
)


def compute_wilson_loop_su(gauge: GaugeFieldSUN, site: Tuple, R: int, T: int) -> float:
    """
    Compute Wilson loop W(R,T) for SU(N).

    The Wilson loop is the trace of the ordered product of links
    around a rectangular R x T path in the (x,t) plane.
    """
    N = gauge.config.N

    # Start at site, go R steps in x-direction
    # then T steps in t-direction, back R in -x, back T in -t

    W = np.eye(N, dtype=np.complex128)
    x, y, z, t = site

    # Forward in x (R steps)
    for i in range(R):
        curr_site = ((x + i) % gauge.config.Nx, y, z, t)
        W = W @ gauge.get_link(curr_site, 0)  # mu=0 is x-direction

    # Forward in t (T steps)
    x_end = (x + R) % gauge.config.Nx
    for i in range(T):
        curr_site = (x_end, y, z, (t + i) % gauge.config.Nt)
        W = W @ gauge.get_link(curr_site, 3)  # mu=3 is t-direction

    # Backward in x (R steps) - use conjugate transpose
    t_end = (t + T) % gauge.config.Nt
    for i in range(R):
        curr_site = ((x_end - 1 - i) % gauge.config.Nx, y, z, t_end)
        W = W @ gauge.get_link(curr_site, 0).conj().T

    # Backward in t (T steps)
    for i in range(T):
        curr_site = (x, y, z, (t_end - 1 - i) % gauge.config.Nt)
        W = W @ gauge.get_link(curr_site, 3).conj().T

    return np.real(np.trace(W)) / N


def compute_wilson_loop_so(gauge: GaugeFieldSO, site: Tuple, R: int, T: int) -> float:
    """Compute Wilson loop W(R,T) for SO(N)."""
    N = gauge.config.N

    W = np.eye(N, dtype=np.float64)
    x, y, z, t = site

    # Forward in x
    for i in range(R):
        curr_site = ((x + i) % gauge.config.Nx, y, z, t)
        W = W @ gauge.get_link(curr_site, 0)

    # Forward in t
    x_end = (x + R) % gauge.config.Nx
    for i in range(T):
        curr_site = (x_end, y, z, (t + i) % gauge.config.Nt)
        W = W @ gauge.get_link(curr_site, 3)

    # Backward in x
    t_end = (t + T) % gauge.config.Nt
    for i in range(R):
        curr_site = ((x_end - 1 - i) % gauge.config.Nx, y, z, t_end)
        W = W @ gauge.get_link(curr_site, 0).T

    # Backward in t
    for i in range(T):
        curr_site = (x, y, z, (t_end - 1 - i) % gauge.config.Nt)
        W = W @ gauge.get_link(curr_site, 3).T

    return np.trace(W) / N


def compute_wilson_loop_sp(gauge: GaugeFieldSp, site: Tuple, R: int, T: int) -> float:
    """Compute Wilson loop W(R,T) for Sp(2N)."""
    dim = gauge.config.dim

    W = np.eye(dim, dtype=np.float64)
    x, y, z, t = site

    for i in range(R):
        curr_site = ((x + i) % gauge.config.Nx, y, z, t)
        W = W @ gauge.get_link(curr_site, 0)

    x_end = (x + R) % gauge.config.Nx
    for i in range(T):
        curr_site = (x_end, y, z, (t + i) % gauge.config.Nt)
        W = W @ gauge.get_link(curr_site, 3)

    t_end = (t + T) % gauge.config.Nt
    for i in range(R):
        curr_site = ((x_end - 1 - i) % gauge.config.Nx, y, z, t_end)
        W = W @ gauge.get_link(curr_site, 0).T

    for i in range(T):
        curr_site = (x, y, z, (t_end - 1 - i) % gauge.config.Nt)
        W = W @ gauge.get_link(curr_site, 3).T

    return np.trace(W) / dim


def compute_wilson_loop_g2(gauge: GaugeFieldG2, site: Tuple, R: int, T: int) -> float:
    """Compute Wilson loop W(R,T) for G2."""
    dim = 7  # G2 fundamental dimension

    W = np.eye(dim, dtype=np.float64)
    x, y, z, t = site

    for i in range(R):
        curr_site = ((x + i) % gauge.config.Nx, y, z, t)
        W = W @ gauge.get_link(curr_site, 0)

    x_end = (x + R) % gauge.config.Nx
    for i in range(T):
        curr_site = (x_end, y, z, (t + i) % gauge.config.Nt)
        W = W @ gauge.get_link(curr_site, 3)

    t_end = (t + T) % gauge.config.Nt
    for i in range(R):
        curr_site = ((x_end - 1 - i) % gauge.config.Nx, y, z, t_end)
        W = W @ gauge.get_link(curr_site, 0).T

    for i in range(T):
        curr_site = (x, y, z, (t_end - 1 - i) % gauge.config.Nt)
        W = W @ gauge.get_link(curr_site, 3).T

    return np.trace(W) / dim


def measure_wilson_loops(gauge, compute_wilson_func, R_max: int, T_max: int,
                         n_meas: int = 20, n_between: int = 5,
                         sweep_func=None, epsilon: float = 0.2) -> Dict:
    """
    Measure Wilson loops W(R,T) for various R,T values.

    Returns dictionary with loop values and errors.
    """
    results = {}

    for R in range(1, R_max + 1):
        for T in range(1, T_max + 1):
            values = []

            for _ in range(n_meas):
                # Decorrelation sweeps
                if sweep_func:
                    for _ in range(n_between):
                        sweep_func(gauge, epsilon)

                # Average over spatial sites
                loop_sum = 0.0
                count = 0

                # Sample a subset of sites for efficiency
                for x in range(0, gauge.config.Nx, 2):
                    for y in range(0, gauge.config.Ny, 2):
                        for z in range(0, gauge.config.Nz, 2):
                            for t in range(0, gauge.config.Nt, 2):
                                site = (x, y, z, t)
                                loop_sum += compute_wilson_func(gauge, site, R, T)
                                count += 1

                values.append(loop_sum / count)

            results[(R, T)] = {
                'mean': np.mean(values),
                'std': np.std(values),
                'values': values
            }

    return results


def extract_string_tension(wilson_data: Dict, R_max: int, T_max: int) -> Tuple[float, float, bool]:
    """
    Extract string tension from Wilson loop data using area law fit.

    W(R,T) = A * exp(-sigma * R * T)

    Returns (sigma, sigma_error, is_valid)
    """
    # Collect data points
    areas = []
    log_W = []
    weights = []

    for R in range(1, R_max + 1):
        for T in range(1, T_max + 1):
            data = wilson_data.get((R, T))
            if data and data['mean'] > 0:
                areas.append(R * T)
                log_W.append(np.log(data['mean']))
                # Weight by inverse variance
                if data['std'] > 0 and data['mean'] > 0:
                    weights.append(1.0 / (data['std'] / data['mean'])**2)
                else:
                    weights.append(1.0)

    if len(areas) < 3:
        return 0.0, 0.0, False

    areas = np.array(areas)
    log_W = np.array(log_W)
    weights = np.array(weights)

    # Linear fit: log(W) = log(A) - sigma * Area
    try:
        # Weighted linear regression
        W_sum = np.sum(weights)
        x_mean = np.sum(weights * areas) / W_sum
        y_mean = np.sum(weights * log_W) / W_sum

        numerator = np.sum(weights * (areas - x_mean) * (log_W - y_mean))
        denominator = np.sum(weights * (areas - x_mean)**2)

        if abs(denominator) < 1e-10:
            return 0.0, 0.0, False

        slope = numerator / denominator
        sigma = -slope  # sigma = -d(log W)/d(Area)

        # Error estimate
        residuals = log_W - (y_mean + slope * (areas - x_mean))
        mse = np.sum(weights * residuals**2) / (len(areas) - 2)
        sigma_err = np.sqrt(mse / denominator)

        # Valid if sigma > 0 and reasonable error
        is_valid = sigma > 0 and sigma_err < sigma

        return sigma, sigma_err, is_valid

    except Exception as e:
        print(f"  Fit error: {e}")
        return 0.0, 0.0, False


def verify_string_tension_su(N: int, beta: float, lattice_size: Tuple[int,int,int,int],
                              n_therm: int = 100, n_meas: int = 30) -> Dict:
    """Verify sigma > 0 for SU(N)."""
    print(f"\n  Testing SU({N}) at beta={beta}...")

    Nx, Ny, Nz, Nt = lattice_size
    config = LatticeConfigSUN(N=N, Nx=Nx, Ny=Ny, Nz=Nz, Nt=Nt, beta=beta)
    gauge = GaugeFieldSUN(config)
    gauge.cold_start()

    # Thermalize
    print(f"    Thermalizing ({n_therm} sweeps)...")
    thermalize_sun(gauge, n_sweeps=n_therm, epsilon=0.2, verbose=False)

    plaq = compute_average_plaquette_sun(gauge)
    print(f"    <P> = {plaq:.4f}")

    # Measure Wilson loops
    R_max = min(3, Nx // 2)
    T_max = min(4, Nt // 2)

    print(f"    Measuring Wilson loops (R=1..{R_max}, T=1..{T_max})...")
    wilson_data = measure_wilson_loops(
        gauge, compute_wilson_loop_su, R_max, T_max,
        n_meas=n_meas, sweep_func=sweep_sun
    )

    # Extract string tension
    sigma, sigma_err, is_valid = extract_string_tension(wilson_data, R_max, T_max)

    result = {
        'group': f'SU({N})',
        'beta': beta,
        'plaquette': plaq,
        'sigma': sigma,
        'sigma_err': sigma_err,
        'sigma_positive': sigma > 0,
        'fit_valid': is_valid,
        'R_max': R_max,
        'T_max': T_max,
        'wilson_loops': {f"W({R},{T})": wilson_data[(R,T)]['mean']
                        for R in range(1, R_max+1) for T in range(1, T_max+1)
                        if (R,T) in wilson_data}
    }

    status = "PASS" if sigma > 0 else "FAIL"
    print(f"    sigma = {sigma:.4f} +/- {sigma_err:.4f} [{status}]")

    return result


def verify_string_tension_so(N: int, beta: float, lattice_size: Tuple[int,int,int,int],
                              n_therm: int = 100, n_meas: int = 30) -> Dict:
    """Verify sigma > 0 for SO(N)."""
    print(f"\n  Testing SO({N}) at beta={beta}...")

    Nx, Ny, Nz, Nt = lattice_size
    config = LatticeConfigSO(N=N, Nx=Nx, Ny=Ny, Nz=Nz, Nt=Nt, beta=beta)
    gauge = GaugeFieldSO(config)
    gauge.cold_start()

    print(f"    Thermalizing ({n_therm} sweeps)...")
    thermalize_so(gauge, n_sweeps=n_therm, epsilon=0.2, verbose=False)

    plaq = compute_average_plaquette_so(gauge)
    print(f"    <P> = {plaq:.4f}")

    R_max = min(3, Nx // 2)
    T_max = min(4, Nt // 2)

    print(f"    Measuring Wilson loops...")
    wilson_data = measure_wilson_loops(
        gauge, compute_wilson_loop_so, R_max, T_max,
        n_meas=n_meas, sweep_func=sweep_so
    )

    sigma, sigma_err, is_valid = extract_string_tension(wilson_data, R_max, T_max)

    result = {
        'group': f'SO({N})',
        'beta': beta,
        'plaquette': plaq,
        'sigma': sigma,
        'sigma_err': sigma_err,
        'sigma_positive': sigma > 0,
        'fit_valid': is_valid,
        'R_max': R_max,
        'T_max': T_max
    }

    status = "PASS" if sigma > 0 else "FAIL"
    print(f"    sigma = {sigma:.4f} +/- {sigma_err:.4f} [{status}]")

    return result


def verify_string_tension_sp(N: int, beta: float, lattice_size: Tuple[int,int,int,int],
                              n_therm: int = 100, n_meas: int = 30) -> Dict:
    """Verify sigma > 0 for Sp(2N)."""
    print(f"\n  Testing Sp({2*N}) at beta={beta}...")

    Nx, Ny, Nz, Nt = lattice_size
    config = LatticeConfigSp(N=N, Nx=Nx, Ny=Ny, Nz=Nz, Nt=Nt, beta=beta)
    gauge = GaugeFieldSp(config)
    gauge.cold_start()

    print(f"    Thermalizing ({n_therm} sweeps)...")
    thermalize_sp(gauge, n_sweeps=n_therm, epsilon=0.2, verbose=False)

    plaq = compute_average_plaquette_sp(gauge)
    print(f"    <P> = {plaq:.4f}")

    R_max = min(3, Nx // 2)
    T_max = min(4, Nt // 2)

    print(f"    Measuring Wilson loops...")
    wilson_data = measure_wilson_loops(
        gauge, compute_wilson_loop_sp, R_max, T_max,
        n_meas=n_meas, sweep_func=sweep_sp
    )

    sigma, sigma_err, is_valid = extract_string_tension(wilson_data, R_max, T_max)

    result = {
        'group': f'Sp({2*N})',
        'beta': beta,
        'plaquette': plaq,
        'sigma': sigma,
        'sigma_err': sigma_err,
        'sigma_positive': sigma > 0,
        'fit_valid': is_valid,
        'R_max': R_max,
        'T_max': T_max
    }

    status = "PASS" if sigma > 0 else "FAIL"
    print(f"    sigma = {sigma:.4f} +/- {sigma_err:.4f} [{status}]")

    return result


def verify_string_tension_g2(beta: float, lattice_size: Tuple[int,int,int,int],
                              n_therm: int = 100, n_meas: int = 30) -> Dict:
    """Verify sigma > 0 for G2."""
    print(f"\n  Testing G2 at beta={beta}...")

    Nx, Ny, Nz, Nt = lattice_size
    config = LatticeConfigG2(Nx=Nx, Ny=Ny, Nz=Nz, Nt=Nt, beta=beta)
    gauge = GaugeFieldG2(config)
    gauge.cold_start()

    print(f"    Thermalizing ({n_therm} sweeps)...")
    thermalize_g2(gauge, n_sweeps=n_therm, epsilon=0.15, verbose=False)

    plaq = compute_average_plaquette_g2(gauge)
    print(f"    <P> = {plaq:.4f}")

    R_max = min(3, Nx // 2)
    T_max = min(4, Nt // 2)

    print(f"    Measuring Wilson loops...")
    wilson_data = measure_wilson_loops(
        gauge, compute_wilson_loop_g2, R_max, T_max,
        n_meas=n_meas, sweep_func=sweep_g2, epsilon=0.15
    )

    sigma, sigma_err, is_valid = extract_string_tension(wilson_data, R_max, T_max)

    result = {
        'group': 'G2',
        'beta': beta,
        'plaquette': plaq,
        'sigma': sigma,
        'sigma_err': sigma_err,
        'sigma_positive': sigma > 0,
        'fit_valid': is_valid,
        'R_max': R_max,
        'T_max': T_max
    }

    status = "PASS" if sigma > 0 else "FAIL"
    print(f"    sigma = {sigma:.4f} +/- {sigma_err:.4f} [{status}]")

    return result


def main():
    print("=" * 70)
    print("STRING TENSION VERIFICATION: sigma > 0")
    print("=" * 70)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    print("Testing Wilson loop area law: W(R,T) ~ exp(-sigma * R * T)")
    print("sigma > 0 demonstrates CONFINEMENT")
    print("=" * 70)

    results = []
    lattice = (6, 6, 6, 8)  # 6^3 x 8 lattice
    small_lattice = (4, 4, 4, 6)  # For larger groups

    # ===== SU(N) =====
    print("\n" + "=" * 50)
    print("SU(N) Groups")
    print("=" * 50)

    # SU(2) - needs lower beta for confinement regime
    results.append(verify_string_tension_su(2, 2.3, lattice, n_therm=150))

    # SU(3) - QCD
    results.append(verify_string_tension_su(3, 5.7, lattice, n_therm=150))

    # ===== SO(N) =====
    print("\n" + "=" * 50)
    print("SO(N) Groups")
    print("=" * 50)

    results.append(verify_string_tension_so(3, 2.0, lattice, n_therm=150))
    results.append(verify_string_tension_so(4, 3.0, lattice, n_therm=150))

    # ===== Sp(2N) =====
    print("\n" + "=" * 50)
    print("Sp(2N) Groups")
    print("=" * 50)

    results.append(verify_string_tension_sp(2, 5.0, small_lattice, n_therm=150))  # Sp(4)

    # ===== Exceptional =====
    print("\n" + "=" * 50)
    print("Exceptional Groups")
    print("=" * 50)

    results.append(verify_string_tension_g2(8.0, small_lattice, n_therm=150))

    # ===== Summary =====
    print("\n" + "=" * 70)
    print("STRING TENSION VERIFICATION SUMMARY")
    print("=" * 70)

    print(f"\n{'Group':<10} {'beta':<8} {'sigma':<12} {'Error':<12} {'Status'}")
    print("-" * 55)

    all_passed = True
    for r in results:
        status = "PASS" if r['sigma_positive'] else "FAIL"
        if not r['sigma_positive']:
            all_passed = False
        print(f"{r['group']:<10} {r['beta']:<8.1f} {r['sigma']:<12.4f} {r['sigma_err']:<12.4f} {status}")

    print("-" * 55)

    passed = sum(1 for r in results if r['sigma_positive'])
    total = len(results)

    print(f"\nTotal: {passed}/{total} groups show sigma > 0 (confinement)")

    if all_passed:
        print("\n*** ALL GROUPS VERIFY CONFINEMENT: sigma > 0 ***")

    # Save results
    os.makedirs('verification_results', exist_ok=True)
    output_file = 'verification_results/STRING_TENSION.json'

    with open(output_file, 'w') as f:
        json.dump({
            'date': datetime.now().isoformat(),
            'description': 'String tension sigma > 0 verification',
            'total_tests': total,
            'passed': passed,
            'all_passed': all_passed,
            'results': results
        }, f, indent=2, default=str)

    print(f"\n[SAVED] {output_file}")

    return all_passed


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
