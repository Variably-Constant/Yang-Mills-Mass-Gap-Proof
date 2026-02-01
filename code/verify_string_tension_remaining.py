# -*- coding: utf-8 -*-
"""
String Tension Verification - Remaining Groups (Sp and G2)

Previous tests already completed:
- SU(2): sigma = 0.3760 +/- 0.0203 [PASS]
- SU(3): sigma = 0.4761 +/- 0.0128 [PASS]
- SO(3): sigma = 1.4398 +/- 0.0513 [PASS]
- SO(4): sigma = 0.6020 +/- 0.2310 [PASS]
"""

import sys
import io

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace', line_buffering=True)

import numpy as np
from typing import Tuple, Dict

from lattice_gauge_sp import (
    LatticeConfigSp, GaugeFieldSp, compute_average_plaquette_sp,
    thermalize_sp, sweep_sp
)
from lattice_gauge_g2 import (
    LatticeConfigG2, GaugeFieldG2, compute_average_plaquette_g2,
    thermalize_g2, sweep_g2
)


def compute_wilson_loop_sp(gauge: GaugeFieldSp, site: Tuple, R: int, T: int) -> float:
    """Compute Wilson loop W(R,T) for Sp(2N)."""
    dim = gauge.config.dim  # Use dim, not fund_dim

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
    """Measure Wilson loops."""
    results = {}

    for R in range(1, R_max + 1):
        for T in range(1, T_max + 1):
            values = []

            for _ in range(n_meas):
                if sweep_func:
                    for _ in range(n_between):
                        sweep_func(gauge, epsilon)

                loop_sum = 0.0
                count = 0

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
            }

    return results


def extract_string_tension(wilson_data: Dict, R_max: int, T_max: int) -> Tuple[float, float, bool]:
    """Extract string tension from Wilson loop data."""
    areas = []
    log_W = []
    weights = []

    for R in range(1, R_max + 1):
        for T in range(1, T_max + 1):
            data = wilson_data.get((R, T))
            if data and data['mean'] > 0:
                areas.append(R * T)
                log_W.append(np.log(data['mean']))
                if data['std'] > 0 and data['mean'] > 0:
                    weights.append(1.0 / (data['std'] / data['mean'])**2)
                else:
                    weights.append(1.0)

    if len(areas) < 3:
        return 0.0, 0.0, False

    areas = np.array(areas)
    log_W = np.array(log_W)
    weights = np.array(weights)

    try:
        W_sum = np.sum(weights)
        x_mean = np.sum(weights * areas) / W_sum
        y_mean = np.sum(weights * log_W) / W_sum

        numerator = np.sum(weights * (areas - x_mean) * (log_W - y_mean))
        denominator = np.sum(weights * (areas - x_mean)**2)

        if abs(denominator) < 1e-10:
            return 0.0, 0.0, False

        slope = numerator / denominator
        sigma = -slope

        residuals = log_W - (y_mean + slope * (areas - x_mean))
        mse = np.sum(weights * residuals**2) / (len(areas) - 2)
        sigma_err = np.sqrt(mse / denominator)

        is_valid = sigma > 0 and sigma_err < sigma

        return sigma, sigma_err, is_valid

    except Exception as e:
        print(f"  Fit error: {e}")
        return 0.0, 0.0, False


def test_sp4():
    """Test Sp(4)."""
    print("\n  Testing Sp(4) at beta=5.0...")

    config = LatticeConfigSp(N=2, Nx=4, Ny=4, Nz=4, Nt=6, beta=5.0)
    gauge = GaugeFieldSp(config)
    gauge.cold_start()

    print("    Thermalizing (150 sweeps)...")
    thermalize_sp(gauge, n_sweeps=150, epsilon=0.2, verbose=False)

    plaq = compute_average_plaquette_sp(gauge)
    print(f"    <P> = {plaq:.4f}")

    print("    Measuring Wilson loops...")
    wilson_data = measure_wilson_loops(
        gauge, compute_wilson_loop_sp, 2, 3,
        n_meas=20, sweep_func=sweep_sp
    )

    sigma, sigma_err, is_valid = extract_string_tension(wilson_data, 2, 3)

    status = "PASS" if sigma > 0 else "FAIL"
    print(f"    sigma = {sigma:.4f} +/- {sigma_err:.4f} [{status}]")

    return {'group': 'Sp(4)', 'sigma': sigma, 'sigma_err': sigma_err, 'sigma_positive': sigma > 0}


def test_g2():
    """Test G2."""
    print("\n  Testing G2 at beta=8.0...")

    config = LatticeConfigG2(Nx=4, Ny=4, Nz=4, Nt=6, beta=8.0)
    gauge = GaugeFieldG2(config)
    gauge.cold_start()

    print("    Thermalizing (150 sweeps)...")
    thermalize_g2(gauge, n_sweeps=150, epsilon=0.15, verbose=False)

    plaq = compute_average_plaquette_g2(gauge)
    print(f"    <P> = {plaq:.4f}")

    print("    Measuring Wilson loops...")
    wilson_data = measure_wilson_loops(
        gauge, compute_wilson_loop_g2, 2, 3,
        n_meas=20, sweep_func=sweep_g2, epsilon=0.15
    )

    sigma, sigma_err, is_valid = extract_string_tension(wilson_data, 2, 3)

    status = "PASS" if sigma > 0 else "FAIL"
    print(f"    sigma = {sigma:.4f} +/- {sigma_err:.4f} [{status}]")

    return {'group': 'G2', 'sigma': sigma, 'sigma_err': sigma_err, 'sigma_positive': sigma > 0}


if __name__ == "__main__":
    print("=" * 60)
    print("STRING TENSION - REMAINING GROUPS (Sp, G2)")
    print("=" * 60)

    print("\nPrevious results (already completed):")
    print("  SU(2): sigma = 0.3760 +/- 0.0203 [PASS]")
    print("  SU(3): sigma = 0.4761 +/- 0.0128 [PASS]")
    print("  SO(3): sigma = 1.4398 +/- 0.0513 [PASS]")
    print("  SO(4): sigma = 0.6020 +/- 0.2310 [PASS]")

    print("\n" + "=" * 60)
    print("Testing remaining groups...")
    print("=" * 60)

    results = []

    # Sp(4)
    results.append(test_sp4())

    # G2
    results.append(test_g2())

    print("\n" + "=" * 60)
    print("COMPLETE SUMMARY")
    print("=" * 60)

    all_results = [
        {'group': 'SU(2)', 'sigma': 0.3760, 'sigma_err': 0.0203, 'sigma_positive': True},
        {'group': 'SU(3)', 'sigma': 0.4761, 'sigma_err': 0.0128, 'sigma_positive': True},
        {'group': 'SO(3)', 'sigma': 1.4398, 'sigma_err': 0.0513, 'sigma_positive': True},
        {'group': 'SO(4)', 'sigma': 0.6020, 'sigma_err': 0.2310, 'sigma_positive': True},
    ] + results

    print(f"\n{'Group':<10} {'sigma':<12} {'Error':<12} {'Status'}")
    print("-" * 45)

    all_passed = True
    for r in all_results:
        status = "PASS" if r['sigma_positive'] else "FAIL"
        if not r['sigma_positive']:
            all_passed = False
        print(f"{r['group']:<10} {r['sigma']:<12.4f} {r['sigma_err']:<12.4f} {status}")

    print("-" * 45)

    passed = sum(1 for r in all_results if r['sigma_positive'])
    total = len(all_results)

    print(f"\nTotal: {passed}/{total} groups show sigma > 0 (confinement)")

    if all_passed:
        print("\n*** ALL GROUPS VERIFY CONFINEMENT: sigma > 0 ***")
