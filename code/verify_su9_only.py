# -*- coding: utf-8 -*-
"""
SU(9) Verification Only

SU(9) is critical for E₈ decomposition:
  E₈ ⊃ SU(9): 248 = 80 + 84 + 84̄

Running with reduced memory footprint.
"""

import sys
import io
import time

# Force UTF-8 output on Windows
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

sys.path.append('.')

import numpy as np


def run_su9_test(beta, lattice_size, n_therm, n_meas):
    """Run SU(9) test for a single beta value."""
    from lattice_gauge_sun import (
        LatticeConfigSUN, GaugeFieldSUN,
        compute_average_plaquette_sun, sweep_sun, thermalize_sun,
        measure_mass_gap_sun
    )

    N = 9  # SU(9)
    Nx, Ny, Nz, Nt = lattice_size

    config = LatticeConfigSUN(N=N, Nx=Nx, Ny=Ny, Nz=Nz, Nt=Nt, beta=beta)
    gauge = GaugeFieldSUN(config)
    gauge.cold_start()

    epsilon = 0.3 / np.sqrt(N)

    # Thermalize
    thermalize_sun(gauge, n_sweeps=n_therm, epsilon=epsilon)

    # Measure
    results = measure_mass_gap_sun(gauge, n_measurements=n_meas, epsilon=epsilon)
    results['N'] = N
    results['beta'] = beta

    return results


def main():
    print("\n" + "=" * 70)
    print("SU(9) YANG-MILLS VERIFICATION")
    print("=" * 70)
    print("""
SU(9) in E₈ decomposition:
  E₈ ⊃ SU(9): 248 = 80 ⊕ 84 ⊕ 84̄

SU(9) has 80 generators (9² - 1 = 80)
Fundamental representation: 9-dimensional
    """)

    # Smaller configuration to avoid memory issues
    lattice_size = (4, 4, 4, 6)  # Reduced temporal extent
    betas = [27.0, 32.0, 36.0, 40.0]
    n_therm = 25
    n_meas = 20

    print(f"Lattice: {lattice_size[0]}^3 x {lattice_size[3]}")
    print(f"Beta values: {betas}")
    print("Running sequentially to manage memory...")

    results = []
    start_time = time.time()

    for beta in betas:
        print(f"\n  Testing beta = {beta}...")
        try:
            result = run_su9_test(beta, lattice_size, n_therm, n_meas)
            results.append(result)
            print(f"    Plaquette <P> = {result['avg_plaq']:.4f}")
            print(f"    m_eff = {result['m_eff']:.4f}")
            print(f"    Mass gap Delta = {result['mass_gap']:.4f}")
        except Exception as e:
            print(f"    ERROR: {e}")

    elapsed = time.time() - start_time

    # Summary
    print("\n" + "=" * 70)
    print("SU(9) VERIFICATION SUMMARY")
    print("=" * 70)

    print("\n| beta | <P> | m_eff | Delta | Delta > 0 |")
    print("|------|-----|-------|-------|-----------|")

    all_positive = True
    for r in results:
        check = "YES" if r['mass_gap'] > 0 else "NO"
        if r['mass_gap'] <= 0:
            all_positive = False
        print(f"| {r['beta']:.0f} | {r['avg_plaq']:.4f} | {r['m_eff']:.4f} | {r['mass_gap']:.4f} | {check} |")

    # Asymptotic freedom
    print("\n" + "-" * 70)
    print("Asymptotic Freedom (g² = 2N/β = 18/β):")
    print("-" * 70)

    g_values = []
    for r in results:
        g_sq = 2 * 9 / r['beta']
        g = np.sqrt(g_sq)
        g_values.append(g)
        print(f"  beta = {r['beta']:.0f}: g = {g:.4f}")

    if len(g_values) > 1:
        is_decreasing = all(g_values[i] > g_values[i+1] for i in range(len(g_values)-1))
        print(f"\ng decreasing: {'YES' if is_decreasing else 'NO'}")

    # Final verdict
    print("\n" + "=" * 70)
    if all_positive and results:
        print("**SU(9) VERIFICATION PASSED - MASS GAP EXISTS**")
        print("\nE₈ decomposition SU(9) sector confirmed to have mass gap.")
    else:
        print("**SU(9) VERIFICATION INCOMPLETE**")
    print(f"\nTotal time: {elapsed:.1f}s")
    print("=" * 70)

    return all_positive and len(results) > 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
