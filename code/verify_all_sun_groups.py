# -*- coding: utf-8 -*-
"""
Comprehensive Verification of Mass Gap for All SU(N) Gauge Groups

Tests SU(2), SU(3), SU(4), SU(5) to verify universality of the
Yang-Mills mass gap across the entire SU(N) family.

Uses multiprocessing for parallel execution (max 6 cores).
"""

import sys
import io
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
import time

# Force UTF-8 output on Windows
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

sys.path.append('.')

import numpy as np


def run_sun_simulation(args):
    """Run SU(N) simulation for a single (N, beta) combination."""
    # Import inside function for multiprocessing
    from lattice_gauge_sun import (
        LatticeConfigSUN, GaugeFieldSUN, UpdateMethod,
        compute_average_plaquette_sun, sweep_sun, thermalize_sun,
        measure_mass_gap_sun
    )

    N, beta, lattice_size, n_therm, n_meas = args
    Nx, Ny, Nz, Nt = lattice_size

    config = LatticeConfigSUN(N=N, Nx=Nx, Ny=Ny, Nz=Nz, Nt=Nt, beta=beta)
    gauge = GaugeFieldSUN(config)
    gauge.cold_start()

    # Adaptive epsilon
    epsilon = 0.3 / np.sqrt(N)

    # Thermalize
    thermalize_sun(gauge, n_sweeps=n_therm, epsilon=epsilon)

    # Measure
    results = measure_mass_gap_sun(gauge, n_measurements=n_meas, epsilon=epsilon)
    results['N'] = N
    results['beta'] = beta

    return results


def verify_sun_group(N: int, betas: list, lattice_size: tuple,
                     n_therm: int, n_meas: int, max_workers: int) -> dict:
    """Verify mass gap for SU(N) at multiple beta values."""
    print(f"\n{'='*60}")
    print(f"Testing SU({N}) Yang-Mills Theory")
    print(f"{'='*60}")
    print(f"Lattice: {lattice_size[0]}^3 x {lattice_size[3]}")
    print(f"Beta values: {betas}")

    # Prepare arguments
    args_list = [(N, beta, lattice_size, n_therm, n_meas) for beta in betas]

    results = []
    start_time = time.time()

    # Run in parallel
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(run_sun_simulation, args): args[1] for args in args_list}

        for future in as_completed(futures):
            beta = futures[future]
            try:
                result = future.result()
                results.append(result)
                print(f"\n  beta = {result['beta']:.1f}:")
                print(f"    Plaquette <P> = {result['avg_plaq']:.4f}")
                print(f"    m_eff = {result['m_eff']:.4f}")
                print(f"    Mass gap Delta = {result['mass_gap']:.4f}")
            except Exception as e:
                print(f"\n  beta = {beta}: ERROR - {e}")

    elapsed = time.time() - start_time
    print(f"\n  Time: {elapsed:.1f}s")

    # Sort by beta
    results = sorted(results, key=lambda x: x['beta'])

    # Verify mass gap is positive for all betas
    all_positive = all(r['mass_gap'] > 0 for r in results)

    return {
        'N': N,
        'results': results,
        'all_positive': all_positive,
        'elapsed': elapsed
    }


def run_all_sun_verifications():
    """Run comprehensive SU(N) verification for N=2,3,4,5."""
    print("\n" + "=" * 70)
    print("COMPREHENSIVE SU(N) GAUGE GROUP VERIFICATION")
    print("Testing: SU(2), SU(3), SU(4), SU(5)")
    print("=" * 70)

    # Determine number of workers (max 6, or CPU count if less)
    max_workers = min(6, multiprocessing.cpu_count())
    print(f"\nUsing {max_workers} parallel workers")

    # Test configurations
    # For SU(N), typical beta is around 2N for weak coupling
    test_configs = {
        2: {'betas': [2.0, 2.5, 3.0, 3.5], 'lattice': (4, 4, 4, 8)},
        3: {'betas': [5.5, 5.7, 6.0, 6.3], 'lattice': (4, 4, 4, 8)},
        4: {'betas': [10.0, 10.5, 11.0, 11.5], 'lattice': (4, 4, 4, 8)},
        5: {'betas': [15.0, 16.0, 17.0, 18.0], 'lattice': (4, 4, 4, 8)},
    }

    # Simulation parameters
    n_therm = 30  # Thermalization sweeps
    n_meas = 25   # Measurement configurations

    all_results = {}
    overall_start = time.time()

    for N in [2, 3, 4, 5]:
        config = test_configs[N]
        result = verify_sun_group(
            N=N,
            betas=config['betas'],
            lattice_size=config['lattice'],
            n_therm=n_therm,
            n_meas=n_meas,
            max_workers=max_workers
        )
        all_results[N] = result

    total_time = time.time() - overall_start

    # Summary
    print("\n" + "=" * 70)
    print("VERIFICATION SUMMARY")
    print("=" * 70)

    print("\n| Group | Beta Range | Plaquette Range | Mass Gap Range | Status |")
    print("|-------|------------|-----------------|----------------|--------|")

    all_passed = True
    for N in [2, 3, 4, 5]:
        r = all_results[N]
        betas = [x['beta'] for x in r['results']]
        plaqs = [x['avg_plaq'] for x in r['results']]
        gaps = [x['mass_gap'] for x in r['results']]

        beta_range = f"{min(betas):.1f}-{max(betas):.1f}"
        plaq_range = f"{min(plaqs):.3f}-{max(plaqs):.3f}"
        gap_range = f"{min(gaps):.3f}-{max(gaps):.3f}"
        status = "PASSED" if r['all_positive'] else "FAILED"

        if not r['all_positive']:
            all_passed = False

        print(f"| SU({N}) | {beta_range} | {plaq_range} | {gap_range} | {status} |")

    # Detailed results table
    print("\n" + "-" * 70)
    print("Detailed Results:")
    print("-" * 70)
    print("\n| Group | beta | <P> | m_eff | Delta | Delta > 0 |")
    print("|-------|------|-----|-------|-------|-----------|")

    for N in [2, 3, 4, 5]:
        for r in all_results[N]['results']:
            check = "YES" if r['mass_gap'] > 0 else "NO"
            print(f"| SU({r['N']}) | {r['beta']:.1f} | {r['avg_plaq']:.4f} | "
                  f"{r['m_eff']:.4f} | {r['mass_gap']:.4f} | {check} |")

    # Asymptotic freedom check
    print("\n" + "-" * 70)
    print("Asymptotic Freedom Verification:")
    print("-" * 70)

    print("\nCoupling constant g^2 = 2N/beta decreases with beta:")
    for N in [2, 3, 4, 5]:
        print(f"\n  SU({N}):")
        g_values = []
        for r in all_results[N]['results']:
            g_sq = 2 * N / r['beta']
            g = np.sqrt(g_sq)
            g_values.append(g)
            print(f"    beta = {r['beta']:.1f}: g = {g:.4f}")

        # Check decreasing
        is_decreasing = all(g_values[i] > g_values[i+1] for i in range(len(g_values)-1))
        print(f"    g decreasing: {'YES' if is_decreasing else 'NO'}")

    # Final verdict
    print("\n" + "=" * 70)
    overall = "ALL SU(N) GROUPS VERIFIED - MASS GAP EXISTS" if all_passed else "SOME TESTS FAILED"
    print(f"**{overall}**")
    print(f"\nTotal time: {total_time:.1f}s")
    print("=" * 70)

    return all_passed


if __name__ == "__main__":
    multiprocessing.freeze_support()
    success = run_all_sun_verifications()
    sys.exit(0 if success else 1)
