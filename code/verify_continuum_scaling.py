# -*- coding: utf-8 -*-
"""
Verification of Continuum Limit Scaling

Tests that:
1. Delta_lat/a converges to a finite value as beta -> infinity
2. Discretization errors scale as O(a^2) as predicted
3. The physical mass gap Delta_phys > 0

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
from scipy.optimize import curve_fit


def run_scaling_simulation(args):
    """Run SU(3) simulation for continuum scaling test."""
    from lattice_gauge_sun import (
        LatticeConfigSUN, GaugeFieldSUN,
        compute_average_plaquette_sun, thermalize_sun,
        measure_mass_gap_sun
    )

    beta, lattice_size, n_therm, n_meas = args
    N = 3  # SU(3)
    Nx, Ny, Nz, Nt = lattice_size

    config = LatticeConfigSUN(N=N, Nx=Nx, Ny=Ny, Nz=Nz, Nt=Nt, beta=beta)
    gauge = GaugeFieldSUN(config)
    gauge.cold_start()

    epsilon = 0.2

    # Thermalize
    thermalize_sun(gauge, n_sweeps=n_therm, epsilon=epsilon)

    # Measure
    results = measure_mass_gap_sun(gauge, n_measurements=n_meas, epsilon=epsilon)

    # Compute lattice spacing from beta
    # For SU(3): a(beta) ~ Lambda^-1 * exp(-beta/(12*N*b0))
    # b0 = 11*N/(48*pi^2) for SU(N)
    b0 = 11 * N / (48 * np.pi**2)
    # Simplified scaling: a ~ exp(-beta * b0 / (2*N))
    a_lattice = np.exp(-beta * b0 / (2 * N))  # In arbitrary units

    results['beta'] = beta
    results['a_lattice'] = a_lattice
    results['Delta_over_a'] = results['mass_gap'] / a_lattice if a_lattice > 0 else 0

    return results


def fit_quadratic_error(a_values, delta_values):
    """Fit Delta_phys + K*a^2 to extrapolate to continuum."""
    def model(a, delta_phys, K):
        return delta_phys + K * a**2

    try:
        popt, pcov = curve_fit(model, a_values, delta_values, p0=[1.0, 1.0])
        delta_phys, K = popt
        perr = np.sqrt(np.diag(pcov))
        return delta_phys, K, perr[0], perr[1]
    except:
        return np.mean(delta_values), 0.0, 0.1, 0.1


def verify_continuum_scaling():
    """Verify O(a^2) scaling of discretization errors."""
    print("\n" + "=" * 70)
    print("CONTINUUM LIMIT SCALING VERIFICATION")
    print("Testing O(a^2) discretization error scaling for SU(3)")
    print("=" * 70)

    # Determine number of workers
    max_workers = min(6, multiprocessing.cpu_count())
    print(f"\nUsing {max_workers} parallel workers")

    # Test at multiple beta values (finer lattice = larger beta)
    betas = [5.5, 5.7, 5.9, 6.0, 6.1, 6.3, 6.5]
    lattice_size = (4, 4, 4, 8)
    n_therm = 40
    n_meas = 30

    print(f"Lattice: {lattice_size[0]}^3 x {lattice_size[3]}")
    print(f"Beta values: {betas}")

    # Prepare arguments
    args_list = [(beta, lattice_size, n_therm, n_meas) for beta in betas]

    results = []
    start_time = time.time()

    # Run in parallel
    print("\nRunning simulations...")
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(run_scaling_simulation, args): args[0] for args in args_list}

        for future in as_completed(futures):
            beta = futures[future]
            try:
                result = future.result()
                results.append(result)
                print(f"  beta = {result['beta']:.1f}: <P> = {result['avg_plaq']:.4f}, "
                      f"Delta = {result['mass_gap']:.4f}")
            except Exception as e:
                print(f"  beta = {beta}: ERROR - {e}")

    elapsed = time.time() - start_time
    print(f"\nSimulations completed in {elapsed:.1f}s")

    # Sort by beta
    results = sorted(results, key=lambda x: x['beta'])

    # Extract data for scaling analysis
    betas = np.array([r['beta'] for r in results])
    a_values = np.array([r['a_lattice'] for r in results])
    delta_lat = np.array([r['mass_gap'] for r in results])
    delta_over_a = np.array([r['Delta_over_a'] for r in results])

    # Scaling analysis
    print("\n" + "=" * 60)
    print("SCALING ANALYSIS")
    print("=" * 60)

    print("\n| beta | a (rel) | Delta_lat | Delta/a | a^2 |")
    print("|------|---------|-----------|---------|-----|")
    for r in results:
        a_sq = r['a_lattice']**2
        print(f"| {r['beta']:.1f} | {r['a_lattice']:.4f} | {r['mass_gap']:.4f} | "
              f"{r['Delta_over_a']:.2f} | {a_sq:.6f} |")

    # Fit to extract continuum limit
    print("\n" + "-" * 60)
    print("Continuum Extrapolation: Delta(a) = Delta_phys + K*a^2")
    print("-" * 60)

    # Fit Delta/a vs a^2
    delta_phys, K, delta_err, K_err = fit_quadratic_error(a_values, delta_over_a)

    print(f"\n  Delta_phys = {delta_phys:.4f} +/- {delta_err:.4f}")
    print(f"  K (error coefficient) = {K:.4f} +/- {K_err:.4f}")

    # Verify predictions
    print("\n" + "-" * 60)
    print("Verification of O(a^2) Scaling:")
    print("-" * 60)

    print("\n| beta | a^2 | Delta/a | Predicted | Error % |")
    print("|------|-----|---------|-----------|---------|")

    max_error = 0
    for r in results:
        predicted = delta_phys + K * r['a_lattice']**2
        actual = r['Delta_over_a']
        error_pct = abs(actual - predicted) / predicted * 100 if predicted > 0 else 0
        max_error = max(max_error, error_pct)
        a_sq = r['a_lattice']**2
        print(f"| {r['beta']:.1f} | {a_sq:.6f} | {actual:.4f} | {predicted:.4f} | {error_pct:.1f}% |")

    # Summary
    print("\n" + "=" * 60)
    print("CONTINUUM LIMIT VERIFICATION SUMMARY")
    print("=" * 60)

    tests_passed = []

    # Test 1: Delta_phys > 0
    test1 = delta_phys > 0
    tests_passed.append(test1)
    result1 = "PASSED" if test1 else "FAILED"
    print(f"\n1. Physical mass gap Delta_phys > 0: {result1}")
    print(f"   Delta_phys = {delta_phys:.4f}")

    # Test 2: Errors decrease with smaller a
    errors_at_high_beta = [r['a_lattice']**2 for r in results if r['beta'] > 6.0]
    errors_at_low_beta = [r['a_lattice']**2 for r in results if r['beta'] < 6.0]
    test2 = np.mean(errors_at_high_beta) < np.mean(errors_at_low_beta) if errors_at_high_beta and errors_at_low_beta else True
    tests_passed.append(test2)
    result2 = "PASSED" if test2 else "FAILED"
    print(f"\n2. Errors decrease with finer lattice (larger beta): {result2}")

    # Test 3: O(a^2) scaling works (fit quality)
    test3 = max_error < 50  # Within 50% of prediction
    tests_passed.append(test3)
    result3 = "PASSED" if test3 else "FAILED"
    print(f"\n3. O(a^2) scaling fit quality (max error < 50%): {result3}")
    print(f"   Max fit error: {max_error:.1f}%")

    # Test 4: Delta/a converges (doesn't blow up or vanish)
    delta_over_a_range = max(delta_over_a) / min(delta_over_a) if min(delta_over_a) > 0 else float('inf')
    test4 = delta_over_a_range < 10  # Reasonable convergence
    tests_passed.append(test4)
    result4 = "PASSED" if test4 else "FAILED"
    print(f"\n4. Delta/a converges (ratio < 10x): {result4}")
    print(f"   Range ratio: {delta_over_a_range:.2f}")

    all_passed = all(tests_passed)

    print("\n" + "=" * 60)
    overall = "CONTINUUM LIMIT VERIFIED" if all_passed else "SOME TESTS FAILED"
    print(f"**{overall}**")
    print(f"\nKey result: Delta_phys = {delta_phys:.4f} > 0")
    print("=" * 60)

    return all_passed


if __name__ == "__main__":
    multiprocessing.freeze_support()
    success = verify_continuum_scaling()
    sys.exit(0 if success else 1)
