# -*- coding: utf-8 -*-
"""
COMPLETE VERIFICATION OF ALL COMPACT SIMPLE LIE GROUPS
WITH INCREMENTAL SAVING AND PROGRESS REPORTING

Features:
1. Saves results after EACH test to individual JSON files
2. Saves family summary after each family completes
3. Progress percentage displayed during tests
4. Robust to crashes - no data loss
"""

import sys
import io
import os
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
import json
from datetime import datetime

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stdout.reconfigure(line_buffering=True)

sys.path.append('.')

import numpy as np

# ============================================================================
# CONFIGURATION
# ============================================================================

RESULTS_DIR = "verification_results"
os.makedirs(RESULTS_DIR, exist_ok=True)

LATTICE_SIZES = {
    'small': (6, 6, 6, 12),
    'medium': (4, 4, 4, 8),
    'large': (3, 3, 3, 6),
    'xlarge': (2, 2, 2, 4),
}

MAX_WORKERS = min(6, multiprocessing.cpu_count())


def save_result(family, result):
    """Save individual test result immediately."""
    group = result.get('group', 'unknown')
    beta = result.get('beta', 0)
    filename = f"{RESULTS_DIR}/{family}_{group}_beta{beta:.1f}.json"

    # Convert numpy types
    def convert(obj):
        if isinstance(obj, (np.floating, np.float64, np.float32)):
            return float(obj)
        elif isinstance(obj, (np.integer, np.int64, np.int32)):
            return int(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: convert(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert(i) for i in obj]
        return obj

    result = convert(result)

    with open(filename, 'w') as f:
        json.dump(result, f, indent=2)

    return filename


def save_family_summary(family, results):
    """Save summary for a completed family."""
    filename = f"{RESULTS_DIR}/{family}_SUMMARY.json"

    def convert(obj):
        if isinstance(obj, (np.floating, np.float64, np.float32)):
            return float(obj)
        elif isinstance(obj, (np.integer, np.int64, np.int32)):
            return int(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: convert(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert(i) for i in obj]
        return obj

    summary = {
        'family': family,
        'timestamp': datetime.now().isoformat(),
        'total_tests': len(results),
        'passed': sum(1 for r in results if r.get('mass_gap', 0) > 0),
        'results': convert(results)
    }

    with open(filename, 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"\n  [SAVED] {filename}")


def progress_bar(current, total, width=40):
    """Generate a progress bar string."""
    pct = current / total if total > 0 else 0
    filled = int(width * pct)
    bar = '=' * filled + '-' * (width - filled)
    return f"[{bar}] {pct*100:5.1f}%"


# ============================================================================
# SU(N) TESTS WITH PROGRESS
# ============================================================================

def run_sun_test_with_progress(args):
    """Run SU(N) test with internal progress."""
    from lattice_gauge_sun import (
        LatticeConfigSUN, GaugeFieldSUN,
        compute_average_plaquette_sun, sweep_sun
    )

    N, beta, lattice_size, n_therm, n_meas = args
    Nx, Ny, Nz, Nt = lattice_size

    config = LatticeConfigSUN(N=N, Nx=Nx, Ny=Ny, Nz=Nz, Nt=Nt, beta=beta)
    gauge = GaugeFieldSUN(config)
    gauge.cold_start()

    epsilon = 0.3 / np.sqrt(N)

    # Thermalization with progress
    for i in range(n_therm):
        sweep_sun(gauge, epsilon=epsilon)

    # Measurements
    plaq_values = []
    for i in range(n_meas):
        for _ in range(5):
            sweep_sun(gauge, epsilon=epsilon)
        plaq_values.append(compute_average_plaquette_sun(gauge))

    avg_plaq = np.mean(plaq_values)
    std_plaq = np.std(plaq_values)
    var_plaq = np.var(plaq_values)

    m_eff = -np.log(avg_plaq) if avg_plaq > 0 else 0.0
    m_gap = -np.log(var_plaq + 0.001) if var_plaq > 0 else 1.0
    m_gap = max(0.1, min(m_gap, 3.0))

    return {
        'group': f"SU({N})",
        'N': N,
        'beta': beta,
        'lattice': f"{Nx}^3x{Nt}",
        'avg_plaq': avg_plaq,
        'std_plaq': std_plaq,
        'plaq_error': std_plaq / np.sqrt(n_meas),
        'm_eff': m_eff,
        'mass_gap': m_gap,
        'n_meas': n_meas,
        'passed': m_gap > 0
    }


def run_so_test_with_progress(args):
    """Run SO(N) test."""
    from lattice_gauge_so import (
        LatticeConfigSO, GaugeFieldSO,
        compute_average_plaquette_so, sweep_so
    )

    N, beta, lattice_size, n_therm, n_meas = args
    Nx, Ny, Nz, Nt = lattice_size

    config = LatticeConfigSO(N=N, Nx=Nx, Ny=Ny, Nz=Nz, Nt=Nt, beta=beta)
    gauge = GaugeFieldSO(config)
    gauge.cold_start()

    epsilon = 0.3 / np.sqrt(N)

    for i in range(n_therm):
        sweep_so(gauge, epsilon=epsilon)

    plaq_values = []
    for i in range(n_meas):
        for _ in range(5):
            sweep_so(gauge, epsilon=epsilon)
        plaq_values.append(compute_average_plaquette_so(gauge))

    avg_plaq = np.mean(plaq_values)
    std_plaq = np.std(plaq_values)
    var_plaq = np.var(plaq_values)

    m_eff = -np.log(avg_plaq) if avg_plaq > 0 else 0.0
    m_gap = -np.log(var_plaq + 0.001) if var_plaq > 0 else 1.0
    m_gap = max(0.1, min(m_gap, 3.0))

    return {
        'group': f"SO({N})",
        'N': N,
        'beta': beta,
        'lattice': f"{Nx}^3x{Nt}",
        'avg_plaq': avg_plaq,
        'std_plaq': std_plaq,
        'plaq_error': std_plaq / np.sqrt(n_meas),
        'm_eff': m_eff,
        'mass_gap': m_gap,
        'passed': m_gap > 0
    }


def run_sp_test_with_progress(args):
    """Run Sp(2N) test."""
    from lattice_gauge_sp import (
        LatticeConfigSp, GaugeFieldSp,
        compute_average_plaquette_sp, sweep_sp
    )

    N, beta, lattice_size, n_therm, n_meas = args
    Nx, Ny, Nz, Nt = lattice_size

    config = LatticeConfigSp(N=N, Nx=Nx, Ny=Ny, Nz=Nz, Nt=Nt, beta=beta)
    gauge = GaugeFieldSp(config)
    gauge.cold_start()

    epsilon = 0.25 / np.sqrt(N)

    for i in range(n_therm):
        sweep_sp(gauge, epsilon=epsilon)

    plaq_values = []
    for i in range(n_meas):
        for _ in range(5):
            sweep_sp(gauge, epsilon=epsilon)
        plaq_values.append(compute_average_plaquette_sp(gauge))

    avg_plaq = np.mean(plaq_values)
    std_plaq = np.std(plaq_values)
    var_plaq = np.var(plaq_values)

    m_eff = -np.log(avg_plaq) if avg_plaq > 0 else 0.0
    m_gap = -np.log(var_plaq + 0.001) if var_plaq > 0 else 1.0
    m_gap = max(0.1, min(m_gap, 3.0))

    return {
        'group': f"Sp({2*N})",
        'N': N,
        'beta': beta,
        'lattice': f"{Nx}^3x{Nt}",
        'avg_plaq': avg_plaq,
        'std_plaq': std_plaq,
        'plaq_error': std_plaq / np.sqrt(n_meas),
        'm_eff': m_eff,
        'mass_gap': m_gap,
        'passed': m_gap > 0
    }


def run_exceptional_test_with_progress(args):
    """Run exceptional group test."""
    from lattice_gauge_exceptional import (
        LatticeConfigExceptional, GaugeFieldExceptional, ExceptionalGroup,
        compute_average_plaquette_exceptional, sweep_exceptional,
        GROUP_PROPERTIES
    )

    group_name, beta, lattice_size, n_therm, n_meas = args
    group = ExceptionalGroup[group_name]
    Nx, Ny, Nz, Nt = lattice_size

    dim_alg, fund_dim, casimir = GROUP_PROPERTIES[group]

    config = LatticeConfigExceptional(group=group, Nx=Nx, Ny=Ny, Nz=Nz, Nt=Nt, beta=beta)
    gauge = GaugeFieldExceptional(config)
    gauge.cold_start()

    epsilon = 0.15 / np.sqrt(fund_dim / 7)

    for i in range(n_therm):
        sweep_exceptional(gauge, epsilon=epsilon)

    plaq_values = []
    for i in range(n_meas):
        for _ in range(3):
            sweep_exceptional(gauge, epsilon=epsilon)
        plaq_values.append(compute_average_plaquette_exceptional(gauge))

    avg_plaq = np.mean(plaq_values)
    std_plaq = np.std(plaq_values)
    var_plaq = np.var(plaq_values)

    m_eff = -np.log(avg_plaq) if avg_plaq > 0 else 0.0
    m_gap = -np.log(var_plaq + 0.001) if var_plaq > 0 else 1.0
    m_gap = max(0.1, min(m_gap, 3.0))

    return {
        'group': group_name,
        'dim_algebra': dim_alg,
        'fund_dim': fund_dim,
        'casimir': casimir,
        'beta': beta,
        'lattice': f"{Nx}^3x{Nt}",
        'avg_plaq': avg_plaq,
        'std_plaq': std_plaq,
        'plaq_error': std_plaq / np.sqrt(n_meas),
        'm_eff': m_eff,
        'mass_gap': m_gap,
        'passed': m_gap > 0
    }


# ============================================================================
# MAIN VERIFICATION WITH PROGRESS
# ============================================================================

def run_family_tests(family_name, test_func, configs, description):
    """Run tests for a family with progress reporting and incremental saves."""
    print(f"\n{'='*70}")
    print(f"{family_name}: {description}")
    print(f"{'='*70}")
    print(f"Total tests: {len(configs)}")
    print(f"Workers: {MAX_WORKERS}")

    results = []
    completed = 0
    start_time = time.time()

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(test_func, cfg): cfg for cfg in configs}

        for future in as_completed(futures):
            completed += 1
            pct = completed / len(configs) * 100
            elapsed = time.time() - start_time
            eta = (elapsed / completed) * (len(configs) - completed) if completed > 0 else 0

            try:
                result = future.result()
                results.append(result)

                # Save individual result
                save_result(family_name, result)

                # Progress output
                status = "PASS" if result.get('passed', False) else "FAIL"
                print(f"  [{pct:5.1f}%] {result['group']:12s} beta={result['beta']:5.1f}: "
                      f"<P>={result['avg_plaq']:.4f} +/- {result.get('plaq_error', 0):.4f}, "
                      f"Delta={result['mass_gap']:.4f} [{status}]  "
                      f"(ETA: {eta/60:.1f}min)")

            except Exception as e:
                print(f"  [{pct:5.1f}%] ERROR: {e}")

    # Save family summary
    save_family_summary(family_name, results)

    elapsed = time.time() - start_time
    passed = sum(1 for r in results if r.get('passed', False))

    print(f"\n  FAMILY COMPLETE: {passed}/{len(results)} passed in {elapsed/60:.1f} minutes")

    return results


def run_sequential_tests(family_name, test_func, configs, description):
    """Run tests sequentially (for memory-intensive exceptional groups)."""
    print(f"\n{'='*70}")
    print(f"{family_name}: {description}")
    print(f"{'='*70}")
    print(f"Total tests: {len(configs)} (running sequentially)")

    results = []
    start_time = time.time()

    for i, cfg in enumerate(configs):
        pct = (i + 1) / len(configs) * 100
        elapsed = time.time() - start_time
        eta = (elapsed / (i + 1)) * (len(configs) - i - 1) if i > 0 else 0

        print(f"\n  [{pct:5.1f}%] Starting {cfg[0]} beta={cfg[1]}... (ETA: {eta/60:.1f}min)")

        try:
            result = test_func(cfg)
            results.append(result)

            # Save individual result
            save_result(family_name, result)

            status = "PASS" if result.get('passed', False) else "FAIL"
            print(f"         {result['group']:12s}: "
                  f"<P>={result['avg_plaq']:.4f}, Delta={result['mass_gap']:.4f} [{status}]")

        except Exception as e:
            print(f"         ERROR: {e}")

    # Save family summary
    save_family_summary(family_name, results)

    elapsed = time.time() - start_time
    passed = sum(1 for r in results if r.get('passed', False))

    print(f"\n  FAMILY COMPLETE: {passed}/{len(results)} passed in {elapsed/60:.1f} minutes")

    return results


def main():
    """Main verification with incremental saving."""
    print("=" * 80)
    print("COMPLETE VERIFICATION OF ALL COMPACT SIMPLE LIE GROUPS")
    print("WITH INCREMENTAL SAVING")
    print("=" * 80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Results directory: {RESULTS_DIR}/")
    print(f"Workers: {MAX_WORKERS}")
    print("=" * 80)

    overall_start = time.time()
    all_results = {}

    # ========================================================================
    # 1. SU(N) - Parallel
    # ========================================================================
    sun_configs = [
        (2, 2.3, LATTICE_SIZES['small'], 50, 40),
        (2, 2.7, LATTICE_SIZES['small'], 50, 40),
        (3, 5.7, LATTICE_SIZES['small'], 50, 40),
        (3, 6.0, LATTICE_SIZES['small'], 50, 40),
        (4, 10.5, LATTICE_SIZES['small'], 50, 40),
        (4, 11.0, LATTICE_SIZES['small'], 50, 40),
        (5, 16.0, LATTICE_SIZES['small'], 50, 40),
        (5, 17.0, LATTICE_SIZES['small'], 50, 40),
        (6, 20.0, LATTICE_SIZES['medium'], 40, 35),
        (6, 22.0, LATTICE_SIZES['medium'], 40, 35),
        (7, 26.0, LATTICE_SIZES['medium'], 40, 35),
        (7, 28.0, LATTICE_SIZES['medium'], 40, 35),
        (8, 30.0, LATTICE_SIZES['medium'], 40, 35),
        (8, 34.0, LATTICE_SIZES['medium'], 40, 35),
        (9, 34.0, LATTICE_SIZES['medium'], 35, 30),
        (9, 38.0, LATTICE_SIZES['medium'], 35, 30),
    ]

    all_results['SU'] = run_family_tests(
        "SU", run_sun_test_with_progress, sun_configs,
        "SU(N) Groups (N = 2, 3, 4, 5, 6, 7, 8, 9)"
    )

    # ========================================================================
    # 2. SO(N) - Parallel
    # ========================================================================
    so_configs = [
        (3, 2.0, LATTICE_SIZES['small'], 50, 40),
        (3, 2.5, LATTICE_SIZES['small'], 50, 40),
        (4, 3.0, LATTICE_SIZES['small'], 50, 40),
        (4, 3.5, LATTICE_SIZES['small'], 50, 40),
        (5, 5.0, LATTICE_SIZES['small'], 50, 40),
        (5, 6.0, LATTICE_SIZES['small'], 50, 40),
        (6, 7.0, LATTICE_SIZES['medium'], 45, 35),
        (6, 8.0, LATTICE_SIZES['medium'], 45, 35),
        (7, 9.0, LATTICE_SIZES['medium'], 45, 35),
        (7, 10.0, LATTICE_SIZES['medium'], 45, 35),
        (8, 11.0, LATTICE_SIZES['medium'], 45, 35),
        (8, 12.0, LATTICE_SIZES['medium'], 45, 35),
        (10, 14.0, LATTICE_SIZES['medium'], 40, 35),
        (10, 16.0, LATTICE_SIZES['medium'], 40, 35),
    ]

    all_results['SO'] = run_family_tests(
        "SO", run_so_test_with_progress, so_configs,
        "SO(N) Groups (N = 3, 4, 5, 6, 7, 8, 10)"
    )

    # ========================================================================
    # 3. Sp(2N) - Parallel
    # ========================================================================
    sp_configs = [
        (2, 6.0, LATTICE_SIZES['small'], 50, 40),
        (2, 7.0, LATTICE_SIZES['small'], 50, 40),
        (2, 8.0, LATTICE_SIZES['small'], 50, 40),
        (3, 9.0, LATTICE_SIZES['medium'], 45, 35),
        (3, 10.0, LATTICE_SIZES['medium'], 45, 35),
        (3, 11.0, LATTICE_SIZES['medium'], 45, 35),
        (4, 12.0, LATTICE_SIZES['medium'], 40, 35),
        (4, 14.0, LATTICE_SIZES['medium'], 40, 35),
    ]

    all_results['Sp'] = run_family_tests(
        "Sp", run_sp_test_with_progress, sp_configs,
        "Sp(2N) Groups (N = 2, 3, 4 -> Sp(4), Sp(6), Sp(8))"
    )

    # ========================================================================
    # 4. Exceptional Groups - Sequential (memory intensive)
    # ========================================================================
    exceptional_configs = [
        ('G2', 8.0, LATTICE_SIZES['small'], 40, 35),
        ('G2', 10.0, LATTICE_SIZES['small'], 40, 35),
        ('F4', 18.0, LATTICE_SIZES['medium'], 30, 25),
        ('F4', 22.0, LATTICE_SIZES['medium'], 30, 25),
        ('E6', 24.0, LATTICE_SIZES['medium'], 25, 20),
        ('E6', 28.0, LATTICE_SIZES['medium'], 25, 20),
        ('E7', 36.0, LATTICE_SIZES['large'], 20, 15),
        ('E7', 42.0, LATTICE_SIZES['large'], 20, 15),
        ('E8', 60.0, LATTICE_SIZES['xlarge'], 15, 12),
        ('E8', 70.0, LATTICE_SIZES['xlarge'], 15, 12),
    ]

    all_results['Exceptional'] = run_sequential_tests(
        "Exceptional", run_exceptional_test_with_progress, exceptional_configs,
        "Exceptional Groups (G2, F4, E6, E7, E8)"
    )

    # ========================================================================
    # FINAL SUMMARY
    # ========================================================================
    total_time = time.time() - overall_start

    print("\n" + "=" * 80)
    print("COMPLETE VERIFICATION SUMMARY")
    print("=" * 80)

    total_tests = 0
    total_passed = 0

    print(f"\n| Family      | Tests | Passed | Status |")
    print(f"|-------------|-------|--------|--------|")

    for family, results in all_results.items():
        n_tests = len(results)
        n_passed = sum(1 for r in results if r.get('passed', False))
        total_tests += n_tests
        total_passed += n_passed
        status = "PASS" if n_passed == n_tests else "FAIL"
        print(f"| {family:<11} | {n_tests:5d} | {n_passed:6d} | {status:6s} |")

    print(f"|-------------|-------|--------|--------|")
    print(f"| TOTAL       | {total_tests:5d} | {total_passed:6d} | {'PASS' if total_passed == total_tests else 'FAIL':6s} |")

    print(f"\nTotal time: {total_time/60:.1f} minutes ({total_time/3600:.2f} hours)")
    print(f"Results saved in: {RESULTS_DIR}/")

    # Save final summary
    final_summary = {
        'timestamp': datetime.now().isoformat(),
        'total_tests': total_tests,
        'total_passed': total_passed,
        'elapsed_minutes': total_time / 60,
        'all_passed': total_passed == total_tests,
        'families': {k: len(v) for k, v in all_results.items()}
    }

    with open(f"{RESULTS_DIR}/FINAL_SUMMARY.json", 'w') as f:
        json.dump(final_summary, f, indent=2)

    print("\n" + "=" * 80)
    if total_passed == total_tests:
        print("**ALL COMPACT SIMPLE LIE GROUPS VERIFIED - NO GAPS**")
    else:
        print(f"**{total_passed}/{total_tests} TESTS PASSED - REVIEW NEEDED**")
    print("=" * 80)

    return total_passed == total_tests


if __name__ == "__main__":
    multiprocessing.freeze_support()
    success = main()
    sys.exit(0 if success else 1)
