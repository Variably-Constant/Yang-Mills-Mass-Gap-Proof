#!/usr/bin/env python3
"""
Main Simulation Runner for Yang-Mills Mass Gap Analysis

This script orchestrates complete lattice gauge theory simulations:
1. Configuration generation (thermalization)
2. Measurement collection
3. Statistical analysis
4. Mass gap extraction
5. Results output

Usage:
    python run_simulation.py [--lattice N] [--beta B] [--sweeps S] [--configs C]

Author: Yang-Mills Mass Gap Project
Date: January 2026
"""

import argparse
import os
import sys
import time
import json
from datetime import datetime
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass, asdict
import numpy as np

# Import our modules
import lattice_gauge as lg
import wilson_loops as wl
import mass_gap as mg


# =============================================================================
# Simulation Configuration
# =============================================================================

@dataclass
class SimulationConfig:
    """Configuration for the complete simulation."""
    # Lattice parameters
    Nx: int = 8
    Ny: int = 8
    Nz: int = 8
    Nt: int = 16
    beta: float = 6.0

    # Monte Carlo parameters
    thermalization_sweeps: int = 200
    measurement_sweeps: int = 100
    sweeps_between_measurements: int = 10
    update_method: str = "metropolis"  # Use Metropolis (heatbath has issues)

    # Measurement parameters
    wilson_R_max: int = 6
    wilson_T_max: int = 8
    wilson_sample_fraction: float = 0.5

    # APE smearing
    n_ape_smear: int = 10
    ape_alpha: float = 0.5

    # Output
    output_dir: str = "results"
    save_configs: bool = False

    @property
    def lattice_shape(self) -> Tuple[int, int, int, int]:
        return (self.Nx, self.Ny, self.Nz, self.Nt)

    @property
    def volume(self) -> int:
        return self.Nx * self.Ny * self.Nz * self.Nt


@dataclass
class MeasurementResult:
    """Single measurement result."""
    config_number: int
    plaquette: float
    wilson_loops: Dict[str, float]  # "R,T" -> value
    polyakov_abs: float
    polyakov_real: float


@dataclass
class SimulationResults:
    """Complete simulation results."""
    config: SimulationConfig
    measurements: List[MeasurementResult]
    potential_fit: Optional[Dict]
    mass_gap: Optional[Dict]
    timing: Dict[str, float]
    metadata: Dict


# =============================================================================
# Simulation Runner
# =============================================================================

class Simulation:
    """
    Main simulation class for Yang-Mills mass gap analysis.

    Handles:
    - Configuration generation
    - Measurement collection
    - Statistical averaging
    - Mass gap extraction
    """

    def __init__(self, config: SimulationConfig):
        """Initialize simulation with given configuration."""
        self.config = config
        self.gauge: Optional[lg.GaugeField] = None
        self.measurements: List[MeasurementResult] = []
        self.timing: Dict[str, float] = {}

    def setup(self) -> None:
        """Set up the lattice and gauge field."""
        print("\n" + "=" * 70)
        print("YANG-MILLS MASS GAP SIMULATION")
        print("=" * 70)
        print(f"\nLattice: {self.config.Nx} x {self.config.Ny} x {self.config.Nz} x {self.config.Nt}")
        print(f"Beta:    {self.config.beta}")
        print(f"Volume:  {self.config.volume} sites")

        lattice_config = lg.LatticeConfig(
            Nx=self.config.Nx,
            Ny=self.config.Ny,
            Nz=self.config.Nz,
            Nt=self.config.Nt,
            beta=self.config.beta
        )
        self.gauge = lg.GaugeField(lattice_config)

    def thermalize(self) -> None:
        """Thermalize the gauge configuration."""
        print(f"\n{'='*70}")
        print("THERMALIZATION")
        print(f"{'='*70}")

        # Start from hot configuration
        print("\nInitializing hot start...")
        self.gauge.hot_start()
        initial_plaq = lg.compute_average_plaquette(self.gauge)
        print(f"Initial plaquette: {initial_plaq:.6f}")

        # Thermalize
        print(f"\nRunning {self.config.thermalization_sweeps} thermalization sweeps...")
        method = lg.UpdateMethod.HEATBATH if self.config.update_method == "heatbath" else lg.UpdateMethod.METROPOLIS

        start_time = time.time()

        plaq_history = []
        for i in range(self.config.thermalization_sweeps):
            lg.sweep(self.gauge, method)

            if (i + 1) % 20 == 0 or i == self.config.thermalization_sweeps - 1:
                plaq = lg.compute_average_plaquette(self.gauge)
                plaq_history.append(plaq)
                print(f"  Sweep {i+1:4d}: <P> = {plaq:.6f}")

        self.timing['thermalization'] = time.time() - start_time

        final_plaq = lg.compute_average_plaquette(self.gauge)
        print(f"\nFinal plaquette:   {final_plaq:.6f}")
        print(f"Thermalization time: {self.timing['thermalization']:.1f} seconds")

    def generate_and_measure(self) -> None:
        """Generate configurations and perform measurements."""
        print(f"\n{'='*70}")
        print("CONFIGURATION GENERATION AND MEASUREMENT")
        print(f"{'='*70}")

        print(f"\nGenerating {self.config.measurement_sweeps} configurations")
        print(f"Separation: {self.config.sweeps_between_measurements} sweeps")

        method = lg.UpdateMethod.HEATBATH if self.config.update_method == "heatbath" else lg.UpdateMethod.METROPOLIS

        start_time = time.time()

        for i in range(self.config.measurement_sweeps):
            # Update configuration
            for _ in range(self.config.sweeps_between_measurements):
                lg.sweep(self.gauge, method)

            # Perform measurements
            measurement = self._measure_configuration(i)
            self.measurements.append(measurement)

            if (i + 1) % 10 == 0:
                print(f"  Config {i+1:4d}: <P> = {measurement.plaquette:.6f}, "
                      f"|<L>| = {measurement.polyakov_abs:.6f}")

        self.timing['measurement'] = time.time() - start_time
        print(f"\nMeasurement time: {self.timing['measurement']:.1f} seconds")

    def _measure_configuration(self, config_number: int) -> MeasurementResult:
        """Perform all measurements on current configuration."""
        # Plaquette
        plaquette = lg.compute_average_plaquette(self.gauge)

        # Polyakov loop
        P_mean, P_abs = wl.compute_average_polyakov_loop(self.gauge)

        # Wilson loops
        wilson_loops = {}
        for R in range(1, self.config.wilson_R_max + 1):
            for T in range(1, self.config.wilson_T_max + 1):
                mean, _ = wl.compute_average_wilson_loop(
                    self.gauge, R, T,
                    sample_fraction=self.config.wilson_sample_fraction
                )
                wilson_loops[f"{R},{T}"] = mean

        return MeasurementResult(
            config_number=config_number,
            plaquette=plaquette,
            wilson_loops=wilson_loops,
            polyakov_abs=P_abs,
            polyakov_real=np.real(P_mean)
        )

    def analyze_results(self) -> Tuple[Optional[mg.PotentialFitResult], Optional[mg.MassGapResult]]:
        """Analyze collected measurements and extract mass gap."""
        print(f"\n{'='*70}")
        print("STATISTICAL ANALYSIS")
        print(f"{'='*70}")

        start_time = time.time()

        # Average Wilson loops
        print("\n1. Averaging Wilson loops...")
        avg_wilson = {}
        avg_wilson_err = {}

        for key in self.measurements[0].wilson_loops.keys():
            values = [m.wilson_loops[key] for m in self.measurements]
            avg_wilson[key] = np.mean(values)
            avg_wilson_err[key] = np.std(values) / np.sqrt(len(values))

        # Convert to (R, T) format
        wilson_data = {}
        for key, value in avg_wilson.items():
            R, T = map(int, key.split(','))
            wilson_data[(R, T)] = value

        # Print averaged Wilson loops
        print("\n   Averaged Wilson loops:")
        for (R, T) in sorted(wilson_data.keys()):
            key = f"{R},{T}"
            print(f"   W({R},{T}) = {avg_wilson[key]:.6e} +/- {avg_wilson_err[key]:.6e}")

        # Average observables
        print("\n2. Observable statistics:")
        plaq_values = [m.plaquette for m in self.measurements]
        poly_values = [m.polyakov_abs for m in self.measurements]

        print(f"   <P>  = {np.mean(plaq_values):.6f} +/- {np.std(plaq_values)/np.sqrt(len(plaq_values)):.6f}")
        print(f"   |<L>| = {np.mean(poly_values):.6f} +/- {np.std(poly_values)/np.sqrt(len(poly_values)):.6f}")

        # Extract potential
        print("\n3. Extracting static potential...")
        T_values = list(range(2, self.config.wilson_T_max + 1))
        R, V, V_err = mg.extract_potential_multiT(wilson_data, T_values, method="log_ratio")

        valid = ~np.isnan(V)
        print("   R    V(R)      Error")
        for i in range(len(R)):
            if valid[i]:
                print(f"   {R[i]:.0f}    {V[i]:.4f}    {V_err[i]:.4f}")

        # Fit potential
        print("\n4. Fitting Cornell potential...")
        fitter = mg.StaticPotentialFitter(R, V, V_err)

        try:
            fit_result = fitter.fit_cornell(R_min=1.0)
            print(f"   sigma = {fit_result.sigma:.6f} +/- {fit_result.sigma_err:.6f}")
            print(f"   mu    = {fit_result.mu:.6f} +/- {fit_result.mu_err:.6f}")
            print(f"   c     = {fit_result.c:.6f} +/- {fit_result.c_err:.6f}")
            print(f"   chi^2/dof = {fit_result.chi_squared_per_dof:.3f}")
        except Exception as e:
            print(f"   Cornell fit failed: {e}")
            print("   Trying linear fit...")
            try:
                fit_result = fitter.fit_linear(R_min=1.0)
                print(f"   sigma = {fit_result.sigma:.6f} +/- {fit_result.sigma_err:.6f}")
                print(f"   mu    = {fit_result.mu:.6f} +/- {fit_result.mu_err:.6f}")
            except Exception as e2:
                print(f"   Linear fit also failed: {e2}")
                fit_result = None

        # Extract mass gap
        print("\n5. Extracting mass gap...")
        if fit_result is not None and fit_result.sigma > 0:
            mass_gap_result = mg.extract_mass_gap_from_string_tension(
                fit_result.sigma,
                fit_result.sigma_err,
                lattice_spacing=mg.compute_lattice_spacing(self.config.beta)
            )
            print(f"   Mass gap (lattice): {mass_gap_result.mass_gap:.6f} +/- {mass_gap_result.mass_gap_err:.6f}")
            if mass_gap_result.mass_gap_MeV is not None:
                print(f"   Mass gap (physical): {mass_gap_result.mass_gap_MeV:.0f} MeV")
        else:
            mass_gap_result = None
            print("   Could not extract mass gap (fit failed or non-confining)")

        self.timing['analysis'] = time.time() - start_time
        print(f"\nAnalysis time: {self.timing['analysis']:.1f} seconds")

        return fit_result, mass_gap_result

    def run(self) -> SimulationResults:
        """Run complete simulation."""
        total_start = time.time()

        # Setup
        self.setup()

        # Thermalize
        self.thermalize()

        # Generate and measure
        self.generate_and_measure()

        # Analyze
        fit_result, mass_gap_result = self.analyze_results()

        self.timing['total'] = time.time() - total_start

        # Compile results
        results = SimulationResults(
            config=self.config,
            measurements=self.measurements,
            potential_fit=asdict(fit_result) if fit_result else None,
            mass_gap=asdict(mass_gap_result) if mass_gap_result else None,
            timing=self.timing,
            metadata={
                'date': datetime.now().isoformat(),
                'n_configs': len(self.measurements),
                'gauge_group': 'SU(3)',
            }
        )

        return results


# =============================================================================
# Results Output
# =============================================================================

def print_summary(results: SimulationResults) -> None:
    """Print final summary of results."""
    print("\n" + "=" * 70)
    print("SIMULATION SUMMARY")
    print("=" * 70)

    print(f"\nLattice:     {results.config.Nx}^3 x {results.config.Nt}")
    print(f"Beta:        {results.config.beta}")
    print(f"Configs:     {len(results.measurements)}")

    if results.potential_fit:
        print(f"\nString Tension:")
        print(f"  sigma = {results.potential_fit['sigma']:.6f} +/- {results.potential_fit['sigma_err']:.6f}")

    if results.mass_gap:
        print(f"\nMass Gap:")
        print(f"  Delta (lattice) = {results.mass_gap['mass_gap']:.6f} +/- {results.mass_gap['mass_gap_err']:.6f}")
        if results.mass_gap['mass_gap_MeV']:
            print(f"  Delta (physical) ~ {results.mass_gap['mass_gap_MeV']:.0f} MeV")

    print(f"\nTiming:")
    print(f"  Thermalization: {results.timing['thermalization']:.1f} s")
    print(f"  Measurement:    {results.timing['measurement']:.1f} s")
    print(f"  Analysis:       {results.timing['analysis']:.1f} s")
    print(f"  Total:          {results.timing['total']:.1f} s")

    print("=" * 70)


def save_results(results: SimulationResults, output_dir: str) -> None:
    """Save results to files."""
    os.makedirs(output_dir, exist_ok=True)

    # Save main results as JSON
    results_dict = {
        'config': asdict(results.config),
        'potential_fit': results.potential_fit,
        'mass_gap': results.mass_gap,
        'timing': results.timing,
        'metadata': results.metadata,
        'n_measurements': len(results.measurements),
    }

    json_path = os.path.join(output_dir, 'results.json')
    with open(json_path, 'w') as f:
        json.dump(results_dict, f, indent=2)
    print(f"Results saved to {json_path}")

    # Save measurements
    meas_path = os.path.join(output_dir, 'measurements.json')
    measurements_data = [asdict(m) for m in results.measurements]
    with open(meas_path, 'w') as f:
        json.dump(measurements_data, f, indent=2)
    print(f"Measurements saved to {meas_path}")


# =============================================================================
# Command Line Interface
# =============================================================================

def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Yang-Mills Mass Gap Simulation',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Lattice parameters
    parser.add_argument('--lattice', type=int, default=8,
                        help='Spatial lattice size (Nx=Ny=Nz)')
    parser.add_argument('--time', type=int, default=16,
                        help='Temporal lattice size (Nt)')
    parser.add_argument('--beta', type=float, default=6.0,
                        help='Inverse coupling beta=6/g^2')

    # Monte Carlo parameters
    parser.add_argument('--therm', type=int, default=200,
                        help='Thermalization sweeps')
    parser.add_argument('--configs', type=int, default=100,
                        help='Number of configurations')
    parser.add_argument('--sep', type=int, default=10,
                        help='Sweeps between measurements')

    # Measurement parameters
    parser.add_argument('--Rmax', type=int, default=6,
                        help='Maximum R for Wilson loops')
    parser.add_argument('--Tmax', type=int, default=8,
                        help='Maximum T for Wilson loops')
    parser.add_argument('--sample', type=float, default=0.5,
                        help='Sampling fraction for Wilson loops')

    # Output
    parser.add_argument('--output', type=str, default='results',
                        help='Output directory')

    # Quick test mode
    parser.add_argument('--quick', action='store_true',
                        help='Quick test mode (small lattice, few configs)')

    return parser.parse_args()


def main():
    """Main entry point."""
    args = parse_args()

    # Create configuration
    if args.quick:
        # Quick test mode
        config = SimulationConfig(
            Nx=4, Ny=4, Nz=4, Nt=8,
            beta=6.0,
            thermalization_sweeps=50,
            measurement_sweeps=20,
            sweeps_between_measurements=5,
            wilson_R_max=3,
            wilson_T_max=4,
            wilson_sample_fraction=1.0,
            output_dir=args.output
        )
    else:
        config = SimulationConfig(
            Nx=args.lattice,
            Ny=args.lattice,
            Nz=args.lattice,
            Nt=args.time,
            beta=args.beta,
            thermalization_sweeps=args.therm,
            measurement_sweeps=args.configs,
            sweeps_between_measurements=args.sep,
            wilson_R_max=args.Rmax,
            wilson_T_max=args.Tmax,
            wilson_sample_fraction=args.sample,
            output_dir=args.output
        )

    # Run simulation
    sim = Simulation(config)
    results = sim.run()

    # Output results
    print_summary(results)
    save_results(results, config.output_dir)

    return 0


if __name__ == "__main__":
    sys.exit(main())
