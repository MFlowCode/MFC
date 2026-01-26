#!/usr/bin/env python3
"""
Visualization script for gas-only chemistry test results.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# =============================================================================
# CONFIGURATION  
# =============================================================================
CASE_DIR = Path(__file__).parent
RESTART_DIR = CASE_DIR / "restart_data"
OUTPUT_DIR = CASE_DIR / "figures"

# Domain parameters
Nx = 200  # Actual output size (m+1)
Lx = 1.0e-2  # 1 cm

# Number of variables: 1 fluid (rho, rho*u, E, alpha) + 10 species = 14
# For model_eqns=2, num_fluids=1: alpha_rho, rho*u, E, alpha, + species
NUM_VARS = 14

# Species names (h2o2.yaml)
SPECIES = ['H2', 'H', 'O', 'O2', 'OH', 'H2O', 'HO2', 'H2O2', 'AR', 'N2']

OUTPUT_DIR.mkdir(exist_ok=True)

def read_restart_data(restart_dir, timestep, num_cells, num_vars):
    """Read restart data for a given timestep"""
    filename = restart_dir / f"lustre_{timestep}.dat"
    
    if not filename.exists():
        print(f"File not found: {filename}")
        return None
    
    data = np.fromfile(filename, dtype=np.float64)
    
    expected_size = num_vars * num_cells
    if len(data) != expected_size:
        # Try to infer num_vars
        num_vars_inferred = len(data) // num_cells
        if len(data) % num_cells == 0:
            print(f"Inferred {num_vars_inferred} variables (expected {num_vars})")
            num_vars = num_vars_inferred
        else:
            print(f"Warning: Data size {len(data)} doesn't match expected {expected_size}")
            return None
    
    data = data.reshape((num_vars, num_cells))
    return data


def read_grid(restart_dir, num_cells):
    """Read grid coordinates"""
    grid_file = restart_dir / "lustre_x_cb.dat"
    
    if grid_file.exists():
        x = np.fromfile(grid_file, dtype=np.float64)
        if len(x) == num_cells + 1:
            x_cc = 0.5 * (x[:-1] + x[1:])
        else:
            x_cc = x[:num_cells]
        return x_cc
    else:
        return np.linspace(0, Lx, num_cells)


def plot_species_profiles(x, data, timestep, output_dir):
    """Plot species mass fraction profiles"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    num_vars = data.shape[0]
    
    # For model_eqns=2, num_fluids=1 with chemistry:
    # Variables: alpha_rho, rho*u, E, alpha, Y_1*rho, Y_2*rho, ..., Y_10*rho
    # So species start at index 4
    
    species_start = 4
    
    # Get total density (alpha_rho for single fluid = rho)
    rho = data[0, :]
    
    # Key species to plot
    key_species = {
        'H2': 0,    # Fuel
        'O2': 3,    # Oxidizer
        'H2O': 5,   # Product
        'OH': 4,    # Radical (flame marker)
    }
    
    colors = {'H2': 'blue', 'O2': 'green', 'H2O': 'red', 'OH': 'orange'}
    
    # Plot 1: Reactants
    ax = axes[0, 0]
    for species, idx in [('H2', 0), ('O2', 3)]:
        if species_start + idx < num_vars:
            rhoY = data[species_start + idx, :]
            Y = rhoY / (rho + 1e-20)
            ax.plot(x * 100, Y, color=colors[species], label=species, linewidth=2)
    ax.set_xlabel('x (cm)')
    ax.set_ylabel('Mass Fraction Y')
    ax.set_title(f'Reactants at t = {timestep} steps')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.05, 1.05)
    
    # Plot 2: Products
    ax = axes[0, 1]
    for species, idx in [('H2O', 5), ('OH', 4)]:
        if species_start + idx < num_vars:
            rhoY = data[species_start + idx, :]
            Y = rhoY / (rho + 1e-20)
            ax.plot(x * 100, Y, color=colors[species], label=species, linewidth=2)
    ax.set_xlabel('x (cm)')
    ax.set_ylabel('Mass Fraction Y')
    ax.set_title(f'Products at t = {timestep} steps')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 3: Density
    ax = axes[1, 0]
    ax.plot(x * 100, rho, 'k-', linewidth=2)
    ax.set_xlabel('x (cm)')
    ax.set_ylabel('Density (kg/m³)')
    ax.set_title(f'Density at t = {timestep} steps')
    ax.grid(True, alpha=0.3)
    
    # Plot 4: Energy/Temperature proxy
    ax = axes[1, 1]
    E = data[2, :]  # Total energy
    # Temperature is roughly proportional to E/rho for ideal gas
    T_proxy = E / (rho + 1e-20) / 1000  # Normalized
    ax.plot(x * 100, T_proxy, 'r-', linewidth=2)
    ax.set_xlabel('x (cm)')
    ax.set_ylabel('E/ρ (kJ/kg)')
    ax.set_title(f'Specific Energy at t = {timestep} steps')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / f'chemistry_t{timestep:04d}.png', dpi=150)
    plt.close()
    
    return True


def plot_time_evolution(restart_dir, x, num_cells, output_dir):
    """Plot time evolution of species"""
    timesteps = [0, 20, 40, 60, 80, 100]
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    colors = plt.cm.viridis(np.linspace(0, 1, len(timesteps)))
    
    for i, ts in enumerate(timesteps):
        data = read_restart_data(restart_dir, ts, num_cells, NUM_VARS)
        if data is None:
            continue
        
        rho = data[0, :]
        species_start = 4
        num_vars = data.shape[0]
        
        # H2
        if species_start + 0 < num_vars:
            Y_H2 = data[species_start + 0, :] / (rho + 1e-20)
            axes[0, 0].plot(x * 100, Y_H2, color=colors[i], label=f't={ts}', linewidth=2)
        
        # O2
        if species_start + 3 < num_vars:
            Y_O2 = data[species_start + 3, :] / (rho + 1e-20)
            axes[0, 1].plot(x * 100, Y_O2, color=colors[i], label=f't={ts}', linewidth=2)
        
        # H2O
        if species_start + 5 < num_vars:
            Y_H2O = data[species_start + 5, :] / (rho + 1e-20)
            axes[1, 0].plot(x * 100, Y_H2O, color=colors[i], label=f't={ts}', linewidth=2)
        
        # Temperature proxy
        E = data[2, :]
        T_proxy = E / (rho + 1e-20) / 1000
        axes[1, 1].plot(x * 100, T_proxy, color=colors[i], label=f't={ts}', linewidth=2)
    
    axes[0, 0].set_xlabel('x (cm)')
    axes[0, 0].set_ylabel('Y_H2')
    axes[0, 0].set_title('H2 (Fuel)')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    axes[0, 1].set_xlabel('x (cm)')
    axes[0, 1].set_ylabel('Y_O2')
    axes[0, 1].set_title('O2 (Oxidizer)')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    axes[1, 0].set_xlabel('x (cm)')
    axes[1, 0].set_ylabel('Y_H2O')
    axes[1, 0].set_title('H2O (Product)')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    axes[1, 1].set_xlabel('x (cm)')
    axes[1, 1].set_ylabel('E/ρ (kJ/kg)')
    axes[1, 1].set_title('Specific Energy')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'chemistry_evolution.png', dpi=150)
    plt.close()
    
    print(f"Saved: {output_dir / 'chemistry_evolution.png'}")


def main():
    print("=" * 60)
    print("Chemistry Test Visualization")
    print("=" * 60)
    
    if not RESTART_DIR.exists():
        print(f"Error: Restart directory not found: {RESTART_DIR}")
        return
    
    x = read_grid(RESTART_DIR, Nx)
    print(f"Grid: {len(x)} cells, x = [{x[0]*100:.4f}, {x[-1]*100:.4f}] cm")
    
    # Check data structure
    data0 = read_restart_data(RESTART_DIR, 0, Nx, NUM_VARS)
    if data0 is not None:
        print(f"Data shape: {data0.shape}")
    
    print("\nGenerating plots...")
    
    for ts in [0, 50, 100]:
        data = read_restart_data(RESTART_DIR, ts, Nx, NUM_VARS)
        if data is not None:
            plot_species_profiles(x, data, ts, OUTPUT_DIR)
            print(f"  Saved plots for timestep {ts}")
    
    plot_time_evolution(RESTART_DIR, x, Nx, OUTPUT_DIR)
    
    print(f"\nAll figures saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
