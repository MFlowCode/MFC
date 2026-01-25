#!/usr/bin/env python3
"""
Visualization script for Phase 1 validation results.
Reads MFC restart data and creates plots of volume fractions, density, etc.
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

# Domain parameters (must match case file)
# Note: MFC output includes m+1 = 200 cells for m=199
Nx = 200  # Actual output size
Lx = 1.0e-3  # 1 mm
NUM_VARS = 11  # Number of variables in restart file

# Create output directory
OUTPUT_DIR.mkdir(exist_ok=True)

# =============================================================================
# DATA READING FUNCTIONS
# =============================================================================

def read_indices(case_dir):
    """Read variable indices from indices.dat"""
    indices_file = case_dir / "indices.dat"
    indices = {}
    
    if indices_file.exists():
        with open(indices_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        idx = int(parts[0])
                        name = parts[1]
                        indices[name] = idx
                    except ValueError:
                        continue
    
    return indices


def read_restart_data(restart_dir, timestep, num_cells):
    """Read restart data for a given timestep"""
    filename = restart_dir / f"lustre_{timestep}.dat"
    
    if not filename.exists():
        print(f"File not found: {filename}")
        return None
    
    # Read binary data
    data = np.fromfile(filename, dtype=np.float64)
    
    # Reshape based on number of variables
    num_vars = len(data) // num_cells
    if len(data) % num_cells != 0:
        print(f"Warning: Data size {len(data)} not divisible by {num_cells}")
        return None
    
    # Reshape to (num_vars, num_cells)
    data = data.reshape((num_vars, num_cells))
    
    return data


def read_grid(restart_dir, num_cells):
    """Read grid coordinates"""
    grid_file = restart_dir / "lustre_x_cb.dat"
    
    if grid_file.exists():
        x = np.fromfile(grid_file, dtype=np.float64)
        # Cell-center coordinates
        if len(x) == num_cells + 1:
            x_cc = 0.5 * (x[:-1] + x[1:])
        else:
            x_cc = x[:num_cells]
        return x_cc
    else:
        # Generate uniform grid
        return np.linspace(0, Lx, num_cells)


# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

def plot_volume_fractions(x, data, timestep, output_dir):
    """Plot volume fractions for all fluids"""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # For 3-fluid model with 6-eqn, the layout is:
    # Assuming indices: alpha_rho_1, alpha_rho_2, alpha_rho_3, rho*u, E, alpha_1, alpha_2, alpha_3, ...
    # This depends on model_eqns. Let's try to infer from data shape.
    
    num_vars = data.shape[0]
    print(f"Timestep {timestep}: {num_vars} variables, {data.shape[1]} cells")
    
    # For 3-fluid 6-eqn model (model_eqns=3):
    # Variables: alpha_rho_1, alpha_rho_2, alpha_rho_3, rho*u, E, alpha_1, alpha_2, alpha_3, int_e_1, int_e_2, int_e_3
    # So alpha starts at index 5 (0-indexed)
    
    if num_vars >= 8:
        # Volume fractions (assuming they start at index 5 for 3-fluid model)
        alpha_start = 5
        alpha_1 = data[alpha_start, :]
        alpha_2 = data[alpha_start + 1, :]
        alpha_3 = data[alpha_start + 2, :]
        
        ax.plot(x * 1000, alpha_1, 'b-', label=r'$\alpha_1$ (Liquid)', linewidth=2)
        ax.plot(x * 1000, alpha_2, 'r--', label=r'$\alpha_2$ (Vapor)', linewidth=2)
        ax.plot(x * 1000, alpha_3, 'g:', label=r'$\alpha_3$ (Air)', linewidth=2)
        ax.plot(x * 1000, alpha_1 + alpha_2 + alpha_3, 'k-', label=r'Sum', linewidth=1, alpha=0.5)
    else:
        # Just plot all variables
        for i in range(min(num_vars, 5)):
            ax.plot(x * 1000, data[i, :], label=f'Var {i}')
    
    ax.set_xlabel('x (mm)', fontsize=12)
    ax.set_ylabel('Volume Fraction', fontsize=12)
    ax.set_title(f'Volume Fractions at t = {timestep} steps', fontsize=14)
    ax.legend(loc='best')
    ax.set_ylim(-0.1, 1.1)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / f'volume_fractions_t{timestep:04d}.png', dpi=150)
    plt.close()


def plot_density(x, data, timestep, output_dir):
    """Plot partial densities"""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    num_vars = data.shape[0]
    
    # Partial densities are the first num_fluids variables
    if num_vars >= 3:
        rho_1 = data[0, :]
        rho_2 = data[1, :]
        rho_3 = data[2, :]
        rho_total = rho_1 + rho_2 + rho_3
        
        ax.semilogy(x * 1000, rho_1 + 1e-10, 'b-', label=r'$\alpha_1 \rho_1$ (Liquid)', linewidth=2)
        ax.semilogy(x * 1000, rho_2 + 1e-10, 'r--', label=r'$\alpha_2 \rho_2$ (Vapor)', linewidth=2)
        ax.semilogy(x * 1000, rho_3 + 1e-10, 'g:', label=r'$\alpha_3 \rho_3$ (Air)', linewidth=2)
        ax.semilogy(x * 1000, rho_total, 'k-', label=r'$\rho_{total}$', linewidth=1)
    
    ax.set_xlabel('x (mm)', fontsize=12)
    ax.set_ylabel('Partial Density (kg/m³)', fontsize=12)
    ax.set_title(f'Partial Densities at t = {timestep} steps', fontsize=14)
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / f'density_t{timestep:04d}.png', dpi=150)
    plt.close()


def plot_time_evolution(restart_dir, x, num_cells, output_dir):
    """Plot time evolution of key quantities"""
    timesteps = [0, 25, 50, 75, 100]
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(timesteps)))
    
    for i, ts in enumerate(timesteps):
        data = read_restart_data(restart_dir, ts, num_cells)
        if data is None:
            continue
        
        # Volume fraction of liquid (alpha_1)
        if data.shape[0] >= 6:
            alpha_1 = data[5, :]
            axes[0, 0].plot(x * 1000, alpha_1, color=colors[i], 
                           label=f't = {ts}', linewidth=2)
        
        # Partial density of liquid
        if data.shape[0] >= 1:
            rho_1 = data[0, :]
            axes[0, 1].plot(x * 1000, rho_1, color=colors[i], 
                           label=f't = {ts}', linewidth=2)
        
        # Partial density of vapor  
        if data.shape[0] >= 2:
            rho_2 = data[1, :]
            axes[1, 0].plot(x * 1000, rho_2, color=colors[i],
                           label=f't = {ts}', linewidth=2)
        
        # Total energy
        if data.shape[0] >= 5:
            E = data[4, :]
            axes[1, 1].plot(x * 1000, E / 1e6, color=colors[i],
                           label=f't = {ts}', linewidth=2)
    
    axes[0, 0].set_xlabel('x (mm)')
    axes[0, 0].set_ylabel(r'$\alpha_1$ (Liquid)')
    axes[0, 0].set_title('Liquid Volume Fraction')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    axes[0, 1].set_xlabel('x (mm)')
    axes[0, 1].set_ylabel(r'$\alpha_1 \rho_1$ (kg/m³)')
    axes[0, 1].set_title('Liquid Partial Density')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    axes[1, 0].set_xlabel('x (mm)')
    axes[1, 0].set_ylabel(r'$\alpha_2 \rho_2$ (kg/m³)')
    axes[1, 0].set_title('Vapor Partial Density')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    axes[1, 1].set_xlabel('x (mm)')
    axes[1, 1].set_ylabel('E (MJ/m³)')
    axes[1, 1].set_title('Total Energy')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'time_evolution.png', dpi=150)
    plt.close()
    
    print(f"Saved: {output_dir / 'time_evolution.png'}")


def plot_interface_position(restart_dir, x, num_cells, output_dir):
    """Track interface position over time"""
    timesteps = list(range(0, 101, 5))
    interface_positions = []
    
    for ts in timesteps:
        data = read_restart_data(restart_dir, ts, num_cells)
        if data is None or data.shape[0] < 6:
            continue
        
        alpha_1 = data[5, :]  # Liquid volume fraction
        
        # Find interface (where alpha_1 = 0.5)
        for i in range(len(alpha_1) - 1):
            if alpha_1[i] > 0.5 and alpha_1[i+1] < 0.5:
                # Linear interpolation
                x_interface = x[i] + (0.5 - alpha_1[i]) / (alpha_1[i+1] - alpha_1[i]) * (x[i+1] - x[i])
                interface_positions.append((ts, x_interface))
                break
    
    if interface_positions:
        ts_arr = np.array([p[0] for p in interface_positions])
        x_arr = np.array([p[1] for p in interface_positions])
        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(ts_arr, x_arr * 1000, 'bo-', linewidth=2, markersize=8)
        ax.set_xlabel('Time Step', fontsize=12)
        ax.set_ylabel('Interface Position (mm)', fontsize=12)
        ax.set_title('Liquid-Gas Interface Position vs Time', fontsize=14)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_dir / 'interface_position.png', dpi=150)
        plt.close()
        
        print(f"Saved: {output_dir / 'interface_position.png'}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 60)
    print("Phase 1 Validation Visualization")
    print("=" * 60)
    
    # Check if data exists
    if not RESTART_DIR.exists():
        print(f"Error: Restart directory not found: {RESTART_DIR}")
        return
    
    # Read grid
    x = read_grid(RESTART_DIR, Nx)
    print(f"Grid: {len(x)} cells, x = [{x[0]*1000:.4f}, {x[-1]*1000:.4f}] mm")
    
    # Read indices
    indices = read_indices(CASE_DIR)
    if indices:
        print(f"Variable indices: {indices}")
    
    # Plot initial and final states
    print("\nGenerating plots...")
    
    for ts in [0, 50, 100]:
        data = read_restart_data(RESTART_DIR, ts, Nx)
        if data is not None:
            plot_volume_fractions(x, data, ts, OUTPUT_DIR)
            plot_density(x, data, ts, OUTPUT_DIR)
            print(f"  Saved plots for timestep {ts}")
    
    # Time evolution
    plot_time_evolution(RESTART_DIR, x, Nx, OUTPUT_DIR)
    
    # Interface tracking
    plot_interface_position(RESTART_DIR, x, Nx, OUTPUT_DIR)
    
    print(f"\nAll figures saved to: {OUTPUT_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
