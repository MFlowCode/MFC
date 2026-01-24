#!/usr/bin/env python3
"""
Visualization script for 2D Burning Droplet simulation results.

This script creates visualizations of:
- Temperature field
- Species mass fractions (fuel, oxidizer, products)
- Velocity field

Usage:
    python viz.py --case_dir <path_to_case_directory>
"""

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable

def read_binary_data(filepath, nx, ny, precision=2):
    """Read MFC binary data file."""
    dtype = np.float64 if precision == 2 else np.float32
    with open(filepath, 'rb') as f:
        data = np.fromfile(f, dtype=dtype)
    return data.reshape((ny + 1, nx + 1))

def find_data_files(case_dir):
    """Find simulation output directories."""
    dirs = []
    for name in os.listdir(case_dir):
        path = os.path.join(case_dir, name)
        if os.path.isdir(path) and name.startswith('D'):
            dirs.append(path)
    return sorted(dirs)

def plot_burning_droplet(data_dir, output_dir=None):
    """Create visualization of burning droplet results."""
    
    # Read grid dimensions from case (you may need to parse case.py)
    # For now, assuming standard dimensions
    nx = 200
    ny = 200
    
    # Read temperature if available
    temp_file = os.path.join(data_dir, 'T.dat')
    if os.path.exists(temp_file):
        T = read_binary_data(temp_file, nx, ny)
    else:
        print(f"Temperature file not found: {temp_file}")
        T = None
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot temperature
    if T is not None:
        ax = axes[0, 0]
        im = ax.imshow(T, origin='lower', cmap='hot', aspect='equal')
        ax.set_title('Temperature (K)')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax)
    
    # Read and plot species mass fractions
    species_names = ['H2', 'O2', 'H2O', 'OH']
    species_files = ['Y_1.dat', 'Y_4.dat', 'Y_6.dat', 'Y_5.dat']
    
    for idx, (name, filename) in enumerate(zip(species_names[1:], species_files[1:])):
        row, col = divmod(idx + 1, 2)
        ax = axes[row, col]
        
        filepath = os.path.join(data_dir, filename)
        if os.path.exists(filepath):
            Y = read_binary_data(filepath, nx, ny)
            im = ax.imshow(Y, origin='lower', cmap='viridis', aspect='equal')
            ax.set_title(f'Mass Fraction: {name}')
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(im, cax=cax)
        else:
            ax.text(0.5, 0.5, f'{filename}\nnot found', 
                   ha='center', va='center', transform=ax.transAxes)
    
    plt.tight_layout()
    
    if output_dir:
        output_file = os.path.join(output_dir, 'burning_droplet.png')
        plt.savefig(output_file, dpi=200, bbox_inches='tight')
        print(f"Saved: {output_file}")
    else:
        plt.show()

def plot_flame_structure(data_dir, output_dir=None):
    """Plot radial profiles through the flame."""
    
    nx = ny = 200
    
    # Read data
    temp_file = os.path.join(data_dir, 'T.dat')
    h2_file = os.path.join(data_dir, 'Y_1.dat')
    o2_file = os.path.join(data_dir, 'Y_4.dat')
    h2o_file = os.path.join(data_dir, 'Y_6.dat')
    
    T = read_binary_data(temp_file, nx, ny) if os.path.exists(temp_file) else None
    Y_H2 = read_binary_data(h2_file, nx, ny) if os.path.exists(h2_file) else None
    Y_O2 = read_binary_data(o2_file, nx, ny) if os.path.exists(o2_file) else None
    Y_H2O = read_binary_data(h2o_file, nx, ny) if os.path.exists(h2o_file) else None
    
    # Extract horizontal profile through center
    center_y = ny // 2
    x = np.linspace(-0.005, 0.005, nx + 1)  # 10 mm domain
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
    
    # Temperature profile
    if T is not None:
        ax1.plot(x * 1000, T[center_y, :], 'r-', linewidth=2, label='Temperature')
        ax1.set_xlabel('x (mm)')
        ax1.set_ylabel('Temperature (K)')
        ax1.set_title('Radial Temperature Profile')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
    
    # Species profiles
    if Y_H2 is not None:
        ax2.plot(x * 1000, Y_H2[center_y, :], 'b-', linewidth=2, label='H$_2$ (Fuel)')
    if Y_O2 is not None:
        ax2.plot(x * 1000, Y_O2[center_y, :], 'g-', linewidth=2, label='O$_2$ (Oxidizer)')
    if Y_H2O is not None:
        ax2.plot(x * 1000, Y_H2O[center_y, :], 'c-', linewidth=2, label='H$_2$O (Product)')
    
    ax2.set_xlabel('x (mm)')
    ax2.set_ylabel('Mass Fraction')
    ax2.set_title('Radial Species Profiles')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if output_dir:
        output_file = os.path.join(output_dir, 'flame_structure.png')
        plt.savefig(output_file, dpi=200, bbox_inches='tight')
        print(f"Saved: {output_file}")
    else:
        plt.show()

def main():
    parser = argparse.ArgumentParser(description='Visualize burning droplet results')
    parser.add_argument('--case_dir', type=str, default='.', 
                       help='Path to case directory')
    parser.add_argument('--output_dir', type=str, default=None,
                       help='Output directory for plots')
    parser.add_argument('--timestep', type=int, default=-1,
                       help='Timestep to visualize (-1 for last)')
    args = parser.parse_args()
    
    # Find data directories
    data_dirs = find_data_files(args.case_dir)
    
    if not data_dirs:
        print(f"No data directories found in {args.case_dir}")
        print("Run the simulation first with:")
        print("  ./mfc.sh run examples/2D_burning_droplet/case.py -t pre_process simulation")
        return
    
    # Select timestep
    if args.timestep == -1:
        data_dir = data_dirs[-1]
    else:
        matching = [d for d in data_dirs if str(args.timestep) in d]
        if matching:
            data_dir = matching[0]
        else:
            print(f"Timestep {args.timestep} not found")
            return
    
    print(f"Visualizing: {data_dir}")
    
    # Create plots
    plot_burning_droplet(data_dir, args.output_dir)
    plot_flame_structure(data_dir, args.output_dir)

if __name__ == "__main__":
    main()
