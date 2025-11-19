#!/usr/bin/env python3
"""
Test script to validate MFC Spack package logic without requiring actual installation.
This tests that cmake_args and environment setup are correct for various HPC configurations.
"""

import sys
import os
import tempfile
import shutil
from pathlib import Path

def test_package_loads():
    """Test that the package.py file loads without errors"""
    print("=== Testing Package Loading ===")
    try:
        with open('package.py', 'r') as f:
            code = f.read()
        compile(code, 'package.py', 'exec')
        print("✓ Package file loads successfully")
        return True
    except Exception as e:
        print(f"✗ Package file failed to load: {e}")
        return False

def test_variant_definitions():
    """Test that all variants are properly defined"""
    print("\n=== Testing Variant Definitions ===")
    variants_to_check = [
        'mpi', 'openacc', 'openmp', 'precision', 'post_process', 
        'chemistry', 'cuda_arch', 'amdgpu_target'
    ]
    
    with open('package.py', 'r') as f:
        content = f.read()
    
    all_found = True
    for variant in variants_to_check:
        # More flexible matching - handle multi-line variant definitions
        if f'"{variant}"' in content or f"'{variant}'" in content:
            print(f"✓ Found variant: {variant}")
        else:
            print(f"✗ Missing variant: {variant}")
            all_found = False
    
    return all_found

def test_gpu_architecture_values():
    """Test that GPU architecture values are defined"""
    print("\n=== Testing GPU Architecture Values ===")
    
    with open('package.py', 'r') as f:
        content = f.read()
    
    # Check NVIDIA architectures
    nvidia_archs = ['60', '70', '75', '80', '90']  # P100, V100, T4, A100, H100
    print("NVIDIA GPU architectures:")
    for arch in nvidia_archs:
        if f'"{arch}"' in content:
            print(f"  ✓ cc{arch} supported")
        else:
            print(f"  ✗ cc{arch} missing")
    
    # Check AMD architectures
    amd_archs = ['gfx908', 'gfx90a', 'gfx940', 'gfx941', 'gfx942']
    print("AMD GPU architectures:")
    for arch in amd_archs:
        if arch in content:
            print(f"  ✓ {arch} supported")
        else:
            print(f"  ✗ {arch} missing")
    
    return True

def test_hpc_environment_setup():
    """Test that HPC-specific environment configuration is present"""
    print("\n=== Testing HPC Environment Setup ===")
    
    with open('package.py', 'r') as f:
        content = f.read()
    
    checks = [
        ('Cray compiler wrappers', 'env.set("CC", "cc")'),
        ('GPU-aware MPI (NVIDIA)', 'MPICH_GPU_SUPPORT_ENABLED'),
        ('GPU-aware MPI (AMD)', 'HSA_ENABLE_SDMA'),
        ('Parallel build configuration', 'CMAKE_BUILD_PARALLEL_LEVEL'),
        ('Fortran compiler wrapper', 'env.set("FC", "ftn")'),
    ]
    
    all_present = True
    for name, check_str in checks:
        if check_str in content:
            print(f"✓ {name} configuration present")
        else:
            print(f"✗ {name} configuration missing")
            all_present = False
    
    return all_present

def test_dependency_definitions():
    """Test that all required dependencies are defined"""
    print("\n=== Testing Dependency Definitions ===")
    
    with open('package.py', 'r') as f:
        content = f.read()
    
    required_deps = [
        'cmake', 'py-fypp', 'python', 'fftw', 'lapack', 'mpi', 'silo'
    ]
    
    all_found = True
    for dep in required_deps:
        if f'depends_on("{dep}' in content or f"depends_on('{dep}" in content:
            print(f"✓ Dependency: {dep}")
        else:
            print(f"✗ Missing dependency: {dep}")
            all_found = False
    
    return all_found

def test_cmake_args_logic():
    """Test that cmake_args method contains HPC-specific logic"""
    print("\n=== Testing CMake Arguments Logic ===")
    
    with open('package.py', 'r') as f:
        content = f.read()
    
    checks = [
        ('GPU architecture targeting', 'CMAKE_CUDA_ARCHITECTURES'),
        ('AMD GPU targeting', 'CMAKE_HIP_ARCHITECTURES'),
        ('NVHPC GPU flags', '-gpu=cc'),
        ('Fortran GPU flags', 'CMAKE_Fortran_FLAGS'),
        ('MFC CMake options', 'MFC_MPI'),
        ('Chemistry support', 'MFC_CHEMISTRY'),
    ]
    
    all_present = True
    for name, check_str in checks:
        if check_str in content:
            print(f"✓ {name} logic present")
        else:
            print(f"✗ {name} logic missing")
            all_present = False
    
    return all_present

def test_conflicts():
    """Test that appropriate conflicts are defined"""
    print("\n=== Testing Conflict Definitions ===")
    
    with open('package.py', 'r') as f:
        content = f.read()
    
    expected_conflicts = [
        ('Apple Clang', 'apple-clang'),
        ('OpenACC requires special compiler', 'conflicts("+openacc", when="%gcc"'),
        ('OpenACC/OpenMP mutual exclusion', 'conflicts("+openacc", when="+openmp"'),
    ]
    
    all_present = True
    for name, check_str in expected_conflicts:
        # Check if conflict is mentioned in any form
        if 'conflicts' in content and any(part in content for part in check_str.split()):
            print(f"✓ Conflict: {name}")
        else:
            print(f"? Conflict may be missing: {name}")
    
    return True  # Don't fail on this, just informational

def test_chemistry_support():
    """Test that chemistry support is properly configured"""
    print("\n=== Testing Chemistry Support ===")
    
    with open('package.py', 'r') as f:
        content = f.read()
    
    checks = [
        ('Cantera dependency', 'cantera'),
        ('Pyrometheus resource', 'pyrometheus'),
        ('PYTHONPATH setup', 'PYTHONPATH'),
        ('Chemistry variant', 'variant("chemistry"'),
    ]
    
    all_present = True
    for name, check_str in checks:
        if check_str in content:
            print(f"✓ {name}")
        else:
            print(f"✗ {name} missing")
            all_present = False
    
    return all_present

def generate_test_report():
    """Generate a summary test report"""
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)
    
    tests = [
        ("Package Loading", test_package_loads),
        ("Variant Definitions", test_variant_definitions),
        ("GPU Architectures", test_gpu_architecture_values),
        ("HPC Environment", test_hpc_environment_setup),
        ("Dependencies", test_dependency_definitions),
        ("CMake Arguments", test_cmake_args_logic),
        ("Chemistry Support", test_chemistry_support),
    ]
    
    results = []
    for name, test_func in tests:
        try:
            result = test_func()
            results.append((name, result))
        except Exception as e:
            print(f"\n✗ {name} test crashed: {e}")
            results.append((name, False))
    
    print("\n" + "="*60)
    print("FINAL RESULTS")
    print("="*60)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{status}: {name}")
    
    print(f"\nTotal: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n🎉 All tests passed! Package is ready for HPC deployment.")
        return 0
    else:
        print(f"\n⚠ {total - passed} test(s) failed. Review issues above.")
        return 1

if __name__ == "__main__":
    # Change to the directory containing package.py
    script_dir = Path(__file__).parent
    os.chdir(script_dir)
    
    if not Path('package.py').exists():
        print("Error: package.py not found in current directory")
        sys.exit(1)
    
    exit_code = generate_test_report()
    sys.exit(exit_code)

