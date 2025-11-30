#!/usr/bin/env python3
"""
Compare local golden files with those from the mflowcode/mfc repository.
Prints the maximum difference for each test case.

Usage:
    python3 compare_goldenfiles.py [--limit N] [--test UUID]
"""

import os
import sys
import argparse
import numpy as np
import urllib.request
import urllib.error
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed


def get_test_case_info(test_path):
    """Extract test case information from case.py file."""
    case_py = test_path / "case.py"
    if not case_py.exists():
        return "Unknown"
    
    try:
        with open(case_py, 'r') as f:
            lines = f.readlines()
            # Look for line 4 (index 3) which has the test case description
            # Format: # 1D -> Example -> hypo_2materials
            if len(lines) >= 4:
                line = lines[3].strip()
                if line.startswith('#'):
                    desc = line[1:].strip()  # Remove leading '#' and whitespace
                    return desc
    except Exception:
        pass
    
    return "Unknown"


def fetch_repo_golden_file(test_uuid):
    """Fetch a golden file from GitHub raw content."""
    url = f"https://raw.githubusercontent.com/mflowcode/mfc/master/tests/{test_uuid}/golden.txt"
    try:
        with urllib.request.urlopen(url, timeout=10) as response:
            return response.read().decode('utf-8')
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None
        return None
    except Exception:
        return None


def parse_golden_data(text):
    """Parse golden file text into numpy array, skipping non-numeric tokens."""
    values = text.split()
    floats = []
    for v in values:
        try:
            floats.append(float(v))
        except ValueError:
            # Skip non-numeric tokens (like file paths)
            continue
    return np.array(floats) if floats else np.array([])


def compare_single_test(test_uuid, local_path, tolerance=1e-12):
    """Compare a single test's golden file."""
    local_test_path = local_path / test_uuid
    local_golden_path = local_test_path / "golden.txt"
    
    if not local_golden_path.exists():
        return test_uuid, "Unknown", "SKIP", "No local golden file"
    
    # Get test case info
    case_info = get_test_case_info(local_test_path)
    
    # Fetch repository version
    repo_content = fetch_repo_golden_file(test_uuid)
    if repo_content is None:
        return test_uuid, case_info, "SKIP", "Not found in repository"
    
    try:
        # Load and compare data
        with open(local_golden_path, 'r') as f:
            local_content = f.read()
        
        repo_data = parse_golden_data(repo_content)
        local_data = parse_golden_data(local_content)
        
        if repo_data.size == 0 or local_data.size == 0:
            return test_uuid, case_info, "SKIP", "Empty or no numeric data"
        
        if repo_data.shape != local_data.shape:
            return test_uuid, case_info, "DIFF", f"Shape mismatch. Repo: {repo_data.shape}, Local: {local_data.shape}"
        
        # Compute differences
        diff = np.abs(repo_data - local_data)
        max_diff = np.max(diff)
        
        if max_diff > tolerance:
            return test_uuid, case_info, "DIFF", f"Max difference = {max_diff:.10e}"
        else:
            return test_uuid, case_info, "OK", f"Max difference = {max_diff:.10e}"
    
    except Exception as e:
        return test_uuid, case_info, "ERROR", str(e)


def compare_golden_files(local_tests_dir, limit=None, specific_test=None, max_workers=10, show_ok=False):
    """Compare local golden files with repository versions."""
    local_path = Path(local_tests_dir)
    
    if not local_path.exists():
        print(f"Error: Local tests directory not found: {local_tests_dir}")
        return
    
    # Get all local test directories
    if specific_test:
        test_dirs = [specific_test]
    else:
        test_dirs = sorted([d for d in os.listdir(local_path) 
                           if os.path.isdir(local_path / d)])
        if limit:
            test_dirs = test_dirs[:limit]
    
    print(f"Comparing {len(test_dirs)} test(s)...")
    print("-" * 80)
    
    differences_found = False
    ok_count = 0
    diff_count = 0
    skip_count = 0
    error_count = 0
    max_diff_overall = 0.0
    max_diff_test = None
    max_diff_case = None
    
    results = []
    
    # Use ThreadPoolExecutor for parallel fetching
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(compare_single_test, test_uuid, local_path): test_uuid 
                  for test_uuid in test_dirs}
        
        for future in as_completed(futures):
            test_uuid, case_info, status, message = future.result()
            results.append((test_uuid, case_info, status, message))
            
            if status == "OK":
                ok_count += 1
                # Extract the max diff value
                try:
                    diff_val = float(message.split('=')[1].strip())
                    if diff_val > max_diff_overall:
                        max_diff_overall = diff_val
                        max_diff_test = test_uuid
                        max_diff_case = case_info
                except:
                    pass
            elif status == "DIFF":
                diff_count += 1
                differences_found = True
                # Extract the max diff value if present
                try:
                    if 'Max difference' in message:
                        diff_val = float(message.split('=')[1].strip())
                        if diff_val > max_diff_overall:
                            max_diff_overall = diff_val
                            max_diff_test = test_uuid
                            max_diff_case = case_info
                except:
                    pass
            elif status == "SKIP":
                skip_count += 1
            elif status == "ERROR":
                error_count += 1
                differences_found = True
    
    # Sort results by status and test UUID
    status_order = {"DIFF": 0, "ERROR": 1, "OK": 2, "SKIP": 3}
    results.sort(key=lambda x: (status_order[x[2]], x[0]))
    
    # Print results
    for test_uuid, case_info, status, message in results:
        if status == "DIFF":
            print(f"  [DIFF]  {test_uuid}: {case_info}")
            print(f"          {message}")
        elif status == "ERROR":
            print(f"  [ERROR] {test_uuid}: {case_info}")
            print(f"          {message}")
        elif status == "OK" and show_ok:
            print(f"  [OK]    {test_uuid}: {case_info}")
            print(f"          {message}")
        # Skip printing SKIP status unless explicitly requested
    
    print("-" * 80)
    print(f"\nSummary:")
    print(f"  OK:     {ok_count}")
    print(f"  DIFF:   {diff_count}")
    print(f"  SKIP:   {skip_count}")
    print(f"  ERROR:  {error_count}")
    
    if max_diff_test:
        print(f"\nMaximum difference: {max_diff_overall:.10e}")
        print(f"  Test: {max_diff_test}")
        print(f"  Case: {max_diff_case}")
    
    if not differences_found:
        print("\n✓ No significant differences found!")
    else:
        print(f"\n✗ Differences found in {diff_count} golden files")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare MFC golden files with repository')
    parser.add_argument('--limit', type=int, help='Limit number of tests to compare')
    parser.add_argument('--test', type=str, help='Compare specific test UUID')
    parser.add_argument('--workers', type=int, default=10, help='Number of parallel workers (default: 10)')
    parser.add_argument('--show-ok', action='store_true', help='Show OK results in detail')
    
    args = parser.parse_args()
    
    script_dir = Path(__file__).parent.absolute()
    tests_dir = script_dir / "tests"
    compare_golden_files(tests_dir, limit=args.limit, specific_test=args.test, 
                        max_workers=args.workers, show_ok=args.show_ok)
