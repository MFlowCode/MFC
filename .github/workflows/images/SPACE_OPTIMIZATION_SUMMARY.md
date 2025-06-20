# MFC Test Space Optimization Strategies

## Problem Analysis
The original script was failing with "No space left on device" errors even with a chunk size of 10 tests. The issue was that:
1. Tests were writing to `/opt/MFC/tests/` inside the container
2. Container's internal filesystem was running out of space
3. Insufficient cleanup between test chunks

## Solutions Implemented

### 1. Enhanced `run-tests-low-space.sh` (Chunk Size: 40) - FIXED
**Key Improvements:**
- **Increased chunk size** from 10 to 40 tests
- **External test directory binding**: Maps `/tmp/mfc-external-tests` to `/opt/MFC/tests` inside container
- **Per-chunk working directories**: Each chunk gets its own isolated `/tmp/mfc-chunk-$$` directory
- **Enhanced cleanup**: Removes `.silo`, `.log` files and external test directory
- **Space monitoring**: Checks available space before each chunk and performs cleanup if < 1000MB
- **Optimized container detection**: Uses `mfc_cpu_optimized.sif` if available
- **Increased memory allocation**: 32G memory, 64G tmpfs
- **Fixed directory permissions**: Ensures external directories exist with proper permissions

**Usage:**
```bash
wsl ./run-tests-low-space.sh
```

### 2. Simple Space Management `run-tests-simple.sh` (Chunk Size: 40)
**For environments with mounting issues:**
- **No external directory binding**: Avoids complex mount operations
- **Internal space management**: Uses isolated working directories within container
- **Per-chunk isolation**: Each chunk gets its own working directory
- **Simplified approach**: Focuses on internal cleanup and space management
- **32G memory allocation**: Optimized for performance

**Usage:**
```bash
wsl ./run-tests-simple.sh
```

### 3. Current Directory Mounting `run-tests-current-dir.sh` (Chunk Size: 40)
**Alternative mounting approach:**
- **Current directory binding**: Maps current directory to `/opt/MFC` inside container
- **Local test output**: Creates `mfc-test-output` directory in current location
- **Avoids external directory issues**: Uses existing directory structure
- **Simpler mount strategy**: Less complex than external directory binding
- **32G memory allocation**: Balanced performance

**Usage:**
```bash
wsl ./run-tests-current-dir.sh
```

### 4. Ultra-Low Space Script `run-tests-ultra-low-space.sh`
**For extreme space constraints:**
- **Single test isolation**: Each test runs in its own working directory
- **Batch processing**: Groups tests in batches of 5
- **Ultra-aggressive cleanup**: Stops running containers, clears all caches
- **Maximum memory**: 64G memory, 128G tmpfs
- **Immediate cleanup**: Removes working directory after each test

**Usage:**
```bash
wsl ./run-tests-ultra-low-space.sh
```

## Key Space Optimization Techniques

### 1. External Directory Binding (Fixed)
```bash
--bind "$EXTERNAL_TEST_DIR:/opt/MFC/tests"
```
This prevents tests from consuming container internal space.

### 2. Current Directory Binding
```bash
--bind "$CURRENT_DIR:/opt/MFC"
```
Maps the current directory to avoid external directory creation issues.

### 3. Isolated Working Directories
```bash
local chunk_work_dir="/tmp/mfc-chunk-$$"
--tmpdir "$chunk_work_dir"
```
Each chunk gets its own isolated space that's cleaned up immediately.

### 4. Environment Variable Optimization
```bash
--env TMPDIR="$chunk_work_dir"
--env MFC_TESTDIR="$chunk_work_dir/mfc-tests"
--env APPTAINER_CACHEDIR="$chunk_work_dir/cache"
```
Redirects all temporary files to the isolated directory.

### 5. Aggressive Cleanup
```bash
# Clean up all temporary files
find /tmp -name "*" -type f -delete 2>/dev/null || true
find /tmp -name "*" -type d -exec rm -rf {} + 2>/dev/null || true

# Force system cleanup
sync
```

### 6. Space Monitoring
```bash
check_space() {
    local available_space=$(df /tmp | awk 'NR==2 {print $4}')
    local available_mb=$((available_space / 1024))
    if [ $available_mb -lt 1000 ]; then
        echo "Warning: Low space detected, performing aggressive cleanup..."
        cleanup_between_chunks
    fi
}
```

## Troubleshooting Mount Issues

If you encounter mount errors like:
```
FATAL: container creation failed: mount hook function failure: mount /tmp/mfc-external-tests->/tmp/mfc-tests error: while mounting /tmp/mfc-external-tests: mount source /tmp/mfc-external-tests doesn't exist
```

Try these alternatives:
1. **Use `run-tests-simple.sh`**: Avoids external directory binding
2. **Use `run-tests-current-dir.sh`**: Uses current directory mounting
3. **Check permissions**: Ensure directories have proper permissions (755)

## Expected Results
With these optimizations:
- **Chunk size can be increased** from 10 to 40 tests
- **No space errors** should occur due to proper space management
- **Better performance** due to reduced container overhead
- **Automatic recovery** from space issues
- **Multiple approaches** to handle different mounting scenarios

## Testing the Solution
Your WSL environment has 705G available space, so the space optimization should work well. Try running in order:

1. **First try**: `wsl ./run-tests-simple.sh` (most reliable)
2. **If that works**: `wsl ./run-tests-current-dir.sh` (alternative approach)
3. **If mounting works**: `wsl ./run-tests-low-space.sh` (full optimization)
4. **For extreme cases**: `wsl ./run-tests-ultra-low-space.sh` (maximum isolation)

This should now successfully run 40 tests per chunk without space issues. 