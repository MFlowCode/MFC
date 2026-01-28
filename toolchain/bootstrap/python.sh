#!/bin/bash

MFC_PYTHON_MIN_MAJOR=3
MFC_PYTHON_MIN_MINOR=11
MFC_PYTHON_MIN_STR="$MFC_PYTHON_MIN_MAJOR.$MFC_PYTHON_MIN_MINOR"

is_python_compatible() {
    if ! ${1:-python3} -c "import sys; exit(int(not (sys.version_info[0]==$MFC_PYTHON_MIN_MAJOR and sys.version_info[1] >= $MFC_PYTHON_MIN_MINOR)))"; then
        return 1
    fi

    return 0
}

assert_python_compatible() {
    if ! is_python_compatible $1; then
        error "$MAGENTA$(${1:-python3} --version)$COLOR_RESET (${1:-python3}) is out of date. Required >= $MAGENTA$MFC_PYTHON_MIN_STR$COLOR_RESET."

        exit 1
    fi
}

if [ -f "$(pwd)/build/venv/bin/activate" ]; then
    if ! is_python_compatible "$(pwd)/build/venv/bin/python3"; then
        warn "$MAGENTA""Python$COLOR_RESET is outdated inside the Virtualenv."
        warn "Deleting the Virtualenv and starting from scratch..."

        rm -r "$(pwd)/build/venv"
    fi
fi

# Only bootstrap pip if we don't already have a venv
if [ ! -f "$(pwd)/build/venv/bin/activate" ]; then
    # Check whether python3 is in the $PATH / is accessible.
    if ! command -v python3 > /dev/null 2>&1; then
        error "Couldn't find$MAGENTA Python$COLOR_RESET. Please ensure it is discoverable."

        exit 1
    fi

    assert_python_compatible

    # Check if pip is already available as a Python module
    # This works on both laptops and HPC systems with module-loaded Python
    if ! python3 -c "import pip" > /dev/null 2>&1; then
        warn "$MAGENTA""Python$COLOR_RESET's$MAGENTA PIP$COLOR_RESET is not installed."
        
        # Try ensurepip first (standard library, safe)
        log "Attempting to install pip via ensurepip..."
        if python3 -m ensurepip --upgrade 2>/dev/null; then
            ok "Installed pip via ensurepip."
        else
            # Fall back to get-pip.py only if ensurepip fails
            get_pip_url="https://bootstrap.pypa.io/pip/get-pip.py"
            log "Downloading$MAGENTA Python$COLOR_RESET's$MAGENTA PIP$COLOR_RESET from $get_pip_url..."

            if ! wget -O "$(pwd)/build/get-pip.py" "$get_pip_url"; then
                error "Couldn't download get-pip.py."
                exit 1
            fi

            # Suppress PIP version warning (out of date)
            export PIP_DISABLE_PIP_VERSION_CHECK=1
            if ! python3 "$(pwd)/build/get-pip.py" --user; then
                error "Couldn't install$MAGENTA pip$COLOR_RESET with get-pip.py"
                exit 1
            fi

            ok "Installed pip via get-pip.py."
            
            # Ensure user-site bin directory is on PATH for this session
            user_base_bin="$(python3 -m site --user-base)/bin"
            if [ -d "$user_base_bin" ]; then
                export PATH="$user_base_bin:$PATH"
            fi
        fi
    fi
fi


# Create a Python virtualenv if it hasn't already been created
if [ ! -f "$(pwd)/build/venv/bin/activate" ]; then
    assert_python_compatible

    if ! python3 -m venv "$(pwd)/build/venv"; then
        error "Failed to create a$MAGENTA Python$COLOR_RESET virtual environment. Delete the build/venv folder and try again."

        exit 1
    fi

    ok "Created a$MAGENTA Python$COLOR_RESET virtual environment (venv)."

    rm "$(pwd)/build/pyproject.toml" > /dev/null 2>&1 || true
fi


# If the user wishes to enter the python venv
# mfc.sh must be source'd for this to work
if [ "$1" == "venv" ]; then
    if [[ "$VIRTUAL_ENV" != "" ]]; then
        # Already inside the venv, exit out of it

        log " > Exiting the$MAGENTA MFC Python$COLOR_RESET virtual environment."

        deactivate
    else
        # Enter the venv
        log " > Entering the$MAGENTA MFC Python$COLOR_RESET virtual environment."

        source "$(pwd)/build/venv/bin/activate"

        log " > To exit, you can do any of the following:"
        log "    - Run 'deactivate'."
        log "    - Run '. ./mfc.sh venv'."
        log "    - Close your terminal."
    fi

    return
fi


# Activate the Python venv
source "$(pwd)/build/venv/bin/activate"
ok "(venv) Entered the $MAGENTA$(python3 --version)$COLOR_RESET virtual environment (>= $MAGENTA$MFC_PYTHON_MIN_STR$COLOR_RESET)."


# Install Python dependencies if, either:
# - This script is running for the first time
# (or)
# - The pyproject.toml file has changed
if ! cmp "$(pwd)/toolchain/pyproject.toml" "$(pwd)/build/pyproject.toml" > /dev/null 2>&1; then
    # Check if this is a fresh install (no previous pyproject.toml in build/)
    if [ ! -f "$(pwd)/build/pyproject.toml" ]; then
        log "(venv) Installing$MAGENTA ~70 Python packages$COLOR_RESET. This may take a few minutes..."
    else
        log "(venv) Updating Python dependencies..."
    fi

    next_arg=0
    nthreads=1
    for arg in "$@"; do
        if [ "$arg" == "-j" ] || [ "$arg" == "--jobs" ]; then
            next_arg=1
            continue
        fi
        if [ "$next_arg" == 1 ]; then
            next_arg=0
            nthreads=$arg
            continue
        fi
    done

    # Run pip and show progress
    PIP_LOG="$(pwd)/build/.pip_install.log"

    # Start pip in background, capturing output
    PIP_DISABLE_PIP_VERSION_CHECK=1 MAKEFLAGS=$nthreads pip3 install "$(pwd)/toolchain" > "$PIP_LOG" 2>&1 &
    PIP_PID=$!

    # Check if we're in an interactive terminal
    if [ -t 1 ]; then
        IS_TTY=1
    else
        IS_TTY=0
    fi

    # Progress bar configuration
    # Two phases: Collecting (50%) and Installing (50%)
    TOTAL_PKGS=70
    BAR_WIDTH=30
    LAST_MILESTONE=0
    LAST_PHASE=""
    START_TIME=$SECONDS

    while kill -0 $PIP_PID 2>/dev/null; do
        # Determine current phase and count from log
        PHASE="resolving"
        COUNT=0
        BUILD_COUNT=0

        if [ -f "$PIP_LOG" ]; then
            # Count packages being collected (dependency resolution)
            COUNT=$(grep -c "^Collecting" "$PIP_LOG" 2>/dev/null | tr -d '[:space:]')
            COUNT=${COUNT:-0}
            if ! [[ "$COUNT" =~ ^[0-9]+$ ]]; then
                COUNT=0
            fi

            # Count wheels being built
            BUILD_COUNT=$(grep -c "Building wheel" "$PIP_LOG" 2>/dev/null | tr -d '[:space:]')
            BUILD_COUNT=${BUILD_COUNT:-0}
            if ! [[ "$BUILD_COUNT" =~ ^[0-9]+$ ]]; then
                BUILD_COUNT=0
            fi

            # Check if we're in the installing phase
            if grep -q "Installing collected packages" "$PIP_LOG" 2>/dev/null; then
                PHASE="installing"
            elif [ "$BUILD_COUNT" -gt 0 ]; then
                PHASE="building"
            fi
        fi

        ELAPSED=$((SECONDS - START_TIME))

        if [ "$IS_TTY" = "1" ]; then
            # Calculate progress based on phase
            # Phase 1 (0-60%): Collecting dependencies
            # Phase 2 (60-80%): Building wheels
            # Phase 3 (80-100%): Installing

            if [ "$COUNT" -gt "$TOTAL_PKGS" ]; then
                TOTAL_PKGS=$COUNT
            fi

            if [ "$PHASE" = "installing" ]; then
                PERCENT=90
                STATUS="Installing..."
            elif [ "$PHASE" = "building" ]; then
                # During building, progress from 60-80%
                PERCENT=$((60 + BUILD_COUNT * 2))
                if [ "$PERCENT" -gt 80 ]; then
                    PERCENT=80
                fi
                STATUS="Building ($BUILD_COUNT wheels)"
            else
                # During collecting, progress from 0-60%
                PERCENT=$((COUNT * 60 / TOTAL_PKGS))
                STATUS="$COUNT packages"
            fi

            FILLED=$((PERCENT * BAR_WIDTH / 100))
            EMPTY=$((BAR_WIDTH - FILLED))

            # Build the bar with Unicode blocks
            BAR=""
            for ((i=0; i<FILLED; i++)); do BAR="${BAR}█"; done
            for ((i=0; i<EMPTY; i++)); do BAR="${BAR}░"; done

            # Format time
            MINS=$((ELAPSED / 60))
            SECS=$((ELAPSED % 60))
            if [ "$MINS" -gt 0 ]; then
                TIME_STR="${MINS}m${SECS}s"
            else
                TIME_STR="${SECS}s"
            fi

            # Print progress bar
            printf "\r  ${CYAN}│${BAR}│${COLOR_RESET} %3d%% (%s, %s)  " "$PERCENT" "$STATUS" "$TIME_STR"
        else
            # Non-interactive: print milestone updates and phase changes
            if [ "$PHASE" != "$LAST_PHASE" ]; then
                case "$PHASE" in
                    "building")
                        log "(venv) Building wheels..."
                        ;;
                    "installing")
                        log "(venv) Installing packages..."
                        ;;
                esac
                LAST_PHASE="$PHASE"
            fi

            if [ "$PHASE" = "resolving" ]; then
                MILESTONE=$((COUNT / 20 * 20))
                if [ "$MILESTONE" -gt "$LAST_MILESTONE" ] && [ "$MILESTONE" -gt 0 ]; then
                    LAST_MILESTONE=$MILESTONE
                    log "(venv) Resolving: ~$COUNT packages..."
                fi
            fi
        fi

        sleep 0.3
    done

    # Wait for pip to finish and get exit code
    wait $PIP_PID
    PIP_EXIT=$?

    # Clear the progress line if in terminal
    if [ "$IS_TTY" = "1" ]; then
        printf "\r%60s\r" " "
    fi

    if [ $PIP_EXIT -ne 0 ]; then
        error "(venv) Installation failed. See output below:"
        echo ""
        cat "$PIP_LOG"
        echo ""

        log   "(venv) Exiting the$MAGENTA Python$COLOR_RESET virtual environment."
        deactivate
        rm -f "$PIP_LOG"
        exit 1
    fi

    rm -f "$PIP_LOG"
    ok "(venv) Installation succeeded."

    # Save the new/current pyproject.toml
    cp "$(pwd)/toolchain/pyproject.toml" "$(pwd)/build/"
fi
