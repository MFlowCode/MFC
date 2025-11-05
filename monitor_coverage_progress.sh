#!/usr/bin/env bash
# Periodic monitoring script for the coverage run

PROCESS_PID=554
LOG_FILE="build/full_coverage_run.log"
STATUS_FILE="build/coverage_progress.txt"

# Check if process is still running
if ! ps -p $PROCESS_PID > /dev/null 2>&1; then
    echo "Process $PROCESS_PID is not running"
    if [ -f "$LOG_FILE" ]; then
        echo "Checking log for completion status..."
        tail -50 "$LOG_FILE"
    fi
    exit 1
fi

# Get current timestamp
echo "=== Coverage Run Progress Check - $(date) ===" | tee -a "$STATUS_FILE"
echo "" | tee -a "$STATUS_FILE"

# Check if log file exists
if [ -f "$LOG_FILE" ]; then
    LINES=$(wc -l < "$LOG_FILE")
    echo "✅ Process running (PID: $PROCESS_PID)" | tee -a "$STATUS_FILE"
    echo "📄 Log file: $LINES lines" | tee -a "$STATUS_FILE"
    echo "" | tee -a "$STATUS_FILE"
    
    # Check what phase we're in
    if grep -q "Running 100% of test suite" "$LOG_FILE" 2>/dev/null; then
        echo "🧪 Phase: TESTING (this takes 2-4 hours)" | tee -a "$STATUS_FILE"
        # Try to get test progress
        tail -50 "$LOG_FILE" | grep -E "^\[|test|Processing" | tail -5 | tee -a "$STATUS_FILE"
    elif grep -q "Building with coverage instrumentation" "$LOG_FILE" 2>/dev/null; then
        echo "🔨 Phase: BUILDING" | tee -a "$STATUS_FILE"
        # Check build progress
        tail -20 "$LOG_FILE" | grep -E "\[.*%\]|Building|Preprocessing" | tail -5 | tee -a "$STATUS_FILE"
    elif grep -q "Generating coverage reports" "$LOG_FILE" 2>/dev/null; then
        echo "📊 Phase: GENERATING REPORTS" | tee -a "$STATUS_FILE"
        tail -10 "$LOG_FILE" | tee -a "$STATUS_FILE"
    else
        echo "🔄 Phase: INITIALIZING" | tee -a "$STATUS_FILE"
        tail -10 "$LOG_FILE" | tee -a "$STATUS_FILE"
    fi
    
    echo "" | tee -a "$STATUS_FILE"
    
    # Check for completion
    if grep -q "FULL COVERAGE RUN COMPLETE" "$LOG_FILE" 2>/dev/null; then
        echo "🎉 COVERAGE RUN COMPLETE!" | tee -a "$STATUS_FILE"
        echo "" | tee -a "$STATUS_FILE"
        echo "=== Final Summary ===" | tee -a "$STATUS_FILE"
        grep -A 20 "COVERAGE SUMMARY" "$LOG_FILE" | tee -a "$STATUS_FILE"
        exit 0
    fi
else
    echo "⏳ Process running (PID: $PROCESS_PID) - log file not created yet" | tee -a "$STATUS_FILE"
    echo "   (Build output is going to terminal)" | tee -a "$STATUS_FILE"
fi

echo "" | tee -a "$STATUS_FILE"
echo "---" | tee -a "$STATUS_FILE"





