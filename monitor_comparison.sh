#!/usr/bin/env bash
# Monitor the coverage comparison runs

PID_FILE="/tmp/baseline_coverage.pid"
LOG_FILE="build/baseline_run.log"

echo "=== Coverage Comparison Monitor ==="
echo ""

if [ -f "$PID_FILE" ]; then
    PID=$(cat "$PID_FILE")
    if ps -p $PID > /dev/null 2>&1; then
        echo "✅ Baseline coverage is RUNNING (PID: $PID)"
        echo ""
        
        if [ -f "$LOG_FILE" ]; then
            LINES=$(wc -l < "$LOG_FILE" 2>/dev/null || echo "0")
            echo "📄 Log file: $LINES lines"
            echo ""
            
            # Check for specific markers
            if grep -q "Building with coverage" "$LOG_FILE" 2>/dev/null; then
                if grep -q "Running 100%" "$LOG_FILE" 2>/dev/null; then
                    echo "🔄 Phase: RUNNING TESTS"
                    grep "Processing:" "$LOG_FILE" 2>/dev/null | tail -1
                elif grep -q "Built target" "$LOG_FILE" 2>/dev/null; then
                    TARGETS=$(grep -c "Built target" "$LOG_FILE" 2>/dev/null)
                    echo "🔄 Phase: BUILDING (${TARGETS} targets complete)"
                else
                    echo "🔄 Phase: BUILDING"
                fi
            else
                echo "🔄 Phase: INITIALIZING"
            fi
            
            echo ""
            echo "=== Recent Log Output (last 10 lines) ==="
            tail -10 "$LOG_FILE" 2>/dev/null || echo "(log not accessible)"
        else
            echo "⏳ Log file not yet created (process initializing)"
        fi
    else
        echo "❌ Baseline coverage process stopped (PID $PID not found)"
        echo ""
        if [ -f "$LOG_FILE" ]; then
            echo "=== Last 20 log lines ==="
            tail -20 "$LOG_FILE"
        fi
    fi
else
    echo "❌ No PID file found - coverage not running"
fi

echo ""
echo "=== Monitor Commands ==="
echo "  Real-time: tail -f build/baseline_run.log"
echo "  Status:    ps -p $(cat $PID_FILE 2>/dev/null || echo 'PID')"
echo ""

