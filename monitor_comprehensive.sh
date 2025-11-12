#!/usr/bin/env bash
# Monitoring script for comprehensive coverage run

PID=6343
LOG_FILE="build/coverage_comparison_full/full_run.log"
STATUS_FILE="build/coverage_comparison_full/monitor_status.txt"

while true; do
    TIMESTAMP=$(date "+%Y-%m-%d %H:%M:%S")
    
    echo "=== Status Check at $TIMESTAMP ===" | tee -a "$STATUS_FILE"
    
    # Check if process is running
    if ! ps -p $PID > /dev/null 2>&1; then
        echo "❌ Process $PID has ended!" | tee -a "$STATUS_FILE"
        echo "" | tee -a "$STATUS_FILE"
        
        # Check for completion
        if [ -f "$LOG_FILE" ]; then
            echo "Log file exists. Checking for completion..." | tee -a "$STATUS_FILE"
            if grep -q "ALL PHASES COMPLETE" "$LOG_FILE"; then
                echo "✅ RUN COMPLETED SUCCESSFULLY!" | tee -a "$STATUS_FILE"
            else
                echo "⚠️ Process ended but may not have completed successfully" | tee -a "$STATUS_FILE"
                echo "Last 20 lines of log:" | tee -a "$STATUS_FILE"
                tail -20 "$LOG_FILE" | tee -a "$STATUS_FILE"
            fi
        else
            echo "❌ Process ended and no log file found" | tee -a "$STATUS_FILE"
        fi
        break
    fi
    
    echo "✅ Process $PID is running" | tee -a "$STATUS_FILE"
    
    # Check for log file
    if [ -f "$LOG_FILE" ]; then
        LINES=$(wc -l < "$LOG_FILE")
        SIZE=$(ls -lh "$LOG_FILE" | awk '{print $5}')
        echo "📄 Log: $LINES lines, $SIZE" | tee -a "$STATUS_FILE"
        
        # Show recent activity
        echo "Recent activity:" | tee -a "$STATUS_FILE"
        tail -5 "$LOG_FILE" | tee -a "$STATUS_FILE"
    else
        echo "⏳ Log file not yet created (buffering)" | tee -a "$STATUS_FILE"
    fi
    
    # Count coverage files
    GCDA_COUNT=$(find build/staging -name "*.gcda" 2>/dev/null | wc -l)
    echo "📊 Coverage files: $GCDA_COUNT .gcda files" | tee -a "$STATUS_FILE"
    
    # Check which phase we're in
    if [ -f "$LOG_FILE" ]; then
        if grep -q "PHASE 1:" "$LOG_FILE" && ! grep -q "Phase 1 Complete" "$LOG_FILE"; then
            echo "🔄 Currently in: PHASE 1 (Baseline)" | tee -a "$STATUS_FILE"
        elif grep -q "Phase 1 Complete" "$LOG_FILE" && ! grep -q "PHASE 2:" "$LOG_FILE"; then
            echo "🔄 Currently in: PHASE 2 (Adding Tests)" | tee -a "$STATUS_FILE"
        elif grep -q "PHASE 3:" "$LOG_FILE" && ! grep -q "Phase 3 Complete" "$LOG_FILE"; then
            echo "🔄 Currently in: PHASE 3 (Expanded)" | tee -a "$STATUS_FILE"
        elif grep -q "Phase 3 Complete" "$LOG_FILE"; then
            echo "🔄 Currently in: PHASE 4 (Analysis)" | tee -a "$STATUS_FILE"
        fi
    fi
    
    echo "---" | tee -a "$STATUS_FILE"
    echo "" | tee -a "$STATUS_FILE"
    
    # Wait 5 minutes
    sleep 300
done

echo "Monitoring ended at $(date)" | tee -a "$STATUS_FILE"





