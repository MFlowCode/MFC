#!/usr/bin/env bash
# Monitor coverage run progress

echo "=== MFC Coverage Run Monitor ==="
echo ""

# Check if process is running
if ps aux | grep -q "[c]overage.sh"; then
    echo "✅ Coverage process is RUNNING"
    echo ""
    
    # Show current stage
    echo "Current stage:"
    tail -5 /tmp/coverage_nohup.log 2>/dev/null || tail -5 build/coverage_run.log 2>/dev/null || echo "  Log not available"
    echo ""
    
    # Show progress
    LINES=$(wc -l < /tmp/coverage_nohup.log 2>/dev/null || echo 0)
    echo "Log size: $LINES lines"
    echo ""
    
    # Check for key milestones
    if grep -q "Built target simulation" /tmp/coverage_nohup.log 2>/dev/null; then
        echo "✅ Build phase: COMPLETE"
    else
        echo "🔄 Build phase: IN PROGRESS"
        BUILD_PCT=$(tail -100 /tmp/coverage_nohup.log 2>/dev/null | grep -oE '\[[0-9]+%\]' | tail -1 || echo "[?%]")
        echo "   $BUILD_PCT"
    fi
    
    if grep -q "Running tests" /tmp/coverage_nohup.log 2>/dev/null; then
        echo "🔄 Test phase: IN PROGRESS"
        TEST_COUNT=$(grep -c "OK" /tmp/coverage_nohup.log 2>/dev/null || echo 0)
        echo "   $TEST_COUNT tests completed"
    elif grep -q "Generating coverage reports" /tmp/coverage_nohup.log 2>/dev/null; then
        echo "✅ Test phase: COMPLETE"
        echo "🔄 Report generation: IN PROGRESS"
    fi
    
    echo ""
    echo "Monitor in real-time:"
    echo "  tail -f /tmp/coverage_nohup.log"
    
else
    echo "⚠️  Coverage process is NOT running"
    echo ""
    
    # Check for results
    if [ -f build/coverage/index.html ]; then
        echo "✅ Coverage report COMPLETE!"
        echo ""
        echo "View results:"
        echo "  open build/coverage/index.html"
        echo "  cat build/coverage/summary.txt"
        echo ""
        
        # Show summary if available
        if [ -f build/coverage/summary.txt ]; then
            echo "Summary:"
            cat build/coverage/summary.txt
        fi
    else
        echo "❌ No coverage report found"
        echo ""
        echo "Check logs:"
        echo "  tail -100 /tmp/coverage_nohup.log"
        echo "  tail -100 build/coverage_run.log"
    fi
fi

echo ""
echo "=== End Monitor ==="






