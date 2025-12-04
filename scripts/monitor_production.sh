#!/bin/bash
# Monitor Production Training

echo "================================"
echo "TRAINING STATUS"
echo "================================"
echo ""

# Check if running
if [ -f training_production.pid ]; then
    PID=$(cat training_production.pid)
    if ps -p $PID > /dev/null; then
        echo "✓ Training is RUNNING (PID: $PID)"
        
        # Show CPU/Memory usage
        ps -p $PID -o pid,pcpu,pmem,etime,cmd | tail -1
    else
        echo "✗ Training is NOT RUNNING"
        echo "  (May have finished or crashed)"
    fi
else
    echo "✗ No PID file found"
    echo "  Training may not have started yet"
fi

echo ""
echo "================================"
echo "LATEST PROGRESS (last 40 lines)"
echo "================================"
echo ""

# Show recent log
if [ -f training_production.log ]; then
    tail -n 40 training_production.log
else
    echo "No log file yet"
fi

echo ""
echo "================================"
echo "EXPERIMENTS"
echo "================================"
echo ""

# Show experiments
if [ -d experiments ]; then
    ls -lth experiments/ | head -5
else
    echo "No experiments directory yet"
fi

echo ""
echo "================================"
echo "COMMANDS"
echo "================================"
echo "Watch live: tail -f training_production.log"
echo "Stop training: kill \$(cat training_production.pid)"
echo "Check metrics: cat experiments/*/metrics.json | jq"
echo ""
