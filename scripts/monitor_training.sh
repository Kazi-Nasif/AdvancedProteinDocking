#!/bin/bash
# Monitor Training Progress

echo "================================"
echo "TRAINING PROGRESS MONITOR"
echo "================================"
echo ""

# Check if training is running
if [ -f training.pid ]; then
    PID=$(cat training.pid)
    if ps -p $PID > /dev/null; then
        echo "✓ Training is RUNNING (PID: $PID)"
    else
        echo "✗ Training is NOT RUNNING"
        echo "  (Process may have finished or crashed)"
    fi
else
    echo "✗ No training.pid file found"
    echo "  (Training may not have been started yet)"
fi

echo ""
echo "================================"
echo "RECENT LOG OUTPUT"
echo "================================"
echo ""

# Show last 30 lines of log
if [ -f training_output.log ]; then
    tail -n 30 training_output.log
else
    echo "No log file found yet"
fi

echo ""
echo "================================"
echo "COMMANDS"
echo "================================"
echo "Monitor live: tail -f training_output.log"
echo "Stop training: kill \$(cat training.pid)"
echo "Check experiments: ls -lh experiments/"
echo ""
