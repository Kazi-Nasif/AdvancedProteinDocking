#!/bin/bash
# Start Training with nohup
# This script runs training in the background

cd "/mnt/bst/bdeng2/knasif/Protein Docking/AdvancedProteinDocking_Step1"

echo "Starting protein docking training..."
echo "Logs will be saved to: nohup.out"
echo ""

# Run training with nohup
nohup python scripts/train.py > training_output.log 2>&1 &

# Get process ID
PID=$!

echo "Training started!"
echo "Process ID: $PID"
echo ""
echo "To monitor progress:"
echo "  tail -f training_output.log"
echo ""
echo "To check if still running:"
echo "  ps -p $PID"
echo ""
echo "To stop training:"
echo "  kill $PID"
echo ""

# Save PID to file
echo $PID > training.pid
echo "PID saved to training.pid"
