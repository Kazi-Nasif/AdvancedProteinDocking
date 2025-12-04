#!/bin/bash
cd "/mnt/bst/bdeng2/knasif/Protein Docking/AdvancedProteinDocking_Step1"

echo "Starting production training..."

# Run with unbuffered output using environment variable
PYTHONUNBUFFERED=1 nohup python scripts/train_production.py > training_production.log 2>&1 &

PID=$!
echo "Training started! PID: $PID"
echo $PID > training_production.pid
echo "Monitor: tail -f training_production.log"
