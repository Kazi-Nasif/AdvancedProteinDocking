#!/bin/bash

echo "================================================================================"
echo "PROTEIN DOCKING VERIFICATION TOOL - DEMONSTRATION"
echo "================================================================================"
echo ""
echo "Testing on 3 representative targets..."
echo ""

# Test 3 targets
targets=("1A2K rigid" "1AHW rigid" "1AVX medium")

for item in "${targets[@]}"; do
    target=$(echo $item | cut -d' ' -f1)
    difficulty=$(echo $item | cut -d' ' -f2)
    
    echo "Testing $target ($difficulty)..."
    python scripts/predict_docking.py \
        --target $target \
        --difficulty $difficulty \
        --evaluate \
        --output-dir demo_results/$target \
        2>&1 | grep -A 20 "RESULTS"
    echo ""
done

echo "================================================================================"
echo "DEMO COMPLETE"
echo "================================================================================"
echo ""
echo "Results saved in demo_results/"
echo "View metrics: cat demo_results/*/metrics.json"
echo ""
