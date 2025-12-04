#!/bin/bash
# Quick Verification Test Script
# Tests the prediction tool on benchmark targets

echo "================================================================================"
echo "PROTEIN DOCKING TOOL - QUICK VERIFICATION TEST"
echo "================================================================================"
echo ""

# Test 1: Single target
echo "TEST 1: Predicting on benchmark target 1A2K..."
echo "--------------------------------------------------------------------------------"
python scripts/predict_docking.py \
    --target 1A2K \
    --difficulty rigid \
    --evaluate \
    --output-dir quick_test/1A2K

if [ $? -eq 0 ]; then
    echo "✓ Test 1 PASSED"
else
    echo "✗ Test 1 FAILED"
    exit 1
fi

echo ""
echo "================================================================================"
echo "TEST RESULTS"
echo "================================================================================"

# Display results
if [ -f "quick_test/1A2K/metrics.json" ]; then
    echo ""
    echo "Metrics for 1A2K:"
    cat quick_test/1A2K/metrics.json
    echo ""
    
    # Extract DockQ score
    dockq=$(python -c "import json; print(json.load(open('quick_test/1A2K/metrics.json'))['dockq'])")
    success=$(python -c "import json; print(json.load(open('quick_test/1A2K/metrics.json'))['success'])")
    
    echo "DockQ Score: $dockq"
    echo "Success: $success"
    
    if [ "$success" == "True" ]; then
        echo ""
        echo "✓✓✓ VERIFICATION PASSED! ✓✓✓"
        echo "The model successfully predicted the docked structure!"
        echo ""
        echo "Results saved to: quick_test/1A2K/"
        echo "  - metrics.json (all metrics)"
        echo "  - predicted_complex.pdb (3D structure)"
    else
        echo ""
        echo "⚠ Prediction did not meet success threshold (DockQ > 0.23)"
        echo "This may indicate an issue with the model or data"
    fi
else
    echo "✗ Results file not found"
    exit 1
fi

echo ""
echo "================================================================================"
echo "NEXT STEPS"
echo "================================================================================"
echo ""
echo "1. View the predicted structure:"
echo "   pymol quick_test/1A2K/predicted_complex.pdb"
echo ""
echo "2. Test more targets:"
echo "   python scripts/predict_docking.py --target 1AHW --difficulty rigid --evaluate"
echo ""
echo "3. Upload your own structures:"
echo "   python scripts/predict_docking.py --receptor your_receptor.pdb --ligand your_ligand.pdb"
echo ""
echo "See PREDICTION_TOOL_GUIDE.md for full documentation"
echo "================================================================================"
