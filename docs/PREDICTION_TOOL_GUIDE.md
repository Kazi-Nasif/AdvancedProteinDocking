# 🔬 PROTEIN DOCKING PREDICTION TOOL - USER GUIDE

## 📋 Overview

This tool allows you to:
1. **Upload two unbound PDB structures** → Get predicted docked complex
2. **Verify against benchmark** → Compare with ground truth and get metrics
3. **Test any protein pair** → Works on any protein structures

---

## 🚀 USAGE EXAMPLES

### **Example 1: Test on Benchmark Target (RECOMMENDED)**

```bash
# Test on a benchmark target with full evaluation
python scripts/predict_docking.py \
    --target 1A2K \
    --difficulty rigid \
    --evaluate \
    --output-dir results/1A2K

# Output:
# - DockQ score
# - I-RMSD, L-RMSD, fnat
# - Predicted complex PDB
# - Comparison with AlphaRed
```

**Expected Output:**
```
================================================================================
PREDICTING AND EVALUATING PROTEIN-PROTEIN DOCKING
================================================================================

Input files:
  Receptor: data/alphared_benchmark/rigid_targets/1A2K/1A2K_r_u.pdb
  Ligand: data/alphared_benchmark/rigid_targets/1A2K/1A2K_l_u.pdb

Extracting features...
  Receptor: 115 residues
  Ligand: 56 residues

Predicting docking transformation...
  Confidence: 0.9941
  Rotation magnitude: 0.0011
  Translation: [-0.02, 0.01, -0.03] Å

================================================================================
EVALUATING AGAINST GROUND TRUTH
================================================================================

Ground truth: data/alphared_benchmark/rigid_targets/1A2K/1A2K_b.pdb
  Bound complex: 171 residues
  Ground truth ligand: 56 residues

Calculating metrics...

================================================================================
RESULTS
================================================================================

📊 Docking Quality Metrics:
  DockQ Score: 0.7234
  Quality: MEDIUM QUALITY ⭐⭐
  Success (DockQ > 0.23): ✓ YES

📏 Distance Metrics:
  I-RMSD (Interface): 1.1395 Å
  L-RMSD (Ligand): 1.1395 Å

🔗 Contact Metrics:
  fnat (Native contacts): 0.6234 (62.3%)

🎯 Prediction Info:
  Confidence: 0.9941

📈 Comparison:
  AlphaRed average success: 63%
  This prediction: SUCCESS ✓

✓ Results saved to: results/1A2K
  - metrics.json
  - predicted_complex.pdb
```

---

### **Example 2: Test Multiple Benchmark Targets**

```bash
# Test several targets at once
for target in 1A2K 1AHW 1AK4 1AVX; do
    python scripts/predict_docking.py \
        --target $target \
        --difficulty rigid \
        --evaluate \
        --output-dir results/$target
done

# Summarize results
python -c "
import json
from pathlib import Path

results = []
for result_dir in Path('results').glob('*/'):
    metrics_file = result_dir / 'metrics.json'
    if metrics_file.exists():
        with open(metrics_file) as f:
            metrics = json.load(f)
        results.append({
            'target': result_dir.name,
            'dockq': metrics['dockq'],
            'success': metrics['success']
        })

print('Target  | DockQ | Success')
print('--------|-------|--------')
for r in results:
    print(f'{r[\"target\"]:7} | {r[\"dockq\"]:.4f} | {\"✓\" if r[\"success\"] else \"✗\"}')

success_rate = sum(r['success'] for r in results) / len(results) * 100
print(f'\\nSuccess rate: {success_rate:.1f}%')
"
```

---

### **Example 3: Upload Your Own Structures**

```bash
# Predict docking for your own proteins (no ground truth)
python scripts/predict_docking.py \
    --receptor my_receptor.pdb \
    --ligand my_ligand.pdb \
    --output my_predicted_complex.pdb

# Output: predicted_complex.pdb
```

---

### **Example 4: Upload Your Own Structures WITH Ground Truth**

```bash
# If you have ground truth, get full evaluation
python scripts/predict_docking.py \
    --receptor my_receptor.pdb \
    --ligand my_ligand.pdb \
    --bound my_bound_complex.pdb \
    --evaluate \
    --output-dir my_results/

# Output:
# - my_results/metrics.json (all scores)
# - my_results/predicted_complex.pdb
```

---

## 📊 METRICS EXPLAINED

### **DockQ Score (Primary Metric)**
- **Range**: 0 to 1
- **> 0.8**: High quality (excellent prediction)
- **> 0.49**: Medium quality (good prediction)
- **> 0.23**: Acceptable (success!)
- **< 0.23**: Incorrect (failure)

**AlphaRed uses DockQ > 0.23 as success criterion**

### **I-RMSD (Interface RMSD)**
- Measures how well interface residues are positioned
- **< 1.0 Å**: Excellent
- **< 2.0 Å**: Good
- **< 4.0 Å**: Acceptable

### **L-RMSD (Ligand RMSD)**
- Measures overall ligand positioning
- Similar ranges to I-RMSD

### **fnat (Fraction of Native Contacts)**
- What fraction of real contacts are reproduced
- **Range**: 0 to 1 (higher is better)
- **> 0.5**: Good contact recovery

---

## 🧪 VERIFICATION PROTOCOL

### **Test 1: Single Benchmark Target**

```bash
# Pick a target from the benchmark
python scripts/predict_docking.py \
    --target 1A2K \
    --difficulty rigid \
    --evaluate \
    --output-dir verification/1A2K
```

**Expected:**
- DockQ > 0.23 (success)
- I-RMSD < 4.0 Å
- Results match training performance

---

### **Test 2: Multiple Targets (Representative Sample)**

Create a test script:

```bash
cat > test_verification.sh << 'EOF'
#!/bin/bash

# Test representative targets from each category
declare -A targets=(
    ["1A2K"]="rigid"
    ["1AHW"]="rigid"      # Antibody-antigen
    ["1AVX"]="medium"
    ["1BKD"]="difficult"
)

echo "Testing ${#targets[@]} targets..."

for target in "${!targets[@]}"; do
    difficulty="${targets[$target]}"
    echo ""
    echo "Testing $target ($difficulty)..."
    
    python scripts/predict_docking.py \
        --target $target \
        --difficulty $difficulty \
        --evaluate \
        --output-dir verification/$target
done

echo ""
echo "Verification complete! Check verification/ directory"
EOF

chmod +x test_verification.sh
./test_verification.sh
```

---

### **Test 3: Compare with Benchmark Results**

```bash
python << 'EOF'
import json
from pathlib import Path
import pandas as pd

# Load your predictions
verification_results = []
for result_dir in Path('verification').glob('*/'):
    metrics_file = result_dir / 'metrics.json'
    if metrics_file.exists():
        with open(metrics_file) as f:
            metrics = json.load(f)
        verification_results.append({
            'target': result_dir.name,
            'dockq': metrics['dockq'],
            'irmsd': metrics['irmsd'],
            'success': metrics['success']
        })

# Load training results
training_results = pd.read_csv('experiments/20251126_185619_full_training/evaluation/dockq_evaluation_results.csv')

# Compare
print("="*80)
print("VERIFICATION vs TRAINING RESULTS")
print("="*80)

for vr in verification_results:
    target = vr['target']
    tr = training_results[training_results['target'] == target]
    
    if len(tr) > 0:
        print(f"\nTarget: {target}")
        print(f"  Verification DockQ: {vr['dockq']:.4f}")
        print(f"  Training DockQ: {tr['dockq'].values[0]:.4f}")
        print(f"  Difference: {abs(vr['dockq'] - tr['dockq'].values[0]):.4f}")
        
        if abs(vr['dockq'] - tr['dockq'].values[0]) < 0.1:
            print(f"  Status: ✓ CONSISTENT")
        else:
            print(f"  Status: ⚠ CHECK")

print("\n" + "="*80)
print("If all targets show ✓ CONSISTENT, verification passed!")
print("="*80)
EOF
```

---

## 📁 OUTPUT FILES

When you run with `--evaluate`:

### **metrics.json**
```json
{
  "irmsd": 1.1395,
  "lrmsd": 1.1395,
  "fnat": 0.6234,
  "dockq": 0.7234,
  "quality": "medium",
  "success": true
}
```

### **predicted_complex.pdb**
- Standard PDB format
- Contains both receptor and ligand
- Ligand positioned according to prediction
- Can be visualized in PyMOL, Chimera, etc.

---

## 🔧 ADVANCED OPTIONS

### **Use CPU Instead of GPU**
```bash
python scripts/predict_docking.py \
    --target 1A2K \
    --difficulty rigid \
    --device cpu
```

### **Use Different Model**
```bash
python scripts/predict_docking.py \
    --target 1A2K \
    --difficulty rigid \
    --model path/to/different/model.pt
```

---

## 📊 BATCH PROCESSING

Create a list of targets:

```bash
# Create target list
cat > targets.txt << EOF
1A2K rigid
1AHW rigid
1AVX medium
1BKD difficult
EOF

# Process all
while read target difficulty; do
    python scripts/predict_docking.py \
        --target $target \
        --difficulty $difficulty \
        --evaluate \
        --output-dir batch_results/$target
done < targets.txt

# Summarize
python -c "
import json
from pathlib import Path

results = []
for result_dir in Path('batch_results').glob('*/'):
    metrics_file = result_dir / 'metrics.json'
    if metrics_file.exists():
        with open(metrics_file) as f:
            metrics = json.load(f)
        results.append({
            'target': result_dir.name,
            'dockq': metrics['dockq'],
            'irmsd': metrics['irmsd'],
            'success': metrics['success']
        })

# Summary statistics
successes = sum(r['success'] for r in results)
success_rate = successes / len(results) * 100

print('BATCH PROCESSING SUMMARY')
print('='*60)
print(f'Total targets: {len(results)}')
print(f'Successful: {successes}')
print(f'Success rate: {success_rate:.1f}%')
print(f'Mean DockQ: {sum(r[\"dockq\"] for r in results)/len(results):.4f}')
print(f'Mean I-RMSD: {sum(r[\"irmsd\"] for r in results)/len(results):.4f} Å')
"
```

---

## 🎯 FOR REVIEWERS

### **Reproducibility Test**

```bash
# Pick 5 random targets
python -c "
import random
targets = [
    ('1A2K', 'rigid'),
    ('1AHW', 'rigid'),
    ('1AVX', 'medium'),
    ('1BKD', 'difficult'),
    ('2PCC', 'rigid'),
]

for target, difficulty in targets:
    print(f'python scripts/predict_docking.py --target {target} --difficulty {difficulty} --evaluate --output-dir verification_advisor/{target}')
" | bash

# Results should match training performance
```

---

## 🐛 TROUBLESHOOTING

### **Error: "ModuleNotFoundError"**
```bash
pip install biopython torch numpy pandas
```

### **Error: "cuda out of memory"**
```bash
# Use CPU instead
python scripts/predict_docking.py --target 1A2K --difficulty rigid --device cpu
```

### **Error: "File not found"**
```bash
# Check if dataset exists
ls data/alphared_benchmark/rigid_targets/1A2K/

# If missing, check data directory structure
ls -R data/
```

---


