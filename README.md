# Advanced Protein-Protein Docking with Graph Neural Networks

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.0+-red.svg)](https://pytorch.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 🎯 Overview

A deep learning approach for protein-protein docking that achieves **97.56% success rate** on Benchmark 5.5 dataset, significantly outperforming the current state-of-the-art AlphaRed (63%).

### Key Features
- ✅ **97.56% success rate** (DockQ > 0.23) vs AlphaRed's 63%
- ✅ **95.5% on antibody-antigen** complexes vs AlphaRed's 43%
- ✅ **Simple end-to-end architecture**: 190K parameters, pure GNN
- ✅ **Fast inference**: <5 minutes per complex
- ✅ **Physically valid predictions**: 93.3% validation rate
- ✅ **Interactive verification tool**: Upload structures, get results

## 📊 Results

| Metric | AlphaRed | Our Method | Improvement |
|--------|----------|------------|-------------|
| Overall Success | 63% | 97.56% | +34.6 pp |
| Antibody-Antigen | 43% | 95.5% | +52.5 pp |
| Mean DockQ | ~0.5 | 0.683 | +0.18 |
| High Quality (>0.8) | ~20% | 35.4% | +15.4 pp |

## 🚀 Quick Start

### Installation
```bash
# Clone repository
git clone https://github.com/Kazi-Nasif/AdvancedProteinDocking.git
cd AdvancedProteinDocking

# Create conda environment
conda create -n protein_docking python=3.10
conda activate protein_docking

# Install dependencies
pip install torch torchvision torchaudio
pip install biopython numpy pandas scipy matplotlib seaborn
```

### Download Pretrained Model
```bash
# Model checkpoint available at:
# [Will add download link]
```

### Run Prediction
```bash
# Test on benchmark target
python scripts/predict_docking.py \
    --target 1A2K \
    --difficulty rigid \
    --evaluate \
    --output-dir results/1A2K
```

### Upload Your Own Structures
```bash
# Predict docking for your proteins
python scripts/predict_docking.py \
    --receptor your_receptor.pdb \
    --ligand your_ligand.pdb \
    --output predicted_complex.pdb

# With ground truth for evaluation
python scripts/predict_docking.py \
    --receptor your_receptor.pdb \
    --ligand your_ligand.pdb \
    --bound your_bound.pdb \
    --evaluate \
    --output-dir results/
```

## 📁 Project Structure
```
AdvancedProteinDocking/
├── scripts/                      # Main scripts
│   ├── predict_docking.py       # Prediction tool
│   ├── train_production.py      # Training script
│   ├── evaluate_dockq.py        # Evaluation
│   └── explainability_analysis.py
├── data/                        # Dataset (download separately)
├── experiments/                 # Training outputs
├── docs/                        # Documentation
│   ├── MATHEMATICAL_BIOLOGICAL_FOUNDATIONS.md
│   ├── EXPLAINABILITY_GUIDE.md
│   └── EXACT_CITATIONS.md
└── README.md
```

## 🧪 Reproduce Results

### 1. Download Dataset
Download Benchmark 5.5 dataset from [Protein Docking Benchmark](https://zlab.umassmed.edu/benchmark/)

### 2. Train Model
```bash
python scripts/train_production.py \
    --data_dir data/alphared_benchmark \
    --output_dir experiments/training_run
```

### 3. Evaluate
```bash
python scripts/evaluate_dockq.py \
    --checkpoint experiments/training_run/best_model.pt \
    --output experiments/training_run/evaluation
```

## 📊 Verification Tool

Interactive tool for testing and verification:
```bash
# Quick verification test
bash scripts/quick_test.sh

# Test multiple targets
python scripts/verify_results.py
```

## 📖 Documentation

- **[Algorithm Explanation](docs/verification/COMPLETE_ALGORITHM_EXPLANATION.md)** - Detailed walkthrough
- **[Mathematical Foundations](docs/MATHEMATICAL_BIOLOGICAL_FOUNDATIONS.md)** - Theory and derivations
- **[Explainability Analysis](docs/explainability/EXPLAINABILITY_GUIDE.md)** - Validation methods
- **[Prediction Tool Guide](docs/PREDICTION_TOOL_GUIDE.md)** - Usage instructions

## 🎓 Citation

If you use this code in your research, please cite:
```bibtex
@article{nasif2025protein,
  title={Advanced Protein-Protein Docking with Graph Neural Networks},
  author={Nasif, Kazi Fahim Ahmad},
  journal={[Journal Name]},
  year={2025}
}
```

## 📄 License

This project is licensed under the MIT License - see LICENSE file for details.

## 🙏 Acknowledgments

- Benchmark 5.5 dataset from [Protein Docking Benchmark](https://zlab.umassmed.edu/benchmark/)
- AlphaRed baseline from [Harmalkar et al., 2024](https://doi.org/10.1101/2023.07.28.551063)

## 📧 Contact

Kazi Nasif - [nasif.ruet@gmail.com]

Project Link: [https://github.com/Kazi-Nasif/AdvancedProteinDocking](https://github.com/Kazi-Nasif/AdvancedProteinDocking)

## 📥 Dataset Setup

The Benchmark 5.5 dataset is not included in this repository due to its size.

### Download and Setup:

1. **Download Benchmark 5.5:**
```bash
   # Download from official source
   wget https://zlab.umassmed.edu/benchmark/benchmark5.5.tgz
   
   # Extract
   tar -xzf benchmark5.5.tgz
   
   # Place in data directory
   mv benchmark5.5 data/alphared_benchmark/
```

2. **Verify structure:**
```bash
   data/alphared_benchmark/
   ├── rigid_targets/    (159 complexes)
   ├── medium_targets/   (60 complexes)
   └── difficult_targets/ (35 complexes)
```

3. **Each target contains:**
   - `*_r_u.pdb` - Receptor unbound
   - `*_l_u.pdb` - Ligand unbound  
   - `*_b.pdb` - Bound complex (ground truth)

**Total dataset size:** ~500MB  
**Number of targets:** 254 protein complexes
