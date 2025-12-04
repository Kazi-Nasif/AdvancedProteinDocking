#!/bin/bash
# Quick Start Script for University Server
# Run this after copying the project to your server

echo "=================================================="
echo "Advanced Protein-Protein Docking Setup"
echo "=================================================="

# Check GPU availability
echo ""
echo "Checking GPU availability..."
nvidia-smi --query-gpu=index,name,memory.total,memory.free --format=csv

# Set up Python environment
echo ""
echo "Setting up Python environment..."
echo "Current Python version:"
python3 --version

# Install dependencies
echo ""
echo "Installing dependencies..."
echo "Run: pip install -r requirements.txt"
echo "Or with conda:"
echo "  conda create -n protein_docking python=3.10"
echo "  conda activate protein_docking"
echo "  pip install -r requirements.txt"

# Key packages to install
echo ""
echo "Key dependencies to install:"
echo "  - PyTorch (with CUDA support for A100)"
echo "  - BioPython"
echo "  - PyTorch Geometric"
echo "  - Other scientific packages (see requirements.txt)"

echo ""
echo "=================================================="
echo "Next Steps:"
echo "=================================================="
echo "1. Copy project to server: scp -r AdvancedProteinDocking user@server:/path/"
echo "2. Install dependencies: pip install -r requirements.txt"
echo "3. Organize your PDB files in data/raw/"
echo "4. Run analysis: python scripts/analyze_benchmark.py"
echo "5. Share PDB directory structure for Step 2"
echo "=================================================="
