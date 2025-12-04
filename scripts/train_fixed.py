#!/usr/bin/env python3 -u
"""Training with unbuffered output"""
import sys
import os

# Force unbuffered output
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 1)
sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', 1)

print("="*80, flush=True)
print("STARTING TRAINING SCRIPT", flush=True)
print("="*80, flush=True)

# Rest of imports
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import json
import time
from pathlib import Path
from datetime import datetime

sys.path.append(str(Path(__file__).parent.parent))

from scripts.step3_data_loader import ProteinDockingLoader
from scripts.step4_feature_extraction import ProteinFeatureExtractor
from scripts.step5_simple_gnn import ProteinDockingModel, prepare_graph_data

print("Imports complete!", flush=True)

# Your DockingTrainer class and other code here...
# (I'll create a simpler version)

def prepare_dataset_simple(difficulty="rigid_targets", max_samples=10):
    """Prepare small dataset quickly"""
    print(f"\nPreparing {max_samples} samples from {difficulty}...", flush=True)
    
    loader = ProteinDockingLoader()
    extractor = ProteinFeatureExtractor()
    
    benchmark_path = Path(loader.benchmark_path) / difficulty
    target_dirs = sorted([d for d in benchmark_path.iterdir() if d.is_dir()])[:max_samples]
    
    dataset = []
    
    for i, target_dir in enumerate(target_dirs):
        target_name = target_dir.name
        print(f"  [{i+1}/{len(target_dirs)}] {target_name}...", end=' ', flush=True)
        
        try:
            features = extractor.extract_full_features(target_name, difficulty)
            receptor_data, ligand_data = prepare_graph_data(features)
            
            target_coords = torch.tensor(
                features['bound']['coords'][len(features['receptor']['coords']):],
                dtype=torch.float32
            )
            
            # Get interface
            interface = loader.get_interface_residues(
                features['raw_data']['bound_complex'],
                features['partners']
            )
            
            ligand_interface_indices = list(range(len(features['ligand']['coords'])))
            
            data_dict = {
                'target_name': target_name,
                'receptor_data': receptor_data,
                'ligand_data': ligand_data,
                'target_coords': target_coords,
                'interface_indices': ligand_interface_indices
            }
            
            dataset.append(data_dict)
            print("✓", flush=True)
            
        except Exception as e:
            print(f"✗ {e}", flush=True)
    
    print(f"\nDataset ready: {len(dataset)} samples\n", flush=True)
    return dataset


def main():
    print("\n" + "="*80, flush=True)
    print("QUICK TRAINING TEST (10 samples, 5 epochs)", flush=True)
    print("="*80, flush=True)
    
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"Device: {device}\n", flush=True)
    
    # Small dataset for testing
    dataset = prepare_dataset_simple("rigid_targets", max_samples=10)
    
    train_data = dataset[:8]
    val_data = dataset[8:]
    
    print(f"Train: {len(train_data)}, Val: {len(val_data)}\n", flush=True)
    
    # Create model
    print("Creating model...", flush=True)
    model = ProteinDockingModel(node_features=26, hidden_dim=64, num_layers=2)
    model = model.to(device)
    optimizer = optim.Adam(model.parameters(), lr=1e-3)
    
    print(f"Parameters: {sum(p.numel() for p in model.parameters()):,}\n", flush=True)
    
    # Quick training loop
    print("="*80, flush=True)
    print("TRAINING", flush=True)
    print("="*80, flush=True)
    
    for epoch in range(1, 6):  # Just 5 epochs
        print(f"\nEpoch {epoch}/5:", flush=True)
        
        model.train()
        epoch_loss = 0
        
        for i, data_dict in enumerate(train_data):
            try:
                # Move to device
                receptor_data = {
                    'node_features': data_dict['receptor_data']['node_features'].to(device),
                    'edge_index': data_dict['receptor_data']['edge_index'].to(device),
                    'coords': data_dict['receptor_data']['coords'].to(device)
                }
                
                ligand_data = {
                    'node_features': data_dict['ligand_data']['node_features'].to(device),
                    'edge_index': data_dict['ligand_data']['edge_index'].to(device),
                    'coords': data_dict['ligand_data']['coords'].to(device)
                }
                
                target_coords = data_dict['target_coords'].to(device)
                
                # Forward
                rotation, translation, confidence = model(receptor_data, ligand_data)
                pred_coords = model.apply_transformation(ligand_data['coords'], rotation, translation)
                
                # Loss (simple RMSD)
                loss = torch.sqrt(torch.mean((pred_coords - target_coords)**2))
                
                # Backward
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()
                
                epoch_loss += loss.item()
                
                if (i + 1) % 2 == 0:
                    print(f"  Sample {i+1}/{len(train_data)}, Loss: {loss.item():.4f}", flush=True)
                    
            except Exception as e:
                print(f"  Error: {e}", flush=True)
        
        avg_loss = epoch_loss / len(train_data)
        print(f"  Epoch {epoch} complete! Avg Loss: {avg_loss:.4f}", flush=True)
    
    print("\n" + "="*80, flush=True)
    print("TRAINING COMPLETE!", flush=True)
    print("="*80, flush=True)
    print("This was a quick test. For full training, increase samples and epochs.", flush=True)


if __name__ == "__main__":
    main()
