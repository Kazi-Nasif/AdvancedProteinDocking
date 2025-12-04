#!/usr/bin/env python3
"""
Comprehensive Evaluation Script
Test trained model on all 254 targets and compare with AlphaRed baseline
"""

import torch
import numpy as np
import pandas as pd
import json
from pathlib import Path
import sys
from tqdm import tqdm

sys.path.append(str(Path(__file__).parent.parent))

from scripts.step3_data_loader import ProteinDockingLoader
from scripts.step4_feature_extraction import ProteinFeatureExtractor
from scripts.step5_simple_gnn import ProteinDockingModel, prepare_graph_data


def load_trained_model(checkpoint_path, device='cuda'):
    """Load the best trained model"""
    print(f"Loading model from: {checkpoint_path}")
    
    model = ProteinDockingModel(
        node_features=26,
        hidden_dim=128,
        num_layers=3
    )
    
    checkpoint = torch.load(checkpoint_path, map_location=device)
    model.load_state_dict(checkpoint['model_state_dict'])
    model = model.to(device)
    model.eval()
    
    print(f"✓ Model loaded from epoch {checkpoint['epoch']}")
    print(f"  Training val loss: {checkpoint['val_loss']:.4f} Å")
    
    return model


def evaluate_single_target(model, features, device='cuda'):
    """Evaluate model on a single target"""
    try:
        # Prepare data
        receptor_data, ligand_data = prepare_graph_data(features)
        
        # Move to device
        receptor_data = {
            'node_features': receptor_data['node_features'].to(device),
            'edge_index': receptor_data['edge_index'].to(device),
            'coords': receptor_data['coords'].to(device)
        }
        
        ligand_data = {
            'node_features': ligand_data['node_features'].to(device),
            'edge_index': ligand_data['edge_index'].to(device),
            'coords': ligand_data['coords'].to(device)
        }
        
        # Get target coordinates
        target_coords = torch.tensor(
            features['bound']['coords'][len(features['receptor']['coords']):],
            dtype=torch.float32
        ).to(device)
        
        # Predict
        with torch.no_grad():
            rotation, translation, confidence = model(receptor_data, ligand_data)
            pred_coords = model.apply_transformation(
                ligand_data['coords'],
                rotation,
                translation
            )
        
        # Calculate RMSD
        diff = pred_coords - target_coords
        rmsd = torch.sqrt(torch.mean(torch.sum(diff**2, dim=1))).item()
        
        # Success criteria
        success = rmsd < 4.0  # CAPRI acceptable quality
        high_quality = rmsd < 1.0  # CAPRI high quality
        medium_quality = rmsd < 2.0  # CAPRI medium quality
        
        return {
            'rmsd': rmsd,
            'confidence': confidence.item(),
            'success': success,
            'high_quality': high_quality,
            'medium_quality': medium_quality,
            'error': None
        }
        
    except Exception as e:
        return {
            'rmsd': None,
            'confidence': None,
            'success': False,
            'high_quality': False,
            'medium_quality': False,
            'error': str(e)
        }


def evaluate_all_targets(model, device='cuda'):
    """Evaluate on all 254 targets"""
    print("\n" + "="*80)
    print("EVALUATING ON ALL 254 TARGETS")
    print("="*80)
    
    loader = ProteinDockingLoader()
    extractor = ProteinFeatureExtractor()
    
    difficulties = ["rigid_targets", "medium_targets", "difficult_targets"]
    
    all_results = []
    
    for difficulty in difficulties:
        print(f"\n{difficulty.upper()}")
        print("-" * 80)
        
        benchmark_path = Path(loader.benchmark_path) / difficulty
        target_dirs = sorted([d for d in benchmark_path.iterdir() if d.is_dir()])
        
        for target_dir in tqdm(target_dirs, desc=f"  Evaluating"):
            target_name = target_dir.name
            
            try:
                # Extract features
                features = extractor.extract_full_features(target_name, difficulty)
                
                # Evaluate
                result = evaluate_single_target(model, features, device)
                
                # Store results
                all_results.append({
                    'target': target_name,
                    'difficulty': difficulty,
                    'rmsd': result['rmsd'],
                    'confidence': result['confidence'],
                    'success': result['success'],
                    'high_quality': result['high_quality'],
                    'medium_quality': result['medium_quality'],
                    'error': result['error']
                })
                
            except Exception as e:
                all_results.append({
                    'target': target_name,
                    'difficulty': difficulty,
                    'rmsd': None,
                    'confidence': None,
                    'success': False,
                    'high_quality': False,
                    'medium_quality': False,
                    'error': str(e)
                })
    
    return pd.DataFrame(all_results)


def analyze_results(df):
    """Analyze and print results"""
    print("\n" + "="*80)
    print("RESULTS ANALYSIS")
    print("="*80)
    
    # Overall statistics
    successful = df[df['success'] == True]
    total = len(df[df['error'].isna()])
    
    print(f"\n{'='*80}")
    print("OVERALL PERFORMANCE")
    print(f"{'='*80}")
    print(f"Total targets evaluated: {total}")
    print(f"Successful (RMSD < 4.0 Å): {len(successful)} ({len(successful)/total*100:.1f}%)")
    print(f"High quality (RMSD < 1.0 Å): {len(df[df['high_quality']==True])} ({len(df[df['high_quality']==True])/total*100:.1f}%)")
    print(f"Medium quality (RMSD < 2.0 Å): {len(df[df['medium_quality']==True])} ({len(df[df['medium_quality']==True])/total*100:.1f}%)")
    print(f"Failed: {len(df[df['error'].notna()])} ({len(df[df['error'].notna()])/len(df)*100:.1f}%)")
    
    # RMSD statistics
    valid_rmsd = df[df['rmsd'].notna()]['rmsd']
    print(f"\nRMSD Statistics:")
    print(f"  Mean: {valid_rmsd.mean():.4f} Å")
    print(f"  Median: {valid_rmsd.median():.4f} Å")
    print(f"  Std: {valid_rmsd.std():.4f} Å")
    print(f"  Min: {valid_rmsd.min():.4f} Å")
    print(f"  Max: {valid_rmsd.max():.4f} Å")
    
    # By difficulty
    print(f"\n{'='*80}")
    print("PERFORMANCE BY DIFFICULTY")
    print(f"{'='*80}")
    
    for difficulty in ["rigid_targets", "medium_targets", "difficult_targets"]:
        subset = df[df['difficulty'] == difficulty]
        subset_valid = subset[subset['error'].isna()]
        subset_success = subset[subset['success'] == True]
        
        if len(subset_valid) > 0:
            print(f"\n{difficulty}:")
            print(f"  Total: {len(subset_valid)}")
            print(f"  Success: {len(subset_success)} ({len(subset_success)/len(subset_valid)*100:.1f}%)")
            print(f"  Mean RMSD: {subset_valid['rmsd'].mean():.4f} Å")
            print(f"  High quality: {len(subset[subset['high_quality']==True])} ({len(subset[subset['high_quality']==True])/len(subset_valid)*100:.1f}%)")
    
    # Comparison with AlphaRed
    print(f"\n{'='*80}")
    print("COMPARISON WITH ALPHARED")
    print(f"{'='*80}")
    print(f"\nAlphaRed (baseline):")
    print(f"  Overall success rate: 63%")
    print(f"  Antibody-antigen: 43%")
    print(f"\nYour Model:")
    print(f"  Overall success rate: {len(successful)/total*100:.1f}%")
    print(f"  Improvement: +{len(successful)/total*100 - 63:.1f} percentage points")
    
    if len(successful)/total*100 > 63:
        print(f"\n🎉 YOU BEAT ALPHARED! 🎉")
    
    return successful, total


def save_results(df, output_dir):
    """Save results to files"""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save full results
    df.to_csv(output_dir / 'evaluation_results.csv', index=False)
    
    # Save summary
    summary = {
        'total_targets': len(df[df['error'].isna()]),
        'successful': len(df[df['success'] == True]),
        'success_rate': len(df[df['success'] == True]) / len(df[df['error'].isna()]) * 100,
        'high_quality': len(df[df['high_quality'] == True]),
        'medium_quality': len(df[df['medium_quality'] == True]),
        'mean_rmsd': float(df[df['rmsd'].notna()]['rmsd'].mean()),
        'median_rmsd': float(df[df['rmsd'].notna()]['rmsd'].median()),
        'by_difficulty': {}
    }
    
    for difficulty in ["rigid_targets", "medium_targets", "difficult_targets"]:
        subset = df[df['difficulty'] == difficulty]
        subset_valid = subset[subset['error'].isna()]
        summary['by_difficulty'][difficulty] = {
            'total': len(subset_valid),
            'successful': len(subset[subset['success'] == True]),
            'success_rate': len(subset[subset['success'] == True]) / len(subset_valid) * 100 if len(subset_valid) > 0 else 0,
            'mean_rmsd': float(subset_valid['rmsd'].mean()) if len(subset_valid) > 0 else None
        }
    
    with open(output_dir / 'evaluation_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\n✓ Results saved to: {output_dir}")
    print(f"  - evaluation_results.csv")
    print(f"  - evaluation_summary.json")


def main():
    # Configuration
    DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'
    CHECKPOINT_PATH = 'experiments/20251126_185619_full_training/best_model.pt'
    OUTPUT_DIR = 'experiments/20251126_185619_full_training/evaluation'
    
    print("="*80)
    print("COMPREHENSIVE MODEL EVALUATION")
    print("="*80)
    print(f"Device: {DEVICE}")
    print(f"Checkpoint: {CHECKPOINT_PATH}")
    
    # Load model
    model = load_trained_model(CHECKPOINT_PATH, DEVICE)
    
    # Evaluate on all targets
    results_df = evaluate_all_targets(model, DEVICE)
    
    # Analyze results
    successful, total = analyze_results(results_df)
    
    # Save results
    save_results(results_df, OUTPUT_DIR)
    
    print("\n" + "="*80)
    print("EVALUATION COMPLETE!")
    print("="*80)


if __name__ == "__main__":
    main()
