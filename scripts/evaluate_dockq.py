#!/usr/bin/env python3
"""
Calculate Complete CAPRI Metrics: DockQ, fnat, L-RMSD, I-RMSD
For direct comparison with AlphaRed and other methods
"""

import torch
import numpy as np
import pandas as pd
import json
from pathlib import Path
import sys
from tqdm import tqdm
from Bio import PDB

sys.path.append(str(Path(__file__).parent.parent))

from scripts.step3_data_loader import ProteinDockingLoader
from scripts.step4_feature_extraction import ProteinFeatureExtractor
from scripts.step5_simple_gnn import ProteinDockingModel, prepare_graph_data


def calculate_fnat(pred_coords_receptor, pred_coords_ligand, 
                   native_coords_receptor, native_coords_ligand, 
                   threshold=5.0):
    """
    Calculate fraction of native contacts (fnat)
    
    Args:
        pred_coords_receptor: Predicted receptor coordinates [N, 3]
        pred_coords_ligand: Predicted ligand coordinates [M, 3]
        native_coords_receptor: Native receptor coordinates [N, 3]
        native_coords_ligand: Native ligand coordinates [M, 3]
        threshold: Contact distance threshold (default 5.0 Å)
    
    Returns:
        fnat: Fraction of native contacts (0 to 1)
    """
    # Get native contacts
    native_contacts = set()
    for i in range(len(native_coords_receptor)):
        for j in range(len(native_coords_ligand)):
            dist = np.linalg.norm(native_coords_receptor[i] - native_coords_ligand[j])
            if dist < threshold:
                native_contacts.add((i, j))
    
    if len(native_contacts) == 0:
        return 0.0
    
    # Get predicted contacts
    pred_contacts = set()
    for i in range(len(pred_coords_receptor)):
        for j in range(len(pred_coords_ligand)):
            dist = np.linalg.norm(pred_coords_receptor[i] - pred_coords_ligand[j])
            if dist < threshold:
                pred_contacts.add((i, j))
    
    # Calculate fnat
    common_contacts = native_contacts.intersection(pred_contacts)
    fnat = len(common_contacts) / len(native_contacts)
    
    return fnat


def calculate_lrmsd(pred_ligand, native_ligand, native_receptor, pred_receptor):
    """
    Calculate Ligand-RMSD after superimposing receptors
    
    Args:
        pred_ligand: Predicted ligand coordinates [M, 3]
        native_ligand: Native ligand coordinates [M, 3]
        native_receptor: Native receptor coordinates [N, 3]
        pred_receptor: Predicted receptor coordinates [N, 3]
    
    Returns:
        L-RMSD in Angstroms
    """
    # Superimpose receptors (align predicted to native)
    # Simple centroid-based alignment
    native_rec_center = native_receptor.mean(axis=0)
    pred_rec_center = pred_receptor.mean(axis=0)
    
    # Translate predicted to align with native
    translation = native_rec_center - pred_rec_center
    aligned_pred_ligand = pred_ligand + translation
    
    # Calculate RMSD of ligand
    diff = aligned_pred_ligand - native_ligand
    lrmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    
    return lrmsd


def calculate_dockq(irmsd, lrmsd, fnat):
    """
    Calculate DockQ score
    
    DockQ formula from Basu & Wallner (2016):
    DockQ = (fnat + 1/(1+(irmsd/1.5)^2) + 1/(1+(lrmsd/8.5)^2)) / 3
    
    Args:
        irmsd: Interface-RMSD
        lrmsd: Ligand-RMSD
        fnat: Fraction of native contacts
    
    Returns:
        DockQ score (0 to 1)
    """
    dockq = (fnat + 
             1.0 / (1.0 + (irmsd / 1.5)**2) + 
             1.0 / (1.0 + (lrmsd / 8.5)**2)) / 3.0
    
    return dockq


def evaluate_with_full_metrics(model, features, device='cuda'):
    """Evaluate with complete CAPRI metrics"""
    try:
        # Prepare data
        receptor_data, ligand_data = prepare_graph_data(features)
        
        # Move to device
        receptor_data = {k: v.to(device) if isinstance(v, torch.Tensor) else v 
                        for k, v in receptor_data.items()}
        ligand_data = {k: v.to(device) if isinstance(v, torch.Tensor) else v 
                      for k, v in ligand_data.items()}
        
        # Get native coordinates
        native_receptor_coords = features['receptor']['coords']
        native_ligand_coords = features['ligand']['coords']
        
        # Get bound complex ligand coordinates (ground truth)
        bound_ligand_coords = features['bound']['coords'][len(native_receptor_coords):]
        
        # Predict
        with torch.no_grad():
            rotation, translation, confidence = model(receptor_data, ligand_data)
            pred_ligand_coords = model.apply_transformation(
                ligand_data['coords'],
                rotation,
                translation
            )
        
        # Convert to numpy
        pred_receptor_np = receptor_data['coords'].cpu().numpy()
        pred_ligand_np = pred_ligand_coords.cpu().numpy()
        bound_ligand_np = bound_ligand_coords
        
        # Calculate I-RMSD (interface residues only - simplified: use all)
        irmsd = np.sqrt(np.mean(np.sum((pred_ligand_np - bound_ligand_np)**2, axis=1)))
        
        # Calculate L-RMSD
        lrmsd = calculate_lrmsd(pred_ligand_np, bound_ligand_np, 
                               native_receptor_coords, pred_receptor_np)
        
        # Calculate fnat
        fnat = calculate_fnat(pred_receptor_np, pred_ligand_np,
                             native_receptor_coords, bound_ligand_np,
                             threshold=5.0)
        
        # Calculate DockQ
        dockq = calculate_dockq(irmsd, lrmsd, fnat)
        
        # CAPRI quality classification
        if dockq >= 0.8:
            quality = 'high'
        elif dockq >= 0.49:
            quality = 'medium'
        elif dockq >= 0.23:
            quality = 'acceptable'
        else:
            quality = 'incorrect'
        
        return {
            'irmsd': irmsd,
            'lrmsd': lrmsd,
            'fnat': fnat,
            'dockq': dockq,
            'quality': quality,
            'success': dockq >= 0.23,  # AlphaRed criterion
            'confidence': confidence.item(),
            'error': None
        }
        
    except Exception as e:
        return {
            'irmsd': None,
            'lrmsd': None,
            'fnat': None,
            'dockq': None,
            'quality': None,
            'success': False,
            'confidence': None,
            'error': str(e)
        }


def evaluate_all_with_dockq(model, device='cuda', max_targets=None):
    """Evaluate all targets with complete CAPRI metrics"""
    print("\n" + "="*80)
    print("COMPREHENSIVE CAPRI EVALUATION (DockQ, fnat, I-RMSD, L-RMSD)")
    print("="*80)
    
    loader = ProteinDockingLoader()
    extractor = ProteinFeatureExtractor()
    
    difficulties = ["rigid_targets", "medium_targets", "difficult_targets"]
    all_results = []
    
    total_processed = 0
    
    for difficulty in difficulties:
        print(f"\n{difficulty.upper()}")
        print("-" * 80)
        
        benchmark_path = Path(loader.benchmark_path) / difficulty
        target_dirs = sorted([d for d in benchmark_path.iterdir() if d.is_dir()])
        
        if max_targets and total_processed >= max_targets:
            break
        
        for target_dir in tqdm(target_dirs, desc=f"  Evaluating"):
            if max_targets and total_processed >= max_targets:
                break
                
            target_name = target_dir.name
            
            try:
                features = extractor.extract_full_features(target_name, difficulty)
                result = evaluate_with_full_metrics(model, features, device)
                
                all_results.append({
                    'target': target_name,
                    'difficulty': difficulty,
                    **result
                })
                
                total_processed += 1
                
            except Exception as e:
                all_results.append({
                    'target': target_name,
                    'difficulty': difficulty,
                    'irmsd': None,
                    'lrmsd': None,
                    'fnat': None,
                    'dockq': None,
                    'quality': None,
                    'success': False,
                    'confidence': None,
                    'error': str(e)
                })
    
    return pd.DataFrame(all_results)


def analyze_dockq_results(df):
    """Analyze results with DockQ metrics"""
    print("\n" + "="*80)
    print("DOCKQ-BASED ANALYSIS (AlphaRed Standard)")
    print("="*80)
    
    valid = df[df['error'].isna()]
    total = len(valid)
    
    # Overall statistics
    print(f"\n{'='*80}")
    print("OVERALL PERFORMANCE")
    print(f"{'='*80}")
    print(f"Total targets: {total}")
    print(f"Success (DockQ > 0.23): {len(valid[valid['success']==True])} ({len(valid[valid['success']==True])/total*100:.1f}%)")
    print(f"High quality (DockQ > 0.8): {len(valid[valid['quality']=='high'])} ({len(valid[valid['quality']=='high'])/total*100:.1f}%)")
    print(f"Medium quality (DockQ > 0.49): {len(valid[valid['quality']=='medium'])} ({len(valid[valid['quality']=='medium'])/total*100:.1f}%)")
    print(f"Acceptable (DockQ > 0.23): {len(valid[valid['quality']=='acceptable'])} ({len(valid[valid['quality']=='acceptable'])/total*100:.1f}%)")
    
    print(f"\nMetric Averages:")
    print(f"  Mean DockQ: {valid['dockq'].mean():.4f}")
    print(f"  Mean I-RMSD: {valid['irmsd'].mean():.4f} Å")
    print(f"  Mean L-RMSD: {valid['lrmsd'].mean():.4f} Å")
    print(f"  Mean fnat: {valid['fnat'].mean():.4f}")
    
    # By difficulty
    print(f"\n{'='*80}")
    print("BY DIFFICULTY (DockQ > 0.23)")
    print(f"{'='*80}")
    
    for difficulty in ["rigid_targets", "medium_targets", "difficult_targets"]:
        subset = valid[valid['difficulty'] == difficulty]
        if len(subset) > 0:
            success = len(subset[subset['success']==True])
            print(f"\n{difficulty}:")
            print(f"  Success: {success}/{len(subset)} ({success/len(subset)*100:.1f}%)")
            print(f"  Mean DockQ: {subset['dockq'].mean():.4f}")
    
    # Comparison
    print(f"\n{'='*80}")
    print("COMPARISON WITH ALPHARED")
    print(f"{'='*80}")
    success_rate = len(valid[valid['success']==True])/total*100
    print(f"AlphaRed: 63% (DockQ > 0.23)")
    print(f"Your Model: {success_rate:.1f}% (DockQ > 0.23)")
    print(f"Improvement: +{success_rate - 63:.1f} pp")
    
    if success_rate > 63:
        print(f"\n🎉 YOU BEAT ALPHARED! 🎉")


def main():
    DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'
    CHECKPOINT_PATH = 'experiments/20251126_185619_full_training/best_model.pt'
    OUTPUT_DIR = 'experiments/20251126_185619_full_training/evaluation'
    
    print("="*80)
    print("FULL CAPRI METRICS EVALUATION")
    print("="*80)
    
    # Load model
    print(f"\nLoading model from: {CHECKPOINT_PATH}")
    model = ProteinDockingModel(node_features=26, hidden_dim=128, num_layers=3)
    checkpoint = torch.load(CHECKPOINT_PATH, map_location=DEVICE)
    model.load_state_dict(checkpoint['model_state_dict'])
    model = model.to(DEVICE)
    model.eval()
    print(f"✓ Model loaded from epoch {checkpoint['epoch']}")
    
    # Evaluate (use max_targets for testing, None for full evaluation)
    results_df = evaluate_all_with_dockq(model, DEVICE, max_targets=None)
    
    # Analyze
    analyze_dockq_results(results_df)
    
    # Save
    output_dir = Path(OUTPUT_DIR)
    results_df.to_csv(output_dir / 'dockq_evaluation_results.csv', index=False)
    
    # Save summary
    valid = results_df[results_df['error'].isna()]
    summary = {
        'total_targets': len(valid),
        'dockq_success_rate': len(valid[valid['success']==True]) / len(valid) * 100,
        'mean_dockq': float(valid['dockq'].mean()),
        'mean_irmsd': float(valid['irmsd'].mean()),
        'mean_lrmsd': float(valid['lrmsd'].mean()),
        'mean_fnat': float(valid['fnat'].mean()),
        'high_quality': len(valid[valid['quality']=='high']),
        'medium_quality': len(valid[valid['quality']=='medium']),
        'acceptable_quality': len(valid[valid['quality']=='acceptable'])
    }
    
    with open(output_dir / 'dockq_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\n✓ Results saved to: {output_dir}")
    print("  - dockq_evaluation_results.csv")
    print("  - dockq_summary.json")


if __name__ == "__main__":
    main()
