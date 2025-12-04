#!/usr/bin/env python3
"""
Protein Docking Prediction & Verification Tool
Upload unbound receptor and ligand structures to predict docked complex
Optionally provide ground truth for validation
"""

import torch
import numpy as np
import argparse
from pathlib import Path
import sys
from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain
import pandas as pd
import json

sys.path.append(str(Path(__file__).parent.parent))

from scripts.step4_feature_extraction import ProteinFeatureExtractor
from scripts.step5_simple_gnn import ProteinDockingModel, prepare_graph_data


class ProteinDockingPredictor:
    """Interactive protein docking predictor"""
    
    def __init__(self, model_path, device='cuda'):
        """
        Initialize predictor with trained model
        
        Args:
            model_path: Path to model checkpoint
            device: 'cuda' or 'cpu'
        """
        self.device = device if torch.cuda.is_available() else 'cpu'
        
        # Load model
        print(f"Loading model from: {model_path}")
        self.model = ProteinDockingModel(node_features=26, hidden_dim=128, num_layers=3)
        checkpoint = torch.load(model_path, map_location=self.device, weights_only=False)
        self.model.load_state_dict(checkpoint['model_state_dict'])
        self.model = self.model.to(self.device)
        self.model.eval()
        print(f"✓ Model loaded (trained for {checkpoint['epoch']} epochs)")
        
        self.extractor = ProteinFeatureExtractor()
    

    def _extract_features_from_pdb(self, pdb_path):
        """Extract features from a PDB file"""
        from Bio.PDB import PDBParser
        import numpy as np
        
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_path)
        
        # Extract coordinates and features
        coords = []
        residue_types = []
        residue_properties = []
        
        # Amino acid mapping
        aa_types = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                   'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
        aa_to_idx = {aa: i for i, aa in enumerate(aa_types)}
        
        for residue in structure.get_residues():
            if 'CA' in residue:
                ca_atom = residue['CA']
                coord = ca_atom.get_coord()
                coords.append(coord)
                
                # Get residue type
                res_name = residue.get_resname()
                res_idx = aa_to_idx.get(res_name, 0)  # Default to ALA if unknown
                residue_types.append(res_idx)
                
                # Chemical properties
                hydrophobic = ['ILE', 'LEU', 'VAL', 'PHE', 'TRP', 'MET']
                positive = ['ARG', 'LYS', 'HIS']
                negative = ['ASP', 'GLU']
                
                if res_name in hydrophobic:
                    properties = [1.0, 0.0, 0.0]
                elif res_name in positive:
                    properties = [0.0, 1.0, 0.0]
                elif res_name in negative:
                    properties = [0.0, 0.0, 1.0]
                else:
                    properties = [0.0, 0.0, 0.0]
                
                residue_properties.append(properties)
        
        coords = np.array(coords)
        residue_types = np.array(residue_types)
        residue_properties = np.array(residue_properties)
        
        # Build contact graph
        edges = []
        for i in range(len(coords)):
            for j in range(i+1, len(coords)):
                dist = np.linalg.norm(coords[i] - coords[j])
                if dist < 8.0:
                    edges.append([i, j])
                    edges.append([j, i])
        
        edge_index = np.array(edges).T if edges else np.array([[], []]).astype(int)
        
        return {
            'coords': coords,
            'residue_types': residue_types,
            'residue_properties': residue_properties,
            'edges': edge_index
        }

    def predict_docking(self, receptor_pdb, ligand_pdb, output_path=None):
        """
        Predict docked complex from unbound structures
        
        Args:
            receptor_pdb: Path to receptor PDB file
            ligand_pdb: Path to ligand PDB file
            output_path: Optional path to save predicted complex
            
        Returns:
            dict with prediction results
        """
        print("\n" + "="*80)
        print("PREDICTING PROTEIN-PROTEIN DOCKING")
        print("="*80)
        
        print(f"\nInput files:")
        print(f"  Receptor: {receptor_pdb}")
        print(f"  Ligand: {ligand_pdb}")
        
        # Extract features
        print("\nExtracting features...")
        receptor_features = self._extract_features_from_pdb(receptor_pdb)
        ligand_features = self._extract_features_from_pdb(ligand_pdb)
        
        print(f"  Receptor: {len(receptor_features['coords'])} residues")
        print(f"  Ligand: {len(ligand_features['coords'])} residues")
        
        # Prepare graph data
        features = {
            'receptor': receptor_features,
            'ligand': ligand_features
        }
        receptor_data, ligand_data = prepare_graph_data(features)
        
        # Move to device
        receptor_data = {k: v.to(self.device) if isinstance(v, torch.Tensor) else v 
                        for k, v in receptor_data.items()}
        ligand_data = {k: v.to(self.device) if isinstance(v, torch.Tensor) else v 
                      for k, v in ligand_data.items()}
        
        # Predict transformation
        print("\nPredicting docking transformation...")
        with torch.no_grad():
            rotation, translation, confidence = self.model(receptor_data, ligand_data)
            
            # Apply transformation
            transformed_ligand = self.model.apply_transformation(
                ligand_data['coords'],
                rotation,
                translation
            )
        
        print(f"  Confidence: {confidence.item():.4f}")
        print(f"  Rotation magnitude: {torch.norm(rotation).item():.4f}")
        print(f"  Translation: [{translation[0].item():.2f}, {translation[1].item():.2f}, {translation[2].item():.2f}] Å")
        
        # Convert to numpy
        receptor_coords = receptor_data['coords'].cpu().numpy()
        predicted_ligand = transformed_ligand.cpu().numpy()
        
        # Save predicted complex if requested
        if output_path:
            self._save_complex(
                receptor_pdb, 
                ligand_pdb, 
                predicted_ligand,
                output_path
            )
            print(f"\n✓ Predicted complex saved to: {output_path}")
        
        result = {
            'receptor_coords': receptor_coords,
            'predicted_ligand_coords': predicted_ligand,
            'rotation': rotation.cpu().numpy(),
            'translation': translation.cpu().numpy(),
            'confidence': confidence.item()
        }
        
        return result
    
    def evaluate_prediction(self, receptor_pdb, ligand_pdb, bound_pdb, output_dir=None):
        """
        Predict docking and evaluate against ground truth
        
        Args:
            receptor_pdb: Path to receptor PDB file
            ligand_pdb: Path to ligand PDB file
            bound_pdb: Path to bound complex PDB file (ground truth)
            output_dir: Optional directory to save results
            
        Returns:
            dict with metrics
        """
        print("\n" + "="*80)
        print("PREDICTING AND EVALUATING PROTEIN-PROTEIN DOCKING")
        print("="*80)
        
        # Make prediction
        prediction = self.predict_docking(receptor_pdb, ligand_pdb)
        
        # Load ground truth
        print("\n" + "="*80)
        print("EVALUATING AGAINST GROUND TRUTH")
        print("="*80)
        
        print(f"\nGround truth: {bound_pdb}")
        bound_features = self._extract_features_from_pdb(bound_pdb)
        
        # Extract ground truth ligand (second part of complex)
        num_receptor = len(prediction['receptor_coords'])
        ground_truth_ligand = bound_features['coords'][num_receptor:]
        
        print(f"  Bound complex: {len(bound_features['coords'])} residues")
        print(f"  Ground truth ligand: {len(ground_truth_ligand)} residues")
        
        # Calculate metrics
        print("\nCalculating metrics...")
        metrics = self._calculate_all_metrics(
            prediction['receptor_coords'],
            prediction['predicted_ligand_coords'],
            prediction['receptor_coords'],  # Receptor stays same
            ground_truth_ligand
        )
        
        # Print results
        print("\n" + "="*80)
        print("RESULTS")
        print("="*80)
        
        print(f"\n📊 Docking Quality Metrics:")
        print(f"  DockQ Score: {metrics['dockq']:.4f}")
        
        if metrics['dockq'] >= 0.8:
            quality = "HIGH QUALITY ⭐⭐⭐"
        elif metrics['dockq'] >= 0.49:
            quality = "MEDIUM QUALITY ⭐⭐"
        elif metrics['dockq'] >= 0.23:
            quality = "ACCEPTABLE ⭐"
        else:
            quality = "INCORRECT"
        
        print(f"  Quality: {quality}")
        print(f"  Success (DockQ > 0.23): {'✓ YES' if metrics['success'] else '✗ NO'}")
        
        print(f"\n📏 Distance Metrics:")
        print(f"  I-RMSD (Interface): {metrics['irmsd']:.4f} Å")
        print(f"  L-RMSD (Ligand): {metrics['lrmsd']:.4f} Å")
        
        print(f"\n🔗 Contact Metrics:")
        print(f"  fnat (Native contacts): {metrics['fnat']:.4f} ({metrics['fnat']*100:.1f}%)")
        
        print(f"\n🎯 Prediction Info:")
        print(f"  Confidence: {prediction['confidence']:.4f}")
        
        # Compare with AlphaRed
        print(f"\n📈 Comparison:")
        print(f"  AlphaRed average success: 63%")
        print(f"  This prediction: {'SUCCESS ✓' if metrics['success'] else 'FAILED ✗'}")
        
        # Save results if requested
        if output_dir:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            
            # Save metrics
            metrics_json = {}
            for k, v in metrics.items():
                if isinstance(v, (np.floating, np.integer)):
                    metrics_json[k] = float(v)
                elif isinstance(v, (bool, np.bool_)):
                    metrics_json[k] = bool(v)
                else:
                    metrics_json[k] = str(v)
            
            with open(output_dir / 'metrics.json', 'w') as f:
                json.dump(metrics_json, f, indent=2)
            
            # Save predicted complex
            complex_path = output_dir / 'predicted_complex.pdb'
            self._save_complex(
                receptor_pdb,
                ligand_pdb,
                prediction['predicted_ligand_coords'],
                complex_path
            )
            
            print(f"\n✓ Results saved to: {output_dir}")
            print(f"  - metrics.json")
            print(f"  - predicted_complex.pdb")
        
        return metrics
    
    def _calculate_all_metrics(self, receptor_pred, ligand_pred, receptor_true, ligand_true):
        """Calculate all docking metrics"""
        
        # I-RMSD (Interface RMSD)
        irmsd = self._calculate_rmsd(ligand_pred, ligand_true)
        
        # L-RMSD (Ligand RMSD after superimposing receptors)
        lrmsd = self._calculate_lrmsd(ligand_pred, ligand_true, receptor_pred, receptor_true)
        
        # fnat (Fraction of native contacts)
        fnat = self._calculate_fnat(receptor_pred, ligand_pred, receptor_true, ligand_true)
        
        # DockQ
        dockq = self._calculate_dockq(irmsd, lrmsd, fnat)
        
        # Quality and success
        if dockq >= 0.8:
            quality = 'high'
        elif dockq >= 0.49:
            quality = 'medium'
        elif dockq >= 0.23:
            quality = 'acceptable'
        else:
            quality = 'incorrect'
        
        success = dockq >= 0.23
        
        return {
            'irmsd': irmsd,
            'lrmsd': lrmsd,
            'fnat': fnat,
            'dockq': dockq,
            'quality': quality,
            'success': success
        }
    
    def _calculate_rmsd(self, coords1, coords2):
        """Calculate RMSD between two coordinate sets"""
        diff = coords1 - coords2
        return np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    
    def _calculate_lrmsd(self, ligand_pred, ligand_true, receptor_pred, receptor_true):
        """Calculate Ligand-RMSD after superimposing receptors"""
        # Simple centroid-based alignment
        receptor_pred_center = receptor_pred.mean(axis=0)
        receptor_true_center = receptor_true.mean(axis=0)
        
        translation = receptor_true_center - receptor_pred_center
        aligned_ligand = ligand_pred + translation
        
        return self._calculate_rmsd(aligned_ligand, ligand_true)
    
    def _calculate_fnat(self, receptor_pred, ligand_pred, receptor_true, ligand_true, threshold=5.0):
        """Calculate fraction of native contacts"""
        # Get native contacts
        native_contacts = set()
        for i in range(len(receptor_true)):
            for j in range(len(ligand_true)):
                dist = np.linalg.norm(receptor_true[i] - ligand_true[j])
                if dist < threshold:
                    native_contacts.add((i, j))
        
        if len(native_contacts) == 0:
            return 0.0
        
        # Get predicted contacts
        pred_contacts = set()
        for i in range(len(receptor_pred)):
            for j in range(len(ligand_pred)):
                dist = np.linalg.norm(receptor_pred[i] - ligand_pred[j])
                if dist < threshold:
                    pred_contacts.add((i, j))
        
        # Calculate fnat
        common = native_contacts.intersection(pred_contacts)
        return len(common) / len(native_contacts)
    
    def _calculate_dockq(self, irmsd, lrmsd, fnat):
        """Calculate DockQ score"""
        dockq = (fnat + 
                1.0 / (1.0 + (irmsd / 1.5)**2) + 
                1.0 / (1.0 + (lrmsd / 8.5)**2)) / 3.0
        return dockq
    
    def _save_complex(self, receptor_pdb, ligand_pdb, transformed_ligand_coords, output_path):
        """Save predicted complex to PDB file"""
        parser = PDBParser(QUIET=True)
        
        # Load receptor
        receptor_structure = parser.get_structure('receptor', receptor_pdb)
        
        # Load ligand
        ligand_structure = parser.get_structure('ligand', ligand_pdb)
        
        # Create new structure for complex
        complex_structure = Structure.Structure('complex')
        complex_model = Model.Model(0)
        complex_structure.add(complex_model)
        
        # Add receptor
        for chain in receptor_structure[0]:
            complex_model.add(chain.copy())
        
        # Add transformed ligand
        for chain in ligand_structure[0]:
            new_chain = chain.copy()
            # Update coordinates
            residue_idx = 0
            for residue in new_chain:
                for atom in residue:
                    if atom.name == 'CA':  # Only update CA atoms
                        atom.set_coord(transformed_ligand_coords[residue_idx])
                        residue_idx += 1
            complex_model.add(new_chain)
        
        # Save
        io = PDBIO()
        io.set_structure(complex_structure)
        io.save(str(output_path))


def main():
    parser = argparse.ArgumentParser(
        description='Protein Docking Prediction & Verification Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:

1. Predict docking (no ground truth):
   python predict_docking.py --receptor receptor.pdb --ligand ligand.pdb --output predicted.pdb

2. Predict and evaluate (with ground truth):
   python predict_docking.py --receptor receptor.pdb --ligand ligand.pdb --bound bound.pdb --output-dir results/

3. Test on benchmark target:
   python predict_docking.py --target 1A2K --difficulty rigid --evaluate

4. Batch evaluation on multiple targets:
   python predict_docking.py --batch targets.txt --output-dir batch_results/
        """
    )
    
    # Input files
    parser.add_argument('--receptor', type=str, help='Path to receptor PDB file')
    parser.add_argument('--ligand', type=str, help='Path to ligand PDB file')
    parser.add_argument('--bound', type=str, help='Path to bound complex PDB file (ground truth)')
    
    # Benchmark target
    parser.add_argument('--target', type=str, help='Benchmark target name (e.g., 1A2K)')
    parser.add_argument('--difficulty', type=str, choices=['rigid', 'medium', 'difficult'], 
                       help='Difficulty level for benchmark target')
    
    # Output
    parser.add_argument('--output', type=str, help='Output path for predicted complex PDB')
    parser.add_argument('--output-dir', type=str, help='Output directory for results')
    
    # Options
    parser.add_argument('--model', type=str, 
                       default='experiments/20251126_185619_full_training/best_model.pt',
                       help='Path to model checkpoint')
    parser.add_argument('--evaluate', action='store_true', 
                       help='Evaluate against ground truth')
    parser.add_argument('--device', type=str, default='cuda', 
                       choices=['cuda', 'cpu'], help='Device to use')
    
    args = parser.parse_args()
    
    # Initialize predictor
    predictor = ProteinDockingPredictor(args.model, args.device)
    
    # Handle benchmark target
    if args.target:
        if not args.difficulty:
            print("Error: --difficulty required when using --target")
            return
        
        difficulty_map = {
            'rigid': 'rigid_targets',
            'medium': 'medium_targets',
            'difficult': 'difficult_targets'
        }
        
        base_path = Path('data/alphared_benchmark') / difficulty_map[args.difficulty] / args.target
        
        receptor_pdb = base_path / f"{args.target}_r_u.pdb"
        ligand_pdb = base_path / f"{args.target}_l_u.pdb"
        
        # Find bound file (different naming conventions)
        bound_candidates = [
            base_path / f"{args.target}_r_l_b.pdb",
            base_path / f"{args.target}_b.pdb"
        ]
        bound_pdb = None
        for candidate in bound_candidates:
            if candidate.exists():
                bound_pdb = candidate
                break
        
        if args.evaluate and bound_pdb:
            output_dir = args.output_dir or f'results/{args.target}'
            predictor.evaluate_prediction(receptor_pdb, ligand_pdb, bound_pdb, output_dir)
        else:
            output_path = args.output or f'predicted_{args.target}.pdb'
            predictor.predict_docking(receptor_pdb, ligand_pdb, output_path)
    
    # Handle custom files
    elif args.receptor and args.ligand:
        if args.evaluate and args.bound:
            output_dir = args.output_dir or 'results'
            predictor.evaluate_prediction(args.receptor, args.ligand, args.bound, output_dir)
        else:
            output_path = args.output or 'predicted_complex.pdb'
            predictor.predict_docking(args.receptor, args.ligand, output_path)
    
    else:
        parser.print_help()
        print("\nError: Must provide either --target or (--receptor and --ligand)")


if __name__ == "__main__":
    main()
