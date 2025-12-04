#!/usr/bin/env python3
"""
Step 4: Feature Extraction for Deep Learning
Extract geometric and chemical features from proteins
"""

import numpy as np
from Bio import PDB
from pathlib import Path
import sys

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))

from scripts.step3_data_loader import ProteinDockingLoader


# Amino acid properties
AA_PROPERTIES = {
    'ALA': [0, 0, 0],  # [hydrophobic, positive, negative]
    'ARG': [0, 1, 0],
    'ASN': [0, 0, 0],
    'ASP': [0, 0, 1],
    'CYS': [0, 0, 0],
    'GLN': [0, 0, 0],
    'GLU': [0, 0, 1],
    'GLY': [0, 0, 0],
    'HIS': [0, 1, 0],
    'ILE': [1, 0, 0],
    'LEU': [1, 0, 0],
    'LYS': [0, 1, 0],
    'MET': [1, 0, 0],
    'PHE': [1, 0, 0],
    'PRO': [0, 0, 0],
    'SER': [0, 0, 0],
    'THR': [0, 0, 0],
    'TRP': [1, 0, 0],
    'TYR': [1, 0, 0],
    'VAL': [1, 0, 0],
}

# One-hot encoding for residue types
AA_TO_INDEX = {aa: i for i, aa in enumerate(sorted(AA_PROPERTIES.keys()))}


class ProteinFeatureExtractor:
    """Extract features from protein structures for deep learning"""
    
    def __init__(self):
        self.loader = ProteinDockingLoader()
    
    def extract_residue_features(self, structure):
        """
        Extract per-residue features
        
        Returns:
            dict with:
                - coords: Nx3 array (CA coordinates)
                - residue_types: N array (residue type indices)
                - residue_properties: Nx3 array (hydrophobic, positive, negative)
                - chain_ids: N array (chain identifiers)
        """
        coords = []
        residue_types = []
        residue_properties = []
        chain_ids = []
        
        model = list(structure.get_models())[0]
        
        for chain in model:
            for residue in chain:
                # Skip non-standard residues
                if residue.id[0] != ' ':
                    continue
                
                # Get CA atom
                if 'CA' not in residue:
                    continue
                
                ca_atom = residue['CA']
                res_name = residue.resname
                
                # Get features
                coords.append(ca_atom.get_coord())
                
                # Residue type (one-hot index)
                if res_name in AA_TO_INDEX:
                    residue_types.append(AA_TO_INDEX[res_name])
                    residue_properties.append(AA_PROPERTIES[res_name])
                else:
                    # Unknown residue - use default
                    residue_types.append(0)
                    residue_properties.append([0, 0, 0])
                
                # Chain ID (convert to number)
                chain_ids.append(ord(chain.id))
        
        return {
            'coords': np.array(coords),
            'residue_types': np.array(residue_types),
            'residue_properties': np.array(residue_properties),
            'chain_ids': np.array(chain_ids),
            'num_residues': len(coords)
        }
    
    def compute_distance_matrix(self, coords):
        """
        Compute pairwise distance matrix
        
        Args:
            coords: Nx3 array
            
        Returns:
            NxN distance matrix
        """
        n = len(coords)
        distances = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i+1, n):
                dist = np.linalg.norm(coords[i] - coords[j])
                distances[i, j] = dist
                distances[j, i] = dist
        
        return distances
    
    def build_contact_graph(self, coords, threshold=8.0):
        """
        Build contact graph (edges between nearby residues)
        
        Args:
            coords: Nx3 array
            threshold: Distance threshold in Angstroms
            
        Returns:
            edge_list: List of (i, j) tuples
            edge_distances: List of distances
        """
        n = len(coords)
        edge_list = []
        edge_distances = []
        
        for i in range(n):
            for j in range(i+1, n):
                dist = np.linalg.norm(coords[i] - coords[j])
                if dist < threshold:
                    edge_list.append((i, j))
                    edge_distances.append(dist)
        
        return edge_list, edge_distances
    
    def extract_full_features(self, target_name, difficulty="rigid_targets"):
        """
        Extract all features for a target
        
        Returns:
            dict with receptor and ligand features
        """
        # Load target
        data = self.loader.load_target(target_name, difficulty)
        
        # Extract features for receptor (unbound)
        print(f"Extracting receptor features...")
        receptor_features = self.extract_residue_features(data['receptor_unbound'])
        receptor_edges, receptor_edge_dists = self.build_contact_graph(
            receptor_features['coords']
        )
        receptor_features['edges'] = receptor_edges
        receptor_features['edge_distances'] = receptor_edge_dists
        
        # Extract features for ligand (unbound)
        print(f"Extracting ligand features...")
        ligand_features = self.extract_residue_features(data['ligand_unbound'])
        ligand_edges, ligand_edge_dists = self.build_contact_graph(
            ligand_features['coords']
        )
        ligand_features['edges'] = ligand_edges
        ligand_features['edge_distances'] = ligand_edge_dists
        
        # Extract features for ground truth (bound)
        print(f"Extracting ground truth features...")
        bound_features = self.extract_residue_features(data['bound_complex'])
        
        return {
            'target_name': target_name,
            'difficulty': difficulty,
            'partners': data['partners'],
            'receptor': receptor_features,
            'ligand': ligand_features,
            'bound': bound_features,
            'raw_data': data
        }


def test_feature_extraction(target_name="1A2K", difficulty="rigid_targets"):
    """Test feature extraction on one target"""
    print("="*80)
    print(f"TESTING FEATURE EXTRACTION ON: {target_name}")
    print("="*80)
    
    extractor = ProteinFeatureExtractor()
    
    print(f"\nExtracting features from {target_name}...")
    features = extractor.extract_full_features(target_name, difficulty)
    
    # Show receptor features
    print(f"\n{'='*80}")
    print("RECEPTOR FEATURES (INPUT)")
    print(f"{'='*80}")
    print(f"  Residues: {features['receptor']['num_residues']}")
    print(f"  Coordinates shape: {features['receptor']['coords'].shape}")
    print(f"  Residue types shape: {features['receptor']['residue_types'].shape}")
    print(f"  Properties shape: {features['receptor']['residue_properties'].shape}")
    print(f"  Contact edges: {len(features['receptor']['edges'])}")
    print(f"  Coordinate range:")
    print(f"    X: [{features['receptor']['coords'][:, 0].min():.2f}, {features['receptor']['coords'][:, 0].max():.2f}]")
    print(f"    Y: [{features['receptor']['coords'][:, 1].min():.2f}, {features['receptor']['coords'][:, 1].max():.2f}]")
    print(f"    Z: [{features['receptor']['coords'][:, 2].min():.2f}, {features['receptor']['coords'][:, 2].max():.2f}]")
    
    # Show ligand features
    print(f"\n{'='*80}")
    print("LIGAND FEATURES (INPUT)")
    print(f"{'='*80}")
    print(f"  Residues: {features['ligand']['num_residues']}")
    print(f"  Coordinates shape: {features['ligand']['coords'].shape}")
    print(f"  Residue types shape: {features['ligand']['residue_types'].shape}")
    print(f"  Properties shape: {features['ligand']['residue_properties'].shape}")
    print(f"  Contact edges: {len(features['ligand']['edges'])}")
    
    # Show ground truth
    print(f"\n{'='*80}")
    print("GROUND TRUTH (BOUND COMPLEX)")
    print(f"{'='*80}")
    print(f"  Total residues: {features['bound']['num_residues']}")
    print(f"  Coordinates shape: {features['bound']['coords'].shape}")
    
    # Show what we'll use for deep learning
    print(f"\n{'='*80}")
    print("FEATURES FOR NEURAL NETWORK")
    print(f"{'='*80}")
    print(f"\nINPUT:")
    print(f"  Receptor:")
    print(f"    - Node features: {features['receptor']['num_residues']} x (3 coords + 20 residue types + 3 properties)")
    print(f"    - Edge features: {len(features['receptor']['edges'])} edges with distances")
    print(f"  Ligand:")
    print(f"    - Node features: {features['ligand']['num_residues']} x (3 coords + 20 residue types + 3 properties)")
    print(f"    - Edge features: {len(features['ligand']['edges'])} edges with distances")
    
    print(f"\nOUTPUT (what we predict):")
    print(f"  - Transformation (rotation + translation) to dock ligand onto receptor")
    print(f"  - Goal: Minimize I-RMSD to ground truth")
    
    # Calculate center of mass for both
    receptor_com = features['receptor']['coords'].mean(axis=0)
    ligand_com = features['ligand']['coords'].mean(axis=0)
    print(f"\nCurrent positions (center of mass):")
    print(f"  Receptor: {receptor_com}")
    print(f"  Ligand: {ligand_com}")
    print(f"  Distance: {np.linalg.norm(receptor_com - ligand_com):.2f} Å")
    
    print(f"\n{'='*80}")
    print("✓ FEATURE EXTRACTION SUCCESSFUL!")
    print(f"{'='*80}")
    
    return features


if __name__ == "__main__":
    # Test on 1A2K
    features = test_feature_extraction("1A2K", "rigid_targets")
    
    print("\n" + "="*80)
    print("NEXT STEPS:")
    print("="*80)
    print("1. ✓ Load protein structures")
    print("2. ✓ Calculate I-RMSD")
    print("3. ✓ Extract features (coordinates, residue types, contacts)")
    print("4. → Next: Build Graph Neural Network")
    print("5. → Then: Train on all 254 targets")
    print("="*80)
