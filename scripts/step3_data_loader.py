#!/usr/bin/env python3
"""
Step 3: Simple Data Loader and I-RMSD Calculator
Load protein structures and calculate Interface-RMSD
"""

import numpy as np
from Bio import PDB
from pathlib import Path
from typing import Tuple, List
import sys

# Paths
PROJECT_ROOT = Path(__file__).parent.parent
ALPHARED_BENCHMARK = "/mnt/bst/bdeng2/knasif/Protein Docking/AlphaRED-main/benchmark"


class ProteinDockingLoader:
    """Simple loader for protein docking data"""
    
    def __init__(self, benchmark_path=ALPHARED_BENCHMARK):
        self.benchmark_path = Path(benchmark_path)
        self.parser = PDB.PDBParser(QUIET=True)
        
    def load_target(self, target_name, difficulty="rigid_targets"):
        """
        Load a single docking target
        
        Returns:
            dict with:
                - receptor_unbound: Structure
                - ligand_unbound: Structure  
                - bound_complex: Structure
                - partners: str (e.g., "AB_C")
        """
        target_dir = self.benchmark_path / difficulty / target_name
        
        if not target_dir.exists():
            raise FileNotFoundError(f"Target not found: {target_dir}")
        
        # Load structures
        receptor_unbound = self.parser.get_structure(
            f"{target_name}_r_u",
            target_dir / f"{target_name}_r_u.pdb"
        )
        
        ligand_unbound = self.parser.get_structure(
            f"{target_name}_l_u",
            target_dir / f"{target_name}_l_u.pdb"
        )
        
        bound_complex = self.parser.get_structure(
            f"{target_name}_b",
            target_dir / f"{target_name}_b.pdb"
        )
        
        # Read partners file
        partners_file = target_dir / "partners"
        with open(partners_file, 'r') as f:
            partners = f.read().strip().split()[1]  # e.g., "AB_C"
        
        return {
            'target_name': target_name,
            'difficulty': difficulty,
            'receptor_unbound': receptor_unbound,
            'ligand_unbound': ligand_unbound,
            'bound_complex': bound_complex,
            'partners': partners
        }
    
    def get_interface_residues(self, structure, partner_info, distance_threshold=8.0):
        """
        Get interface residues (residues within threshold distance of partner)
        
        Args:
            structure: Bio.PDB Structure
            partner_info: str like "AB_C" (receptor_ligand)
            distance_threshold: Distance in Angstroms (default 8.0)
            
        Returns:
            dict: {
                'receptor_interface': list of residues,
                'ligand_interface': list of residues
            }
        """
        # Parse partner info
        parts = partner_info.split('_')
        receptor_chains = list(parts[0])  # e.g., ['A', 'B']
        ligand_chains = list(parts[1])    # e.g., ['C']
        
        model = list(structure.get_models())[0]
        
        # Get atoms for receptor and ligand
        receptor_atoms = []
        ligand_atoms = []
        
        for chain in model:
            if chain.id in receptor_chains:
                receptor_atoms.extend(list(chain.get_atoms()))
            elif chain.id in ligand_chains:
                ligand_atoms.extend(list(chain.get_atoms()))
        
        # Find interface residues
        receptor_interface_residues = set()
        ligand_interface_residues = set()
        
        for r_atom in receptor_atoms:
            for l_atom in ligand_atoms:
                distance = r_atom - l_atom  # Biopython calculates distance
                if distance < distance_threshold:
                    receptor_interface_residues.add(r_atom.get_parent())
                    ligand_interface_residues.add(l_atom.get_parent())
        
        return {
            'receptor_interface': list(receptor_interface_residues),
            'ligand_interface': list(ligand_interface_residues)
        }
    
    def calculate_rmsd(self, coords1, coords2):
        """
        Calculate RMSD between two sets of coordinates
        
        Args:
            coords1: Nx3 numpy array
            coords2: Nx3 numpy array
            
        Returns:
            float: RMSD in Angstroms
        """
        if len(coords1) != len(coords2):
            raise ValueError("Coordinate arrays must have same length")
        
        diff = coords1 - coords2
        return np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    
    def calculate_interface_rmsd(self, structure1, structure2, partner_info):
        """
        Calculate Interface-RMSD between two structures
        
        Args:
            structure1: First structure (e.g., predicted)
            structure2: Second structure (e.g., ground truth)
            partner_info: str like "AB_C"
            
        Returns:
            float: Interface-RMSD in Angstroms
        """
        # Get interface residues from structure2 (ground truth)
        interface = self.get_interface_residues(structure2, partner_info)
        
        # Get CA atoms from interface residues
        model1 = list(structure1.get_models())[0]
        model2 = list(structure2.get_models())[0]
        
        # Collect CA coordinates from interface
        coords1 = []
        coords2 = []
        
        # Get residue IDs from ground truth
        interface_residue_ids = set()
        for res in interface['receptor_interface'] + interface['ligand_interface']:
            interface_residue_ids.add((res.get_parent().id, res.id))
        
        # Match residues between structures
        for chain in model2:
            if chain.id not in [c.id for c in model1]:
                continue
                
            chain1 = model1[chain.id]
            
            for res2 in chain:
                res_id = (chain.id, res2.id)
                if res_id not in interface_residue_ids:
                    continue
                
                # Try to find matching residue in structure1
                try:
                    res1 = chain1[res2.id]
                    
                    # Get CA atoms
                    if 'CA' in res1 and 'CA' in res2:
                        coords1.append(res1['CA'].get_coord())
                        coords2.append(res2['CA'].get_coord())
                except:
                    continue
        
        if len(coords1) == 0:
            return float('inf')
        
        coords1 = np.array(coords1)
        coords2 = np.array(coords2)
        
        # Calculate RMSD
        rmsd = self.calculate_rmsd(coords1, coords2)
        
        return rmsd


def test_single_target(target_name="1A2K", difficulty="rigid_targets"):
    """Test the loader on a single target"""
    print("="*80)
    print(f"TESTING DATA LOADER ON: {target_name}")
    print("="*80)
    
    loader = ProteinDockingLoader()
    
    # Load target
    print(f"\n1. Loading target {target_name}...")
    data = loader.load_target(target_name, difficulty)
    print(f"   ✓ Loaded successfully!")
    print(f"   Partners: {data['partners']}")
    
    # Get interface residues
    print(f"\n2. Finding interface residues...")
    interface = loader.get_interface_residues(
        data['bound_complex'], 
        data['partners']
    )
    print(f"   ✓ Receptor interface: {len(interface['receptor_interface'])} residues")
    print(f"   ✓ Ligand interface: {len(interface['ligand_interface'])} residues")
    
    # Calculate I-RMSD (bound to itself = should be ~0)
    print(f"\n3. Calculating I-RMSD (bound to itself, should be ~0)...")
    irmsd = loader.calculate_interface_rmsd(
        data['bound_complex'],
        data['bound_complex'],
        data['partners']
    )
    print(f"   ✓ I-RMSD: {irmsd:.4f} Å")
    
    # Show what we'll need for deep learning
    print(f"\n4. Data for Deep Learning:")
    print(f"   INPUT:")
    print(f"   - Receptor unbound: {len(list(data['receptor_unbound'].get_atoms()))} atoms")
    print(f"   - Ligand unbound: {len(list(data['ligand_unbound'].get_atoms()))} atoms")
    print(f"   GROUND TRUTH:")
    print(f"   - Bound complex: {len(list(data['bound_complex'].get_atoms()))} atoms")
    print(f"   TASK:")
    print(f"   - Predict how to dock receptor + ligand → minimize I-RMSD")
    
    print("\n" + "="*80)
    print("✓ DATA LOADER TEST SUCCESSFUL!")
    print("="*80)
    
    return data, interface, irmsd


if __name__ == "__main__":
    # Test on 1A2K
    data, interface, irmsd = test_single_target("1A2K", "rigid_targets")
    
    print("\n" + "="*80)
    print("NEXT STEPS:")
    print("="*80)
    print("1. ✓ We can load protein structures")
    print("2. ✓ We can identify interface residues")
    print("3. ✓ We can calculate I-RMSD")
    print("4. → Next: Extract features for deep learning")
    print("5. → Then: Build the neural network")
    print("="*80)
