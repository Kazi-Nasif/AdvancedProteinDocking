import sys
sys.path.append('/mnt/bst/bdeng2/knasif/Protein Docking/AdvancedProteinDocking_Step1')

print("Starting quick test...")

from scripts.step3_data_loader import ProteinDockingLoader
from scripts.step4_feature_extraction import ProteinFeatureExtractor
from scripts.step5_simple_gnn import prepare_graph_data

print("Imports successful!")

print("Loading one target...")
extractor = ProteinFeatureExtractor()
features = extractor.extract_full_features("1A2K", "rigid_targets")

print("Preparing graph data...")
receptor_data, ligand_data = prepare_graph_data(features)

print("SUCCESS! Data loading works.")
print(f"Receptor: {receptor_data['node_features'].shape}")
print(f"Ligand: {ligand_data['node_features'].shape}")
