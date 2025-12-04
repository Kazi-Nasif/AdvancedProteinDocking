#!/usr/bin/env python3
"""
Step 5: Simple Graph Neural Network for Protein Docking
Pure PyTorch implementation - no external docking tools
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from pathlib import Path
import sys

# Add parent to path
sys.path.append(str(Path(__file__).parent.parent))

from scripts.step4_feature_extraction import ProteinFeatureExtractor


class ProteinGraphConv(nn.Module):
    """Simple Graph Convolution Layer for proteins"""
    
    def __init__(self, in_features, out_features):
        super().__init__()
        self.fc = nn.Linear(in_features, out_features)
        self.norm = nn.LayerNorm(out_features)
        
    def forward(self, x, edge_index):
        """
        Args:
            x: Node features [N, in_features]
            edge_index: Edge connections [2, num_edges]
        Returns:
            Updated node features [N, out_features]
        """
        # Transform features
        x_transformed = self.fc(x)
        
        # Simple message passing: average neighbors
        x_aggregated = x_transformed.clone()
        
        for i in range(edge_index.shape[1]):
            src, dst = edge_index[0, i], edge_index[1, i]
            x_aggregated[dst] += x_transformed[src]
        
        # Normalize and activate
        x_out = self.norm(x_aggregated)
        x_out = F.relu(x_out)
        
        return x_out


class ProteinEncoder(nn.Module):
    """Encode protein structure into embedding"""
    
    def __init__(self, node_features=26, hidden_dim=128, num_layers=3):
        """
        Args:
            node_features: Input features per node (3 coords + 20 types + 3 props)
            hidden_dim: Hidden dimension size
            num_layers: Number of graph conv layers
        """
        super().__init__()
        
        # Input projection
        self.input_proj = nn.Linear(node_features, hidden_dim)
        
        # Graph convolution layers
        self.conv_layers = nn.ModuleList([
            ProteinGraphConv(hidden_dim, hidden_dim)
            for _ in range(num_layers)
        ])
        
        # Global pooling (graph-level representation)
        self.pool = nn.Linear(hidden_dim, hidden_dim)
        
    def forward(self, node_features, edge_index):
        """
        Args:
            node_features: [N, node_features]
            edge_index: [2, num_edges]
        Returns:
            node_embeddings: [N, hidden_dim]
            graph_embedding: [hidden_dim]
        """
        # Project input
        x = self.input_proj(node_features)
        x = F.relu(x)
        
        # Apply graph convolutions
        for conv in self.conv_layers:
            x = conv(x, edge_index) + x  # Residual connection
        
        # Global pooling (mean)
        graph_embedding = torch.mean(x, dim=0)
        graph_embedding = self.pool(graph_embedding)
        
        return x, graph_embedding


class DockingHead(nn.Module):
    """Predict docking transformation"""
    
    def __init__(self, hidden_dim=128):
        super().__init__()
        
        # Combine receptor and ligand embeddings
        self.combine = nn.Sequential(
            nn.Linear(hidden_dim * 2, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU()
        )
        
        # Predict rotation (as rotation vector - 3D)
        self.rotation_head = nn.Linear(hidden_dim, 3)
        
        # Predict translation (3D vector)
        self.translation_head = nn.Linear(hidden_dim, 3)
        
        # Confidence score
        self.confidence_head = nn.Linear(hidden_dim, 1)
        
    def forward(self, receptor_embedding, ligand_embedding):
        """
        Args:
            receptor_embedding: [hidden_dim]
            ligand_embedding: [hidden_dim]
        Returns:
            rotation: [3] - rotation vector
            translation: [3] - translation vector
            confidence: [1] - confidence score
        """
        # Combine embeddings
        combined = torch.cat([receptor_embedding, ligand_embedding], dim=0)
        features = self.combine(combined)
        
        # Predict transformation
        rotation = self.rotation_head(features)
        translation = self.translation_head(features)
        confidence = torch.sigmoid(self.confidence_head(features))
        
        return rotation, translation, confidence


class ProteinDockingModel(nn.Module):
    """Complete protein docking model"""
    
    def __init__(self, node_features=26, hidden_dim=128, num_layers=3):
        super().__init__()
        
        # Separate encoders for receptor and ligand
        self.receptor_encoder = ProteinEncoder(node_features, hidden_dim, num_layers)
        self.ligand_encoder = ProteinEncoder(node_features, hidden_dim, num_layers)
        
        # Docking head
        self.docking_head = DockingHead(hidden_dim)
        
    def forward(self, receptor_data, ligand_data):
        """
        Args:
            receptor_data: dict with 'node_features' and 'edge_index'
            ligand_data: dict with 'node_features' and 'edge_index'
        Returns:
            rotation: [3]
            translation: [3]
            confidence: [1]
        """
        # Encode receptor
        _, receptor_embedding = self.receptor_encoder(
            receptor_data['node_features'],
            receptor_data['edge_index']
        )
        
        # Encode ligand
        _, ligand_embedding = self.ligand_encoder(
            ligand_data['node_features'],
            ligand_data['edge_index']
        )
        
        # Predict docking transformation
        rotation, translation, confidence = self.docking_head(
            receptor_embedding, ligand_embedding
        )
        
        return rotation, translation, confidence
    
    def apply_transformation(self, coords, rotation, translation):
        """
        Apply rotation and translation to coordinates
        
        Args:
            coords: [N, 3] coordinates
            rotation: [3] rotation vector (axis-angle representation)
            translation: [3] translation vector
        Returns:
            transformed_coords: [N, 3]
        """
        # Convert rotation vector to rotation matrix (simplified)
        # In practice, use proper axis-angle to matrix conversion
        angle = torch.norm(rotation)
        if angle < 1e-6:
            rotation_matrix = torch.eye(3, device=coords.device)
        else:
            axis = rotation / angle
            # Rodrigues' rotation formula (simplified)
            K = torch.zeros(3, 3, device=coords.device)
            K[0, 1] = -axis[2]
            K[0, 2] = axis[1]
            K[1, 0] = axis[2]
            K[1, 2] = -axis[0]
            K[2, 0] = -axis[1]
            K[2, 1] = axis[0]
            
            rotation_matrix = (torch.eye(3, device=coords.device) + 
                             torch.sin(angle) * K + 
                             (1 - torch.cos(angle)) * torch.mm(K, K))
        
        # Apply rotation and translation
        coords_rotated = torch.mm(coords, rotation_matrix.t())
        coords_transformed = coords_rotated + translation.unsqueeze(0)
        
        return coords_transformed


def prepare_graph_data(features):
    """
    Convert extracted features to PyTorch tensors
    
    Args:
        features: Output from ProteinFeatureExtractor
    Returns:
        receptor_data, ligand_data as dicts
    """
    # Prepare receptor data
    receptor_coords = torch.tensor(features['receptor']['coords'], dtype=torch.float32)
    receptor_types = torch.tensor(features['receptor']['residue_types'], dtype=torch.long)
    receptor_props = torch.tensor(features['receptor']['residue_properties'], dtype=torch.float32)
    
    # One-hot encode residue types
    receptor_types_onehot = F.one_hot(receptor_types, num_classes=20).float()
    
    # Combine features: [coords, one-hot types, properties]
    receptor_features = torch.cat([
        receptor_coords,
        receptor_types_onehot,
        receptor_props
    ], dim=1)
    
    # Prepare edge index
    receptor_edges = features['receptor']['edges']
    receptor_edge_index = torch.tensor(
        [[e[0] for e in receptor_edges] + [e[1] for e in receptor_edges],
         [e[1] for e in receptor_edges] + [e[0] for e in receptor_edges]],
        dtype=torch.long
    )
    
    receptor_data = {
        'node_features': receptor_features,
        'edge_index': receptor_edge_index,
        'coords': receptor_coords
    }
    
    # Prepare ligand data (same process)
    ligand_coords = torch.tensor(features['ligand']['coords'], dtype=torch.float32)
    ligand_types = torch.tensor(features['ligand']['residue_types'], dtype=torch.long)
    ligand_props = torch.tensor(features['ligand']['residue_properties'], dtype=torch.float32)
    
    ligand_types_onehot = F.one_hot(ligand_types, num_classes=20).float()
    
    ligand_features = torch.cat([
        ligand_coords,
        ligand_types_onehot,
        ligand_props
    ], dim=1)
    
    ligand_edges = features['ligand']['edges']
    ligand_edge_index = torch.tensor(
        [[e[0] for e in ligand_edges] + [e[1] for e in ligand_edges],
         [e[1] for e in ligand_edges] + [e[0] for e in ligand_edges]],
        dtype=torch.long
    )
    
    ligand_data = {
        'node_features': ligand_features,
        'edge_index': ligand_edge_index,
        'coords': ligand_coords
    }
    
    return receptor_data, ligand_data


def test_model(target_name="1A2K", difficulty="rigid_targets"):
    """Test the model on one target"""
    print("="*80)
    print(f"TESTING PROTEIN DOCKING MODEL ON: {target_name}")
    print("="*80)
    
    # Extract features
    print("\n1. Extracting features...")
    extractor = ProteinFeatureExtractor()
    features = extractor.extract_full_features(target_name, difficulty)
    print("   ✓ Features extracted")
    
    # Prepare data
    print("\n2. Preparing graph data...")
    receptor_data, ligand_data = prepare_graph_data(features)
    print(f"   ✓ Receptor: {receptor_data['node_features'].shape}")
    print(f"   ✓ Ligand: {ligand_data['node_features'].shape}")
    
    # Create model
    print("\n3. Creating model...")
    model = ProteinDockingModel(
        node_features=26,  # 3 coords + 20 types + 3 props
        hidden_dim=128,
        num_layers=3
    )
    print(f"   ✓ Model created")
    print(f"   Parameters: {sum(p.numel() for p in model.parameters()):,}")
    
    # Forward pass
    print("\n4. Running forward pass...")
    with torch.no_grad():
        rotation, translation, confidence = model(receptor_data, ligand_data)
    
    print(f"   ✓ Predicted transformation:")
    print(f"     Rotation: {rotation.numpy()}")
    print(f"     Translation: {translation.numpy()}")
    print(f"     Confidence: {confidence.item():.4f}")
    
    # Apply transformation to ligand
    print("\n5. Applying transformation...")
    ligand_coords_transformed = model.apply_transformation(
        ligand_data['coords'],
        rotation,
        translation
    )
    
    # Calculate distance moved
    original_center = ligand_data['coords'].mean(dim=0)
    new_center = ligand_coords_transformed.mean(dim=0)
    distance_moved = torch.norm(new_center - original_center).item()
    
    print(f"   ✓ Ligand moved: {distance_moved:.2f} Å")
    print(f"   Original center: {original_center.numpy()}")
    print(f"   New center: {new_center.numpy()}")
    
    print("\n" + "="*80)
    print("✓ MODEL TEST SUCCESSFUL!")
    print("="*80)
    
    return model, features, receptor_data, ligand_data


if __name__ == "__main__":
    # Test on 1A2K
    model, features, receptor_data, ligand_data = test_model("1A2K", "rigid_targets")
    
    print("\n" + "="*80)
    print("NEXT STEPS:")
    print("="*80)
    print("1. ✓ Load and extract features")
    print("2. ✓ Build Graph Neural Network")
    print("3. ✓ Forward pass works!")
    print("4. → Next: Training loop with I-RMSD loss")
    print("5. → Then: Train on all 254 targets")
    print("="*80)
    print("\nModel Architecture:")
    print("  - Separate encoders for receptor and ligand")
    print("  - 3 layers of graph convolutions")
    print("  - Predicts: rotation + translation + confidence")
    print("  - Pure PyTorch - no external tools!")
    print("="*80)
