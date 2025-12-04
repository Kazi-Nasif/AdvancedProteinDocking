#!/usr/bin/env python3
"""
Step 6: Training Script for Protein Docking Model
Train on all 254 targets from Benchmark 5.5
"""

import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import json
import time
from pathlib import Path
from datetime import datetime
import sys

# Add parent to path
sys.path.append(str(Path(__file__).parent.parent))

from scripts.step3_data_loader import ProteinDockingLoader
from scripts.step4_feature_extraction import ProteinFeatureExtractor
from scripts.step5_simple_gnn import ProteinDockingModel, prepare_graph_data


class DockingTrainer:
    """Trainer for protein docking model"""
    
    def __init__(self, 
                 model,
                 device='cuda',
                 learning_rate=1e-3,
                 save_dir='experiments'):
        self.model = model.to(device)
        self.device = device
        self.optimizer = optim.Adam(model.parameters(), lr=learning_rate)
        
        # Create save directory
        self.save_dir = Path(save_dir) / datetime.now().strftime('%Y%m%d_%H%M%S')
        self.save_dir.mkdir(parents=True, exist_ok=True)
        
        # Logging
        self.log_file = self.save_dir / 'training.log'
        self.metrics_file = self.save_dir / 'metrics.json'
        
        # Initialize metrics
        self.metrics = {
            'train_losses': [],
            'val_losses': [],
            'epoch_times': [],
            'best_val_loss': float('inf'),
            'best_epoch': 0
        }
        
        self.log(f"Model created with {sum(p.numel() for p in model.parameters()):,} parameters")
        self.log(f"Device: {device}")
        self.log(f"Save directory: {self.save_dir}")
    
    def log(self, message):
        """Log message to file and console"""
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        log_msg = f"[{timestamp}] {message}"
        print(log_msg)
        with open(self.log_file, 'a') as f:
            f.write(log_msg + '\n')
    
    def save_checkpoint(self, epoch, is_best=False):
        """Save model checkpoint"""
        checkpoint = {
            'epoch': epoch,
            'model_state_dict': self.model.state_dict(),
            'optimizer_state_dict': self.optimizer.state_dict(),
            'metrics': self.metrics
        }
        
        # Save regular checkpoint
        checkpoint_path = self.save_dir / f'checkpoint_epoch_{epoch}.pt'
        torch.save(checkpoint, checkpoint_path)
        
        # Save best model
        if is_best:
            best_path = self.save_dir / 'best_model.pt'
            torch.save(checkpoint, best_path)
            self.log(f"Saved best model at epoch {epoch}")
    
    def calculate_irmsd_loss(self, pred_coords, target_coords, interface_indices):
        """
        Calculate I-RMSD loss between predicted and target coordinates
        
        Args:
            pred_coords: [N, 3] predicted coordinates
            target_coords: [N, 3] target coordinates
            interface_indices: indices of interface residues
        
        Returns:
            loss: I-RMSD loss
        """
        if len(interface_indices) == 0:
            return torch.tensor(0.0, device=pred_coords.device)
        
        # Get interface coordinates
        pred_interface = pred_coords[interface_indices]
        target_interface = target_coords[interface_indices]
        
        # Calculate RMSD
        diff = pred_interface - target_interface
        rmsd = torch.sqrt(torch.mean(torch.sum(diff**2, dim=1)))
        
        return rmsd
    
    def train_epoch(self, train_data):
        """Train for one epoch"""
        self.model.train()
        total_loss = 0.0
        num_samples = 0
        
        for i, data_dict in enumerate(train_data):
            try:
                # Move data to device
                receptor_data = {
                    'node_features': data_dict['receptor_data']['node_features'].to(self.device),
                    'edge_index': data_dict['receptor_data']['edge_index'].to(self.device),
                    'coords': data_dict['receptor_data']['coords'].to(self.device)
                }
                
                ligand_data = {
                    'node_features': data_dict['ligand_data']['node_features'].to(self.device),
                    'edge_index': data_dict['ligand_data']['edge_index'].to(self.device),
                    'coords': data_dict['ligand_data']['coords'].to(self.device)
                }
                
                target_coords = data_dict['target_coords'].to(self.device)
                interface_indices = data_dict['interface_indices']
                
                # Forward pass
                rotation, translation, confidence = self.model(receptor_data, ligand_data)
                
                # Apply transformation
                pred_coords = self.model.apply_transformation(
                    ligand_data['coords'],
                    rotation,
                    translation
                )
                
                # Calculate loss (I-RMSD)
                loss = self.calculate_irmsd_loss(pred_coords, target_coords, interface_indices)
                
                # Add confidence loss (encourage high confidence for good predictions)
                confidence_loss = -torch.log(confidence + 1e-6) if loss < 4.0 else torch.log(confidence + 1e-6)
                total_loss_value = loss + 0.1 * confidence_loss
                
                # Backward pass
                self.optimizer.zero_grad()
                total_loss_value.backward()
                torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
                self.optimizer.step()
                
                total_loss += loss.item()
                num_samples += 1
                
                if (i + 1) % 10 == 0:
                    self.log(f"  Sample {i+1}/{len(train_data)}, Loss: {loss.item():.4f}")
                
            except Exception as e:
                self.log(f"  Error processing sample {i}: {e}")
                continue
        
        avg_loss = total_loss / max(num_samples, 1)
        return avg_loss
    
    def validate(self, val_data):
        """Validate model"""
        self.model.eval()
        total_loss = 0.0
        num_samples = 0
        
        with torch.no_grad():
            for i, data_dict in enumerate(val_data):
                try:
                    # Move data to device
                    receptor_data = {
                        'node_features': data_dict['receptor_data']['node_features'].to(self.device),
                        'edge_index': data_dict['receptor_data']['edge_index'].to(self.device),
                        'coords': data_dict['receptor_data']['coords'].to(self.device)
                    }
                    
                    ligand_data = {
                        'node_features': data_dict['ligand_data']['node_features'].to(self.device),
                        'edge_index': data_dict['ligand_data']['edge_index'].to(self.device),
                        'coords': data_dict['ligand_data']['coords'].to(self.device)
                    }
                    
                    target_coords = data_dict['target_coords'].to(self.device)
                    interface_indices = data_dict['interface_indices']
                    
                    # Forward pass
                    rotation, translation, confidence = self.model(receptor_data, ligand_data)
                    
                    # Apply transformation
                    pred_coords = self.model.apply_transformation(
                        ligand_data['coords'],
                        rotation,
                        translation
                    )
                    
                    # Calculate loss
                    loss = self.calculate_irmsd_loss(pred_coords, target_coords, interface_indices)
                    
                    total_loss += loss.item()
                    num_samples += 1
                    
                except Exception as e:
                    self.log(f"  Error validating sample {i}: {e}")
                    continue
        
        avg_loss = total_loss / max(num_samples, 1)
        return avg_loss
    
    def train(self, train_data, val_data, num_epochs=50, save_every=5):
        """Main training loop"""
        self.log("="*80)
        self.log("STARTING TRAINING")
        self.log("="*80)
        self.log(f"Training samples: {len(train_data)}")
        self.log(f"Validation samples: {len(val_data)}")
        self.log(f"Number of epochs: {num_epochs}")
        self.log("="*80)
        
        for epoch in range(1, num_epochs + 1):
            epoch_start = time.time()
            
            self.log(f"\nEpoch {epoch}/{num_epochs}")
            self.log("-" * 80)
            
            # Train
            train_loss = self.train_epoch(train_data)
            self.log(f"Training Loss: {train_loss:.4f}")
            
            # Validate
            val_loss = self.validate(val_data)
            self.log(f"Validation Loss: {val_loss:.4f}")
            
            # Update metrics
            epoch_time = time.time() - epoch_start
            self.metrics['train_losses'].append(train_loss)
            self.metrics['val_losses'].append(val_loss)
            self.metrics['epoch_times'].append(epoch_time)
            
            # Save best model
            is_best = val_loss < self.metrics['best_val_loss']
            if is_best:
                self.metrics['best_val_loss'] = val_loss
                self.metrics['best_epoch'] = epoch
            
            # Save checkpoint
            if epoch % save_every == 0 or is_best:
                self.save_checkpoint(epoch, is_best)
            
            # Save metrics
            with open(self.metrics_file, 'w') as f:
                json.dump(self.metrics, f, indent=2)
            
            self.log(f"Epoch time: {epoch_time:.2f}s")
            self.log(f"Best validation loss: {self.metrics['best_val_loss']:.4f} (epoch {self.metrics['best_epoch']})")
        
        self.log("\n" + "="*80)
        self.log("TRAINING COMPLETE")
        self.log("="*80)
        self.log(f"Best model: epoch {self.metrics['best_epoch']}, val_loss: {self.metrics['best_val_loss']:.4f}")
        self.log(f"Results saved to: {self.save_dir}")


def prepare_dataset(difficulty="rigid_targets", max_samples=None):
    """Prepare dataset from benchmark"""
    print(f"Preparing dataset from {difficulty}...")
    
    loader = ProteinDockingLoader()
    extractor = ProteinFeatureExtractor()
    
    benchmark_path = Path(loader.benchmark_path) / difficulty
    target_dirs = sorted([d for d in benchmark_path.iterdir() if d.is_dir()])
    
    if max_samples:
        target_dirs = target_dirs[:max_samples]
    
    dataset = []
    
    for i, target_dir in enumerate(target_dirs):
        target_name = target_dir.name
        
        try:
            print(f"Processing {i+1}/{len(target_dirs)}: {target_name}...", end=' ')
            
            # Extract features
            features = extractor.extract_full_features(target_name, difficulty)
            
            # Prepare graph data
            receptor_data, ligand_data = prepare_graph_data(features)
            
            # Get target coordinates (from bound complex)
            target_coords = torch.tensor(
                features['bound']['coords'][len(features['receptor']['coords']):],
                dtype=torch.float32
            )
            
            # Get interface indices
            interface = loader.get_interface_residues(
                features['raw_data']['bound_complex'],
                features['partners']
            )
            ligand_interface_indices = [
                i for i, res in enumerate(features['raw_data']['ligand_unbound'].get_residues())
                if res in interface['ligand_interface']
            ]
            
            data_dict = {
                'target_name': target_name,
                'receptor_data': receptor_data,
                'ligand_data': ligand_data,
                'target_coords': target_coords,
                'interface_indices': ligand_interface_indices
            }
            
            dataset.append(data_dict)
            print("✓")
            
        except Exception as e:
            print(f"✗ Error: {e}")
            continue
    
    print(f"Dataset prepared: {len(dataset)} samples")
    return dataset


def main():
    """Main training function"""
    print("="*80)
    print("PROTEIN DOCKING MODEL TRAINING")
    print("="*80)
    
    # Settings
    DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'
    NUM_EPOCHS = 50
    LEARNING_RATE = 1e-3
    
    print(f"\nDevice: {DEVICE}")
    print(f"Epochs: {NUM_EPOCHS}")
    print(f"Learning rate: {LEARNING_RATE}")
    
    # Prepare dataset (start with rigid targets only for faster initial training)
    print("\n" + "="*80)
    print("PREPARING DATASET")
    print("="*80)
    
    # Load rigid targets first (easier to learn)
    rigid_dataset = prepare_dataset("rigid_targets", max_samples=50)  # Start small
    
    # Split train/val (80/20)
    split_idx = int(0.8 * len(rigid_dataset))
    train_data = rigid_dataset[:split_idx]
    val_data = rigid_dataset[split_idx:]
    
    print(f"\nTrain: {len(train_data)} samples")
    print(f"Val: {len(val_data)} samples")
    
    # Create model
    print("\n" + "="*80)
    print("CREATING MODEL")
    print("="*80)
    
    model = ProteinDockingModel(
        node_features=26,
        hidden_dim=128,
        num_layers=3
    )
    
    print(f"Parameters: {sum(p.numel() for p in model.parameters()):,}")
    
    # Create trainer
    trainer = DockingTrainer(
        model=model,
        device=DEVICE,
        learning_rate=LEARNING_RATE,
        save_dir='experiments'
    )
    
    # Train
    trainer.train(
        train_data=train_data,
        val_data=val_data,
        num_epochs=NUM_EPOCHS,
        save_every=5
    )
    
    print("\n" + "="*80)
    print("DONE!")
    print("="*80)


if __name__ == "__main__":
    main()
