#!/usr/bin/env python3
"""
Benchmark 5.5 Dataset Analysis
Analyzes the protein-protein docking benchmark dataset
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
from pathlib import Path
import sys

# Add project root to path
PROJECT_ROOT = Path("/home/claude/AdvancedProteinDocking")
sys.path.append(str(PROJECT_ROOT / 'src'))

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 10

class BenchmarkAnalyzer:
    """Analyzer for Benchmark 5.5 dataset"""
    
    def __init__(self, excel_path='/mnt/project/Table_BM5_5.xlsx'):
        """Initialize analyzer"""
        self.excel_path = excel_path
        self.df = None
        self.results_dir = PROJECT_ROOT / 'results'
        self.results_dir.mkdir(exist_ok=True)
        
    def load_data(self):
        """Load the benchmark dataset"""
        print("Loading Benchmark 5.5 dataset...")
        self.df = pd.read_excel(self.excel_path)
        print(f"✓ Loaded {len(self.df)} protein complexes")
        return self.df
    
    def basic_statistics(self):
        """Compute basic statistics"""
        print("\n" + "="*80)
        print("BASIC DATASET STATISTICS")
        print("="*80)
        
        print(f"\nTotal number of complexes: {len(self.df)}")
        print(f"\nColumns: {self.df.columns.tolist()}")
        
        # Category distribution
        print("\n--- Category Distribution ---")
        category_counts = self.df['Category'].value_counts()
        print(category_counts)
        
        # Difficulty distribution
        print("\n--- Difficulty Distribution ---")
        difficulty_counts = self.df['Difficulty'].value_counts()
        print(difficulty_counts)
        
        # I-RMSD statistics
        print("\n--- Interface-RMSD (I-RMSD) Statistics ---")
        print(f"Mean: {self.df['I-RMSD'].mean():.3f} Å")
        print(f"Median: {self.df['I-RMSD'].median():.3f} Å")
        print(f"Std: {self.df['I-RMSD'].std():.3f} Å")
        print(f"Min: {self.df['I-RMSD'].min():.3f} Å")
        print(f"Max: {self.df['I-RMSD'].max():.3f} Å")
        print(f"< 4.0 Å (near-native): {(self.df['I-RMSD'] < 4.0).sum()} ({(self.df['I-RMSD'] < 4.0).sum()/len(self.df)*100:.1f}%)")
        
        # ΔASA statistics
        print("\n--- ΔASA Statistics ---")
        print(f"Mean: {self.df['ΔASA'].mean():.1f} Ų")
        print(f"Median: {self.df['ΔASA'].median():.1f} Ų")
        print(f"Std: {self.df['ΔASA'].std():.1f} Ų")
        print(f"Min: {self.df['ΔASA'].min():.1f} Ų")
        print(f"Max: {self.df['ΔASA'].max():.1f} Ų")
        
        # Save statistics
        stats = {
            'total_complexes': len(self.df),
            'categories': category_counts.to_dict(),
            'difficulties': difficulty_counts.to_dict(),
            'irmsd': {
                'mean': float(self.df['I-RMSD'].mean()),
                'median': float(self.df['I-RMSD'].median()),
                'std': float(self.df['I-RMSD'].std()),
                'min': float(self.df['I-RMSD'].min()),
                'max': float(self.df['I-RMSD'].max()),
                'near_native_count': int((self.df['I-RMSD'] < 4.0).sum()),
                'near_native_percent': float((self.df['I-RMSD'] < 4.0).sum()/len(self.df)*100)
            },
            'dasa': {
                'mean': float(self.df['ΔASA'].mean()),
                'median': float(self.df['ΔASA'].median()),
                'std': float(self.df['ΔASA'].std()),
                'min': float(self.df['ΔASA'].min()),
                'max': float(self.df['ΔASA'].max())
            }
        }
        
        stats_path = self.results_dir / 'benchmark_statistics.json'
        with open(stats_path, 'w') as f:
            json.dump(stats, f, indent=2)
        print(f"\n✓ Statistics saved to: {stats_path}")
        
        return stats
    
    def analyze_categories(self):
        """Analyze different categories"""
        print("\n" + "="*80)
        print("CATEGORY ANALYSIS")
        print("="*80)
        
        # Unique categories
        categories = self.df['Category'].unique()
        print(f"\nUnique categories: {list(categories)}")
        
        # Category meanings (common in protein docking)
        category_meanings = {
            'AA': 'Antibody-Antigen',
            'AB': 'Antibody-Bound',
            'AS': 'Antibody-Substrate',
            'A': 'Antibody',
            'E': 'Enzyme',
            'ER': 'Enzyme-Receptor',
            'EI': 'Enzyme-Inhibitor',
            'ES': 'Enzyme-Substrate',
            'O': 'Other',
            'OG': 'Other-G protein',
            'OR': 'Other-Receptor',
            'OX': 'Other-Other'
        }
        
        print("\nCategory Analysis:")
        for cat in sorted(categories):
            cat_data = self.df[self.df['Category'] == cat]
            meaning = category_meanings.get(cat, 'Unknown')
            print(f"\n{cat} ({meaning}): {len(cat_data)} complexes")
            print(f"  I-RMSD: {cat_data['I-RMSD'].mean():.3f} ± {cat_data['I-RMSD'].std():.3f} Å")
            print(f"  ΔASA: {cat_data['ΔASA'].mean():.1f} ± {cat_data['ΔASA'].std():.1f} Ų")
            
            # Difficulty breakdown
            diff_counts = cat_data['Difficulty'].value_counts()
            for diff, count in diff_counts.items():
                print(f"  - {diff}: {count}")
        
        # Antibody complexes analysis (important for our goal)
        antibody_cats = ['AA', 'AB', 'AS', 'A']
        antibody_data = self.df[self.df['Category'].isin(antibody_cats)]
        
        print("\n" + "-"*80)
        print(f"ANTIBODY COMPLEXES (AA, AB, AS, A): {len(antibody_data)} total")
        print(f"  I-RMSD: {antibody_data['I-RMSD'].mean():.3f} ± {antibody_data['I-RMSD'].std():.3f} Å")
        print(f"  Success rate if I-RMSD < 4.0 Å: {(antibody_data['I-RMSD'] < 4.0).sum()/len(antibody_data)*100:.1f}%")
        print(f"  (AlphaRed achieved 43% success on antibody-antigen targets)")
        print("-"*80)
    
    def analyze_difficulty(self):
        """Analyze difficulty levels"""
        print("\n" + "="*80)
        print("DIFFICULTY ANALYSIS")
        print("="*80)
        
        for difficulty in ['Rigid-body', 'Medium', 'Difficult']:
            diff_data = self.df[self.df['Difficulty'] == difficulty]
            if len(diff_data) > 0:
                print(f"\n{difficulty}: {len(diff_data)} complexes")
                print(f"  I-RMSD: {diff_data['I-RMSD'].mean():.3f} ± {diff_data['I-RMSD'].std():.3f} Å")
                print(f"  ΔASA: {diff_data['ΔASA'].mean():.1f} ± {diff_data['ΔASA'].std():.1f} Ų")
                print(f"  Near-native (I-RMSD < 4.0): {(diff_data['I-RMSD'] < 4.0).sum()} ({(diff_data['I-RMSD'] < 4.0).sum()/len(diff_data)*100:.1f}%)")
    
    def analyze_pdb_structure(self):
        """Analyze PDB ID structure"""
        print("\n" + "="*80)
        print("PDB STRUCTURE ANALYSIS")
        print("="*80)
        
        print("\nSample PDB IDs:")
        print(self.df[['Complex', 'PDB ID 1', 'PDB ID 2']].head(10))
        
        # Parse complex names to understand structure
        print("\nComplex name structure:")
        for idx in range(min(5, len(self.df))):
            complex_name = self.df.iloc[idx]['Complex']
            pdb1 = self.df.iloc[idx]['PDB ID 1']
            pdb2 = self.df.iloc[idx]['PDB ID 2']
            print(f"  {complex_name} -> PDB1: {pdb1}, PDB2: {pdb2}")
    
    def create_visualizations(self):
        """Create comprehensive visualizations"""
        print("\n" + "="*80)
        print("CREATING VISUALIZATIONS")
        print("="*80)
        
        fig_dir = self.results_dir / 'figures'
        fig_dir.mkdir(exist_ok=True)
        
        # 1. Distribution of I-RMSD
        plt.figure(figsize=(12, 5))
        plt.subplot(1, 2, 1)
        plt.hist(self.df['I-RMSD'], bins=50, edgecolor='black', alpha=0.7)
        plt.axvline(4.0, color='red', linestyle='--', linewidth=2, label='Near-native threshold (4.0 Å)')
        plt.xlabel('Interface-RMSD (Å)')
        plt.ylabel('Frequency')
        plt.title('Distribution of Interface-RMSD')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 2. Distribution of ΔASA
        plt.subplot(1, 2, 2)
        plt.hist(self.df['ΔASA'], bins=50, edgecolor='black', alpha=0.7, color='green')
        plt.xlabel('ΔASA (Ų)')
        plt.ylabel('Frequency')
        plt.title('Distribution of ΔASA')
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(fig_dir / '01_distributions.png', dpi=300, bbox_inches='tight')
        print(f"✓ Saved: {fig_dir / '01_distributions.png'}")
        plt.close()
        
        # 3. I-RMSD by Category
        plt.figure(figsize=(14, 6))
        categories_sorted = self.df.groupby('Category')['I-RMSD'].median().sort_values().index
        sns.boxplot(data=self.df, x='Category', y='I-RMSD', order=categories_sorted)
        plt.axhline(4.0, color='red', linestyle='--', linewidth=2, label='Near-native threshold')
        plt.xticks(rotation=45)
        plt.ylabel('Interface-RMSD (Å)')
        plt.title('Interface-RMSD Distribution by Category')
        plt.legend()
        plt.grid(True, alpha=0.3, axis='y')
        plt.tight_layout()
        plt.savefig(fig_dir / '02_irmsd_by_category.png', dpi=300, bbox_inches='tight')
        print(f"✓ Saved: {fig_dir / '02_irmsd_by_category.png'}")
        plt.close()
        
        # 4. I-RMSD by Difficulty
        plt.figure(figsize=(10, 6))
        difficulty_order = ['Rigid-body', 'Medium', 'Difficult']
        sns.boxplot(data=self.df, x='Difficulty', y='I-RMSD', order=difficulty_order)
        plt.axhline(4.0, color='red', linestyle='--', linewidth=2, label='Near-native threshold')
        plt.ylabel('Interface-RMSD (Å)')
        plt.title('Interface-RMSD Distribution by Difficulty')
        plt.legend()
        plt.grid(True, alpha=0.3, axis='y')
        plt.tight_layout()
        plt.savefig(fig_dir / '03_irmsd_by_difficulty.png', dpi=300, bbox_inches='tight')
        print(f"✓ Saved: {fig_dir / '03_irmsd_by_difficulty.png'}")
        plt.close()
        
        # 5. ΔASA vs I-RMSD
        plt.figure(figsize=(10, 8))
        scatter = plt.scatter(self.df['ΔASA'], self.df['I-RMSD'], 
                            c=self.df['Difficulty'].map({'Rigid-body': 0, 'Medium': 1, 'Difficult': 2}),
                            cmap='viridis', alpha=0.6, s=50)
        plt.axhline(4.0, color='red', linestyle='--', linewidth=2, alpha=0.5)
        plt.xlabel('ΔASA (Ų)')
        plt.ylabel('Interface-RMSD (Å)')
        plt.title('ΔASA vs Interface-RMSD')
        cbar = plt.colorbar(scatter, ticks=[0, 1, 2])
        cbar.set_label('Difficulty')
        cbar.ax.set_yticklabels(['Rigid', 'Medium', 'Difficult'])
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(fig_dir / '04_dasa_vs_irmsd.png', dpi=300, bbox_inches='tight')
        print(f"✓ Saved: {fig_dir / '04_dasa_vs_irmsd.png'}")
        plt.close()
        
        # 6. Category and Difficulty Heatmap
        plt.figure(figsize=(12, 6))
        pivot = pd.crosstab(self.df['Category'], self.df['Difficulty'])
        sns.heatmap(pivot, annot=True, fmt='d', cmap='YlOrRd', cbar_kws={'label': 'Count'})
        plt.title('Number of Complexes by Category and Difficulty')
        plt.tight_layout()
        plt.savefig(fig_dir / '05_category_difficulty_heatmap.png', dpi=300, bbox_inches='tight')
        print(f"✓ Saved: {fig_dir / '05_category_difficulty_heatmap.png'}")
        plt.close()
        
        # 7. Success rate estimation (I-RMSD < 4.0)
        plt.figure(figsize=(14, 6))
        
        plt.subplot(1, 2, 1)
        category_success = self.df.groupby('Category').apply(
            lambda x: (x['I-RMSD'] < 4.0).sum() / len(x) * 100
        ).sort_values(ascending=False)
        category_success.plot(kind='bar', color='steelblue')
        plt.axhline(63, color='red', linestyle='--', linewidth=2, label='AlphaRed overall (63%)')
        plt.ylabel('Success Rate (%)')
        plt.xlabel('Category')
        plt.title('Potential Success Rate by Category\n(I-RMSD < 4.0 Å)')
        plt.xticks(rotation=45)
        plt.legend()
        plt.grid(True, alpha=0.3, axis='y')
        
        plt.subplot(1, 2, 2)
        difficulty_success = self.df.groupby('Difficulty').apply(
            lambda x: (x['I-RMSD'] < 4.0).sum() / len(x) * 100
        ).reindex(['Rigid-body', 'Medium', 'Difficult'])
        difficulty_success.plot(kind='bar', color='coral')
        plt.axhline(63, color='red', linestyle='--', linewidth=2, label='AlphaRed overall (63%)')
        plt.ylabel('Success Rate (%)')
        plt.xlabel('Difficulty')
        plt.title('Potential Success Rate by Difficulty\n(I-RMSD < 4.0 Å)')
        plt.xticks(rotation=45)
        plt.legend()
        plt.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        plt.savefig(fig_dir / '06_success_rates.png', dpi=300, bbox_inches='tight')
        print(f"✓ Saved: {fig_dir / '06_success_rates.png'}")
        plt.close()
        
        print(f"\n✓ All visualizations saved to: {fig_dir}")
    
    def generate_target_list(self):
        """Generate a list of all targets for download"""
        print("\n" + "="*80)
        print("GENERATING TARGET LIST")
        print("="*80)
        
        target_list = []
        for idx, row in self.df.iterrows():
            target_info = {
                'complex_id': row['Complex'],
                'category': row['Category'],
                'difficulty': row['Difficulty'],
                'pdb_id_1': row['PDB ID 1'],
                'protein_1': row['Protein 1'],
                'pdb_id_2': row['PDB ID 2'],
                'protein_2': row['Protein 2'],
                'irmsd': row['I-RMSD'],
                'dasa': row['ΔASA']
            }
            target_list.append(target_info)
        
        # Save as JSON
        target_list_path = PROJECT_ROOT / 'data' / 'benchmark' / 'target_list.json'
        with open(target_list_path, 'w') as f:
            json.dump(target_list, f, indent=2)
        print(f"✓ Target list saved to: {target_list_path}")
        
        # Save as CSV for easy viewing
        target_csv_path = PROJECT_ROOT / 'data' / 'benchmark' / 'target_list.csv'
        pd.DataFrame(target_list).to_csv(target_csv_path, index=False)
        print(f"✓ Target list CSV saved to: {target_csv_path}")
        
        return target_list
    
    def create_report(self):
        """Create a comprehensive text report"""
        print("\n" + "="*80)
        print("CREATING COMPREHENSIVE REPORT")
        print("="*80)
        
        report = []
        report.append("="*80)
        report.append("BENCHMARK 5.5 DATASET - COMPREHENSIVE ANALYSIS REPORT")
        report.append("="*80)
        report.append("")
        
        # Dataset Overview
        report.append("1. DATASET OVERVIEW")
        report.append("-" * 80)
        report.append(f"Total Complexes: {len(self.df)}")
        report.append(f"Source: Docking Benchmark 5.5 (DB5.5)")
        report.append(f"Data Type: Protein-Protein Complexes")
        report.append(f"Structure: Bound and Unbound PDB files")
        report.append("")
        
        # Key Metrics
        report.append("2. KEY METRICS")
        report.append("-" * 80)
        report.append("Interface-RMSD (I-RMSD):")
        report.append(f"  - Range: {self.df['I-RMSD'].min():.2f} - {self.df['I-RMSD'].max():.2f} Å")
        report.append(f"  - Mean ± Std: {self.df['I-RMSD'].mean():.2f} ± {self.df['I-RMSD'].std():.2f} Å")
        report.append(f"  - Median: {self.df['I-RMSD'].median():.2f} Å")
        report.append(f"  - Near-native (< 4.0 Å): {(self.df['I-RMSD'] < 4.0).sum()} ({(self.df['I-RMSD'] < 4.0).sum()/len(self.df)*100:.1f}%)")
        report.append("")
        report.append("ΔASA (Delta Accessible Surface Area):")
        report.append(f"  - Range: {self.df['ΔASA'].min():.0f} - {self.df['ΔASA'].max():.0f} Ų")
        report.append(f"  - Mean ± Std: {self.df['ΔASA'].mean():.0f} ± {self.df['ΔASA'].std():.0f} Ų")
        report.append(f"  - Median: {self.df['ΔASA'].median():.0f} Ų")
        report.append("")
        
        # Category Distribution
        report.append("3. CATEGORY DISTRIBUTION")
        report.append("-" * 80)
        cat_counts = self.df['Category'].value_counts()
        for cat, count in cat_counts.items():
            report.append(f"{cat}: {count} complexes ({count/len(self.df)*100:.1f}%)")
        report.append("")
        
        # Difficulty Distribution
        report.append("4. DIFFICULTY DISTRIBUTION")
        report.append("-" * 80)
        diff_counts = self.df['Difficulty'].value_counts()
        for diff in ['Rigid-body', 'Medium', 'Difficult']:
            if diff in diff_counts:
                count = diff_counts[diff]
                report.append(f"{diff}: {count} complexes ({count/len(self.df)*100:.1f}%)")
        report.append("")
        
        # Baseline Performance
        report.append("5. BASELINE PERFORMANCE (AlphaRed)")
        report.append("-" * 80)
        report.append("Overall Success Rate: 63% (DockQ > 0.23)")
        report.append("Antibody-Antigen Success Rate: 43%")
        report.append("Interface-pLDDT Threshold: 85")
        report.append("")
        
        # Our Goal
        report.append("6. PROJECT GOALS")
        report.append("-" * 80)
        report.append("Target: Beat AlphaRed's 63% success rate")
        report.append("Focus: Improve antibody-antigen docking (>43%)")
        report.append("Validation: Benchmark 5.5 dataset")
        report.append("Metrics: I-RMSD, ΔASA, DockQ, fnat")
        report.append("Target Conferences: IJCAI, ICML, NeurIPS")
        report.append("")
        
        report.append("="*80)
        report.append("END OF REPORT")
        report.append("="*80)
        
        # Save report
        report_path = self.results_dir / 'benchmark_analysis_report.txt'
        with open(report_path, 'w') as f:
            f.write('\n'.join(report))
        
        print(f"✓ Report saved to: {report_path}")
        
        # Print to console
        print("\n" + '\n'.join(report))
    
    def run_full_analysis(self):
        """Run complete analysis pipeline"""
        print("\n" + "="*80)
        print("STARTING COMPREHENSIVE BENCHMARK ANALYSIS")
        print("="*80)
        
        # Load data
        self.load_data()
        
        # Run analyses
        self.basic_statistics()
        self.analyze_categories()
        self.analyze_difficulty()
        self.analyze_pdb_structure()
        
        # Create visualizations
        self.create_visualizations()
        
        # Generate target list
        self.generate_target_list()
        
        # Create report
        self.create_report()
        
        print("\n" + "="*80)
        print("ANALYSIS COMPLETE!")
        print("="*80)
        print(f"Results saved to: {self.results_dir}")
        print("="*80)

def main():
    """Main function"""
    analyzer = BenchmarkAnalyzer()
    analyzer.run_full_analysis()

if __name__ == "__main__":
    main()
