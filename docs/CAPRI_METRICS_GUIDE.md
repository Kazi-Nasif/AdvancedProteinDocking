# 📊 COMPLETE CAPRI METRICS EVALUATION

## What AlphaRed and Other Papers Actually Use:

### **PRIMARY METRIC: DockQ**

AlphaRed and most modern papers use **DockQ** as the main evaluation metric.

**DockQ Formula:**
```
DockQ = (fnat + 1/(1+(I-RMSD/1.5)²) + 1/(1+(L-RMSD/8.5)²)) / 3
```

**DockQ combines 3 metrics:**
1. **I-RMSD** (Interface-RMSD) - What you already calculated! ✓
2. **fnat** (Fraction of native contacts)
3. **L-RMSD** (Ligand-RMSD)

---

## 📈 Standard Evaluation Metrics:

### 1. **DockQ Score** (Main)
- **Range**: 0 to 1
- **> 0.8**: High quality (CAPRI)
- **> 0.49**: Medium quality
- **> 0.23**: Acceptable ← **AlphaRed's success threshold**
- **< 0.23**: Incorrect

### 2. **I-RMSD (Interface-RMSD)** ✓ You have this!
- RMSD of interface residues
- **< 1.0 Å**: High quality
- **< 2.0 Å**: Medium quality  
- **< 4.0 Å**: Acceptable
- Your results: Mean 2.20 Å ✓

### 3. **fnat (Fraction of Native Contacts)**
- Fraction of native residue-residue contacts reproduced
- Contact = Cβ atoms within 5 Å
- **Range**: 0 to 1 (higher is better)

### 4. **L-RMSD (Ligand-RMSD)**
- RMSD of ligand after superimposing receptors
- Measures rigid-body positioning accuracy
- **< 1.0 Å**: Excellent
- **< 5.0 Å**: Good
- **< 10 Å**: Acceptable

### 5. **ΔASA (Delta ASA)**
- **NOT used for evaluation!**
- Only for characterizing interface size
- Papers report it in dataset description, not results

---


