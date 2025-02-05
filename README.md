# Leveraging Conformational Diversity for Enhanced Structure-Based Virtual Screening: Insights from Molecular Dynamics Simulations of HIV-1 Protease–Ligand Complexes

This repository contains molecular dynamics (MD) simulation data and computational scripts for analyzing the impact of conformational diversity on structure-based virtual screening (SBVS). The study focuses on HIV-1 protease–ligand complexes, exploring protein and ligand flexibility using MD simulations and evaluating their effects on docking accuracy.

## 📁 Directory Structure

## 📂 Directory Structure

```
📂 1_MD/ - Molecular dynamics simulations
├── 📂 1_pro-lig/ - MD simulations of protein-ligand complexes (subset included)
├── 📂 2_protein/ - MD simulations of apo protein (PDB ID: 3IXO)
└── 📂 3_ligand/ - MD simulations of ligands from protein-ligand complexes (subset included)

📂 2_super_CA/ - Superimposition of protein and ligand conformations
└── 📂 superimpose/ - Aligns 1_pro-lig, 2_protein, and 3_ligand

📂 3_cluster/ - Clustering protein and ligand conformations
├── 📂 cluster_protein/ - Clustering apo protein conformations
└── 📂 cluster_ligand/ - Clustering ligand conformations

📂 4_FF/ - Force field parameters for protein-ligand interaction calculations

📂 5_Energy/ - Computation of electrostatic and van der Waals interaction energies
├── 📂 1_rec_lig/ - Energy calculations using original protein-ligand complexes
├── 📂 2_MD_pro/ - Energy calculations using MD-generated protein conformations
├── 📂 3_MD_lig/ - Energy calculations using MD-generated ligand conformations
└── 📂 4_MD_pro_lig/ - Energy calculations using MD-generated protein-ligand complexes


## 📌 Project Overview

### 🔬 Background
- Structure-based virtual screening (SBVS) often assumes a fixed protein conformation, limiting its accuracy in identifying potential ligands.
- This study investigates how incorporating multiple protein and ligand conformations improves SBVS accuracy.
- Molecular dynamics (MD) simulations were performed to generate a diverse set of conformations for HIV-1 protease and its ligands.
- Force field-based energy calculations were conducted to evaluate binding free energies.

### 🛠 Methods
1. **Molecular Dynamics Simulations**: Conducted MD simulations for protein-ligand complexes, apo proteins, and ligands.
2. **Superimposition**: Aligned protein and ligand structures to ensure structural consistency.
3. **Clustering**: Used RMSD-based clustering to reduce redundant conformations while preserving structural diversity.
4. **Force Field Preparation**: Generated force field parameters for interaction energy calculations.
5. **Energy Calculations**: Computed electrostatic (Eele) and van der Waals (Evdw) interaction energies.

## 🚀 Usage
To analyze the provided datasets or run additional simulations:

1. Clone the repository:
   ```bash
   git clone https://github.com/your-repo-name.git
   cd your-repo-name
