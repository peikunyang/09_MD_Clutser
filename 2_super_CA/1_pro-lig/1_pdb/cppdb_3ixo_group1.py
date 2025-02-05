import os
import MDAnalysis as mda

u = mda.Universe("../../../1_MD/2_protein/3ixo_group1/step4_equilibration.pdb")
selection = u.select_atoms("not resname HOH CLA POT")
selection.write("3ixo_group1/step4_equilibration.pdb")

