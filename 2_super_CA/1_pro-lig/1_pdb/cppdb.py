import os
import MDAnalysis as mda

pdb=sorted(os.listdir("../../../1_MD/1_pro-lig/pdb"))
for i in range (len(pdb)):
  u = mda.Universe("../../../1_MD/1_pro-lig/pdb/%s/step4_equilibration.pdb"%(pdb[i]))
  selection = u.select_atoms("not resname HOH CLA POT")
  selection.write("pdb/%s.pdb"%(pdb[i]))


