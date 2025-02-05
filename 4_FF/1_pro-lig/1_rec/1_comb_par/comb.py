import os

pdbs=sorted(os.listdir("../par_res"))
chains = ['A', 'B']
residues = [8, 23, 25, 27, 28, 29, 30, 32, 47, 48, 49, 50, 81, 82, 84]

fw=open("exe_exe","w")
for pdb in pdbs:
  fw.write("cat ")
  for chain in chains:
    for residue in residues:
      fw.write("../par_res/%s/%s_%s "%(pdb,chain,residue))
  fw.write("> ../par/%s\n"%(pdb))
fw.close()


