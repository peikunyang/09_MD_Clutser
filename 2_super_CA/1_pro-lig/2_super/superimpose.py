import os
import numpy as np
import torch
from Bio.PDB import PDBParser, Superimposer, PDBIO

def get_residues(structure, chain_id, residue_ids):
  residues = []
  for model in structure:
    for i in range (len(chain_id)):
      chain = model[chain_id[i]]
      for residue in chain:
        if residue.id[1] in residue_ids:
          residues.append(residue)
  return residues

def extract_ca_atoms(residues):
  ca_atoms = []
  for residue in residues:
    for atom in residue:
      if atom.name == 'CA':
        ca_atoms.append(atom)
  return ca_atoms

def Super(fw,pdb,res):
  structure1 = parser.get_structure('protein1','../1_pdb/3ixo_group1/step4_equilibration.pdb')
  structure2 = parser.get_structure('protein2','../1_pdb/pdb/%s'%(pdb))
  chain_id =["A","B"]
  residues1 = get_residues(structure1, chain_id, res)
  residues2 = get_residues(structure2, chain_id, res)
  ca_atoms1 = extract_ca_atoms(residues1)
  ca_atoms2 = extract_ca_atoms(residues2)
  super_imposer = Superimposer()
  super_imposer.set_atoms(ca_atoms1, ca_atoms2)
  super_imposer.apply(structure2.get_atoms())
  fw.write("%4s %6.3f\n"%(pdb[:4],super_imposer.rms))
  io = PDBIO()
  io.set_structure(structure2)
  io.save('pdb/%s'%(pdb))

def Main():
  global parser
  parser = PDBParser(QUIET=True)
  res=[8,23,25,27,28,29,30,32,47,48,49,50,81,82,84]
  pdb=sorted(os.listdir("../1_pdb/pdb"))
  fw=open("rmsd","w")
  for i in range (len(pdb)): # len(pdb)
    Super(fw,pdb[i],res)
  fw.close()

Main()



