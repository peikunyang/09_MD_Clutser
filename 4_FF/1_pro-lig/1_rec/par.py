import os
import numpy as np
import parmed

def Parameter(psf,params,pdb,chain_id,residue_id):
  Atom=[]
  Par=[]
  for atom in psf.atoms:
    if atom.residue.chain==chain_id and atom.residue.number==residue_id:
      atom_type=atom.type.strip()
      atom_charge=atom.charge
      residue_name=atom.residue.name
      if not atom_type:
        print(f"Warning: Missing atom type for atom '{atom.name}' in residue {atom.residue.name} {atom.residue.number} chain {atom.residue.chain}")
        continue
      try:
        vdw_param=params.atom_types[atom_type]
        vdw_radius=vdw_param.rmin
        vdw_well_depth=vdw_param.epsilon
      except KeyError:
        print(f"Warning: Missing VdW parameters for atom type '{atom_type}' in residue {atom.residue.name} {atom.residue.number} chain {atom.residue.chain}")
        vdw_radius='N_A'
        vdw_well_depth='N_A'
      pdb_atoms=pdb.atoms
      pdb_atom=None
      for pdb_atom in pdb_atoms:
        if pdb_atom.name==atom.name and pdb_atom.residue.idx==atom.residue.idx:
          break
      if pdb_atom is None:
        print(f"Error: Could not find matching atom '{atom.name}' in PDB for residue {atom.residue.name} {atom.residue.number}")
        continue
      x,y,z=pdb_atom.xx,pdb_atom.xy,pdb_atom.xz
      Atom.extend([residue_name,atom.name,atom_type])
      Par.extend([x,y,z,atom_charge,vdw_radius,vdw_well_depth])
  Atomn=np.array(Atom).reshape(-1,3)
  Parn=np.array(Par).reshape(-1,6)
  return Atomn,Parn

def Out_Par(pdb,chain_id,residue_id,Atom,Par):
  fw=open("par_res/%s/%s_%s"%(pdb,chain_id[3],(residue_id%100)),'w')
  for i in range (Atom.shape[0]):
    fw.write("%-4s %-4s %-7s %8.4f %8.4f %8.4f %7.4f %9.6f %9.6f\n"%(Atom[i][0],Atom[i][1],Atom[i][2],Par[i][0],Par[i][1],Par[i][2],Par[i][3],Par[i][4],Par[i][5]))
  fw.close()

def Res_ID(pdb):
  residue_id=[8,23,25,27,28,29,30,32,47,48,49,50,81,82,84]
  if ((pdb=="1d4h")|(pdb=="1d4i")|(pdb=="1d4j")|(pdb=="1ebw")|(pdb=="1ebz")|(pdb=="1ec0")|(pdb=="1ec1")|(pdb=="1ec2")|(pdb=="1ec3")|(pdb=="1m0b")|
      (pdb=="1zsf")|(pdb=="1zsr")|(pdb=="2bqv")|(pdb=="2cej")|(pdb=="2cem")|(pdb=="2cen")|(pdb=="2uxz")|(pdb=="2uy0")|(pdb=="2wkz")):
    residue_id=[108,123,125,127,128,129,130,132,147,148,149,150,181,182,184]
  elif ((pdb=="1c70")|(pdb=="1npv")|(pdb=="1npw")|(pdb=="2bb9")|(pdb=="2bbb")):
    residue_id=[208,223,225,227,228,229,230,232,247,248,249,250,281,282,284]
  return residue_id

def Main():
  param_files=['../../par/par_all36m_prot.prm','../../par/top_all36_prot.rtf','../../par/par_all36_cgenff.prm','../../par/top_all36_cgenff.rtf','../../par/toppar_water_ions.str']
  pdb_name=sorted(os.listdir('../../../2_MD/pro-lig/MD/pdb/'))
  chain_id=['PROA','PROB']
  for i in range (len(pdb_name)): # len(pdb_name)
    os.mkdir("par_res/%s"%(pdb_name[i]))
    residue_id=[]
    residue_id.append([8,23,25,27,28,29,30,32,47,48,49,50,81,82,84])
    residue_id.append(Res_ID(pdb_name[i]))
    psf=parmed.charmm.CharmmPsfFile('../../../2_MD/pro-lig/MD/pdb/%s/step3_input.psf'%(pdb_name[i]))
    pdb=parmed.load_file('../../../2_MD/pro-lig/MD/pdb/%s/step4_equilibration.pdb'%(pdb_name[i]))
    params=parmed.charmm.CharmmParameterSet(*param_files)
    for j in range (len(chain_id)): # len(chain_id)
      for k in range (len(residue_id[j])): # len(residue_id)
        Atom,Par=Parameter(psf,params,pdb,chain_id[j],residue_id[j][k])
        Out_Par(pdb_name[i],chain_id[j],residue_id[j][k],Atom,Par)

Main()


