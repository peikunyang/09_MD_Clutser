import os
import numpy as np
import parmed

def Read_Lig():
  fr=open("ligand_info","r")
  pdbs=[]
  ligs=[]
  for line in fr:
    lx=line.split()
    pdbs.append(lx[0])
    ligs.append(lx[2])
  fr.close()
  return pdbs,ligs

def Parameter(psf, params, pdb, residue_name):
    Atom = []
    Par = []
    for atom in psf.atoms:
        if atom.residue.chain[:3] == 'HET' and atom.residue.name == residue_name:
            atom_type = atom.type.strip()
            atom_charge = atom.charge
            if not atom_type:
                print(f"Warning: Missing atom type for atom '{atom.name}' in residue {atom.residue.name} {atom.residue.number} chain {atom.residue.chain}")
                continue
            try:
                vdw_param = params.atom_types[atom_type]
                vdw_radius = vdw_param.rmin
                vdw_well_depth = vdw_param.epsilon
            except KeyError:
                print(f"Warning: Missing VdW parameters for atom type '{atom_type}' in residue {atom.residue.name} {atom.residue.number} chain {atom.residue.chain}")
                vdw_radius = 'N/A'
                vdw_well_depth = 'N/A'
            pdb_atom = None
            for pdb_atom in pdb.atoms:
                if pdb_atom.name == atom.name and pdb_atom.residue.idx == atom.residue.idx:
                    break
            if pdb_atom is None:
                print(f"Error: Could not find matching atom '{atom.name}' in PDB for residue {atom.residue.name} {atom.residue.number}")
                continue
            x, y, z = pdb_atom.xx, pdb_atom.xy, pdb_atom.xz
            Atom.extend([residue_name, atom.name, atom_type])
            Par.extend([x, y, z, atom_charge, vdw_radius, vdw_well_depth])
    Atomn = np.array(Atom).reshape(-1, 3)
    Parn = np.array(Par).reshape(-1, 6)
    return Atomn, Parn

def Out_Par(lig, Atom, Par):
    output_dir = "par"
    os.makedirs(output_dir, exist_ok=True)
    with open(f"{output_dir}/{lig[:4]}", 'w') as fw:
        for i in range(Atom.shape[0]):
            fw.write(f"{Atom[i][0]:4s} {Atom[i][1]:4s} {Atom[i][2]:7s} {Par[i][0]:8.4f} {Par[i][1]:8.4f} {Par[i][2]:8.4f} {Par[i][3]:7.4f} {Par[i][4]:9.6f} {Par[i][5]:9.6f}\n")

def Search(pdb,pdb_t,lig_t):
    idx = pdb_t.index(pdb)
    residue_name = lig_t[idx]
    return residue_name

def Main():
    pdb_t,lig_t = Read_Lig()
    param_files = ['../../0_par/par_all36m_prot.prm','../../0_par/top_all36_prot.rtf','../../0_par/par_all36_cgenff.prm','../../0_par/top_all36_cgenff.rtf',
        '../../0_par/toppar_water_ions.str']
    lig_files = sorted(os.listdir('../../../2_superimpose/1_pro-lig/3_seperate/2_lig/pdb'))
    params = parmed.charmm.CharmmParameterSet(*param_files)
    for lig in lig_files:
        psf_file = f'../../../1_MD/1_pro-lig/pdb/{lig[:4]}/step3_input.psf'
        pdb_file = f'../../../1_MD/1_pro-lig/pdb/{lig[:4]}/step4_equilibration.pdb'
        psf = parmed.charmm.CharmmPsfFile(psf_file)
        pdb = parmed.load_file(pdb_file)
        residue_name = Search(lig[:4],pdb_t,lig_t)
        Atom, Par = Parameter(psf, params, pdb, residue_name)
        Out_Par(lig, Atom, Par)

Main()

