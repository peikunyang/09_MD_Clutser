import os
import numpy as np
import parmed

def Parameter(psf, params, pdb):
    Atom = []
    Par = []
    for atom in psf.atoms:
        if atom.residue.chain[:3] == 'HET':
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
            Atom.extend([atom.residue.name, atom.name, atom_type])
            Par.extend([x, y, z, atom_charge, vdw_radius, vdw_well_depth])
    Atomn = np.array(Atom).reshape(-1, 3)
    Parn = np.array(Par).reshape(-1, 6)
    return Atomn, Parn

def Out_Par(lig, Atom, Par):
    output_dir = "par"
    os.makedirs(output_dir, exist_ok=True)
    with open(f"{output_dir}/{lig}", 'w') as fw:
        for i in range(Atom.shape[0]):
            fw.write(f"{Atom[i][0]:4s} {Atom[i][1]:4s} {Atom[i][2]:7s} {Par[i][0]:8.4f} {Par[i][1]:8.4f} {Par[i][2]:8.4f} {Par[i][3]:7.4f} {Par[i][4]:9.6f} {Par[i][5]:9.6f}\n")

def Main():
    param_files = [
        '../par/par_all36m_prot.prm',
        '../par/top_all36_prot.rtf',
        '../par/par_all36_cgenff.prm',
        '../par/top_all36_cgenff.rtf',
        '../par/toppar_water_ions.str'
    ]
    lig_files = sorted(os.listdir('../../3_superimpose/lig/super_MD1/pdb'))
    params = parmed.charmm.CharmmParameterSet(*param_files)
    for lig in lig_files:
        psf_file = f'../../2_MD/ligand/MD1/pdb/{lig[:4]}/step3_input.psf'
        pdb_file = f'../../2_MD/ligand/MD1/pdb/{lig[:4]}/step4_equilibration.pdb'
        if not os.path.exists(psf_file) or not os.path.exists(pdb_file):
          print(f"Error: PSF file or PDB file does not exist for {lig}")
          continue
        psf = parmed.charmm.CharmmPsfFile(psf_file)
        pdb = parmed.load_file(pdb_file)
        Atom, Par = Parameter(psf, params, pdb)
        Out_Par(lig, Atom, Par)

Main()

