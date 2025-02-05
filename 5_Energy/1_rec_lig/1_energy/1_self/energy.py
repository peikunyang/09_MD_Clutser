import os
import numpy as np

def Read_Par(pdb, rec):
    with open("../../../../4_FF/1_pro-lig/%s/par/%s" % (rec, pdb), 'r') as fr:
        charges = []
        vdw_params = []
        for line in fr:
            lx = line.split()
            charges.append(float(lx[6]))
            vdw_params.extend([float(lx[7]), float(lx[8])])
    charges_array = np.array(charges).reshape(-1)
    vdw_params_array = np.array(vdw_params).reshape(-1, 2)
    return charges_array, vdw_params_array

def Read_Crd(pdb, rec):
    with open("../../../../2_super_CA/1_pro-lig/3_seperate/%s/pdb/%s" % (rec, pdb), 'r') as fr:
        coords = []
        for line in fr:
            lx = line.split()
            coords.extend([float(lx[6]), float(lx[7]), float(lx[8])])
    coords_array = np.array(coords).reshape(-1, 3)
    return coords_array

def E_ele(coords_ligand, coords_receptor, charges_ligand, charges_receptor):
    e = 0
    for i in range(coords_ligand.shape[0]):
        for j in range(coords_receptor.shape[0]):
            e += 332 * charges_ligand[i] * charges_receptor[j] / np.linalg.norm(coords_ligand[i] - coords_receptor[j])
    return e

def E_vdw(coords_ligand, coords_receptor, vdw_ligand, vdw_receptor):
    e = 0
    for i in range(coords_ligand.shape[0]):
        for j in range(coords_receptor.shape[0]):
            r = np.linalg.norm(coords_ligand[i] - coords_receptor[j])
            rij = vdw_ligand[i][0] + vdw_receptor[j][0]
            sij = np.sqrt(vdw_ligand[i][1] * vdw_receptor[j][1])
            e += sij * ((rij / r) ** 12 - 2 * (rij / r) ** 6)
    return e

def Main():
    ligands = sorted(os.listdir('../../../../2_super_CA/1_pro-lig/3_seperate/2_lig/pdb/'))
    with open("Energy", "w") as fw:
        for ligand in ligands:
            ligand_charges, ligand_vdw = Read_Par(ligand[:4], "2_lig")
            ligand_coords = Read_Crd(ligand, "2_lig")
            receptor_charges, receptor_vdw = Read_Par(ligand[:4], "1_rec")
            receptor_coords = Read_Crd('%s.pdb' % (ligand[:4]), "1_rec")

            e_ele_val = E_ele(ligand_coords, receptor_coords, ligand_charges, receptor_charges)
            e_vdw_val = E_vdw(ligand_coords, receptor_coords, ligand_vdw, receptor_vdw)
            e_net_val = e_ele_val + e_vdw_val

            fw.write("%4s %9.3f %9.3f %9.3f\n" % (ligand[:4], e_net_val, e_ele_val, e_vdw_val))

Main()

