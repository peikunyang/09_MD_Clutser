import os
import torch

Div = 1.0

def read_par(file_path):
    chg, vdw = [], []
    with open(file_path, 'r') as fr:
        for line in fr:
            lx = line.split()
            chg.append(float(lx[6]))
            vdw.extend([float(lx[7]), float(lx[8])])
    chgn = torch.tensor(chg, dtype=torch.float32)
    vdwn = torch.tensor(vdw, dtype=torch.float32).reshape(-1, 2)
    return chgn, vdwn

def read_trj(file_path):
    crd = []
    count = 0
    with open(file_path, 'r') as fr:
        for line in fr:
            if line.startswith('MODEL'):
                count += 1
            else:
                lx = line.split()
                crd.extend([float(lx[0]), float(lx[1]), float(lx[2])])
    crdn = torch.tensor(crd, dtype=torch.float32).reshape(count, -1, 3)
    return crdn

def read_single_crd(file_path):
    crd = []
    with open(file_path, 'r') as fr:
        for line in fr:
            if line.startswith('ATOM'):
                lx = line.split()
                crd.extend([float(lx[6]), float(lx[7]), float(lx[8])])
    if len(crd) == 0:
        raise ValueError(f"No valid coordinate data found in {file_path}.")
    crdn = torch.tensor(crd, dtype=torch.float32).reshape(-1, 3)
    return crdn

def e_ele(crdl, crdr, chgl, chgr):
    dist = torch.norm(crdl[:, None, :] - crdr[None, :, :], dim=2)
    e = 332 * torch.sum(chgl[:, None] * chgr[None, :] / dist)
    return e.item()

def e_vdw(crdl, crdr, vdwl, vdwr):
    dist = torch.norm(crdl[:, None, :] - crdr[None, :, :], dim=2)
    rij = vdwl[:, 0][:, None] + vdwr[:, 0][None, :]
    sij = torch.sqrt(vdwl[:, 1][:, None] * vdwr[:, 1][None, :])
    e = torch.sum(sij * ((rij / dist) ** 12 - 2 * (rij / dist) ** 6))
    return e.item()

def compute_energy(pdb_lig, pdb_rec):
    chgl, vdwl = read_par(os.path.join("../../../../4_FF/3_ligand/par", pdb_lig))
    crdl = read_trj(os.path.join("../../../../3_cluster/3_ligand/div_%3.1f/pdb"%(Div), pdb_lig))
    chgr, vdwr = read_par(os.path.join("../../../../4_FF/1_pro-lig/1_rec/par", pdb_rec[:4]))
    crdr = read_single_crd(os.path.join("../../../../2_super_CA/1_pro-lig/3_seperate/1_rec/pdb", pdb_rec))
    energies = []
    for i in range(crdl.shape[0]):
        e_ele_value = e_ele(crdl[i], crdr, chgl, chgr)
        e_vdw_value = e_vdw(crdl[i], crdr, vdwl, vdwr)
        e_net = e_ele_value + e_vdw_value
        energies.append((i, min(e_ele_value, 999999), min(e_vdw_value, 999999), min(e_net, 999999)))

    # 取出 e_net 最小的構型
    best_energy = min(energies, key=lambda x: x[3])
    i, e_ele_value, e_vdw_value, e_net = best_energy
    return pdb_rec, i, e_ele_value, e_vdw_value, e_net  # 回傳最小的 e_net 和對應的構型編號

def main():
    pdb_lig_dir = '../../../../3_cluster/3_ligand/div_%3.1f/pdb'%(Div)
    pdb_rec_dir = '../../../../2_super_CA/1_pro-lig/3_seperate/1_rec/pdb/'
    lig_files = sorted(os.listdir(pdb_lig_dir))
    results = []
    for pdb_lig in lig_files:
        pdb_rec = pdb_lig
        result = compute_energy(pdb_lig, pdb_rec)
        results.append(result)

    with open("Energy", "w") as fw:
        for i in range(0, len(results), 3):
            line_results = results[i:i+3]
            for pdb_rec, best_model, e_ele, e_vdw, e_net in line_results:
                fw.write(f"{pdb_rec[:4]:4s} Model {best_model:3d} {e_net:9.2f} {e_ele:9.2f} {e_vdw:9.2f}     ")
            fw.write("\n")

main()

