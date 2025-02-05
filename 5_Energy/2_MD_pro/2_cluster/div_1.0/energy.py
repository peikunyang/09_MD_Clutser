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

def read_crd(file_path):
    crd = []
    with open(file_path, 'r') as fr:
        for line in fr:
            lx = line.split()
            crd.extend([float(lx[6]), float(lx[7]), float(lx[8])])
    crdn = torch.tensor(crd, dtype=torch.float32).reshape(-1, 3)
    return crdn

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

def compute_energy(crdr, chgr, vdwr, crdl, chgl, vdwl):
    energies = []
    for i in range(crdr.shape[0]):
        e_ele_value = e_ele(crdl, crdr[i], chgl, chgr)
        e_vdw_value = e_vdw(crdl, crdr[i], vdwl, vdwr)
        e_net = e_ele_value + e_vdw_value
        energies.append((i, min(e_ele_value, 99999), min(e_vdw_value, 99999), min(e_net, 99999)))
    energies_sorted = sorted(energies, key=lambda x: x[3])[:4]
    indices = [e[0] for e in energies_sorted]
    e_elen = [e[1] for e in energies_sorted]
    e_vdwn = [e[2] for e in energies_sorted]
    e_netn = [e[3] for e in energies_sorted]

    return e_netn, e_elen, e_vdwn, indices

def main():
    pdb_lig_list = sorted(os.listdir('../../../../2_super_CA/1_pro-lig/3_seperate/2_lig/pdb'))
    crdr = read_trj(f"../../../../3_cluster/1_protein/div_%3.1f/cluster_%3.1f"%(Div,Div))
    chgr, vdwr = read_par(f"../../../../5_FF/2_protein/par/3ixo_group1")
    with open("Energy", "w") as fw:
      for pdb_lig in pdb_lig_list:
        crdl = read_crd(f"../../../../2_super_CA/1_pro-lig/3_seperate/2_lig/pdb/{pdb_lig}")
        chgl, vdwl = read_par(f"../../../../5_FF/1_pro-lig/2_lig/par/{pdb_lig[:4]}")
        e_netn, e_elen, e_vdwn, indices = compute_energy(crdr, chgr, vdwr, crdl, chgl, vdwl)
        fw.write(f"{pdb_lig[:4]:4s}")
        for i in range(4):
          fw.write(f"  {indices[i]:6d} {e_netn[i]:8.2f} {e_elen[i]:8.2f} {e_vdwn[i]:8.2f}    ")
        fw.write("\n")

main()

