import os
import torch

def read_par(file_path, device):
    chg, vdw = [], []
    with open(file_path, 'r') as fr:
        for line in fr:
            lx = line.split()
            if len(lx) >= 9:  # 確保有足夠的元素
                chg.append(float(lx[6]))
                vdw.extend([float(lx[7]), float(lx[8])])
    chgn = torch.tensor(chg, dtype=torch.float32).to(device)
    vdwn = torch.tensor(vdw, dtype=torch.float32).reshape(-1, 2).to(device)
    return chgn, vdwn

def read_crd(file_path, device):
    crd = []
    with open(file_path, 'r') as fr:
        for line in fr:
            lx = line.split()
            if len(lx) >= 9:  # 確保有足夠的元素
                crd.extend([float(lx[6]), float(lx[7]), float(lx[8])])
    crdn = torch.tensor(crd, dtype=torch.float32).reshape(-1, 3).to(device)
    return crdn

def read_trj(file_path, device):
    crd = []
    count = 0
    with open(file_path, 'r') as fr:
        for line in fr:
            if line.startswith('MODEL'):
                count += 1
            if line.startswith('ATOM'):
                lx = line.split()
                crd.extend([float(lx[6]), float(lx[7]), float(lx[8])])
    crdn = torch.tensor(crd, dtype=torch.float32).reshape(count, -1, 3).to(device)
    return crdn

def e_ele(crdl, crdr, chgl, chgr, epsilon=1e-6):
    with torch.no_grad():  # 避免追蹤計算圖，節省記憶體
        dist = torch.norm(crdl[:, None, :] - crdr[None, :, :], dim=2) + epsilon
        e = 332 * torch.sum(chgl[:, None] * chgr[None, :] / dist)
    return e.item()

def e_vdw(crdl, crdr, vdwl, vdwr, epsilon=1e-6):
    with torch.no_grad():  # 避免追蹤計算圖，節省記憶體
        dist = torch.norm(crdl[:, None, :] - crdr[None, :, :], dim=2) + epsilon
        rij = vdwl[:, 0][:, None] + vdwr[:, 0][None, :]
        sij = torch.sqrt(vdwl[:, 1][:, None] * vdwr[:, 1][None, :])
        e = torch.sum(sij * ((rij / dist) ** 12 - 2 * (rij / dist) ** 6))
    return e.item()

def compute_energy(pdb_lig, pdb_rec, chgr, vdwr, crdr, device):
    chgl, vdwl = read_par(f"../../../5_FF/1_pro-lig/2_lig/par/{pdb_lig[:4]}", device)
    crdl = read_crd(f"../../../2_super_CA/1_pro-lig/3_seperate/2_lig/pdb/{pdb_lig}", device)
    energies = []

    for i in range(crdr.shape[0]):
        e_ele_value = e_ele(crdl, crdr[i], chgl, chgr)
        e_vdw_value = e_vdw(crdl, crdr[i], vdwl, vdwr)
        e_net = e_ele_value + e_vdw_value
        energies.append((i, min(e_ele_value, 99999), min(e_vdw_value, 99999), min(e_net, 99999)))

        # 清理 GPU 記憶體，避免累積
        torch.cuda.empty_cache()

    energies_sorted = sorted(energies, key=lambda x: x[3])[:4]

    indices = [e[0] for e in energies_sorted]
    e_elen = [e[1] for e in energies_sorted]
    e_vdwn = [e[2] for e in energies_sorted]
    e_netn = [e[3] for e in energies_sorted]

    return pdb_rec, pdb_lig, e_netn, e_elen, e_vdwn, indices

def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    pdb_rec_dir = '../../../2_super_CA/2_protein/4_comb/pdb'
    rec_files = sorted(os.listdir(pdb_rec_dir))
    pdb_lig_list = sorted(os.listdir('../../../2_super_CA/1_pro-lig/3_seperate/2_lig/pdb'))

    for pdb_rec in rec_files:
        chgr, vdwr = read_par(f"../../../5_FF/2_protein/par/{pdb_rec[:11]}", device)
        crdr = read_trj(f"../../../2_super_CA/2_protein/4_comb/pdb/{pdb_rec}", device)

        for pdb_lig in pdb_lig_list:
            print(pdb_lig, flush=True)
            if os.path.isfile(os.path.join(pdb_rec_dir, pdb_rec)):
                result = compute_energy(pdb_lig, pdb_rec, chgr, vdwr, crdr, device)

                with open("energy", "a") as fw:  # 改為 "a" 模式
                    pdb_rec, pdb_lig, e_netn, e_elen, e_vdwn, indices = result
                    fw.write(f"{pdb_rec[:-4]:14s} {pdb_lig[:-4]:<4s}")
                    for i in range(1):
                        fw.write(f"  {indices[i]:6d} {e_netn[i]:8.2f} {e_elen[i]:8.2f} {e_vdwn[i]:8.2f}    ")
                    fw.write("\n")

            # 清理 GPU 記憶體，避免累積
            torch.cuda.empty_cache()

if __name__ == "__main__":
    main()

