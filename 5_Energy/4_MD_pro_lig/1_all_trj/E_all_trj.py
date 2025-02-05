import os
import torch

epilson = 1e-9
div_rec = 0.0
div_lig = 0.0

def read_par(file_path):
    chg, vdw = [], []
    with open(file_path, 'r') as fr:
        for line in fr:
            lx = line.split()
            chg.append(float(lx[6]))
            vdw.extend([float(lx[7]), float(lx[8])])
    chgn = torch.tensor(chg, dtype=torch.float32).to('cuda')
    vdwn = torch.tensor(vdw, dtype=torch.float32).reshape(-1, 2).to('cuda')
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
    crdn = torch.tensor(crd, dtype=torch.float32).reshape(count, -1, 3).to('cuda')
    return crdn

def e_ele(crdl, crdr, chgl, chgr, dist):
    e = 332 * torch.sum(chgl[None, :, None, None] * chgr[None, None, None, :] / dist, dim=(1, 3))
    return e

def e_vdw(crdl, crdr, vdwl, vdwr, dist):
    rij = vdwl[:, 0][None, :, None, None] + vdwr[:, 0][None, None, None, :]
    sij = torch.sqrt(vdwl[:, 1][None, :, None, None] * vdwr[:, 1][None, None, None, :])
    e = torch.sum(sij * ((rij / dist) ** 12 - 2 * (rij / dist) ** 6), dim=(1, 3))
    return e

def compute_energy(lig, chgr, vdwr, crdr):
    chgl, vdwl = read_par('../../../4_FF/3_ligand/par/%s' % (lig))
    crdl = read_trj('../../../3_cluster/3_ligand/div_%3.1f/pdb/%s' % (div_lig, lig))
    
    batch_size = 50  # 設定批次大小以降低 GPU 記憶體使用量
    num_batches_crdl = (crdl.shape[0] + batch_size - 1) // batch_size
    num_batches_crdr = (crdr.shape[0] + batch_size - 1) // batch_size
    e_net_min = float('inf')
    global_min_idx_crdl, global_min_idx_crdr = -1, -1
    e_ele_min, e_vdw_min = 99999, 99999

    for batch_idx_crdl in range(num_batches_crdl):
        crdl_batch = crdl[batch_idx_crdl * batch_size:(batch_idx_crdl + 1) * batch_size, :, :]
        
        for batch_idx_crdr in range(num_batches_crdr):
            crdr_batch = crdr[batch_idx_crdr * batch_size:(batch_idx_crdr + 1) * batch_size, :, :]
            
            dist = torch.norm(crdl_batch[:, :, None, None, :] - crdr_batch[None, None, :, :, :], dim=4) + epilson
            e_ele_value = e_ele(crdl_batch, crdr_batch, chgl, chgr, dist)
            e_vdw_value = e_vdw(crdl_batch, crdr_batch, vdwl, vdwr, dist)
            e_net = e_ele_value + e_vdw_value
            
            # 找到 e_net 中的極小值及其索引
            min_val, min_idx = torch.min(e_net.view(-1), dim=0)
            if min_val.item() < e_net_min:
                e_net_min = min_val.item()
                local_min_idx_crdl, local_min_idx_crdr = torch.unravel_index(min_idx, e_net.shape)
                global_min_idx_crdl = batch_idx_crdl * batch_size + local_min_idx_crdl.item()
                global_min_idx_crdr = batch_idx_crdr * batch_size + local_min_idx_crdr.item()
                e_ele_min = e_ele_value[local_min_idx_crdl, local_min_idx_crdr].item()
                e_vdw_min = e_vdw_value[local_min_idx_crdl, local_min_idx_crdr].item()
            
            # 清除不再需要的變數以釋放 GPU 記憶體
            del crdr_batch, dist, e_ele_value, e_vdw_value, e_net
            torch.cuda.empty_cache()
        
        del crdl_batch
        torch.cuda.empty_cache()

    # 輸出結果到檔案，改為追加模式
    with open("Energy", "a") as fw:
        fw.write(f"{lig[:4]} {global_min_idx_crdr:6d} {global_min_idx_crdl:6d}  {e_net_min:9.2f} {e_ele_min:9.2f} {e_vdw_min:9.2f}\n")

def main():
    chgr, vdwr = read_par('../../../4_FF/2_protein/par/3ixo_group1')
    crdr = read_trj('../../../3_cluster/2_protein/div_%3.1f/cluster_%3.1f' % (div_rec, div_rec))
    lig_files = sorted(os.listdir('../../../3_cluster/3_ligand/div_%3.1f/pdb/'%(div_lig)))
    for lig in lig_files:
        compute_energy(lig, chgr, vdwr, crdr)

main()

