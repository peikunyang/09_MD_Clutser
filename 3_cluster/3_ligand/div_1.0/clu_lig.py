import os
import torch
import time

threshold = 1.0
max_num_trj = 10000  # 控制讀取的模型數量

device = torch.device("cpu")  # 固定為 CPU

def read_trj(lig):
    models = []
    current_model = []
    count = 0
    with open(f'../../../2_super_CA/3_ligand/pdb/%s'%(lig), 'r') as fr:
        for line in fr:
            if line.startswith('MODEL'):
                if current_model:
                    models.append(current_model)
                current_model = []
                count += 1
            elif line.startswith('ATOM') and count <= max_num_trj:
                lx = line.split()
                current_model.extend([float(lx[6]), float(lx[7]), float(lx[8])])
        if current_model:  # 添加最後一個模型
            models.append(current_model)
    return models

def calculate_rmsd(crd1, crd2):
    differences = crd1 - crd2
    squared_differences = differences ** 2
    mean_squared_difference = torch.mean(squared_differences)
    rmsd = torch.sqrt(mean_squared_difference)
    return rmsd

def cluster(clu_model, model, threshold):
    if not clu_model:
        clu_model.append(model[0])
    for i in range(len(model)):  # 直接跳過第0個，因為它已加入
        crd1 = torch.tensor(model[i], device=device).reshape(-1, 3)  # 使用 CPU 進行運算
        add_to_cluster = True
        for j in range(len(clu_model)):
            crd2 = torch.tensor(clu_model[j], device=device).reshape(-1, 3)  # 使用 CPU 進行運算
            rmsd = calculate_rmsd(crd1, crd2)
            if rmsd < threshold:
                add_to_cluster = False
                break
        if add_to_cluster:
            clu_model.append(model[i])
        del crd1  # 清理記憶體

def Out_crd(clu_model, ligand):
    with open(f'pdb/%s' % (ligand), 'w') as fw:  # 用 `with` 確保文件正確關閉
        for i in range(len(clu_model)):
            fw.write('MODEL %6d\n' % i)
            for j in range(0, len(clu_model[i]), 3):  # 每3個值寫入一次
                fw.write('%7.3f %7.3f %7.3f\n' % (clu_model[i][j], clu_model[i][j + 1], clu_model[i][j + 2]))

def main():
    start_time = time.time()  # 開始計時
    ligands = sorted(os.listdir('../../../2_super_CA/3_ligand/pdb'))
    for ligand in ligands:
        clu_model = []
        models = read_trj(ligand)
        cluster(clu_model, models, threshold)
        Out_crd(clu_model, ligand)
    end_time = time.time()  # 結束計時
    total_time = end_time - start_time
    print(f"Total execution time: {total_time:.4f} seconds", flush=True)  # 顯示總執行時間

main()

