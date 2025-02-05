import torch
import time

threshold = 1.0
max_num_trj = 100000
batch_size = 100000  # 每次處理最多 10,000 個軌跡

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

def read_trj(i):
    models = []
    current_model = []
    count = 0
    with open(f'../../../2_super_CA/2_protein/4_comb/pdb/3ixo_group1_{i}.pdb', 'r') as fr:
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

def calculate_rmsd_tensor(crd1, crd2):
    # crd1 和 crd2 是兩個 Tensor，形狀為 (N, 3)
    differences = crd1 - crd2
    squared_differences = differences ** 2
    mean_squared_difference = torch.mean(squared_differences, dim=(1, 2))  # 沿每個模型計算平均
    rmsd = torch.sqrt(mean_squared_difference)
    return rmsd

def cluster(clu_model, model, threshold):
    if not clu_model:
        clu_model.append(model[0])

    # 將 model 和 clu_model 轉換成張量，並移動到 GPU 上
    model_tensor = torch.tensor(model, device=device).reshape(len(model), -1, 3)
    clu_tensor = torch.tensor(clu_model, device=device).reshape(len(clu_model), -1, 3)

    # 將 model_tensor 分成多個 batch 來處理
    for i in range(0, model_tensor.shape[0], batch_size):  # 從 0 開始
        batch_end = min(i + batch_size, model_tensor.shape[0])  # 確保不超出範圍
        crd1_batch = model_tensor[i:batch_end]  # 取出一個批次的 crd1

        for j in range(crd1_batch.shape[0]):  # 對每個 batch 中的模型進行處理
            crd1 = crd1_batch[j].unsqueeze(0).expand(clu_tensor.shape[0], -1, -1)  # 廣播 crd1

            # 平行計算 RMSD
            rmsd_tensor = calculate_rmsd_tensor(clu_tensor, crd1)
            
            # 判斷是否所有 RMSD 都超過閾值
            if torch.all(rmsd_tensor >= threshold):
                clu_model.append(model[i + j])
                # 將新的模型加入 clu_tensor，並更新 clu_tensor
                clu_tensor = torch.cat([clu_tensor, crd1[:1]], dim=0)

def Out_crd(clu_model):
    with open(f'cluster_{threshold:.1f}', 'w') as fw:  # 用 `with` 確保文件正確關閉
        for i in range(len(clu_model)):
            fw.write('MODEL %6d\n' % i)
            for j in range(0, len(clu_model[i]), 3):  # 每3個值寫入一次
                fw.write('%7.3f %7.3f %7.3f\n' % (clu_model[i][j], clu_model[i][j+1], clu_model[i][j+2]))

def main():
    start_time = time.time()  # 開始計時
    clu_model = []  # 將 clu_model 的初始化放在迴圈外部，累積所有檔案的模型
    for i in range(1, 2):  # 讀取多個檔案
        models = read_trj(i)
        cluster(clu_model, models, threshold)
        print(f"Total clustered models after file {i}: {len(clu_model)}", flush=True)
    Out_crd(clu_model)
    end_time = time.time()  # 結束計時
    total_time = end_time - start_time
    print(f"Total execution time: {total_time:.4f} seconds")  # 顯示總執行時間

main()

