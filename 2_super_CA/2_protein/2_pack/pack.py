import mdtraj as md
import glob
import os

# 找到所有的 .dcd 檔案
dcd_files = glob.glob('../1_dcd/dcd/*.dcd')

# 逐個處理每個 .dcd 檔案
for dcd_file in dcd_files:
    # 加載軌跡文件
    traj = md.load(dcd_file, top='../../../1_MD/2_protein_1/3ixo_group1/step4_equilibration.pdb')
    
    # 使用MDTraj的wrap功能來重新封裝分子
    traj.image_molecules(inplace=True)
    
    # 獲取檔案名稱，並構建輸出檔案名稱
    file_name = os.path.basename(dcd_file)
    output_file = f'dcd/processed_{file_name}'
    
    # 保存處理後的軌跡
    traj.save(output_file)

