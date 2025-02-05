import MDAnalysis as mda
from MDAnalysis.coordinates.DCD import DCDWriter
import os

input_dcd = '../../../1_MD/2_protein_1/3ixo_group1/step5.dcd'
output_prefix = 'dcd'
n_parts = 10

# 確保輸出目錄存在
if not os.path.exists(output_prefix):
    os.makedirs(output_prefix)

# 讀取DCD檔案
try:
    u = mda.Universe('../../../1_MD/2_protein_1/3ixo_group1/step5.pdb', input_dcd)
except FileNotFoundError as e:
    print(f"Error: {e}")
    exit(1)

# 確認已成功讀取文件
print(f"Total frames in the input DCD file: {len(u.trajectory)}")

# 計算每部分包含的frames數
n_frames = len(u.trajectory)
frames_per_part = n_frames // n_parts

# 檢查是否可以均勻分配frames
if n_frames % n_parts != 0:
    print(f"Warning: The total number of frames {n_frames} is not evenly divisible by {n_parts}. "
          f"The last part will contain {n_frames % frames_per_part} extra frames.")

# 拆分並寫入新的DCD檔案
for i in range(n_parts):
    start_frame = i * frames_per_part
    end_frame = (i + 1) * frames_per_part if i < n_parts - 1 else n_frames
    output_dcd = f'{output_prefix}/{output_prefix}_{i + 1}.dcd'

    print(f"Writing frames {start_frame} to {end_frame} into {output_dcd}")

    with DCDWriter(output_dcd, n_atoms=u.atoms.n_atoms) as writer:
        for ts in u.trajectory[start_frame:end_frame]:
            writer.write(u.atoms)

print("Splitting complete.")

