import os
import re

def create_directory_if_not_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def merge_pdb_files(input_directories, output_file):
    file_written = False  # 追蹤是否有成功寫入任何內容
    with open(output_file, 'w') as outfile:
        for directory in input_directories:
            pdb_file_path = os.path.join(directory, 'superimposed.pdb')
            if os.path.exists(pdb_file_path):
                print(f"Processing {pdb_file_path}...")
                with open(pdb_file_path, 'r') as infile:
                    for line in infile:
                        if not line.startswith('TER'):
                            outfile.write(line)
                            file_written = True
            else:
                print(f"Warning: {pdb_file_path} not found!")
    
    if not file_written:
        print(f"No content written to {output_file}, check if source files are empty or missing.")

def extract_number(s):
    # 提取目錄名稱中的數字部分
    match = re.search(r'\d+', s)
    return int(match.group()) if match else float('inf')

def main():
    # 更新資料來源的根目錄
    input_base_directory = '../3_super/rmsd'
    output_directory = './pdb'
    
    # 確保 'pdb' 資料夾存在
    create_directory_if_not_exists(output_directory)
    
    # 生成 'processed_dcd_1' 到 'processed_dcd_100' 的完整路徑
    subdirectories = [os.path.join(input_base_directory, f"processed_dcd_{i}") for i in range(1, 11)]
    
    # 每 10 個子目錄合併到一個檔案中
    for i in range(0, len(subdirectories), 10):
        batch_subdirectories = subdirectories[i:i + 10]
        output_file = os.path.join(output_directory, f'3ixo_group1_{i//10 + 1}.pdb')
        merge_pdb_files(batch_subdirectories, output_file)

if __name__ == "__main__":
    main()

