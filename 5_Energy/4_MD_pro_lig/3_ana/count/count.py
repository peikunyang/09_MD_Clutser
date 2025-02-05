#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# rec 列表中的每個數值代表要處理檔案的後綴
rec = [0.4, 0.6, 0.8, 1.0, 1.2, 1.4]

def count_negative_values(input_file):
    """
    讀取 input_file 檔案，計算每一個欄位中小於 0 的數值個數，
    回傳 (header, count_list)：
        header: 原始檔案第一行的欄位名稱列表
        count_list: 每個欄位中小於 0 的數值個數列表
    若檔案不存在則回傳 (None, None)
    """
    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"找不到檔案: {input_file}")
        return None, None

    # 取得第一行作為 header
    header = lines[0].strip().split()
    count_list = [0] * len(header)

    # 從第二行開始計算每個欄位中小於 0 的數值
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue  # 忽略空行
        tokens = line.split()
        for i, token in enumerate(tokens):
            try:
                num = float(token)
                if num < 0:
                    count_list[i] += 1
            except ValueError:
                # 若無法轉換成數字則跳過
                continue
    return header, count_list

def main():
    # 用來存放各個檔案的結果，key 為 rec 值，value 為 count_list
    results = {}
    header_global = None

    for div in rec:
        # 根據 rec 數值格式化 input 檔案名稱，例如: ../Energy/Energy_0.5
        input_file = f"../Energy/Energy_{div:3.1f}"
        head, count_list = count_negative_values(input_file)
        if head is None:
            # 若檔案不存在則略過
            continue
        # 假設各檔案的 header 都相同，將第一個成功讀取的 header 記錄下來
        if header_global is None:
            header_global = head
        else:
            if header_global != head:
                print(f"警告: {input_file} 的 header 與其他檔案不一致")
        results[div] = count_list

    # 若沒有讀到任何結果則結束程式
    if not results or header_global is None:
        print("沒有讀取到有效的檔案結果")
        return

    # 輸出合併結果到一個檔案中
    output_file = "count_neg_E.txt"
    # 構造輸出的表格，第一列為 header：第一欄顯示「Column」，後續欄位為各 rec 值
    with open(output_file, "w", encoding="utf-8") as f_out:
        # 依 rec 值排序（例如 0.5, 0.6, ...）
        rec_sorted = sorted(results.keys())
        header_line = "\t" + "\t".join(f"rec_{div:3.1f}" for div in rec_sorted)
        f_out.write(header_line + "\n")
        
        # 對於每個原始檔案的欄位，輸出該欄位在各個檔案中的計數結果
        for i, col_name in enumerate(header_global):
            line = col_name[4:]
            for div in rec_sorted:
                count = results[div][i]
                line += "\t" + str(count)
            f_out.write(line + "\n")
    print(f"所有檔案處理完成，結果合併輸出至 {output_file}")

if __name__ == '__main__':
    main()

