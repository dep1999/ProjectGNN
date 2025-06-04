import pandas as pd
import os

# 自动切换工作目录到脚本所在位置
os.chdir(os.path.dirname(os.path.abspath(__file__)))
# 设置数据文件所在的目录路径
data_dir = r"d:\Googledownload\depmap data"


# ==== 设置 ====
# 表达矩阵路径
expr_path = "CCLE_expression.csv"

# Gene symbol 列表（用户定义）
symbols_of_interest = ["KRAS", "TP53", "DUSP6", "ERK2", "AKT1"]

# ==== 加载表达矩阵（预览） ====
expr_df = pd.read_csv(expr_path, nrows=5)
ensembl_ids_in_expr = expr_df.iloc[:, 0].tolist()  # 第一列是 Ensembl ID

# ==== 加载 symbol ↔ Ensembl ID 映射表 ====
# 建议使用 BioMart 或 HGNC 导出的 TSV
# 这里只是一个示例映射，可以换成你自己的大表
symbol_map = {
    "KRAS": "ENSG00000133703",
    "TP53": "ENSG00000141510",
    "DUSP6": "ENSG00000139318",
    "ERK2": "ENSG00000100030",
    "AKT1": "ENSG00000142208"
}

# ==== 映射 Symbol → Ensembl ID 并过滤有效基因 ====
mapped_ids = [symbol_map[s] for s in symbols_of_interest if s in symbol_map]
valid_ids = [eid for eid in mapped_ids if eid in ensembl_ids_in_expr]

print("用户请求的 symbol：", symbols_of_interest)
print("映射后的 Ensembl ID：", mapped_ids)
print("表达矩阵中存在的有效 Ensembl ID：", valid_ids)

# ==== 可保存结果或用于后续提取 ====
with open("mapped_ensembl_ids.txt", "w") as f:
    for eid in valid_ids:
        f.write(eid + "\n")

print("✅ 有效 Ensembl ID 已保存至 mapped_ensembl_ids.txt")
