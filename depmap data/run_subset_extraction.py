
import pandas as pd
import os

# 设置数据文件所在的目录路径
data_dir = r"d:\Googledownload\depmap data"

# ==== 设置 ====
# 表达矩阵文件路径（CCLE_expression.csv）
expr_path = os.path.join(data_dir, "CCLE_expression.csv")

# 感兴趣的基因（可替换为你自己的列表）
genes_of_interest = ["KRAS", "TP53", "DUSP6", "ERK2", "AKT1"]

# 感兴趣的 DepMap_ID（细胞系）
samples_of_interest = ["ACH-000001", "ACH-000002", "ACH-000003"]

# ==== 读取 ====
# 加载表达矩阵（行为基因，列为 DepMap_ID）
expr_df = pd.read_csv(expr_path)
expr_df = expr_df.rename(columns={expr_df.columns[0]: "gene"})  # 第一列为基因名

# 筛选基因
expr_df = expr_df[expr_df["gene"].isin(genes_of_interest)]

# 设置 gene 为索引
expr_df.set_index("gene", inplace=True)

# 筛选样本（列）
valid_cols = [col for col in samples_of_interest if col in expr_df.columns]
expr_subset = expr_df[valid_cols]

# ==== 输出 ====
output_path = os.path.join(data_dir, "gene_expression_subset.csv")
expr_subset.to_csv(output_path)
print(f"已保存 {len(genes_of_interest)} 个基因在 {len(valid_cols)} 个细胞系中的表达矩阵")
