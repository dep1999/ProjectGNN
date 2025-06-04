import pandas as pd
import os

# ==== 设置部分 ====
# 设置数据文件所在的目录路径
data_dir = r"d:\Googledownload\depmap data"

# 表达矩阵文件路径（使用绝对路径）
expr_path = os.path.join(data_dir, "CCLE_expression.csv")

# 感兴趣的基因（可替换为你自己的列表）
genes_of_interest = ["KRAS", "TP53", "DUSP6", "ERK2", "AKT1"]

# 感兴趣的 DepMap_ID（细胞系 ID 列表）
samples_of_interest = ["ACH-000001", "ACH-000002", "ACH-000003"]

# ==== 加载表达矩阵 ====
print("正在加载表达矩阵...")
expr_df = pd.read_csv(expr_path)
expr_df = expr_df.rename(columns={expr_df.columns[0]: "gene"})  # 第一列为基因名

# ==== 检查存在性 ====
missing_genes = [g for g in genes_of_interest if g not in expr_df["gene"].values]
missing_samples = [s for s in samples_of_interest if s not in expr_df.columns]

print(f"已识别 {len(genes_of_interest) - len(missing_genes)} 个有效基因，{len(samples_of_interest) - len(missing_samples)} 个有效样本")
if missing_genes:
    print("⚠️ 缺失的基因：", missing_genes)
if missing_samples:
    print("⚠️ 缺失的样本列：", missing_samples)

# ==== 子集筛选 ====
expr_subset = expr_df[expr_df["gene"].isin(genes_of_interest)]
expr_subset = expr_subset[["gene"] + [s for s in samples_of_interest if s in expr_df.columns]]

# ==== 保存输出 ====
output_path = os.path.join(data_dir, "gene_expression_subset.csv")
expr_subset.to_csv(output_path, index=False)
print(f"✅ 已保存子集表达矩阵至: {output_path}")
