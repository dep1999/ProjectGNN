import pandas as pd
from scipy.stats import ttest_ind
import numpy as np
import os
import time

# 添加计时功能
start_time = time.time()

# 设置数据文件所在的目录路径
data_dir = r"d:\Googledownload\depmap data"

# 读取表达矩阵（CCLE_expression.csv）
print("正在读取表达矩阵...")
expr_df = pd.read_csv(os.path.join(data_dir, "CCLE_expression.csv"), index_col=0)
expr_df = expr_df.T  # 行是 DepMap_ID
print(f"表达矩阵形状: {expr_df.shape}")

# 读取 metadata 和 mutation 信息
print("正在读取样本信息和突变数据...")
sample_info = pd.read_csv(os.path.join(data_dir, "sample_info.csv"))
mut_df = pd.read_csv(os.path.join(data_dir, "CCLE_mutations.csv"), low_memory=False)
print(f"突变数据形状: {mut_df.shape}")

# 查看可能的突变注释字段，优先使用 Protein_Change
mut_kras_g12c = mut_df[
    (mut_df["Hugo_Symbol"] == "KRAS") &
    (mut_df["Protein_Change"].str.contains("G12C", na=False))
]
kras_g12c_ids = mut_kras_g12c["DepMap_ID"].unique().tolist()

# 与表达数据对齐
kras_g12c_ids = list(set(kras_g12c_ids).intersection(set(expr_df.index)))
other_ids = list(set(expr_df.index) - set(kras_g12c_ids))

print(f"表达数据中找到 {len(kras_g12c_ids)} 个 KRAS G12C 阳性样本")
print(f"对照组样本数量: {len(other_ids)}")

# 差异表达分析
print("开始差异表达分析...")
results = []
total_genes = len(expr_df.columns)

# 添加进度显示
for i, gene in enumerate(expr_df.columns):
    if i % 1000 == 0:
        elapsed = time.time() - start_time
        print(f"处理中: {i}/{total_genes} 基因 ({i/total_genes*100:.1f}%), 已用时间: {elapsed:.1f}秒")
    
    expr_kras = expr_df.loc[kras_g12c_ids, gene].dropna()
    expr_other = expr_df.loc[other_ids, gene].dropna()

    if len(expr_kras) >= 3 and len(expr_other) >= 3:
        stat, pval = ttest_ind(expr_kras, expr_other, equal_var=False)
        lfc = np.log2(expr_kras.mean() + 1) - np.log2(expr_other.mean() + 1)
        results.append((gene, lfc, pval))
    else:
        results.append((gene, np.nan, np.nan))

# 保存结果表格
df = pd.DataFrame(results, columns=["gene", "RNAseq_LFC", "RNAseq_pval"])
df.to_csv(os.path.join(data_dir, "KRAS_G12C_RNAseq_differential_expression.csv"), index=False)

total_time = time.time() - start_time
print(f"差异分析完成，总用时: {total_time:.1f}秒")
print(f"结果已保存为 {os.path.join(data_dir, 'KRAS_G12C_RNAseq_differential_expression.csv')}")
