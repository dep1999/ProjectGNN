#功能：
#自动识别 KRAS G12C 阳性 vs 其他样本
#针对每个基因计算：log2 fold change，p 值，全样本均值，全样本标准差
#输出结果表格：gene，log2 fold change，p 值，全样本均值，全样本标准差
#最终输出成一张结构标准的特征表，用于推荐系统
import pandas as pd
from scipy.stats import ttest_ind
import numpy as np
import os

#设置数据文件所在的目录路径
data_dir = r"d:\Googledownload\depmap data"
# ==== 配置 ====
expr_path = os.path.join(data_dir,"CCLE_expression.csv")
mutation_path = os.path.join(data_dir,"CCLE_mutations.csv")
sample_info_path = os.path.join(data_dir, "sample_info.csv")

# ==== 加载数据 ====
print("📥 加载表达矩阵...")
expr_df = pd.read_csv(expr_path, index_col=0)  # 行是 DepMap_ID，列是 "KRAS (3845)" 等
print(f"表达矩阵维度: {expr_df.shape}")

# 清洗列名，提取 gene symbol
expr_df.columns = expr_df.columns.str.extract(r"^([A-Z0-9\-]+)")[0]

# 加载突变和样本信息
mut_df = pd.read_csv(mutation_path, low_memory=False)
sample_info = pd.read_csv(sample_info_path)

# 找出 KRAS G12C 阳性样本
mut_kras_g12c = mut_df[
    (mut_df["Hugo_Symbol"] == "KRAS") &
    (mut_df["Protein_Change"].str.contains("G12C", na=False))
]
g12c_ids = set(mut_kras_g12c["DepMap_ID"].unique())

# 表达数据中存在的 G12C 样本
g12c_in_expr = list(g12c_ids.intersection(set(expr_df.index)))
non_g12c_ids = list(set(expr_df.index) - set(g12c_in_expr))

print(f"✅ 表达数据中找到 KRAS G12C 阳性样本数: {len(g12c_in_expr)}")
print(f"✅ 非 G12C 样本数: {len(non_g12c_ids)}")

# ==== 差异分析 ====
results = []
for gene in expr_df.columns:
    g12c_vals = expr_df.loc[g12c_in_expr, gene].dropna()
    other_vals = expr_df.loc[non_g12c_ids, gene].dropna()
    
    if len(g12c_vals) >= 3 and len(other_vals) >= 3:
        stat, pval = ttest_ind(g12c_vals, other_vals, equal_var=False)
        lfc = np.log2(g12c_vals.mean() + 1) - np.log2(other_vals.mean() + 1)
    else:
        lfc, pval = np.nan, np.nan
    
    mean_expr = expr_df[gene].mean()
    std_expr = expr_df[gene].std()
    
    results.append((gene, lfc, pval, mean_expr, std_expr))

# ==== 输出结果 ====
df_out = pd.DataFrame(results, columns=["gene", "RNAseq_LFC", "RNAseq_pval", "RNAseq_mean", "RNAseq_std"])
output_path = os.path.join(data_dir, "RNAseq_expression_features_all_genes.csv")
df_out.to_csv(output_path, index=False)
print(f"✅ 特征提取完成，结果已保存至 {output_path}")
