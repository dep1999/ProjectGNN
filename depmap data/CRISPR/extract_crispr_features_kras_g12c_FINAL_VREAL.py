import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import os
import sys

# 确保脚本在其所在文件夹中运行
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
print(f"📂 工作目录已设置为: {os.getcwd()}")

# ==== 文件路径 ====
effect_path = "CRISPR_gene_effect.csv"
dependency_path = "CRISPR_gene_dependency.csv"
mutation_path = "CCLE_mutations.csv"
sample_info_path = "sample_info.csv"

# ==== 加载数据 ====
print("📥 加载数据中...")
effect_df = pd.read_csv(effect_path)
dependency_df = pd.read_csv(dependency_path)
mut_df = pd.read_csv(mutation_path, low_memory=False)
sample_info = pd.read_csv(sample_info_path)

# ==== 设置 DepMap_ID 为行索引（不转置） ====
effect_df = effect_df.set_index(effect_df.columns[0])
dependency_df = dependency_df.set_index(dependency_df.columns[0])

# ==== 获取 KRAS G12C 阳性 DepMap ID ====
g12c_samples = mut_df[
    (mut_df["Hugo_Symbol"] == "KRAS") &
    (mut_df["Protein_Change"].str.contains("G12C", na=False))
]["DepMap_ID"].unique().tolist()

print(f"✅ 找到 KRAS G12C 阳性细胞系数：{len(g12c_samples)}")

# 保留在 effect 中存在的样本
common_samples = list(set(g12c_samples).intersection(set(effect_df.index)))
other_samples = list(set(effect_df.index) - set(g12c_samples))

# ==== 计算 CRISPR KO 特征 ====
print("⚙️ 正在计算 KO 特征...")
gene_stats = []
for gene in effect_df.columns:
    try:
        g12c_vals = effect_df.loc[common_samples, gene].dropna()
        other_vals = effect_df.loc[other_samples, gene].dropna()
    except KeyError:
        continue

    if len(g12c_vals) >= 3 and len(other_vals) >= 3:
        lfc = g12c_vals.mean() - other_vals.mean()
        pval = ttest_ind(g12c_vals, other_vals, equal_var=False).pvalue
    else:
        lfc, pval = np.nan, np.nan

    gene_stats.append({
        "gene": gene,
        "crispr_mean_effect": g12c_vals.mean() if not g12c_vals.empty else np.nan,
        "crispr_std_effect": g12c_vals.std() if not g12c_vals.empty else np.nan,
        "crispr_consistency": 1 - g12c_vals.std() if not g12c_vals.empty else np.nan,
        "crispr_effect_LFC": lfc,
        "crispr_effect_pval": pval
    })

ko_df = pd.DataFrame(gene_stats)

# ==== 添加 dependency 分数 ====
print("📊 添加依赖性分数...")
dependency_means = dependency_df.loc[common_samples].mean()
ko_df["crispr_dependency_score"] = ko_df["gene"].map(dependency_means)

# ==== 保存结果 ====
ko_df.to_csv("CRISPR_features_KRAS_G12C_FINAL_VREAL.csv", index=False)
print("✅ 最终特征提取完成，结果已保存至 CRISPR_features_KRAS_G12C_FINAL_VREAL.csv")
