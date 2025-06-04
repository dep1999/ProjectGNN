import pandas as pd

# 读取 metadata 文件，筛选 NSCLC 背景样本
sample_info = pd.read_csv("sample_info.csv")
nsclc_samples = sample_info[sample_info["primary_disease"].str.contains("NSCLC", na=False)]
nsclc_ids = nsclc_samples["DepMap_ID"].tolist()

# 读取 KRAS G12C 突变信息
mut_df = pd.read_csv("CCLE_mutations.csv", low_memory=False)
mut_kras_g12c = mut_df[(mut_df["Hugo_Symbol"] == "KRAS") & (mut_df["Variant_Classification"].str.contains("G12C", na=False))]
kras_g12c_ids = mut_kras_g12c["DepMap_ID"].unique().tolist()

# 交集样本 = NSCLC 且含 KRAS G12C 的细胞系
target_ids = list(set(nsclc_ids).intersection(set(kras_g12c_ids)))

print(f"共识别出 {len(target_ids)} 个 KRAS G12C 阳性 NSCLC 细胞系")

# 读取 CRISPR KO 表格
crispr_df = pd.read_csv("CRISPR_gene_effect.csv")
crispr_df.set_index("DepMap_ID", inplace=True)

# 计算所有基因在 KRAS G12C 阳性细胞系中的 KO 平均值与一致性指标（标准差）
crispr_kras = crispr_df.loc[target_ids].T  # 行变为基因，列为细胞系
crispr_summary = pd.DataFrame({
    "gene": crispr_kras.index,
    "crispr_mean_effect": crispr_kras.mean(axis=1),
    "crispr_std_effect": crispr_kras.std(axis=1),
    "crispr_consistency": 1 - crispr_kras.std(axis=1)  # 标准差越小 → 一致性越高
}).reset_index(drop=True)

# 保存结果
crispr_summary.to_csv("KRAS_G12C_crispr_summary.csv", index=False)
print("已保存 KO 平均值与一致性文件：KRAS_G12C_crispr_summary.csv")
