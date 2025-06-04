import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import sys
# 确保脚本在其所在文件夹中运行
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
print(f"📂 工作目录已设置为: {os.getcwd()}")

# ==== 加载数据 ====
data_path = "CRISPR_features_KRAS_G12C_FINAL_VREAL.csv"
df = pd.read_csv(data_path)

# ==== 数据预处理 ====
df = df.dropna()
df["gene_symbol"] = df["gene"].str.extract(r"^([A-Z0-9\-]+)")

# ==== 1️⃣ SCATTER PLOT: LFC vs Dependency Score ====
plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=df, x="crispr_effect_LFC", y="crispr_dependency_score",
    hue="crispr_consistency", size="crispr_mean_effect", palette="viridis", legend=False
)
plt.title("Scatter: KO Effect LFC vs Dependency Score")
plt.xlabel("crispr_effect_LFC")
plt.ylabel("crispr_dependency_score")
plt.tight_layout()
plt.savefig("crispr_scatter_lfc_dependency.png")
plt.close()

# ==== 2️⃣ HEATMAP: Top 20 Genes by KO LFC ====
top_lfc = df.sort_values(by="crispr_effect_LFC", ascending=False).head(20)
heatmap_data = top_lfc[["crispr_mean_effect", "crispr_consistency", "crispr_dependency_score"]]
heatmap_data.index = top_lfc["gene_symbol"]

plt.figure(figsize=(10, 8))
sns.heatmap(heatmap_data, cmap="RdBu_r", annot=True, fmt=".2f", cbar=True)
plt.title("Heatmap of Top 20 LFC Genes (KO Features)")
plt.tight_layout()
plt.savefig("crispr_heatmap_topLFC.png")
plt.close()

print("✅ 补充图生成完成：")
print(" - scatter: crispr_scatter_lfc_dependency.png")
print(" - heatmap: crispr_heatmap_topLFC.png")
