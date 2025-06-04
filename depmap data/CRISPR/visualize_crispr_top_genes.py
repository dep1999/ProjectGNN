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

# ==== 排序并提取前20名一致性最高的基因 ====
top_consistent = df.sort_values(by="crispr_consistency", ascending=False).head(20)

# ==== 可视化：Top 20 一致性基因的 KO 分数均值 ====
plt.figure(figsize=(10, 6))
sns.barplot(x="crispr_mean_effect", y="gene", data=top_consistent, palette="viridis")
plt.title("Top 20 CRISPR Consistent Genes in KRAS G12C+ Cells")
plt.xlabel("Mean KO Effect Score")
plt.ylabel("Gene")
plt.tight_layout()
plt.savefig("top20_crispr_consistency_genes.png")
plt.close()

# ==== 富集分析预处理：提取基因 symbol（去除括号内 Entrez ID） ====
top_consistent["gene_symbol"] = top_consistent["gene"].str.extract(r"^([A-Z0-9\-]+)")
top_consistent["gene_symbol"].to_csv("top20_crispr_genes_for_enrichment.txt", index=False, header=False)

print("✅ 可视化图像已保存为: top20_crispr_consistency_genes.png")
print("✅ 富集分析基因列表已保存为: top20_crispr_genes_for_enrichment.txt")
