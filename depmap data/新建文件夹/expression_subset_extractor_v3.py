import pandas as pd
import os

#设置
#  设置数据文件所在目录
data_dir = r"D:\Googledownload\depmap data" 
expr_path = os.path.join(data_dir, "CCLE_expression.csv")  
#构建完整表达矩阵文件路径

# 感兴趣的基因名（symbol），与列名中的格式如 "KRAS (3845)" 匹配
genes_of_interest = ["KRAS", "TP53", "DUSP6", "ERK2", "AKT1"]

# 感兴趣的细胞系 ID（行名）
samples_of_interest = ["ACH-001113", "ACH-001289", "ACH-001339"]

# ==== 加载并转置表达矩阵 ====
print("📥 正在加载表达矩阵...")
expr_df = pd.read_csv(expr_path, index_col=0)  # 行是 DepMap_ID，列是 gene(symbol + EntrezID)

print(f"✅ 原始维度: {expr_df.shape}")
print("🔁 正在转置表达矩阵...")
expr_df = expr_df.T  # 转置：行→gene，列→DepMap_ID

# ==== 清洗 gene 名（提取符号前半部分） ====
expr_df.index = expr_df.index.str.extract(r"^([A-Z0-9\-]+)")[0]  # 提取基因symbol部分

# ==== 筛选目标基因和样本 ====
valid_genes = [g for g in genes_of_interest if g in expr_df.index]
valid_samples = [s for s in samples_of_interest if s in expr_df.columns]

if not valid_genes:
    print("⚠️ 无有效基因匹配！")
if not valid_samples:
    print("⚠️ 无有效细胞系 ID 匹配！")

expr_subset = expr_df.loc[valid_genes, valid_samples]

# ==== 保存 ====
output_path = os.path.join(data_dir, "gene_expression_subset_fixed.csv")
expr_subset.to_csv(output_path)
print(f"✅ 子集表达矩阵已保存至: {output_path}")
