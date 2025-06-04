# -*- coding: utf-8 -*-
"""
脚本功能：
- 从 HGNC 提供的 protein-coding_gene.txt 文件中读取基因名（symbol）与其别名（alias_symbol）
- 构建一个包含所有基因关键词的集合（标准名 + 所有别名）
- 将结果保存为一个文本文件，每行一个关键词，便于共现分析等后续使用
"""

import pandas as pd
import os
import sys

# 确保脚本在其所在文件夹中运行
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
print(f"📂 工作目录已设置为: {os.getcwd()}")


# ------------------------------------------------------
# 步骤 1：读取 HGNC 文件
# 文件为 TSV（tab-separated values，制表符分隔文本文件）
# ------------------------------------------------------
input_file = "protein-coding_gene.txt"  # 替换为你本地的文件路径

# pandas.read_csv 是用于读取表格型数据的函数
# 参数 sep='\t' 表示使用制表符分隔字段
# dtype=str 保证所有数据以字符串形式读取（避免数值转化错误）
df = pd.read_csv(input_file, sep="\t", dtype=str)

# ------------------------------------------------------
# 步骤 2：处理 alias_symbol 字段
# 有些别名字段为空，因此先填充空字符串
# ------------------------------------------------------
df["alias_symbol"] = df["alias_symbol"].fillna("")

# ------------------------------------------------------
# 步骤 3：构建关键词集合（set 不重复）
# 包含所有 symbol 和所有 alias_symbol 中的每个别名
# ------------------------------------------------------
gene_keywords = set()

# 添加标准 symbol
for symbol in df["symbol"]:
    gene_keywords.add(symbol.strip())

# 添加别名 alias（可能是用逗号分隔的一串别名）
for alias_str in df["alias_symbol"]:
    # 分割 alias 字符串，注意可能存在多个 alias 用逗号隔开
    aliases = alias_str.split(",")
    for alias in aliases:
        alias = alias.strip()
        if alias:
            gene_keywords.add(alias)

# ------------------------------------------------------
# 步骤 4：保存结果到本地文件
# 每行一个基因关键词，便于下游 NLP 分词、共现分析使用
# ------------------------------------------------------
output_file = "gene_keywords_list.txt"

with open(output_file, "w", encoding="utf-8") as f:
    for keyword in sorted(gene_keywords):
        f.write(keyword + "\n")

print(f"✅ 共提取关键词数：{len(gene_keywords)}")
print(f"📄 关键词列表已保存为：{output_file}")
