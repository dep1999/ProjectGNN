# -*- coding: utf-8 -*-
"""
该脚本用于处理 HGNC 转换后的 CSV 文件，提取主要字段并验证数据有效性。
提取字段：
- symbol（标准基因符号）
- name（标准基因名称）
- alias_symbol（别名符号）
- alias_name（别名描述）
"""

import pandas as pd
import os
import sys

# 确保脚本在当前目录下运行
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
print(f"📂 工作目录已设置为: {os.getcwd()}")

# Step 1: 加载 CSV 文件（替换为你的路径）
input_file = "protein-coding_gene.csv"  # 你自己的文件路径
output_file = "processed_gene_info.csv"

# Step 2: 读取 CSV 文件
try:
    # 首先尝试UTF-8编码
    df = pd.read_csv(input_file, dtype=str)
except UnicodeDecodeError:
    # 如果UTF-8失败，尝试其他常见编码
    try:
        df = pd.read_csv(input_file, dtype=str, encoding='latin1')
        print("ℹ️ 使用latin1编码读取文件成功")
    except Exception as e:
        try:
            df = pd.read_csv(input_file, dtype=str, encoding='gbk')
            print("ℹ️ 使用GBK编码读取文件成功")
        except Exception as e:
            try:
                df = pd.read_csv(input_file, dtype=str, encoding='cp1252')
                print("ℹ️ 使用cp1252编码读取文件成功")
            except Exception as e:
                print(f"❌ 无法读取文件，请检查文件编码: {str(e)}")
                sys.exit(1)

df = df.fillna("")  # 空值填为""

# Step 3: 提取目标列
if not all(col in df.columns for col in ["symbol", "name", "alias_symbol", "alias_name"]):
    raise ValueError("❌ 缺少必要字段，请确认CSV文件包含 'symbol', 'name', 'alias_symbol', 'alias_name'")

extracted_df = df[["symbol", "name", "alias_symbol", "alias_name"]].copy()

# Step 4: 验证数据正确性
def validate_row(row):
    issues = []
    if not row["symbol"]:
        issues.append("缺失symbol")
    if "," in row["alias_symbol"] and "," in row["alias_name"]:
        num_aliases = len(row["alias_symbol"].split(","))
        num_alias_names = len(row["alias_name"].split(","))
        if num_aliases != num_alias_names:
            issues.append("alias数量不一致")
    return "; ".join(issues)

# 添加验证结果列
extracted_df["Validation"] = extracted_df.apply(validate_row, axis=1)

# Step 5: 导出
extracted_df.to_csv(output_file, index=False)
print(f"✅ 处理完成，结果保存为：{output_file}")
print("📋 包含字段：symbol, name, alias_symbol, alias_name, Validation")
