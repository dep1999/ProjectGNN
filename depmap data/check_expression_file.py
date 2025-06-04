import pandas as pd
import os

# 设置数据文件所在的目录路径
data_dir = r"d:\Googledownload\depmap data"

# 使用绝对路径读取文件
expr_path = os.path.join(data_dir, "CCLE_expression.csv")

# 读取前10行查看文件结构
try:
    expr_df = pd.read_csv(expr_path, nrows=10)
    print("文件读取成功！")
    print("前10个基因ID:")
    print(expr_df.iloc[:, 0].tolist())
    print("\n列名:")
    print(expr_df.columns.tolist()[:5])  # 只显示前5个列名
except Exception as e:
    print(f"读取文件时出错: {e}")