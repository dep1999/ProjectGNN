
import pandas as pd
import os

# ���������ļ����ڵ�Ŀ¼·��
data_dir = r"d:\Googledownload\depmap data"

# ==== ���� ====
# �������ļ�·����CCLE_expression.csv��
expr_path = os.path.join(data_dir, "CCLE_expression.csv")

# ����Ȥ�Ļ��򣨿��滻Ϊ���Լ����б�
genes_of_interest = ["KRAS", "TP53", "DUSP6", "ERK2", "AKT1"]

# ����Ȥ�� DepMap_ID��ϸ��ϵ��
samples_of_interest = ["ACH-000001", "ACH-000002", "ACH-000003"]

# ==== ��ȡ ====
# ���ر�������Ϊ������Ϊ DepMap_ID��
expr_df = pd.read_csv(expr_path)
expr_df = expr_df.rename(columns={expr_df.columns[0]: "gene"})  # ��һ��Ϊ������

# ɸѡ����
expr_df = expr_df[expr_df["gene"].isin(genes_of_interest)]

# ���� gene Ϊ����
expr_df.set_index("gene", inplace=True)

# ɸѡ�������У�
valid_cols = [col for col in samples_of_interest if col in expr_df.columns]
expr_subset = expr_df[valid_cols]

# ==== ��� ====
output_path = os.path.join(data_dir, "gene_expression_subset.csv")
expr_subset.to_csv(output_path)
print(f"�ѱ��� {len(genes_of_interest)} �������� {len(valid_cols)} ��ϸ��ϵ�еı�����")
