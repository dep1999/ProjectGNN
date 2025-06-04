# KRAS G12C突变型癌症耐药基因预测系统

## 项目概述

本项目构建了一个基于图神经网络（GNN）的多源知识图谱增强推荐系统，用于KRAS G12C突变非小细胞肺癌（NSCLC）的耐药因子预测与推荐。系统整合了多种数据源，包括CRISPR KO数据、基因表达数据、文献共现数据和通路关系数据，构建了一个综合性的知识图谱，并应用GNN模型进行耐药基因预测。

### 主要功能

- **多源数据整合**：整合DepMap数据集（CRISPR KO、突变信息、RNA-seq表达数据）
- **知识图谱构建**：构建PPI网络、文献共现网络和通路关系网络
- **图嵌入学习**：应用GCN和GAT进行图嵌入学习
- **多目标优化**：设计多目标优化框架进行耐药基因预测
- **可视化分析**：提供知识图谱和预测结果的可视化功能

## 项目结构

```
.
├── data_processing/             # 数据处理模块
│   └── knowledge_graph_builder.py  # 知识图谱构建
├── models/                      # 模型模块
│   ├── gnn_model.py            # GNN模型定义
│   └── train_evaluate.py       # 模型训练和评估
├── depmap data/                # 原始数据目录
├── results/                    # 结果输出目录
├── main.py                     # 主程序
├── README.md                   # 项目说明
└── requirements.txt            # 依赖包列表
```

## 安装说明

### 环境要求

- Python 3.7+
- PyTorch 1.9+
- PyTorch Geometric 2.0+

### 安装步骤

1. 克隆项目代码

```bash
git clone <repository-url>
cd <repository-directory>
```

2. 安装依赖包

```bash
pip install -r requirements.txt
```

## 使用方法

### 数据准备

将DepMap数据集放置在`depmap data`目录下，包括：

- CRISPR_gene_effect.csv：CRISPR KO数据
- CCLE_expression.csv：RNA-seq表达数据
- CCLE_mutations.csv：突变信息数据
- sample_info.csv：样本信息数据

### 运行程序

```bash
python main.py --model_type gcn --hidden_channels 64 --epochs 200
```

### 命令行参数

- `--data_dir`：数据目录，默认为`depmap data`
- `--model_type`：模型类型，可选`gcn`、`gat`或`multi`，默认为`gcn`
- `--hidden_channels`：隐藏层维度，默认为64
- `--epochs`：训练轮数，默认为200
- `--lr`：学习率，默认为0.001
- `--weight_decay`：权重衰减，默认为5e-4
- `--patience`：早停耐心值，默认为20
- `--seed`：随机种子，默认为42
- `--output_dir`：输出目录，默认为`results`
- `--known_resistant_genes`：已知耐药基因列表文件，默认为None
- `--top_k`：预测结果展示的基因数量，默认为20

### 输出结果

程序运行后，将在`results`目录下生成以下文件：

- `knowledge_graph.png`：知识图谱可视化
- `labeled_knowledge_graph.png`：带标签的知识图谱可视化
- `prediction_results.png`：预测结果可视化
- `prediction_results.csv`：预测结果数据
- `training_curve_xxx.png`：训练曲线
- `xxx_model.pt`：训练好的模型
- `graph_data/`：知识图谱数据
  - `nodes.csv`：节点数据
  - `edges.csv`：边数据

## 模型说明

### GCN模型

图卷积神经网络（GCN）模型通过聚合邻居节点的特征来学习节点表示，适合捕捉全局结构特征。

### GAT模型

图注意力网络（GAT）模型引入了注意力机制，能够自适应地为不同邻居分配不同的权重，适合识别邻居中的重要节点。

### 多目标优化模型

多目标优化模型同时优化多个目标，包括耐药性预测、结构新颖性、靶向可行性等，能够提供更全面的耐药基因推荐。

## 参考文献

1. Kipf, T. N., & Welling, M. (2016). Semi-supervised classification with graph convolutional networks. arXiv preprint arXiv:1609.02907.
2. Veličković, P., Cucurull, G., Casanova, A., Romero, A., Lio, P., & Bengio, Y. (2017). Graph attention networks. arXiv preprint arXiv:1710.10903.
3. Tsherniak, A., Vazquez, F., Montgomery, P. G., Weir, B. A., Kryukov, G., Cowley, G. S., ... & Meyers, R. M. (2017). Defining a cancer dependency map. Cell, 170(3), 564-576.