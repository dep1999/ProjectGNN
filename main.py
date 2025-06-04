import os
import sys
import torch
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from sklearn.model_selection import train_test_split
from torch_geometric.utils import to_networkx
import argparse

# 添加项目根目录到系统路径
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# 导入自定义模块
from data_processing.knowledge_graph_builder import KnowledgeGraphBuilder
from models.gnn_model import GCNModel, GATModel, MultiObjectiveGNN
from models.train_evaluate import GNNTrainer

def parse_args():
    """
    解析命令行参数
    """
    parser = argparse.ArgumentParser(description='KRAS G12C突变型癌症耐药基因预测系统')
    parser.add_argument('--data_dir', type=str, default='depmap data', help='数据目录')
    parser.add_argument('--model_type', type=str, default='gcn', choices=['gcn', 'gat', 'multi'], help='模型类型')
    parser.add_argument('--hidden_channels', type=int, default=64, help='隐藏层维度')
    parser.add_argument('--epochs', type=int, default=200, help='训练轮数')
    parser.add_argument('--lr', type=float, default=0.001, help='学习率')
    parser.add_argument('--weight_decay', type=float, default=5e-4, help='权重衰减')
    parser.add_argument('--patience', type=int, default=20, help='早停耐心值')
    parser.add_argument('--seed', type=int, default=42, help='随机种子')
    parser.add_argument('--output_dir', type=str, default='results', help='输出目录')
    parser.add_argument('--known_resistant_genes', type=str, default=None, help='已知耐药基因列表文件')
    parser.add_argument('--top_k', type=int, default=20, help='预测结果展示的基因数量')
    
    return parser.parse_args()

def set_seed(seed):
    """
    设置随机种子
    """
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    torch.backends.cudnn.deterministic = True

def load_known_resistant_genes(file_path):
    """
    加载已知耐药基因列表
    """
    if file_path and os.path.exists(file_path):
        with open(file_path, 'r') as f:
            genes = [line.strip() for line in f.readlines()]
        print(f"加载了 {len(genes)} 个已知耐药基因")
        return genes
    else:
        print("未提供已知耐药基因列表，将使用CRISPR KO数据中的前20个基因作为已知耐药基因")
        return None

def create_labels(data, known_genes=None):
    """
    创建标签
    如果提供了已知耐药基因列表，则将这些基因的标签设为1，其他设为0
    否则，使用CRISPR KO数据中的前20个基因作为已知耐药基因
    """
    num_nodes = data.num_nodes
    labels = torch.zeros(num_nodes)
    
    if known_genes:
        # 使用提供的已知耐药基因列表
        for gene in known_genes:
            if gene in data.node_map:
                labels[data.node_map[gene]] = 1.0
    else:
        # 使用CRISPR KO数据中的前20个基因作为已知耐药基因
        try:
            crispr_file = os.path.join('depmap data', 'CRISPR', 'crispr_top20_candidates_summary.csv')
            if os.path.exists(crispr_file):
                crispr_df = pd.read_csv(crispr_file)
                top_genes = crispr_df['gene'].tolist()
                for gene in top_genes:
                    if gene in data.node_map:
                        labels[data.node_map[gene]] = 1.0
                print(f"使用CRISPR KO数据中的前 {len(top_genes)} 个基因作为已知耐药基因")
            else:
                # 随机选择20个基因作为已知耐药基因（仅用于演示）
                indices = np.random.choice(num_nodes, 20, replace=False)
                labels[indices] = 1.0
                print("未找到CRISPR KO数据，随机选择20个基因作为已知耐药基因")
        except Exception as e:
            print(f"创建标签时出错: {e}")
            # 随机选择20个基因作为已知耐药基因（仅用于演示）
            indices = np.random.choice(num_nodes, 20, replace=False)
            labels[indices] = 1.0
            print("随机选择20个基因作为已知耐药基因")
    
    print(f"创建了 {labels.sum().item()} 个正样本标签")
    return labels

def visualize_graph(G, node_labels=None, output_path=None):
    """
    可视化知识图谱
    """
    plt.figure(figsize=(12, 12))
    
    # 设置节点颜色
    if node_labels is not None:
        node_colors = []
        for node in G.nodes():
            if node_labels[node] > 0:
                node_colors.append('red')  # 已知耐药基因
            else:
                node_colors.append('skyblue')  # 其他基因
    else:
        node_colors = 'skyblue'
    
    # 设置节点大小
    node_sizes = []
    for node in G.nodes():
        degree = G.degree(node)
        node_sizes.append(100 + degree * 10)  # 节点大小与度成正比
    
    # 设置边颜色
    edge_colors = []
    for u, v, data in G.edges(data=True):
        if data.get('type') == 'PPI':
            edge_colors.append('gray')
        elif data.get('type') == 'Literature':
            edge_colors.append('green')
        elif data.get('type') == 'Pathway':
            edge_colors.append('blue')
        else:
            edge_colors.append('black')
    
    # 绘制图谱
    pos = nx.spring_layout(G, seed=42)  # 使用弹簧布局
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, alpha=0.8)
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=0.5, alpha=0.5)
    
    # 添加标签
    labels = {}
    for node in G.nodes():
        if G.nodes[node].get('name'):
            labels[node] = G.nodes[node]['name']
    nx.draw_networkx_labels(G, pos, labels, font_size=8, font_color='black')
    
    plt.title('KRAS G12C突变型癌症耐药基因知识图谱')
    plt.axis('off')
    
    # 保存图像
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"图谱可视化已保存到 {output_path}")
    
    plt.close()

def visualize_predictions(predictions, data, top_k=20, output_path=None):
    """
    可视化预测结果
    """
    # 获取预测分数
    if isinstance(predictions, list):
        # 多目标模型的预测结果
        scores = predictions[0].squeeze().cpu().numpy()
    else:
        # 单目标模型的预测结果
        scores = predictions.squeeze().cpu().numpy()
    
    # 创建预测结果DataFrame
    result_df = pd.DataFrame({
        'gene': [data.reverse_node_map[i] for i in range(len(scores))],
        'score': scores
    })
    
    # 按预测分数排序
    result_df = result_df.sort_values('score', ascending=False).reset_index(drop=True)
    
    # 保存完整预测结果
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    result_df.to_csv(output_path.replace('.png', '.csv'), index=False)
    print(f"预测结果已保存到 {output_path.replace('.png', '.csv')}")
    
    # 可视化前top_k个预测结果
    plt.figure(figsize=(12, 8))
    top_k_df = result_df.head(top_k)
    plt.barh(top_k_df['gene'][::-1], top_k_df['score'][::-1], color='skyblue')
    plt.xlabel('预测分数')
    plt.ylabel('基因')
    plt.title(f'KRAS G12C突变型癌症耐药基因预测结果 (Top {top_k})')
    plt.grid(axis='x', linestyle='--', alpha=0.7)
    
    # 保存图像
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"预测结果可视化已保存到 {output_path}")
    
    plt.close()
    
    return result_df

def main():
    """
    主函数
    """
    # 解析命令行参数
    args = parse_args()
    
    # 设置随机种子
    set_seed(args.seed)
    
    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 加载已知耐药基因列表
    known_genes = load_known_resistant_genes(args.known_resistant_genes)
    
    # 创建知识图谱构建器
    print("\n=== 构建知识图谱 ===")
    kg_builder = KnowledgeGraphBuilder(data_dir=args.data_dir)
    
    # 加载数据
    crispr_df = kg_builder.load_crispr_data()
    expression_df = kg_builder.load_expression_data()
    literature_df = kg_builder.load_literature_data()
    
    # 创建节点特征
    merged_df = kg_builder.create_node_features(crispr_df, expression_df, literature_df)
    
    # 构建PPI网络
    kg_builder.build_ppi_network()
    
    # 添加文献共现边
    kg_builder.add_literature_edges(literature_df)
    
    # 添加通路关系边
    kg_builder.add_pathway_edges()
    
    # 创建PyG数据对象
    data = kg_builder.create_pyg_data()
    
    # 保存图谱数据
    node_df, edge_df = kg_builder.save_graph(output_dir=os.path.join(args.output_dir, 'graph_data'))
    
    # 可视化知识图谱
    G = kg_builder.to_networkx()
    visualize_graph(G, output_path=os.path.join(args.output_dir, 'knowledge_graph.png'))
    
    # 创建标签
    labels = create_labels(data, known_genes)
    
    # 可视化带标签的知识图谱
    visualize_graph(G, node_labels=labels, output_path=os.path.join(args.output_dir, 'labeled_knowledge_graph.png'))
    
    # 划分训练集和测试集
    indices = list(range(data.num_nodes))
    train_indices, test_indices = train_test_split(indices, test_size=0.2, random_state=args.seed)
    train_mask = torch.zeros(data.num_nodes, dtype=torch.bool)
    test_mask = torch.zeros(data.num_nodes, dtype=torch.bool)
    train_mask[train_indices] = True
    test_mask[test_indices] = True
    
    # 初始化训练器
    print("\n=== 训练GNN模型 ===")
    trainer = GNNTrainer(
        model_type=args.model_type,
        hidden_channels=args.hidden_channels,
        num_classes=1,
        lr=args.lr,
        weight_decay=args.weight_decay
    )
    trainer.initialize_model(data.num_features)
    
    # 训练模型
    trainer.train(
        data,
        labels,
        train_mask,
        test_mask,
        epochs=args.epochs,
        patience=args.patience
    )
    
    # 预测
    print("\n=== 预测耐药基因 ===")
    predictions, embeddings = trainer.predict(data)
    
    # 计算测试集指标
    metrics = trainer.calculate_metrics(labels[test_mask], predictions[test_mask])
    print(f"测试集指标: {metrics}")
    
    # 保存模型
    model_path = os.path.join(args.output_dir, f'{args.model_type}_model.pt')
    trainer.save_model(model_path)
    
    # 可视化预测结果
    result_df = visualize_predictions(
        predictions,
        data,
        top_k=args.top_k,
        output_path=os.path.join(args.output_dir, 'prediction_results.png')
    )
    
    # 输出预测结果
    print("\n=== 预测结果 ===")
    print(f"Top {args.top_k} 耐药基因预测结果:")
    for i, (gene, score) in enumerate(zip(result_df['gene'].head(args.top_k), result_df['score'].head(args.top_k))):
        print(f"{i+1}. {gene}: {score:.4f}")
    
    print("\n=== 完成 ===")
    print(f"所有结果已保存到 {os.path.abspath(args.output_dir)} 目录")

if __name__ == "__main__":
    main()