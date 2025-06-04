import pandas as pd
import numpy as np
import os
import torch
from torch_geometric.data import Data
import networkx as nx
from sklearn.preprocessing import StandardScaler

class KnowledgeGraphBuilder:
    """
    知识图谱构建类，用于KRAS G12C突变型癌症耐药基因预测
    整合多源数据构建知识图谱，包括PPI网络、文献共现和通路关系
    """
    def __init__(self, data_dir="../depmap data"):
        self.data_dir = data_dir
        self.node_features = None
        self.edge_index = None
        self.edge_attr = None
        self.node_map = {}
        self.reverse_node_map = {}
        self.node_types = {}
        self.edge_types = {}
        
    def load_crispr_data(self, file_path="CRISPR/CRISPR_features_KRAS_G12C_FINAL_VREAL.csv"):
        """
        加载CRISPR KO数据
        """
        file_path = os.path.join(self.data_dir, file_path)
        crispr_df = pd.read_csv(file_path)
        print(f"加载了 {len(crispr_df)} 个基因的CRISPR KO数据")
        return crispr_df
    
    def load_expression_data(self, file_path="RNAseq全基因表达特征-来自（expression_feature_extractor_full.py）/RNAseq_expression_features_all_genes.csv"):
        """
        加载RNA-seq表达数据
        """
        file_path = os.path.join(self.data_dir, file_path)
        expression_df = pd.read_csv(file_path)
        print(f"加载了 {len(expression_df)} 个基因的RNA-seq表达数据")
        return expression_df
    
    def load_literature_data(self, file_path="../Papers from PubMed/KRAS_literature_collection.csv"):
        """
        加载文献共现数据
        """
        file_path = os.path.join(os.path.dirname(self.data_dir), file_path)
        if os.path.exists(file_path):
            literature_df = pd.read_csv(file_path)
            print(f"加载了 {len(literature_df)} 条文献共现数据")
            return literature_df
        else:
            print(f"文件不存在: {file_path}")
            # 创建一个空的DataFrame作为替代
            return pd.DataFrame(columns=['gene', 'lit_kras_freq', 'lit_nsclc_freq', 'lit_total'])
    
    def create_node_features(self, crispr_df, expression_df, literature_df=None):
        """
        创建节点特征矩阵
        整合CRISPR KO、RNA-seq表达和文献共现数据
        """
        # 合并CRISPR和表达数据
        merged_df = pd.merge(crispr_df, expression_df, on='gene', how='outer')
        
        # 如果有文献数据，也合并进来
        if literature_df is not None and not literature_df.empty:
            merged_df = pd.merge(merged_df, literature_df, on='gene', how='outer')
        
        # 填充缺失值
        merged_df = merged_df.fillna(0)
        
        # 创建节点映射
        for i, gene in enumerate(merged_df['gene']):
            self.node_map[gene] = i
            self.reverse_node_map[i] = gene
            self.node_types[i] = 'Gene'
        
        # 选择特征列
        feature_cols = [col for col in merged_df.columns if col != 'gene']
        
        # 标准化特征
        scaler = StandardScaler()
        features = scaler.fit_transform(merged_df[feature_cols])
        
        self.node_features = torch.FloatTensor(features)
        print(f"创建了 {len(self.node_map)} 个节点的特征矩阵，特征维度为 {self.node_features.shape[1]}")
        
        return merged_df
    
    def build_ppi_network(self, ppi_file=None, threshold=0.7):
        """
        构建PPI网络
        如果没有提供PPI文件，则根据CRISPR KO相关性构建
        """
        if ppi_file and os.path.exists(ppi_file):
            # 从文件加载PPI网络
            ppi_df = pd.read_csv(ppi_file)
            edges = []
            edge_weights = []
            edge_types = []
            
            for _, row in ppi_df.iterrows():
                source = row['source']
                target = row['target']
                weight = row['weight']
                edge_type = row['edge_type']
                
                if source in self.node_map and target in self.node_map:
                    edges.append([self.node_map[source], self.node_map[target]])
                    edge_weights.append(weight)
                    edge_types.append(edge_type)
            
            self.edge_index = torch.LongTensor(edges).t().contiguous()
            self.edge_attr = torch.FloatTensor(edge_weights).view(-1, 1)
            
            # 存储边类型
            for i, edge_type in enumerate(edge_types):
                self.edge_types[i] = edge_type
                
            print(f"从文件加载了 {len(edges)} 条PPI网络边")
        else:
            # 根据CRISPR KO相关性构建PPI网络
            print("未提供PPI文件，将根据特征相关性构建网络")
            
            # 计算节点特征的相关性矩阵
            features_np = self.node_features.numpy()
            corr_matrix = np.corrcoef(features_np)
            
            # 根据相关性阈值构建边
            edges = []
            edge_weights = []
            
            for i in range(len(corr_matrix)):
                for j in range(i+1, len(corr_matrix)):
                    if abs(corr_matrix[i, j]) > threshold:
                        edges.append([i, j])
                        # 双向边
                        edges.append([j, i])
                        edge_weights.append(abs(corr_matrix[i, j]))
                        edge_weights.append(abs(corr_matrix[i, j]))
                        # 存储边类型
                        self.edge_types[len(self.edge_types)] = 'PPI'
                        self.edge_types[len(self.edge_types)] = 'PPI'
            
            self.edge_index = torch.LongTensor(edges).t().contiguous()
            self.edge_attr = torch.FloatTensor(edge_weights).view(-1, 1)
            
            print(f"根据特征相关性构建了 {len(edges)} 条网络边")
    
    def add_literature_edges(self, literature_df, threshold=5):
        """
        添加文献共现边
        """
        if literature_df is None or literature_df.empty:
            print("没有文献共现数据，跳过添加文献共现边")
            return
        
        # 获取当前边的数量
        current_edges = self.edge_index.shape[1] if self.edge_index is not None else 0
        current_edge_types = len(self.edge_types)
        
        # 创建文献共现边
        new_edges = []
        new_edge_weights = []
        
        for i, row_i in literature_df.iterrows():
            gene_i = row_i['gene']
            if gene_i not in self.node_map:
                continue
                
            for j, row_j in literature_df.iterrows():
                if i >= j:
                    continue
                    
                gene_j = row_j['gene']
                if gene_j not in self.node_map:
                    continue
                
                # 计算共现频率
                cooccurrence = min(row_i['lit_kras_freq'], row_j['lit_kras_freq'])
                
                if cooccurrence > threshold:
                    new_edges.append([self.node_map[gene_i], self.node_map[gene_j]])
                    new_edges.append([self.node_map[gene_j], self.node_map[gene_i]])
                    new_edge_weights.append(cooccurrence / 100.0)  # 归一化权重
                    new_edge_weights.append(cooccurrence / 100.0)
                    # 存储边类型
                    self.edge_types[current_edge_types + len(new_edges) - 2] = 'Literature'
                    self.edge_types[current_edge_types + len(new_edges) - 1] = 'Literature'
        
        if new_edges:
            new_edge_index = torch.LongTensor(new_edges).t().contiguous()
            new_edge_attr = torch.FloatTensor(new_edge_weights).view(-1, 1)
            
            # 合并现有边和新边
            if self.edge_index is not None:
                self.edge_index = torch.cat([self.edge_index, new_edge_index], dim=1)
                self.edge_attr = torch.cat([self.edge_attr, new_edge_attr], dim=0)
            else:
                self.edge_index = new_edge_index
                self.edge_attr = new_edge_attr
            
            print(f"添加了 {len(new_edges)} 条文献共现边")
        else:
            print("没有添加文献共现边")
    
    def add_pathway_edges(self, pathway_file=None):
        """
        添加通路关系边
        """
        if pathway_file and os.path.exists(pathway_file):
            # 从文件加载通路关系
            pathway_df = pd.read_csv(pathway_file)
            
            # 获取当前边的数量
            current_edges = self.edge_index.shape[1] if self.edge_index is not None else 0
            current_edge_types = len(self.edge_types)
            
            # 创建通路关系边
            new_edges = []
            new_edge_weights = []
            
            for _, row in pathway_df.iterrows():
                source = row['source']
                target = row['target']
                weight = row['weight']
                edge_type = row['edge_type']
                
                if source in self.node_map and target in self.node_map and edge_type == 'Pathway_Member':
                    new_edges.append([self.node_map[source], self.node_map[target]])
                    new_edges.append([self.node_map[target], self.node_map[source]])
                    new_edge_weights.append(weight)
                    new_edge_weights.append(weight)
                    # 存储边类型
                    self.edge_types[current_edge_types + len(new_edges) - 2] = 'Pathway'
                    self.edge_types[current_edge_types + len(new_edges) - 1] = 'Pathway'
            
            if new_edges:
                new_edge_index = torch.LongTensor(new_edges).t().contiguous()
                new_edge_attr = torch.FloatTensor(new_edge_weights).view(-1, 1)
                
                # 合并现有边和新边
                if self.edge_index is not None:
                    self.edge_index = torch.cat([self.edge_index, new_edge_index], dim=1)
                    self.edge_attr = torch.cat([self.edge_attr, new_edge_attr], dim=0)
                else:
                    self.edge_index = new_edge_index
                    self.edge_attr = new_edge_attr
                
                print(f"添加了 {len(new_edges)} 条通路关系边")
            else:
                print("没有添加通路关系边")
        else:
            print("未提供通路关系文件，跳过添加通路关系边")
    
    def create_pyg_data(self):
        """
        创建PyTorch Geometric数据对象
        """
        if self.node_features is None or self.edge_index is None:
            raise ValueError("节点特征或边索引为空，请先构建知识图谱")
        
        data = Data(x=self.node_features, edge_index=self.edge_index, edge_attr=self.edge_attr)
        
        # 添加节点映射和类型信息
        data.node_map = self.node_map
        data.reverse_node_map = self.reverse_node_map
        data.node_types = self.node_types
        data.edge_types = self.edge_types
        
        return data
    
    def save_graph(self, output_dir="../graph_data"):
        """
        保存知识图谱数据
        """
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # 保存节点特征
        node_df = pd.DataFrame({
            'node_id': list(self.reverse_node_map.keys()),
            'name': list(self.reverse_node_map.values()),
            'type': [self.node_types[i] for i in self.reverse_node_map.keys()]
        })
        node_df.to_csv(os.path.join(output_dir, 'nodes.csv'), index=False)
        
        # 保存边
        edge_list = []
        for i in range(self.edge_index.shape[1]):
            source = self.reverse_node_map[self.edge_index[0, i].item()]
            target = self.reverse_node_map[self.edge_index[1, i].item()]
            weight = self.edge_attr[i].item()
            edge_type = self.edge_types.get(i, 'Unknown')
            edge_list.append([source, target, edge_type, weight])
        
        edge_df = pd.DataFrame(edge_list, columns=['source', 'target', 'edge_type', 'weight'])
        edge_df.to_csv(os.path.join(output_dir, 'edges.csv'), index=False)
        
        print(f"知识图谱数据已保存到 {output_dir}")
        
        return node_df, edge_df
    
    def to_networkx(self):
        """
        转换为NetworkX图对象，用于可视化
        """
        G = nx.DiGraph()
        
        # 添加节点
        for node_id, name in self.reverse_node_map.items():
            G.add_node(node_id, name=name, type=self.node_types.get(node_id, 'Unknown'))
        
        # 添加边
        for i in range(self.edge_index.shape[1]):
            source = self.edge_index[0, i].item()
            target = self.edge_index[1, i].item()
            weight = self.edge_attr[i].item()
            edge_type = self.edge_types.get(i, 'Unknown')
            G.add_edge(source, target, weight=weight, type=edge_type)
        
        return G

# 示例用法
if __name__ == "__main__":
    # 创建知识图谱构建器
    kg_builder = KnowledgeGraphBuilder(data_dir="../depmap data")
    
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
    kg_builder.save_graph()
    
    print("知识图谱构建完成！")
    print(f"节点数量: {data.num_nodes}")
    print(f"边数量: {data.num_edges}")