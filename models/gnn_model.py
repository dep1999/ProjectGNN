import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, GATConv, global_mean_pool
from torch_geometric.data import Data, Batch

class GCNModel(nn.Module):
    """
    图卷积神经网络模型，用于KRAS G12C突变型癌症耐药基因预测
    """
    def __init__(self, num_node_features, hidden_channels, num_classes):
        super(GCNModel, self).__init__()
        self.conv1 = GCNConv(num_node_features, hidden_channels)
        self.conv2 = GCNConv(hidden_channels, hidden_channels)
        self.conv3 = GCNConv(hidden_channels, hidden_channels)
        self.lin = nn.Linear(hidden_channels, num_classes)
        self.dropout = nn.Dropout(0.3)
        
    def forward(self, x, edge_index, edge_weight=None, batch=None):
        # 第一层图卷积
        x = self.conv1(x, edge_index, edge_weight)
        x = F.relu(x)
        x = self.dropout(x)
        
        # 第二层图卷积
        x = self.conv2(x, edge_index, edge_weight)
        x = F.relu(x)
        x = self.dropout(x)
        
        # 第三层图卷积
        x = self.conv3(x, edge_index, edge_weight)
        x = F.relu(x)
        
        # 如果是图分类任务，需要池化
        if batch is not None:
            x = global_mean_pool(x, batch)
            
        # 输出层
        x = self.lin(x)
        
        return x
    
    def get_embeddings(self, x, edge_index, edge_weight=None):
        """获取节点嵌入向量"""
        x = self.conv1(x, edge_index, edge_weight)
        x = F.relu(x)
        
        x = self.conv2(x, edge_index, edge_weight)
        x = F.relu(x)
        
        x = self.conv3(x, edge_index, edge_weight)
        
        return x

class GATModel(nn.Module):
    """
    图注意力网络模型，用于KRAS G12C突变型癌症耐药基因预测
    """
    def __init__(self, num_node_features, hidden_channels, num_classes, heads=4):
        super(GATModel, self).__init__()
        self.conv1 = GATConv(num_node_features, hidden_channels, heads=heads, dropout=0.3)
        # 第二层输入维度为 hidden_channels * heads
        self.conv2 = GATConv(hidden_channels * heads, hidden_channels, heads=heads, dropout=0.3)
        # 第三层输入维度为 hidden_channels * heads，输出维度为 hidden_channels
        self.conv3 = GATConv(hidden_channels * heads, hidden_channels, heads=1, dropout=0.3)
        self.lin = nn.Linear(hidden_channels, num_classes)
        self.dropout = nn.Dropout(0.3)
        
    def forward(self, x, edge_index, edge_weight=None, batch=None):
        # 第一层图注意力
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = self.dropout(x)
        
        # 第二层图注意力
        x = self.conv2(x, edge_index)
        x = F.relu(x)
        x = self.dropout(x)
        
        # 第三层图注意力
        x = self.conv3(x, edge_index)
        x = F.relu(x)
        
        # 如果是图分类任务，需要池化
        if batch is not None:
            x = global_mean_pool(x, batch)
            
        # 输出层
        x = self.lin(x)
        
        return x
    
    def get_embeddings(self, x, edge_index):
        """获取节点嵌入向量"""
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        
        x = self.conv2(x, edge_index)
        x = F.relu(x)
        
        x = self.conv3(x, edge_index)
        
        return x

class MultiObjectiveGNN(nn.Module):
    """
    多目标优化的GNN模型，用于KRAS G12C突变型癌症耐药基因预测
    同时优化多个目标：耐药性预测、结构新颖性、靶向可行性等
    """
    def __init__(self, num_node_features, hidden_channels, num_tasks=3):
        super(MultiObjectiveGNN, self).__init__()
        # 基础GNN模型，可以是GCN或GAT
        self.gnn = GCNModel(num_node_features, hidden_channels, hidden_channels)
        
        # 多任务输出层
        self.task_heads = nn.ModuleList()
        for _ in range(num_tasks):
            self.task_heads.append(nn.Sequential(
                nn.Linear(hidden_channels, hidden_channels),
                nn.ReLU(),
                nn.Dropout(0.3),
                nn.Linear(hidden_channels, 1)
            ))
        
    def forward(self, x, edge_index, edge_weight=None, batch=None):
        # 获取节点嵌入
        embeddings = self.gnn.get_embeddings(x, edge_index, edge_weight)
        
        # 如果是图分类任务，需要池化
        if batch is not None:
            embeddings = global_mean_pool(embeddings, batch)
        
        # 多任务输出
        outputs = []
        for task_head in self.task_heads:
            outputs.append(task_head(embeddings))
        
        return outputs, embeddings