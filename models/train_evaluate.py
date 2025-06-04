import os
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import pandas as pd
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, f1_score, precision_recall_curve, auc
import matplotlib.pyplot as plt
import sys

# 添加项目根目录到系统路径
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# 导入自定义模型
from models.gnn_model import GCNModel, GATModel, MultiObjectiveGNN

class GNNTrainer:
    """
    GNN模型训练和评估类
    用于KRAS G12C突变型癌症耐药基因预测
    """
    def __init__(self, model_type='gcn', hidden_channels=64, num_classes=1, 
                 lr=0.001, weight_decay=5e-4, device=None):
        self.model_type = model_type
        self.hidden_channels = hidden_channels
        self.num_classes = num_classes
        self.lr = lr
        self.weight_decay = weight_decay
        
        # 设置设备
        if device is None:
            self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        else:
            self.device = device
            
        self.model = None
        self.optimizer = None
        self.criterion = None
        
    def initialize_model(self, num_node_features):
        """
        初始化模型
        """
        if self.model_type.lower() == 'gcn':
            self.model = GCNModel(num_node_features, self.hidden_channels, self.num_classes)
        elif self.model_type.lower() == 'gat':
            self.model = GATModel(num_node_features, self.hidden_channels, self.num_classes)
        elif self.model_type.lower() == 'multi':
            self.model = MultiObjectiveGNN(num_node_features, self.hidden_channels, num_tasks=self.num_classes)
        else:
            raise ValueError(f"不支持的模型类型: {self.model_type}")
        
        self.model = self.model.to(self.device)
        self.optimizer = optim.Adam(self.model.parameters(), lr=self.lr, weight_decay=self.weight_decay)
        
        # 二分类任务使用BCE损失，多分类任务使用交叉熵损失
        if self.num_classes == 1:
            self.criterion = nn.BCEWithLogitsLoss()
        else:
            self.criterion = nn.CrossEntropyLoss()
            
        print(f"初始化 {self.model_type} 模型，特征维度: {num_node_features}, 隐藏层维度: {self.hidden_channels}")
        
    def train(self, data, labels, train_mask, val_mask=None, epochs=200, patience=20, verbose=True):
        """
        训练模型
        
        参数:
        - data: PyG数据对象
        - labels: 节点标签
        - train_mask: 训练集掩码
        - val_mask: 验证集掩码
        - epochs: 训练轮数
        - patience: 早停耐心值
        - verbose: 是否打印训练过程
        """
        if self.model is None:
            raise ValueError("模型未初始化，请先调用initialize_model方法")
        
        # 将数据移动到设备
        data = data.to(self.device)
        labels = labels.to(self.device)
        
        # 早停设置
        best_val_loss = float('inf')
        best_epoch = 0
        no_improve = 0
        
        # 训练历史记录
        train_losses = []
        val_losses = []
        
        for epoch in range(epochs):
            # 训练模式
            self.model.train()
            self.optimizer.zero_grad()
            
            # 前向传播
            if self.model_type.lower() == 'multi':
                outputs, _ = self.model(data.x, data.edge_index, data.edge_attr)
                loss = 0
                for i, output in enumerate(outputs):
                    if i < len(labels):
                        loss += self.criterion(output[train_mask], labels[i][train_mask].view(-1, 1))
            else:
                output = self.model(data.x, data.edge_index, data.edge_attr)
                if self.num_classes == 1:
                    loss = self.criterion(output[train_mask], labels[train_mask].view(-1, 1))
                else:
                    loss = self.criterion(output[train_mask], labels[train_mask])
            
            # 反向传播
            loss.backward()
            self.optimizer.step()
            
            train_losses.append(loss.item())
            
            # 验证
            if val_mask is not None:
                val_loss = self.evaluate(data, labels, val_mask)
                val_losses.append(val_loss)
                
                # 早停检查
                if val_loss < best_val_loss:
                    best_val_loss = val_loss
                    best_epoch = epoch
                    no_improve = 0
                else:
                    no_improve += 1
                    if no_improve >= patience:
                        if verbose:
                            print(f"早停: {patience} 轮未改善，最佳轮次: {best_epoch+1}")
                        break
                
                if verbose and (epoch + 1) % 10 == 0:
                    print(f"轮次: {epoch+1}/{epochs}, 训练损失: {loss.item():.4f}, 验证损失: {val_loss:.4f}")
            else:
                if verbose and (epoch + 1) % 10 == 0:
                    print(f"轮次: {epoch+1}/{epochs}, 训练损失: {loss.item():.4f}")
        
        # 绘制训练曲线
        if val_mask is not None:
            self.plot_training_curve(train_losses, val_losses)
        
        return train_losses, val_losses
    
    def evaluate(self, data, labels, mask):
        """
        评估模型
        
        参数:
        - data: PyG数据对象
        - labels: 节点标签
        - mask: 评估集掩码
        
        返回:
        - loss: 损失值
        """
        self.model.eval()
        with torch.no_grad():
            if self.model_type.lower() == 'multi':
                outputs, _ = self.model(data.x, data.edge_index, data.edge_attr)
                loss = 0
                for i, output in enumerate(outputs):
                    if i < len(labels):
                        loss += self.criterion(output[mask], labels[i][mask].view(-1, 1))
            else:
                output = self.model(data.x, data.edge_index, data.edge_attr)
                if self.num_classes == 1:
                    loss = self.criterion(output[mask], labels[mask].view(-1, 1))
                else:
                    loss = self.criterion(output[mask], labels[mask])
        
        return loss.item()
    
    def predict(self, data):
        """
        预测
        
        参数:
        - data: PyG数据对象
        
        返回:
        - predictions: 预测结果
        - embeddings: 节点嵌入
        """
        self.model.eval()
        with torch.no_grad():
            if self.model_type.lower() == 'multi':
                predictions, embeddings = self.model(data.x, data.edge_index, data.edge_attr)
            else:
                predictions = self.model(data.x, data.edge_index, data.edge_attr)
                # 获取节点嵌入
                if hasattr(self.model, 'get_embeddings'):
                    embeddings = self.model.get_embeddings(data.x, data.edge_index, data.edge_attr)
                else:
                    embeddings = None
        
        return predictions, embeddings
    
    def calculate_metrics(self, y_true, y_pred, threshold=0.5):
        """
        计算评估指标
        
        参数:
        - y_true: 真实标签
        - y_pred: 预测分数
        - threshold: 二分类阈值
        
        返回:
        - metrics: 评估指标字典
        """
        # 转换为numpy数组
        if isinstance(y_true, torch.Tensor):
            y_true = y_true.cpu().numpy()
        if isinstance(y_pred, torch.Tensor):
            y_pred = y_pred.cpu().numpy()
        
        # 二分类预测
        if self.num_classes == 1:
            y_pred_binary = (y_pred > threshold).astype(int)
            
            # 计算AUC
            try:
                roc_auc = roc_auc_score(y_true, y_pred)
            except:
                roc_auc = 0.5
            
            # 计算F1分数
            f1 = f1_score(y_true, y_pred_binary)
            
            # 计算PR曲线下面积
            precision, recall, _ = precision_recall_curve(y_true, y_pred)
            pr_auc = auc(recall, precision)
            
            metrics = {
                'roc_auc': roc_auc,
                'f1_score': f1,
                'pr_auc': pr_auc
            }
        else:
            # 多分类预测
            y_pred_class = np.argmax(y_pred, axis=1)
            
            # 计算F1分数（宏平均）
            f1 = f1_score(y_true, y_pred_class, average='macro')
            
            metrics = {
                'f1_score': f1
            }
        
        return metrics
    
    def plot_training_curve(self, train_losses, val_losses=None):
        """
        绘制训练曲线
        
        参数:
        - train_losses: 训练损失列表
        - val_losses: 验证损失列表
        """
        plt.figure(figsize=(10, 6))
        plt.plot(train_losses, label='训练损失')
        if val_losses:
            plt.plot(val_losses, label='验证损失')
        plt.xlabel('轮次')
        plt.ylabel('损失')
        plt.title('训练曲线')
        plt.legend()
        plt.grid(True)
        
        # 保存图像
        os.makedirs('../results', exist_ok=True)
        plt.savefig(f'../results/training_curve_{self.model_type}.png')
        plt.close()
    
    def save_model(self, path):
        """
        保存模型
        
        参数:
        - path: 保存路径
        """
        os.makedirs(os.path.dirname(path), exist_ok=True)
        torch.save({
            'model_state_dict': self.model.state_dict(),
            'optimizer_state_dict': self.optimizer.state_dict(),
            'model_type': self.model_type,
            'hidden_channels': self.hidden_channels,
            'num_classes': self.num_classes
        }, path)
        print(f"模型已保存到 {path}")
    
    def load_model(self, path, num_node_features):
        """
        加载模型
        
        参数:
        - path: 模型路径
        - num_node_features: 节点特征维度
        """
        checkpoint = torch.load(path, map_location=self.device)
        
        # 更新模型参数
        self.model_type = checkpoint['model_type']
        self.hidden_channels = checkpoint['hidden_channels']
        self.num_classes = checkpoint['num_classes']
        
        # 初始化模型
        self.initialize_model(num_node_features)
        
        # 加载模型参数
        self.model.load_state_dict(checkpoint['model_state_dict'])
        self.optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        
        print(f"模型已从 {path} 加载")

# 示例用法
if __name__ == "__main__":
    # 假设我们已经有了图数据
    from data_processing.knowledge_graph_builder import KnowledgeGraphBuilder
    
    # 创建知识图谱
    kg_builder = KnowledgeGraphBuilder()
    # ... 加载数据和构建图谱 ...
    
    # 创建PyG数据对象
    data = kg_builder.create_pyg_data()
    
    # 创建标签（示例）
    num_nodes = data.num_nodes
    labels = torch.zeros(num_nodes)
    # ... 设置已知耐药基因的标签为1 ...
    
    # 划分训练集和测试集
    indices = list(range(num_nodes))
    train_indices, test_indices = train_test_split(indices, test_size=0.2, random_state=42)
    train_mask = torch.zeros(num_nodes, dtype=torch.bool)
    test_mask = torch.zeros(num_nodes, dtype=torch.bool)
    train_mask[train_indices] = True
    test_mask[test_indices] = True
    
    # 初始化训练器
    trainer = GNNTrainer(model_type='gcn', hidden_channels=64, num_classes=1)
    trainer.initialize_model(data.num_features)
    
    # 训练模型
    trainer.train(data, labels, train_mask, test_mask, epochs=200, patience=20)
    
    # 预测
    predictions, embeddings = trainer.predict(data)
    
    # 计算测试集指标
    metrics = trainer.calculate_metrics(labels[test_mask], predictions[test_mask])
    print(f"测试集指标: {metrics}")
    
    # 保存模型
    trainer.save_model('../models/gcn_model.pt')