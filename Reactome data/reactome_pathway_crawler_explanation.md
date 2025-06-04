# Reactome通路爬取脚本详细解释

## 脚本概述

这个Python脚本用于从Reactome数据库中爬取KRAS相关通路数据，特别是MAPK/RAS信号通路的信息。Reactome是一个开放的、经过同行评审的生物学通路数据库，提供了详细的分子相互作用和生物学过程信息。脚本通过Reactome的ContentService API获取数据，并将结果保存为JSON和CSV文件格式。

## 导入的库

```python
import os                # 用于操作文件系统和路径
import json              # 用于处理JSON格式数据
import pandas as pd      # 用于数据处理和分析
import requests          # 用于发送HTTP请求
from tqdm import tqdm    # 用于显示进度条
import time              # 用于时间相关操作
import logging           # 用于日志记录
```

## 全局配置

```python
# 配置日志
logging.basicConfig(
    level=logging.INFO,                                      # 设置日志级别为INFO
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',  # 设置日志格式
    handlers=[
        logging.FileHandler("reactome_crawler.log"),         # 将日志输出到文件
        logging.StreamHandler()                             # 将日志输出到控制台
    ]
)

logger = logging.getLogger("Reactome爬虫")                   # 创建日志记录器

# 创建数据存储目录
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pathway_data")  # 设置输出目录为脚本所在目录下的pathway_data文件夹
os.makedirs(OUTPUT_DIR, exist_ok=True)                      # 创建输出目录，如果已存在则不报错

# Reactome API基础URL
BASE_URL = "https://reactome.org/ContentService/data"         # Reactome API的基础URL

# 关键基因和通路ID
TARGET_GENES = ["KRAS", "NRAS", "HRAS", "BRAF", "MEK1", "MEK2", "ERK1", "ERK2", "PIK3CA", "AKT1"]  # 目标基因列表
MAPK_PATHWAY_ID = "R-HSA-5673001"  # MAPK信号通路ID
RAS_PATHWAY_ID = "R-HSA-162582"    # RAS信号通路ID
```

## 函数详解

### 1. get_gene_info

```python
def get_gene_info(gene_symbol):
    """获取基因的Reactome信息"""
```

**功能**：获取指定基因在Reactome数据库中的信息。

**参数**：
- `gene_symbol`：基因符号（如"KRAS"）

**返回值**：包含基因信息的JSON数据，如果请求失败则返回None。

**工作流程**：
1. 构建API请求URL
2. 发送GET请求
3. 检查响应状态码，如果成功则返回JSON数据
4. 如果发生异常，记录错误并返回None

### 2. get_pathway_info

```python
def get_pathway_info(pathway_id):
    """获取通路详细信息"""
```

**功能**：获取指定通路ID的详细信息。

**参数**：
- `pathway_id`：通路ID（如"R-HSA-5673001"）

**返回值**：包含通路详细信息的JSON数据，如果请求失败则返回None。

**工作流程**：
1. 构建API请求URL，访问通路包含的事件
2. 发送GET请求
3. 检查响应状态码，如果成功则返回JSON数据
4. 如果发生异常，记录错误并返回None

### 3. get_pathway_participants

```python
def get_pathway_participants(pathway_id):
    """获取通路参与者信息"""
```

**功能**：获取指定通路的参与分子信息。

**参数**：
- `pathway_id`：通路ID

**返回值**：包含通路参与分子信息的JSON数据，如果请求失败则返回None。

**工作流程**：
1. 构建API请求URL，访问通路参与分子
2. 发送GET请求
3. 检查响应状态码，如果成功则返回JSON数据
4. 如果发生异常，记录错误并返回None

### 4. get_entity_interactions

```python
def get_entity_interactions(entity_id):
    """获取实体的相互作用信息"""
```

**功能**：获取指定实体ID的相互作用信息。

**参数**：
- `entity_id`：实体ID

**返回值**：包含实体相互作用信息的JSON数据，如果请求失败则返回None。

**工作流程**：
1. 构建API请求URL，访问实体的相互作用详情
2. 发送GET请求
3. 检查响应状态码，如果成功则返回JSON数据
4. 如果发生异常，记录错误并返回None

### 5. get_pathway_hierarchy

```python
def get_pathway_hierarchy(pathway_id):
    """获取通路层次结构"""
```

**功能**：获取指定通路的层次结构信息。

**参数**：
- `pathway_id`：通路ID

**返回值**：包含通路层次结构的JSON数据，如果请求失败则返回None。

**工作流程**：
1. 构建API请求URL，访问通路层次结构
2. 发送GET请求
3. 检查响应状态码，如果成功则返回JSON数据
4. 如果发生异常，记录错误并返回None

### 6. extract_pathway_relations

```python
def extract_pathway_relations(pathway_data):
    """从通路数据中提取关系信息"""
```

**功能**：从通路数据中提取实体间的关系信息。

**参数**：
- `pathway_data`：通路数据（通常是get_pathway_info函数的返回值）

**返回值**：包含关系信息的列表，每个关系是一个字典，包含源实体、目标实体、关系类型和事件ID。

**工作流程**：
1. 初始化一个空的关系列表
2. 检查通路数据是否有效
3. 如果通路数据是列表类型，遍历其中的每个事件
4. 对于每个事件，检查是否包含输入和输出实体
5. 如果包含，则为每对输入和输出实体创建一个关系
6. 返回关系列表

### 7. crawl_gene_pathways

```python
def crawl_gene_pathways():
    """爬取目标基因相关的通路信息"""
```

**功能**：爬取目标基因列表中所有基因相关的通路信息。

**返回值**：包含所有基因通路数据的字典，键为基因符号，值为该基因的通路数据。

**工作流程**：
1. 初始化一个空字典，用于存储所有基因的数据
2. 遍历目标基因列表
3. 对于每个基因，获取其Reactome信息
4. 提取基因参与的通路
5. 获取每个通路的详细信息和参与者
6. 提取通路关系
7. 将数据保存为JSON文件
8. 返回所有基因数据的字典

### 8. crawl_specific_pathways

```python
def crawl_specific_pathways():
    """爬取特定通路的详细信息"""
```

**功能**：爬取特定通路（MAPK和RAS信号通路）的详细信息。

**返回值**：包含特定通路数据的字典，键为通路ID，值为该通路的详细数据。

**工作流程**：
1. 定义要爬取的通路ID列表
2. 初始化一个空字典，用于存储通路数据
3. 遍历通路ID列表
4. 获取每个通路的详细信息、参与者和层次结构
5. 提取通路关系
6. 将数据保存为JSON和CSV文件
7. 返回通路数据字典

### 9. generate_network_files

```python
def generate_network_files():
    """生成网络分析文件"""
```

**功能**：生成用于网络分析的节点和边文件。

**工作流程**：
1. 查找所有关系CSV文件
2. 读取并合并所有关系数据
3. 去除重复关系
4. 保存合并后的关系数据为CSV文件
5. 提取所有节点（源实体和目标实体）
6. 保存节点数据为CSV文件
7. 记录生成的节点和关系数量

### 10. main

```python
def main():
    """主函数"""
```

**功能**：脚本的主函数，协调整个数据爬取过程。

**工作流程**：
1. 爬取目标基因相关的通路信息
2. 爬取特定通路的详细信息
3. 生成网络分析文件
4. 记录爬取完成的信息

## 执行流程

当脚本被直接运行时，会执行main函数：

```python
if __name__ == "__main__":
    main()
```

## 数据处理流程

1. **获取基因通路**：首先获取目标基因（KRAS、NRAS等）相关的通路信息。
2. **获取特定通路**：获取MAPK和RAS信号通路的详细信息。
3. **提取关系**：从通路数据中提取实体间的关系信息。
4. **生成网络文件**：合并所有关系数据，生成用于网络分析的节点和边文件。

## 输出文件

脚本会在`pathway_data`目录下生成以下文件：

1. `{基因}_pathway_data.json`：包含每个目标基因的通路数据。
2. `pathway_{通路ID}_data.json`：包含特定通路的详细数据。
3. `pathway_{通路ID}_relations.csv`：包含特定通路的关系数据。
4. `combined_pathway_relations.csv`：包含所有通路的合并关系数据。
5. `pathway_nodes.csv`：包含所有节点（实体）的数据。

此外，还会生成一个日志文件`reactome_crawler.log`，记录脚本执行过程中的信息和错误。