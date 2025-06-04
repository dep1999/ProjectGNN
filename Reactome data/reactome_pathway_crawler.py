#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Reactome数据爬取脚本
用于获取KRAS相关通路数据，特别是MAPK/RAS信号通路的信息
"""

import os
import json
import pandas as pd
import requests
from tqdm import tqdm
import time
import logging

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("reactome_crawler.log"),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger("Reactome爬虫")

# 创建数据存储目录
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pathway_data")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Reactome API基础URL
BASE_URL = "https://reactome.org/ContentService/data"

# 关键基因和通路ID
TARGET_GENES = ["KRAS", "NRAS", "HRAS", "BRAF", "MEK1", "MEK2", "ERK1", "ERK2", "PIK3CA", "AKT1"]
MAPK_PATHWAY_ID = "R-HSA-5673001"  # MAPK信号通路ID
RAS_PATHWAY_ID = "R-HSA-162582"    # RAS信号通路ID

def get_gene_info(gene_symbol):
    """获取基因的Reactome信息"""
    url = f"{BASE_URL}/query/{gene_symbol}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        logger.error(f"获取{gene_symbol}信息失败: {e}")
        return None

def get_pathway_info(pathway_id):
    """获取通路详细信息"""
    url = f"{BASE_URL}/pathway/{pathway_id}/containedEvents"
    try:
        response = requests.get(url)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        logger.error(f"获取通路{pathway_id}信息失败: {e}")
        return None

def get_pathway_participants(pathway_id):
    """获取通路参与者信息"""
    url = f"{BASE_URL}/pathway/{pathway_id}/participatingMolecules"
    try:
        response = requests.get(url)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        logger.error(f"获取通路{pathway_id}参与者信息失败: {e}")
        return None

def get_entity_interactions(entity_id):
    """获取实体的相互作用信息"""
    url = f"{BASE_URL}/interactors/static/molecule/{entity_id}/details"
    try:
        response = requests.get(url)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        logger.error(f"获取实体{entity_id}相互作用信息失败: {e}")
        return None

def get_pathway_hierarchy(pathway_id):
    """获取通路层次结构"""
    url = f"{BASE_URL}/pathways/hierarchy/{pathway_id}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        logger.error(f"获取通路{pathway_id}层次结构失败: {e}")
        return None

def extract_pathway_relations(pathway_data):
    """从通路数据中提取关系信息"""
    relations = []
    
    if not pathway_data:
        return relations
    
    # 处理不同类型的通路数据
    if isinstance(pathway_data, list):
        for event in pathway_data:
            if 'inputs' in event and 'outputs' in event:
                for input_entity in event.get('inputs', []):
                    for output_entity in event.get('outputs', []):
                        relations.append({
                            'source': input_entity.get('displayName', ''),
                            'target': output_entity.get('displayName', ''),
                            'relation_type': event.get('displayName', 'interacts with'),
                            'event_id': event.get('stId', '')
                        })
    
    return relations

def crawl_gene_pathways():
    """爬取目标基因相关的通路信息"""
    all_gene_data = {}
    
    for gene in tqdm(TARGET_GENES, desc="爬取基因通路信息"):
        logger.info(f"开始爬取{gene}相关通路信息")
        gene_info = get_gene_info(gene)
        
        if not gene_info:
            continue
            
        gene_data = {
            'info': gene_info,
            'pathways': [],
            'interactions': []
        }
        
        # 获取基因参与的通路
        for item in gene_info:
            if item.get('exactType') == 'Pathway':
                pathway_id = item.get('stId')
                pathway_info = get_pathway_info(pathway_id)
                participants = get_pathway_participants(pathway_id)
                
                pathway_data = {
                    'id': pathway_id,
                    'name': item.get('displayName', ''),
                    'info': pathway_info,
                    'participants': participants
                }
                
                gene_data['pathways'].append(pathway_data)
                
                # 提取通路关系
                relations = extract_pathway_relations(pathway_info)
                gene_data['interactions'].extend(relations)
        
        all_gene_data[gene] = gene_data
        
        # 保存单个基因数据
        with open(os.path.join(OUTPUT_DIR, f"{gene}_pathway_data.json"), 'w', encoding='utf-8') as f:
            json.dump(gene_data, f, ensure_ascii=False, indent=2)
            
        logger.info(f"完成{gene}相关通路信息爬取，共获取{len(gene_data['pathways'])}个通路")
        
        # 避免请求过于频繁
        time.sleep(1)
    
    return all_gene_data

def crawl_specific_pathways():
    """爬取特定通路的详细信息"""
    pathway_ids = [MAPK_PATHWAY_ID, RAS_PATHWAY_ID]
    pathway_data = {}
    
    for pathway_id in tqdm(pathway_ids, desc="爬取特定通路信息"):
        logger.info(f"开始爬取通路{pathway_id}详细信息")
        
        # 获取通路详细信息
        info = get_pathway_info(pathway_id)
        participants = get_pathway_participants(pathway_id)
        hierarchy = get_pathway_hierarchy(pathway_id)
        
        # 提取通路关系
        relations = extract_pathway_relations(info)
        
        pathway_data[pathway_id] = {
            'info': info,
            'participants': participants,
            'hierarchy': hierarchy,
            'relations': relations
        }
        
        # 保存通路数据
        with open(os.path.join(OUTPUT_DIR, f"pathway_{pathway_id}_data.json"), 'w', encoding='utf-8') as f:
            json.dump(pathway_data[pathway_id], f, ensure_ascii=False, indent=2)
        
        # 将关系数据转换为CSV格式
        if relations:
            relations_df = pd.DataFrame(relations)
            relations_df.to_csv(os.path.join(OUTPUT_DIR, f"pathway_{pathway_id}_relations.csv"), index=False)
        
        logger.info(f"完成通路{pathway_id}详细信息爬取")
        
        # 避免请求过于频繁
        time.sleep(1)
    
    return pathway_data

def generate_network_files():
    """生成网络分析文件"""
    # 合并所有关系数据
    all_relations = []
    relation_files = [f for f in os.listdir(OUTPUT_DIR) if f.endswith('_relations.csv')]
    
    for file in relation_files:
        file_path = os.path.join(OUTPUT_DIR, file)
        try:
            df = pd.read_csv(file_path)
            all_relations.append(df)
        except Exception as e:
            logger.error(f"读取文件{file}失败: {e}")
    
    if all_relations:
        # 合并所有关系
        combined_relations = pd.concat(all_relations, ignore_index=True)
        
        # 去除重复关系
        combined_relations = combined_relations.drop_duplicates()
        
        # 保存合并后的关系数据
        combined_relations.to_csv(os.path.join(OUTPUT_DIR, "combined_pathway_relations.csv"), index=False)
        
        # 生成节点文件
        nodes = set(combined_relations['source'].tolist() + combined_relations['target'].tolist())
        nodes_df = pd.DataFrame({'node_name': list(nodes)})
        nodes_df.to_csv(os.path.join(OUTPUT_DIR, "pathway_nodes.csv"), index=False)
        
        logger.info(f"生成网络文件完成，共有{len(nodes)}个节点和{len(combined_relations)}条关系")

def main():
    """主函数"""
    logger.info("开始Reactome数据爬取")
    
    # 爬取基因相关通路
    gene_data = crawl_gene_pathways()
    
    # 爬取特定通路
    pathway_data = crawl_specific_pathways()
    
    # 生成网络文件
    generate_network_files()
    
    logger.info("Reactome数据爬取完成")

if __name__ == "__main__":
    main()