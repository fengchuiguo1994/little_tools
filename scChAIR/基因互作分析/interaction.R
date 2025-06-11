####### 2.1 代谢网络的构建
# 加载igraph包
library(igraph)
# 假设代谢反应数据包含反应物和产物
reactions <- data.frame(from = c("glucose", "pyruvate", "glucose"), to = c("pyruvate", "lactate", "acetyl-CoA"))
# 创建代谢网络图
metabolic_network <- graph_from_data_frame(reactions, directed = TRUE)
# 绘制代谢网络
plot(metabolic_network, vertex.size=30, vertex.label.cex=1.5, edge.arrow.size=0.5, main="Metabolic Network")

########### 2.2 代谢网络分析
# 计算每个节点的度中心性（即每个代谢物参与的反应数量）
degree_centrality <- degree(metabolic_network)
print(degree_centrality)
# 查找网络中度中心性最大的节点（重要代谢物）
important_metabolite <- names(degree_centrality[which.max(degree_centrality)])
print(important_metabolite)

############ 3.1 基因调控网络的构建
# 假设基因调控数据包含基因和转录因子的关系
regulatory_relations <- data.frame(gene = c("GeneA", "GeneB", "GeneC"), regulator = c("TF1", "TF2", "TF1"))
# 创建基因调控网络图
gene_regulatory_network <- graph_from_data_frame(regulatory_relations, directed = TRUE)
# 绘制基因调控网络
plot(gene_regulatory_network, vertex.size=30, vertex.label.cex=1.5, edge.arrow.size=0.5, main="Gene Regulatory Network", edge.color="blue")

########### 3.2 基因调控网络分析，我们可以计算基因调控网络的拓扑特性，如节点的度中心性、网络的模块化等。
# 计算基因节点的度中心性
gene_degree_centrality <- degree(gene_regulatory_network, v = V(gene_regulatory_network)$name)
print(gene_degree_centrality)
# 查找度中心性最大的基因（重要基因）
important_gene <- names(gene_degree_centrality[which.max(gene_degree_centrality)])
print(important_gene)

############# 4.1 使用ggraph绘制复杂网络。ggraph包是基于ggplot2的，用于网络图的可视化，支持更为复杂的网络布局和样式。# 加载ggraph包
library(ggraph)
# 假设我们要绘制的是一个组合的代谢网络和基因调控网络
combined_network <- graph.union(metabolic_network, gene_regulatory_network)
# 使用ggraph绘制网络
ggraph(combined_network, layout = 'fr') +
  geom_edge_link(aes(color = "gray"), alpha = 0.5) +
  geom_node_point(aes(color = 'lightblue', size = degree(combined_network))) +
  geom_node_text(aes(label = name), vjust = 1.5) +
  theme_void() +
  labs(title = "Combined Metabolic and Gene Regulatory Network")

############## 5.1 模块化分析。模块化分析帮助我们识别网络中的子网络（模块），这些模块通常代表了生物学功能相关的基因或代谢物。# 使用igraph的社区发现算法进行模块化分析
community <- cluster_walktrap(gene_regulatory_network)
plot(community, gene_regulatory_network, main="Community Detection in Gene Regulatory Network")

################ 5.2 基因集富集分析。基因集富集分析（GSEA）有助于识别在特定条件下显著变化的基因集。我们可以使用clusterProfiler包来进行这种分析。# 加载clusterProfiler包
library(clusterProfiler)
# 假设gene_list是基因的差异表达数据，包含基因名称和相应的log2FoldChange值
gene_list <- read.csv("gene_expression_results.csv")
gene_list <- sort(gene_list$log2FoldChange, decreasing = TRUE)
# 进行GSEA分析
gsea_result <- gseKEGG(geneList = gene_list, organism = "hsa")
# 查看结果
summary(gsea_result)
