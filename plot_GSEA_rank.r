# GSEA分析GSEA的输入数据为差异基因排序列表，排序指标为log2FC值：

## 1.读取数据

rm(list=ls())
# 加载R包
library(ggplot2)
library(tibble)
library(ggrepel)
library(tidyverse)
library(dplyr)
library(patchwork)

##### 01、加载数据

# 加载：KRAS vs NS siRNA（log2FC）
group1 <- read.csv("./data/science.adk0775_data_s1.csv" )
head(group1)

# 提取FDR值和FC值
group1 <- group1[,c("external_gene_name","logFC","FDR" )]
group1 <- group1[group1$external_gene_name!="", ]
group1 <- distinct(group1,external_gene_name,.keep_all = T)
rownames(group1) <- group1$external_gene_name
head(group1)

# 增加一列上下调：log2FC > 0.5, adj. p < 0.05
group1$g1 <- "normal"
group1$g1[group1$logFC >0.5 & group1$FDR < 0.05 ] <- "up"
group1$g1[group1$logFC < -0.5 & group1$FDR < 0.05 ] <- "down"
table(group1$g1)
# down normal     up 
# 675  12876   1043 

# 获得差异基因
DEGs <- group1[group1$g1!="normal",]
head(DEGs)

# 得到排序列表
DEGs <- DEGs[order(DEGs$logFC, decreasing = T),]
genelist <- DEGs$logFC
names(genelist) <- DEGs$external_gene_name
head(genelist)
# NPPC   SMTNL2    KCNK3     CHPF PCDHGA12  EXOC3L1 
# 2.860472 2.712068 2.551034 2.536558 2.499657 2.446140

tail(genelist)

## 2、读取 50 Hallmark gene sets 通路并富集：

## === HALLMARK通路富集
geneset <- read.gmt("data/h.all.v2024.1.Hs.symbols.gmt")
table(geneset$term)
geneset$term <- gsub(pattern = "HALLMARK_","", geneset$term)

# 运行,输出全部结果
egmt <- GSEA(genelist, TERM2GENE=geneset, pvalueCutoff = 1, minGSSize = 1, maxGSSize = 500000)
colnames(egmt@result)
head(egmt[, 1:6])

## 3、使用 ggplot2 定制化绘图
# 绘图
data <- egmt[,c("ID", "NES","setSize","pvalue")]
data$setSize_1 <- data$setSize/10
head(data)

# 给Y轴的通路名设置为因子，排序
data <- data[order(data$NES, decreasing = T),]
data$ID <- factor(data$ID, levels = data$ID)
data$xlab <- 1:49
head(data)
summary(data$NES)
summary(data$setSize_1)

# 图中标出的通路名字
label <- c( "KRAS_SIGNALING_DN", "INTERFERON_ALPHA_RESPONSE",
  "UV_RESPONSE_DN", "EMT", "TNFA_SIGNALING_VIA_NFKB", "MYC_TARGETS_V2",
  "KRAS_SIGNALING_UP", "MYC_TARGETS_V1", "G2M_CHECKPOINT", "E2F_TARGETS" )

# 没有EMT通路
data_label <- data[data$ID %in% label,]
data_label
# 图中通路的颜色
data_label$col <- c("black", "grey80","grey80","grey80","grey80","#ffb882","grey80","#ff5eff","#ff5eff")


p <- ggplot(data = data, aes(x = xlab, y = NES)) +
  geom_point(aes(size = setSize_1, alpha = -log10(pvalue)), shape = 21, stroke = 0.7,fill = "#0000ff", colour = "black") + # stroke：设置点的边框宽度。
  scale_size_continuous(range = c(0.2, 8)) + 
  xlab(label = "Hallmark gene sets") + 
  ylab(label = "Normalized enrichment score (NES)") + 
  theme_classic(base_size = 15) + 
  scale_x_continuous(breaks = seq(0, 50, by = 10), labels = seq(0, 50, by = 10) ) + # 设置x轴的刻度线和刻度标签
  scale_y_continuous( breaks = seq(-4, 2.3, by = 1), labels = seq(-4, 2.3, by = 1) ) + 
  guides(size = guide_legend(title = "Detection\n(proportion)"),
         alpha = guide_legend(title = "Significance\n(-log10 p-val.)") ) +  # 修改图例标题
  theme(
    axis.line = element_line(color = "black", size = 0.6), # 加粗x轴和y轴的线条
    axis.text = element_text(face = "bold"), # 加粗x轴和y轴的标签
    axis.title = element_text( size = 13)    # 加粗x轴和y轴的标题
  )

p


# 添加label：vjust（垂直调整）或hjust（水平调整）
p3 <- p + 
  geom_text_repel(data= data_label, aes(x = xlab, y = NES,label = ID), size = 3, color = data_label$col, 
                  force=20,                # 将重叠的文本标签相互推开的强度。force 参数的值越大，标签之间的排斥力度也越大，这会导致标签在图中更分散地排列
                  point.padding = 0.5,     # 设置文本标签与对应点之间的最小距离
                  min.segment.length = 0,  # 长度大于0就可以添加引线
                  hjust = 1.2,             # 文本标签的右侧与指定位置对齐
                  segment.color="grey20",
                  segment.size=0.3,        # 设置引导线的粗细
                  segment.alpha=0.8,       # 文本标签中连接线段的透明度
                  nudge_y=-0.1             # 在y轴方向上微调标签位置
                  ) 
p3

# 保存
ggsave(filename = "Figure1C.pdf", plot = p3, width = 6.2, height = 5)








