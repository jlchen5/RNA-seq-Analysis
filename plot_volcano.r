##### import R packages
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
rm(list=ls())

##### read data
dataset <- read.csv("sample.csv",header = T)

##### set pvalue and logFC
cut_off_pvalue = 0.05
cut_off_logFC = 1

##### save‘up’, ‘Down’ and ‘Stable’ to change column;change column for color setting
dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_logFC, 
                        ifelse(dataset$log2FoldChange> cut_off_logFC ,'Up','Down'),
                        'Stable')
##### plot
p <- ggplot(
  dataset, 
  aes(x = log2FoldChange, 
      y = -log10(pvalue), 
      colour=change)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  
  # line
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  
  # axis
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_bw()+
  
  # theme
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()
  )

##### mark genes
dataset$label = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= 2, as.character(dataset$gene),"")
p + geom_text_repel(data = dataset, aes(x = dataset$log2FoldChange, 
                                      y = -log10(dataset$pvalue), 
                                      label = label),
                  size = 3,box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"), 
                  segment.color = "black", 
                  show.legend = FALSE)
