##### import your data
forheatmap <- read.csv(file = 'forheatmap.csv',header = T)

##### change your rownames
rownames(forheatmap) <- forheatmap[,1]
forheatmap <- forheatmap[,-1]

##### plot heatmap
pheatmap::pheatmap(forheatmap,scale = 'row',cluster_col = F,show_rownames=T,angle_col = 0) 
