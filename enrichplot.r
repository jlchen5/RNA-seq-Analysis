library(clusterProfiler)
library(dplyr)

data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

eego <- enrichGO(gene          = gene,
                 universe      = names(geneList),
                 OrgDb         = org.Hs.eg.db,
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

df <- data.frame(eego) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 6) %>%
  arrange(desc(pvalue))

ratio <- lapply(df$GeneRatio,function(x){as.numeric(eval(parse(text = x)))}) %>% unlist()
df$ratio <- ratio

df$Description <- factor(df$Description,levels = df$Description)


ggplot(df) +
  ggforce::geom_link(aes(x = 0,y = Description,
                         xend = -log10(pvalue),yend = Description,
                         alpha = stat(index),
                         color = ONTOLOGY,
                         size = after_stat(index)),
  n = 500,
  # color = "#FF0033",
  show.legend = F) +
  geom_point(aes(x = -log10(pvalue),y = Description),
             color = "black",
             fill = "white",size = 6,shape = 21) +
  geom_line(aes(x = ratio*100,y = Description,group = 1),
            orientation = "y",linewidth = 1,color = "#FFCC00") +
  scale_x_continuous(sec.axis = sec_axis(~./100,
                                         labels = scales::label_percent(),
                                         name = "Percent of geneRatio")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        axis.text = element_text(color = "black")) +
  ylab("") + xlab("-log10 Pvalue") +
  facet_wrap(~ONTOLOGY,scales = "free",ncol = 1) +
  scale_color_brewer(palette = "Set1")

