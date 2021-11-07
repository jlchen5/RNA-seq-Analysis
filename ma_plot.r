##### Load packages #####
rm(list = ls())

library(tidyverse)
library(magrittr)
library(glue)
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggsci)
library(scales)
library(latex2exp)

##### Load data #####

degdata <- fread("miR_DEG.csv")
colnames(degdata)[1] <- "gene"

degdata[is.na(padj), padj := 1][]
degdata[, baseMean := log2(baseMean)][]

##### Plot #####

degdata[, type := "ns"][]
degdata[log2FoldChange > 1 & padj < 0.05, type := "up"][log2FoldChange < -1 & padj < 0.05, type := "down"][]

labelGene <- degdata[order(pvalue, decreasing = T)][type == "up"][1:10]

pal_nejm()(8) %>% show_col()
typeColor <- structure(
  c(pal_nejm()(2), "gray80"),
  names = c("up", "down", "ns")
)

ggplot(degdata, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = type, size = pvalue), show.legend = F) +
  geom_text_repel(
    data = labelGene, aes(label = gene),
    size = 3, fontface = 3,
    nudge_x = .5, nudge_y = .5) +
  scale_radius(range = c(.1, 2)) +
  scale_color_manual(values = typeColor) +
  labs(
    x = TeX("$log_{2}(base\\,Mean)$"),
    y = TeX("$log_{2}(Fold\\,Change)$")) +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line())
