library(DESeq2)
library(ggpubr)
library(ggplot2)

dataframe <- read.csv("merge.csv", header = TRUE, row.names = 1)
countdata <- dataframe[-(1:5), ]
row_names <- row.names(countdata)
name_replace <- gsub("\\.\\w+", "", row.names(countdata))
row.names(countdata) <- name_replace
countdata <- countdata[rowSums(countdata) > 0, ]

coldata <-
  read.table(
    "phenotype.csv",
    row.names = 1,
    header = TRUE,
    sep = ","
  )

countdata <- countdata[row.names(coldata)]
coldata$group <- as.factor(coldata$group)

dds <-
  DESeqDataSetFromMatrix(countData = countdata,
                         colData = coldata,
                         design = ~ group)

vsd <- assay(vst(dds, blind = FALSE))
dds2 <- DESeq(dds, parallel = T)
suppressMessages(dds2)
res <-
  results(
    dds2,
    contrast = c('group', 'treatment', 'control'),
    pAdjustMethod = 'fdr',
    alpha = 0.05
  )
plotMA(res, alpha = 0.05, ylim = c(-3, 3))

res[which(res$padj %in% NA), 'sig'] <- 'no diff'
res[which(res$log2FoldChange >= 1 &
            res$padj < 0.05), 'sig'] <-
  'up (p.adj < 0.05, log2FC >= 1)'
res[which(res$log2FoldChange <= -1 &
            res$padj < 0.05), 'sig'] <-
  'down (p.adj < 0.05, log2FC <= -1)'
res[which(abs(res$log2FoldChange) < 1 |
            res$padj >= 0.05), 'sig'] <- 'no diff'

res[which(res$padj %in% NA), 'sigg'] <- 'no diff'
res[which(res$log2FoldChange >= 1 &
            res$padj < 0.05), 'sigg'] <- 'up'
res[which(res$log2FoldChange <= -1 &
            res$padj < 0.05), 'sigg'] <- 'down'
res[which(abs(res$log2FoldChange) < 1 |
            res$padj >= 0.05), 'sigg'] <- 'no diff'

ggres <- as(res, "data.frame")
ggres$label <- ""

shared <- intersect(mhv, hc)
mhv_sp <- setdiff(mhv, hc)
hc_sp <- setdiff(hc, mhv)

name1 <- data.frame(id = mhv, name = mhv_name)
name2 <- data.frame(id = hc, name = hc_name)
name_all <- rbind(name1, name2)
rownames(name_all) <- name_all[, 1]

shared_list <- name1[which(mhv %in% shared),]
mhv_list <- name1[which(mhv %in% mhv_sp),]
hc_list <- name2[which(hc %in% hc_sp),]


list <- row.names(ggres[which(ggres$sigg != 'no diff' &
                                row.names(ggres) %in% mhv_sp), ])

ggres[which(ggres$sigg != 'no diff' &
              row.names(ggres) %in% mhv_sp), 'label'] <-
  mhv_list[match(list, mhv_list[, 1]), 'name']
ggres[which(ggres$sigg != 'no diff' &
              row.names(ggres) %in% mhv_sp), 'tom'] <- 'appear'

list <- row.names(ggres[which(ggres$sigg != 'no diff' &
                                row.names(ggres) %in% hc_sp), ])

ggres[which(ggres$sigg != 'no diff' &
              row.names(ggres) %in% hc_sp), 'label'] <-
  hc_list[match(list, hc_list[, 1]), 'name']
ggres[which(ggres$sigg != 'no diff' &
              row.names(ggres) %in% hc_sp), 'tom'] <- 'disappear'



ggplot() +
  geom_point(
    data = ggres[which(is.na(ggres$tom)), ],
    mapping = aes(log2FoldChange,-log(padj, 10), color = "total"),
    alpha = 0.3,
    size = 2,
    show.legend = F,
    shape = 20
  ) +
  geom_point(
    data = ggres[which(!is.na(ggres$tom)), ],
    mapping = aes(log2FoldChange,-log(padj, 10), color = tom),
    alpha = 1,
    size = 3,
    shape = 20
  ) +
  geom_text(
    data = ggres[which(!is.na(ggres$tom)), ],
    mapping = aes(
      log2FoldChange,
      -log(padj, 10),
      label = label,
      color = tom
    ),
    hjust = 0,
    vjust = 0
  ) +
  scale_color_manual(values = c('#ff7500', '#008aff', 'gray30', 'gray30', 'red2')) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(color = 'black', fill = 'transparent'),
    legend.position = c(0.26, 0.92),
    legend.title = element_blank(),
    legend.key = element_rect(fill = 'transparent'),
    legend.background = element_rect(fill = 'transparent')
  ) +
  geom_vline(xintercept = c(-1, 1),
             color = 'black',
             size = 0.25) +
  geom_hline(yintercept = -log(0.05, 10),
             color = 'black',
             size = 0.25) +
  labs(x = 'log2 Fold Change', y = '-log10 p-value', color = NA) +
  xlim(-5, 5)
