#!/usr/bin/env Rscript
chr.size = data.frame(
  chr1 = 0,
  chr2 = 249250621,
  chr3 = 492449994,
  chr4 = 690472424,
  chr5 = 881626700,
  chr6 = 1062541960,
  chr7 = 1233657027,
  chr8 = 1392795690,
  chr9 = 1539159712,
  chr10 = 1680373143,
  chr11 = 1815907890,
  chr12 = 1950914406,
  chr13 = 2084766301,
  chr14 = 2199936179,
  chr15 = 2307285719,
  chr16 = 2409817111,
  chr17 = 2500171864,
  chr18 = 2581367074,
  chr19 = 2659444322,
  chr20 = 2718573305,
  chr21 = 2781598825,
  chr22 = 2829728720,
  chr23 = 2881033286,
  stringsAsFactors = F
)

chr.size[2, ] <- c(paste0("chr", seq(1, 22, 1)), "end")
chr.size[1, ] <- as.numeric(chr.size[1, ])

chr.label <- data.frame()
for (i in seq(1, 22, 1)) {
  chr.label[1, paste0("chr", i)] <-
    (as.numeric(chr.size[1, paste0("chr", i)]) + as.numeric(chr.size[1, paste0("chr", i +
                                                                                 1)])) / 2
}
chr.label[2, ] <- paste0("chr", seq(1, 22, 1))

dfraw <-
  read.csv(
    "BvPC_beta/uniloci.result.info.tsv",
    sep = "\t",
    header = F,
    row.names = 1
  )
dfraw <- na.omit(dfraw)
dfraw$X = NA
dfraw$Color = NA
dfraw$V10 <- as.numeric(dfraw$V10)
for (i in paste0("chr", seq(1, 22, 1))) {
  dfraw[dfraw$V9 == i, ]$X = dfraw[dfraw$V9 == i, ]$V10 + as.numeric(chr.size[1, i])
  dfraw[dfraw$V9 == i, ]$Color = chr.size[2, i]
}
dfraw$V4 <- as.numeric(dfraw$V4)
dfraw$X <- as.numeric(dfraw$X)
dfraw$Color <- as.factor(dfraw$Color)
dfraw <- na.omit(dfraw)
dfraw$label = NA
dfraw[dfraw$V4 %in% tail(sort(dfraw[dfraw$V2 > 0, ]$V4), 5), ]$label <-
  paste(rownames(dfraw[dfraw$V4 %in% tail(sort(dfraw[dfraw$V2 > 0, ]$V4), 5), ]), dfraw[dfraw$V4 %in% tail(sort(dfraw[dfraw$V2 >
                                                                                                                        0, ]$V4), 5), ]$V11, sep = "\n")
dfraw[dfraw$V4 %in% tail(sort(dfraw[dfraw$V2 < 0, ]$V4), 5), ]$label <-
  paste(rownames(dfraw[dfraw$V4 %in% tail(sort(dfraw[dfraw$V2 < 0, ]$V4), 5), ]), dfraw[dfraw$V4 %in% tail(sort(dfraw[dfraw$V2 <
                                                                                                                        0, ]$V4), 5), ]$V11, sep = "\n")

ggplot() +
  geom_point(data = dfraw[dfraw$V2 > 0, ],
             aes(x = X, y = V4, color = Color),
             show.legend = F) +
  geom_point(data = dfraw[dfraw$V2 < 0, ],
             aes(x = X, y = -V4, color = Color),
             show.legend = F) +
  geom_hline(yintercept = c(-3, 3),
             linetype = 2,
             alpha = 0.5) +
  geom_hline(yintercept = 0,
             linetype = 1,
             alpha = 0.3) +
  geom_vline(
    xintercept = as.numeric(chr.size[1, ]),
    linetype = 1,
    alpha = 0.1
  ) +
  geom_label_repel(
    data = dfraw[!is.na(dfraw$label) & dfraw$V2 > 0, ],
    aes(x = X, y = V4, label = label),
    box.padding   = 1,
    max.overlaps = Inf,
    point.padding = 0.3,
    force         = 100,
    segment.color = '#000000'
  ) +
  geom_label_repel(
    data = dfraw[!is.na(dfraw$label) & dfraw$V2 < 0, ],
    aes(x = X, y = -V4, label = label),
    box.padding   = 1,
    max.overlaps = Inf,
    point.padding = 0.3,
    force         = 100,
    segment.color = '#000000'
  ) +
  ylim(-7, 7) +
  ylab(label = expression(-log[10](P))) +
  xlab(label = "Chromosome") +
  scale_x_continuous(
    breaks = as.numeric(chr.label[1, ]),
    label = chr.label[2, ],
    expand = c(0.01, 0.01)
  ) +
  scale_colour_manual(
    values = c(
      "#0073C2FF",
      "#EFC000FF",
      "#868686FF",
      "#CD534CFF",
      "#7AA6DCFF",
      "#003C67FF",
      "#8F7700FF",
      "#3B3B3BFF",
      "#A73030FF",
      "#4A6990FF",
      "#0073C2FF",
      "#EFC000FF",
      "#868686FF",
      "#CD534CFF",
      "#7AA6DCFF",
      "#003C67FF",
      "#8F7700FF",
      "#3B3B3BFF",
      "#A73030FF",
      "#4A6990FF",
      "#0073C2FF",
      "#EFC000FF"
    )
  ) +
  theme(
    plot.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "#000000"),
    axis.ticks = element_line(colour = "#000000"),
    axis.ticks.length.y = unit(2, units = "mm"),
    axis.ticks.length.x = unit(1, units = "mm"),
    axis.title = element_text(colour = "#000000"),
    axis.text = element_text(colour = "#000000")
  )