ozone <-
  data.frame(
    group1 = c(
      rep(1, 6),
      rep(2, 5),
      rep(3, 4),
      rep(4, 3),
      rep(5, 2),
      6,
      2:7,
      3:7,
      4:7,
      5:7,
      6:7,
      7
    ),
    group2 = c(
      2:7,
      3:7,
      4:7,
      5:7,
      6:7,
      7,
      rep(1, 6),
      rep(2, 5),
      rep(3, 4),
      rep(4, 3),
      rep(5, 2),
      6
    ),
    P = as.numeric(rawdat[1, 5:46])
  )
ozone.dists <-
  as.matrix(dist(cbind(ozone$group1, ozone$group2)))ozone.dists.inv <-
  1 / ozone.distsdiag(ozone.dists.inv) <- 0
Moran.I(ozone$P, ozone.dists.inv)

euldis_dens <-
  read.csv("Moran.tsv",
           sep = "\t",
           header = F)
euldis_dens_part <-
  euldis_dens[euldis_dens$V50 < 0.05,]tmpp <-
  do.call(c, lapply(5:46, function(x)
    sum(euldis_dens_part[, x])))tmp <-
  matrix(nrow = 7, ncol = 7)rownames(tmp) <-
  c("TC", "TE", "P5", "P10", "P15", "P20", "PN")colnames(tmp) <-
  c("TC", "TE", "P5", "P10", "P15", "P20", "PN")tmp[1, 2:7] <-
  tmpp[1:6]tmp[2, 3:7] <-
  tmpp[7:11]tmp[3, 4:7] <-
  tmpp[12:15]tmp[4, 5:7] <-
  tmpp[16:18]tmp[5, 6:7] <-
  tmpp[19:20]tmp[6, 7] <-
  tmpp[21]tmp[2:7, 1] <-
  tmpp[22:27]tmp[3:7, 2] <-
  tmpp[28:32]tmp[4:7, 3] <-
  tmpp[33:36]tmp[5:7, 4] <-
  tmpp[37:39]tmp[6:7, 5] <- tmpp[40:41]tmp[7, 6] <- tmpp[42]
temp <- reshape2:::melt.matrix(tmp, na.rm = TRUE)
p1 <-
  ggplot(temp, aes(Var2, reorder(Var1, desc(Var1)), fill = value)) +
  geom_tile(color = "#ffffff", aes(width = 1, height = 1), size = 1) +
  geom_text(
    aes(label = value),
    color = "#111111",
    size = 3,
    fontface = "bold"
  ) +
  coord_fixed() +
  scale_fill_gradientn(colors = pal_material("red")(10)[c(1:3, rep(4:10, each = 6))], name = "Count") +
  scale_x_discrete(position = "top") +
  labs(title = "Moran's I P-value < 0.05\n(n=79,153)") +
  theme(
    plot.title = element_text(
      colour = "#000000",
      face = "bold",
      size = 14,
      hjust = 0.5
    ),
    axis.text = element_text(
      colour = "#000000",
      face = "bold",
      size = 10
    ),
    legend.title = element_text(
      colour = "#000000",
      face = "bold",
      size = 6
    ),
    legend.text = element_text(
      colour = "#000000",
      face = "bold",
      size = 6
    ),
    axis.text.x.top = element_text(vjust = -1.2),
    axis.text.y = element_text(hjust = 1.2),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    plot.background = element_blank()
  )
euldis_dens_part <-
  euldis_dens[euldis_dens$V50 >= 0.05,]tmpp <-
  do.call(c, lapply(5:46, function(x)
    sum(euldis_dens_part[, x])))tmp <-
  matrix(nrow = 7, ncol = 7)rownames(tmp) <-
  c("TC", "TE", "P5", "P10", "P15", "P20", "PN")colnames(tmp) <-
  c("TC", "TE", "P5", "P10", "P15", "P20", "PN")tmp[1, 2:7] <-
  tmpp[1:6]tmp[2, 3:7] <-
  tmpp[7:11]tmp[3, 4:7] <-
  tmpp[12:15]tmp[4, 5:7] <-
  tmpp[16:18]tmp[5, 6:7] <-
  tmpp[19:20]tmp[6, 7] <-
  tmpp[21]tmp[2:7, 1] <-
  tmpp[22:27]tmp[3:7, 2] <-
  tmpp[28:32]tmp[4:7, 3] <-
  tmpp[33:36]tmp[5:7, 4] <-
  tmpp[37:39]tmp[6:7, 5] <- tmpp[40:41]tmp[7, 6] <- tmpp[42]
temp <- reshape2:::melt.matrix(tmp, na.rm = TRUE)
p2 <-
  ggplot(temp, aes(Var2, reorder(Var1, desc(Var1)), fill = value)) + geom_tile(color = "#ffffff",      aes(width = 1, height = 1),      size = 1) + geom_text(
    aes(label = value),
    color = "#111111",
    size = 3,
    fontface = "bold"
  ) + coord_fixed() + scale_fill_gradientn(colors = pal_material("red")(10)[c(1:3, rep(4:10, each = 6))], name = "Count") + scale_x_discrete(position = "top") + labs(title = "Moran's I P-value >= 0.05\n(n=33,417)") + theme(
    plot.title = element_text(
      colour = "#000000",
      face = "bold",
      size = 14,
      hjust = 0.5
    ),
    axis.text = element_text(
      colour = "#000000",
      face = "bold",
      size = 10
    ),
    legend.title = element_text(
      colour = "#000000",
      face = "bold",
      size = 6
    ),
    legend.text = element_text(
      colour = "#000000",
      face = "bold",
      size = 6
    ),
    axis.text.x.top = element_text(vjust = -1.2),
    axis.text.y = element_text(hjust = 1.2),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    plot.background = element_blank()
  )
euldis_dens_part <-
  euldis_dens[euldis_dens$V50 < 0.001,]
  tmpp <-
  do.call(c, lapply(5:46, function(x)
    sum(euldis_dens_part[, x])))tmp <-
  matrix(nrow = 7, ncol = 7)
  rownames(tmp) <-
  c("TC", "TE", "P5", "P10", "P15", "P20", "PN")
  colnames(tmp) <-
  c("TC", "TE", "P5", "P10", "P15", "P20", "PN")
  tmp[1, 2:7] <- tmpp[1:6]
  tmp[2, 3:7] <- tmpp[7:11]
  tmp[3, 4:7] <- tmpp[12:15]
  tmp[4, 5:7] <- tmpp[16:18]
  tmp[5, 6:7] <- tmpp[19:20]
  tmp[6, 7] <- tmpp[21]
  tmp[2:7, 1] <- tmpp[22:27]
  tmp[3:7, 2] <- tmpp[28:32]
  tmp[4:7, 3] <- tmpp[33:36]
  tmp[5:7, 4] <- tmpp[37:39]
  tmp[6:7, 5] <- tmpp[40:41]
  tmp[7, 6] <- tmpp[42]
temp <- reshape2:::melt.matrix(tmp, na.rm = TRUE)
p3 <-
  ggplot(temp, aes(Var2, reorder(Var1, desc(Var1)), fill = value)) + geom_tile(color = "#ffffff",      aes(width = 1, height = 1),      size = 1) + geom_text(
    aes(label = value),
    color = "#111111",
    size = 3,
    fontface = "bold"
  ) +
  coord_fixed() +
  scale_fill_gradientn(colors = pal_material("red")(10)[c(1:3, rep(4:10, each = 6))], name = "Count") +
  scale_x_discrete(position = "top") +
  labs(title = "Moran's I P-value < 0.001\n(n=61,600)") +
  theme(
    plot.title = element_text(
      colour = "#000000",
      face = "bold",
      size = 14,
      hjust = 0.5
    ),
    axis.text = element_text(
      colour = "#000000",
      face = "bold",
      size = 10
    ),
    legend.title = element_text(
      colour = "#000000",
      face = "bold",
      size = 6
    ),
    legend.text = element_text(
      colour = "#000000",
      face = "bold",
      size = 6
    ),
    axis.text.x.top = element_text(vjust = -1.2),
    axis.text.y = element_text(hjust = 1.2),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    plot.background = element_blank()
  )
euldis_dens_part <-
  euldis_dens[euldis_dens$V50 >= 0.001,]
  tmpp <-
  do.call(c, lapply(5:46, function(x)
    sum(euldis_dens_part[, x])))
    tmp <-
  matrix(nrow = 7, ncol = 7)
  rownames(tmp) <-
  c("TC", "TE", "P5", "P10", "P15", "P20", "PN")
  colnames(tmp) <-
  c("TC", "TE", "P5", "P10", "P15", "P20", "PN")
  tmp[1, 2:7] <-
  tmpp[1:6]
  tmp[2, 3:7] <-
  tmpp[7:11]
  tmp[3, 4:7] <-
  tmpp[12:15]
  tmp[4, 5:7] <-
  tmpp[16:18]
  tmp[5, 6:7] <-
  tmpp[19:20]
  tmp[6, 7] <-
  tmpp[21]
  tmp[2:7, 1] <-
  tmpp[22:27]
  tmp[3:7, 2] <-
  tmpp[28:32]
  tmp[4:7, 3] <-
  tmpp[33:36]
  tmp[5:7, 4] <-
  tmpp[37:39]
  tmp[6:7, 5] <- tmpp[40:41]
  tmp[7, 6] <- tmpp[42]
temp <- reshape2:::melt.matrix(tmp, na.rm = TRUE)
p4 <-
  ggplot(temp, aes(Var2, reorder(Var1, desc(Var1)), fill = value)) + geom_tile(color = "#ffffff",      aes(width = 1, height = 1),      size = 1) + geom_text(
    aes(label = value),
    color = "#111111",
    size = 3,
    fontface = "bold"
  ) +
  coord_fixed() + scale_fill_gradientn(colors = pal_material("red")(10)[c(1:3, rep(4:10, each = 6))], name = "Count") + scale_x_discrete(position = "top") + labs(title = "Moran's I P-value >= 0.001\n(n=50,970)") + theme(
    plot.title = element_text(
      colour = "#000000",
      face = "bold",
      size = 14,
      hjust = 0.5
    ),
    axis.text = element_text(
      colour = "#000000",
      face = "bold",
      size = 10
    ),
    legend.title = element_text(
      colour = "#000000",
      face = "bold",
      size = 6
    ),
    legend.text = element_text(
      colour = "#000000",
      face = "bold",
      size = 6
    ),
    axis.text.x.top = element_text(vjust = -1.2),
    axis.text.y = element_text(hjust = 1.2),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    plot.background = element_blank()
  )
grid.arrange(p1, p2, p3, p4, nrow = 2)
head(euldis_dens)
ggplot() + geom_density(stat = "density", aes(x = euldis_dens$V47)) + labs(title = "Moran's I observed distribution") + theme(
  plot.title = element_text(
    colour = "#000000",
    face = "bold",
    size = 14,
    hjust = 0.5
  ),
  axis.text = element_text(
    colour = "#000000",
    face = "bold",
    size = 10
  ),
  legend.title = element_text(
    colour = "#000000",
    face = "bold",
    size = 6
  ),
  legend.text = element_text(
    colour = "#000000",
    face = "bold",
    size = 6
  ),
  axis.text.x.top = element_text(vjust = -1.2),
  axis.text.y = element_text(hjust = 1.2),
  axis.title = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_line(colour = "#000000",    size = 0.4668623),
  axis.line = element_line(colour = "#000000",   size = 0.4668623),
  plot.background = element_blank()
)euldis_dens_part <-
  euldis_dens[euldis_dens$V47 >= 0.09528642,]tmpp <-
  do.call(c, lapply(5:46, function(x)
    sum(euldis_dens_part[, x])))tmp <-
  matrix(nrow = 7, ncol = 7)rownames(tmp) <-
  c("TC", "TE", "P5", "P10", "P15", "P20", "PN")colnames(tmp) <-
  c("TC", "TE", "P5", "P10", "P15", "P20", "PN")tmp[1, 2:7] <-
  tmpp[1:6]tmp[2, 3:7] <-
  tmpp[7:11]tmp[3, 4:7] <-
  tmpp[12:15]tmp[4, 5:7] <-
  tmpp[16:18]tmp[5, 6:7] <-
  tmpp[19:20]tmp[6, 7] <-
  tmpp[21]tmp[2:7, 1] <-
  tmpp[22:27]tmp[3:7, 2] <-
  tmpp[28:32]tmp[4:7, 3] <-
  tmpp[33:36]tmp[5:7, 4] <-
  tmpp[37:39]tmp[6:7, 5] <-
  tmpp[40:41]tmp[7, 6] <-
  tmpp[42]temp <- reshape2:::melt.matrix(tmp, na.rm = TRUE)

ggplot(temp, aes(Var2, reorder(Var1, desc(Var1)), fill = value)) +
  geom_tile(color = "#ffffff",      aes(width = 1, height = 1),      size = 1) +
  geom_text(
    aes(label = value),
    color = "#111111",
    size = 3,
    fontface = "bold"
  ) +
  coord_fixed() + scale_fill_gradientn(colors = pal_material("red")(10)[c(1:3, rep(4:10, each = 6))], name = "Count") +
  scale_x_discrete(position = "top") +
  labs(title = "Moran's I observed Component 1\n(n=32,305)") +
  theme(
    plot.title = element_text(
      colour = "#000000",
      face = "bold",
      size = 14,
      hjust = 0.5
    ),
    axis.text = element_text(
      colour = "#000000",
      face = "bold",
      size = 10
    ),
    legend.title = element_text(
      colour = "#000000",
      face = "bold",
      size = 6
    ),
    legend.text = element_text(
      colour = "#000000",
      face = "bold",
      size = 6
    ),
    axis.text.x.top = element_text(vjust = -1.2),
    axis.text.y = element_text(hjust = 1.2),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    plot.background = element_blank()
  )
y <- as.numeric(euldis_dens$V47)
BIC <-
  mclustBIC(y, modelNames = c("V"), G = 3)mod4 <-
  Mclust(y, x = BIC)summary(mod4, parameters = TRUE)table(mod4$classification)
min(y[which(mod4$classification == 1)])min(y[which(mod4$classification == 2)])min(y[which(mod4$classification == 3)])
dtf$Class <- 0
dtf[dtf$V4 >= 0,]$Class = 1
dtf[dtf$V4 >= 7,]$Class = 2
dtf[dtf$V4 >= 16,]$Class = 3
dtf[dtf$V4 >= 25,]$Class = 4
dtf[dtf$V4 >= 36,]$Class = 5
dtf$Class <-
  as.factor(dtf$Class)
library(mixtools)
mixmdl = normalmixEM(euldis_dens$V47, k = 3, mu = c(0, 0.1, 0.22))
plot(mixmdl, which = 2)
title(sub = "Moran's I observed distribution")
lines(density(euldis_dens$V47), lty = 1, lwd = 1)
write.table(
  euldis_dens_part,
  file = "MoranC1.tsv",
  quote = F,
  sep = "\t",
  col.names = F,
  row.names = F
)


rawdat <-
  read.csv(
    "RRBS_sites_merge.tsv",
    header = F,
    sep = "\t",
    check.names = F
  )
  patient_list <-
  sprintf("%02d", c(1, 5:15))
  c12 <- c(
    "dodgerblue2",
    "#FF7F00",
    "#8F7700",
    "gold1",
    "skyblue2",
    "#FB9A99",
    "palegreen2",
    "#CAB2D6",
    "#FDBF6F",
    "gray70",
    "khaki2",
    "maroon"
  )
  myPlots <-
  list()
  for (j in 1:25) {
    f = sample(c(1:nrow(euldis_dens_part)), 1)
    k = rownames(rawdat[rawdat$V1 == euldis_dens_part[f, 1] & rawdat$V2 == euldis_dens_part[f, 2],])
    sitedat <-tibble(beta = as.numeric(rawdat[k, 3:ncol(rawdat)]),
    patient = rep(patient_list, each = 7),
    loc_no = rep((1:7), 12),
    is_complete = 1
    )
    filter <- unique(sitedat[is.na(sitedat$beta), ]$patient)
    sitedat[which(sitedat$patient %in% filter),]$is_complete = 0
    res.spr <- cor.test(sitedat$loc_no,round(sitedat$beta, 2),conf.level = 0.95,method = "kendall",exact = FALSE)
    TP <- res.spr$p.value
    SIG = 0
    if ((!is.na(TP)) & (TP < 0.05)) {
    SIG = 1
    }
    TR <- as.numeric(res.spr$estimate)
    SSIG <- rep(0.3, 12)
    for (patientno in 1:12) {
    if (!patient_list[patientno] %in% filter) {
    singledat <- sitedat[sitedat$patient == patient_list[patientno], ]
    res.single.spr <- cor.test(singledat$loc_no,round(singledat$beta, 2),conf.level = 0.95,method = "kendall",exact = FALSE)
    if ((!is.na(res.single.spr$p.value)) & (res.single.spr$p.value < 0.05)) {
    SSIG[patientno] = 0.9
    }
    }
    }
    stat.test <- sitedat %>%  pairwise_t_test(beta ~ loc_no, paired = F,p.adjust.method = "none") %>% add_xy_position(x = "loc_no")
    myPlots[[j]] <- ggplot() +
    geom_point( aes_string( x = sitedat$loc_no, y = round(sitedat$beta, 2), alpha = as.factor(sitedat$is_complete)), size = 0.85, show.legend = F ) +
    scale_alpha_manual(values = c(0.3, 0.8)) +
    geom_line(aes_string( x = sitedat[sitedat$is_complete == 1, ]$loc_no, y = round(sitedat[sitedat$is_complete == 1, ]$beta, 2), size = as.factor(sitedat[sitedat$is_complete == 1, ]$patient), color = as.factor(sitedat[sitedat$is_complete == 1, ]$patient))) +
                                                                           scale_color_manual(values = c12) +  scale_x_continuous(
                                                                             breaks = c(1:7),
                                                                             limits = c(0.8, 7.2),
                                                                             expand = c(0, 0),
                                                                             labels = c("TC", "TE", "P5", "P10", "P15", "P20", "PN")
                                                                           ) +  stat_pvalue_manual(stat.test, hide.ns = T, label = "p.adj.signif") +  scale_size_manual(values = SSIG, guide = 'none') +  xlab(label = "Location") +  ylab(label = "Methylation ratio (%)") +
                                                                           labs(
                                                                             title = paste(str_to_title(rawdat[k, 1]), comma(rawdat[k, 2]), sep = " "),
                                                                             subtitle = paste(
                                                                               "Kendall, p = ",
                                                                               format(TP, digits = 4),
                                                                               ", tau = ",
                                                                               format(TR, digits = 2),
                                                                               sep = ""
                                                                             ),
                                                                             color = "Patient"
                                                                           ) +  theme(
                                                                             legend.background = element_blank(),
                                                                             legend.spacing = unit(0.1, units = "mm"),
                                                                             legend.key.size = unit(4, 0.5, units = "mm"),
                                                                             plot.title = element_text(
                                                                               colour = "#000000",
                                                                               face = "bold",
                                                                               size = 12,
                                                                               vjust = -1,
                                                                             ),
                                                                             panel.grid.major = element_blank(),
                                                                             panel.grid.minor = element_blank(),
                                                                             plot.background = element_blank(),
                                                                             panel.background = element_blank(),
                                                                             axis.line = element_line(colour = "#000000"),
                                                                             axis.ticks = element_line(colour = "#000000"),
                                                                             axis.ticks.length.y = unit(2, units = "mm"),
                                                                             axis.ticks.length.x = unit(1, units = "mm"),
                                                                             legend.box.background = element_blank(),
                                                                             legend.key = element_blank(),
                                                                             axis.title = element_text(
                                                                               colour = "#000000",
                                                                               face = "bold",
                                                                               size = 12
                                                                             ),
                                                                             axis.text = element_text(
                                                                               colour = "#000000",
                                                                               face = "bold",
                                                                               size = 10
                                                                             )
                                                                           )
  }
do.call(grid.arrange, c(myPlots, list(ncol = 5)))
