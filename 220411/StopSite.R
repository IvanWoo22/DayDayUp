library(ggplot2)
library(ggseqlogo)
library(gridExtra)

setwd("~/fat/")

# control <- read.csv("NJU6429/Control.tsv", sep = "\t", header = F)
# treatment <-
#   read.csv("NJU6428/Treatment.tsv", sep = "\t", header = F)

control_length <-
  read.csv("NJU6429/Control_length.tsv",
           sep = "\t",
           header = F)
treatment_length <-
  read.csv("NJU6428/Treatment_length.tsv",
           sep = "\t",
           header = F)

control_fa <-
  read.csv("NJU6429/Control.fa",
           sep = "\t",
           header = F)
treatment_fa <-
  read.csv("NJU6428/Treatment.fa",
           sep = "\t",
           header = F)

# ggplot() +
#   geom_col(
#     aes(
#       x = control$V1 - 1,
#       y = control$V2,
#       group = "Control"
#     ),
#     fill = "blue",
#     alpha = 0.5
#   ) +
#   geom_col(
#     aes(
#       x = treatment$V1 - 1,
#       y = treatment$V2,
#       group = "Treatment"
#     ),
#     fill = "red",
#     alpha = 0.5
#   ) +
#   xlab(label = "Stop site") +
#   ylab(label = "Reads count") +
#   scale_y_continuous(limits = c(0, NA),
#                      expand = c(0, 0)) +
#   scale_x_continuous(
#     limits = c(9, 61),
#     expand = c(0, 0),
#     breaks = seq(10, 60, 2),
#     labels = seq(10, 60, 2)
#   ) +
#   theme(
#     legend.background = element_blank(),
#     legend.spacing = unit(0.1, units = "mm"),
#     legend.key.size = unit(3.2, 0.2, units = "mm"),
#     plot.title = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     plot.background = element_blank(),
#     panel.background = element_blank(),
#     axis.line = element_line(colour = "#000000"),
#     axis.ticks = element_line(colour = "#000000"),
#     axis.ticks.length.y = unit(2, units = "mm"),
#     axis.ticks.length.x = unit(1, units = "mm"),
#     axis.title = element_text(
#       colour = "#000000",
#       face = "bold",
#       size = 12
#     ),
#     axis.text = element_text(
#       colour = "#000000",
#       face = "bold",
#       size = 10
#     )
#   )

control_length$Group = "Control"
treatment_length$Group = "Treatment"
input <- rbind(control_length, treatment_length)
input$Group <- as.factor(input$Group)

p0 <- ggplot(input) +
  geom_density(aes(x = V1 - 1,
                   fill = Group,
                   color = Group), alpha = 0.1) +
  xlab(label = "Stop site") +
  ylab(label = "Density") +
  scale_y_continuous(limits = c(0, NA),
                     expand = c(0, 0)) +
  scale_x_continuous(
    limits = c(9, 57),
    expand = c(0, 0),
    breaks = seq(10, 56, 2),
    labels = seq(10, 56, 2)
  ) +
  theme(
    legend.position = c(0.8, 0.8),
    legend.background = element_blank(),
    legend.spacing = unit(0.1, units = "mm"),
    legend.key.size = unit(3.2, 0.2, units = "mm"),
    plot.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "#000000"),
    axis.ticks = element_line(colour = "#000000"),
    axis.ticks.length.y = unit(2, units = "mm"),
    axis.ticks.length.x = unit(1, units = "mm"),
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

set.seed(0411)

control_sample <-
  control_fa[sample(1:nrow(control_fa), size = 1000000),]
treatment_sample <-
  treatment_fa[sample(1:nrow(treatment_fa), size = 1000000),]

p1 <-
  ggseqlogo(control_sample, method = "prob") +
  geom_vline(xintercept = 16.5) +
  geom_vline(xintercept = 17.5) +
  theme(
    axis.title = element_text(
      colour = "#000000",
      face = "bold",
      size = 12
    ),
    axis.text = element_text(
      colour = "#000000",
      face = "bold",
      size = 8
    )
  )

p2 <-
  ggseqlogo(control_sample, method = "bits") +
  geom_vline(xintercept = 16.5) +
  geom_vline(xintercept = 17.5) +
  theme(
    axis.title = element_text(
      colour = "#000000",
      face = "bold",
      size = 12
    ),
    axis.text = element_text(
      colour = "#000000",
      face = "bold",
      size = 8
    )
  )

p3 <-
  ggseqlogo(treatment_sample, method = "prob") +
  geom_vline(xintercept = 16.5) +
  geom_vline(xintercept = 17.5) +
  theme(
    axis.title = element_text(
      colour = "#000000",
      face = "bold",
      size = 12
    ),
    axis.text = element_text(
      colour = "#000000",
      face = "bold",
      size = 8
    )
  )

p4 <-
  ggseqlogo(treatment_sample, method = "bits") +
  geom_vline(xintercept = 16.5) +
  geom_vline(xintercept = 17.5) +
  theme(
    axis.title = element_text(
      colour = "#000000",
      face = "bold",
      size = 12
    ),
    axis.text = element_text(
      colour = "#000000",
      face = "bold",
      size = 8
    )
  )

lay <- rbind(c(1, 1), c(1, 1), c(1, 1), c(2, 3), c(4, 5))
grid.arrange(p0, p1, p2, p3, p4, layout_matrix = lay)

front_motif <-
  read.csv(
    "front_motif.tsv",
    sep = "\t",
    header = F,
    row.names = 1
  )
ggplot() +
  geom_density(
    aes(x = front_motif$V2),
    alpha = 0.1,
    fill = "#5050FF",
    color = "#5050FF"
  ) +
  geom_density(
    aes(x = front_motif$V3),
    alpha = 0.1,
    fill = "#CE3D32",
    color = "#CE3D32"
  ) +
  scale_x_continuous(
    limits = c(0, 51),
    expand = c(0, 0),
    breaks = seq(1, 50, 1)
  )

front_motif_loss <-
  read.csv("front_motif_loss.tsv",
           sep = "\t",
           header = F)

ggseqlogo(front_motif_loss, method = "bits") +
  theme(
    axis.title = element_text(
      colour = "#000000",
      face = "bold",
      size = 12
    ),
    axis.text = element_text(
      colour = "#000000",
      face = "bold",
      size = 8
    )
  )
ggseqlogo(front_motif_loss, method = "prob") +
  theme(
    axis.title = element_text(
      colour = "#000000",
      face = "bold",
      size = 12
    ),
    axis.text = element_text(
      colour = "#000000",
      face = "bold",
      size = 8
    )
  )
