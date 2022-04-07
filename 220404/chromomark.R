#!/usr/bin/env Rscript
require(RIdeogram)
data(human_karyotype, package="RIdeogram")
human_karyotype <- human_karyotype[1:22,]
NJU9221_raw <-
  read.csv("fat/NJU9221_gap_input.tsv",
           sep = "\t",
           header = F)
NJU9226_raw <-
  read.csv("fat/NJU9226_gap_input.tsv",
           sep = "\t",
           header = F)
NJU9233_raw <-
  read.csv("fat/NJU9233_gap_input.tsv",
           sep = "\t",
           header = F)
gap <-
  read.csv("fat/gap_input.tsv",
           sep = "\t",
           header = F)
colnames(gap) <- c("Chr", "Start", "End")
gap$Value = 1

NJU9221_mark = data.frame(
  Type = "NJU9221",
  Shape = "triangle",
  Chr = NJU9221_raw$V1,
  Start = NJU9221_raw$V2,
  End = NJU9221_raw$V3,
  color = "3B4992"
)

NJU9226_mark = data.frame(
  Type = "NJU9226",
  Shape = "triangle",
  Chr = NJU9226_raw$V1,
  Start = NJU9226_raw$V2,
  End = NJU9226_raw$V3,
  color = "EE0000"
)

NJU9233_mark = data.frame(
  Type = "NJU9233",
  Shape = "circle",
  Chr = NJU9233_raw$V1,
  Start = NJU9233_raw$V2,
  End = NJU9233_raw$V3,
  color = "008B45"
)
mark <- rbind(NJU9221_mark, NJU9226_mark, NJU9233_mark)

ideogram(
  karyotype = human_karyotype,
  overlaid = gap,
  label = mark,
  label_type = "marker",
  colorset1 = c("white", "red", "blue")
)
convertSVG("chromosome.svg", device = "pdf")