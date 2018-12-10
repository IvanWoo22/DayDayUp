# RNAseq with R
## Prepare
### Install R packages. 
Open R and copy the following commands into R console:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("RNAseq123", version = "3.8")
install.packages(“R.utils”)
```
### Download data. 
```
url <- “https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file” 
utils::download.file(url, destfile="GSE63310_RAW.tar", mode="wb") 
utils::untar("GSE63310_RAW.tar", exdir = ".")
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt","GSM1545538_purep53.txt", 
  "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt", 
  "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt","GSM1545545_JMS9-P8c.txt")
for(i in paste(files, ".gz", sep="")) R.utils::gunzip(i, overwrite=TRUE)
```
Save file downloaded at: http://bioinf.wehi.edu.au/software/MSigDB/mouse_c2_v5p1.rdata
### Check on R console.
```
library(limma)
library(edgeR)
library(Mus.musculus)
library(RColorBrewer)
library(Glimma)
library(gplots)
read.delim(files[1], nrow=5)
```
## Then try like 
