rm(list = ls())

rawdat1 <- t(read.table("10Xr_18s_score_output.tsv"))
mdat1 <- matrix(as.numeric(rawdat1[2:3,]),2)
colnames(mdat1) <- rawdat1[1,]


rawdat2 <- t(read.table("6Xr_18s_score_output.tsv"))
mdat2 <- matrix(as.numeric(rawdat2[2:3,]),2)
colnames(mdat2) <- rawdat2[1,]


par(mfrow=c(2,1),mar=c(2,0,1,0))
barplot(height = mdat1[2,],
        width = 1,
        yaxt = "n",
        xaxt = "n",
        space = 0,
        col = c(1))
axis(1,seq(0.5,1940.5,1),rawdat1[1,],mgp=c(2,0,0),cex.axis=0.17)
axis(1,seq(0.5,1940.5,1),c(1:1941),mgp=c(2,0.25,0),cex.axis=0.17)


barplot(height = mdat2[2,],
        width = 1,
        yaxt = "n",
        xaxt = "n",
        space = 0,
        col = c(1))
axis(1,seq(0.5,1940.5,1),rawdat2[1,],mgp=c(2,0,0),cex.axis=0.17)
axis(1,seq(0.5,1940.5,1),c(1:1941),mgp=c(2,0.25,0),cex.axis=0.17)
