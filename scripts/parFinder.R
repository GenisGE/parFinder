
col12 <- c("#FFB7FC","#ef8a62","#fddbc7","#d1e5f0","#FF106A","#FFB7FC")[c(1,5)]
names(col12) <- c("homomorphic", "heteromorphic")

args <- commandArgs(trailingOnly=TRUE)

infile <- args[1]
sexassign <- args[2]
outprefix <- args[3]

# whether use median or mean of window to classify them
useMedian <- FALSE

# read in input files
dat <- read.table(infile, h=T, sep=" ")
samplesex <- read.table(sexassign,h=T)
sex <- samplesex$Inferred_sex
names(sex) <- gsub(".idxstats", "", samplesex$Sample)

# subset columns with depth information
depmat <- dat[,-(1:4)]

# first scaffold should be a well-behaved autosomal scaffold, for normalization
normwins <- dat$scaffold == dat$scaffold[1]
# remaining scaffolds are sex-linked scaffolds
sexwins <- dat$scaffold != dat$scaffold[1]

# normalize by either median or mean of autosomal windows
if(useMedian){
    denom <- apply(depmat[normwins,]/dat$len[normwins], 2, median)
} else {
    denom <- apply(depmat[normwins,]/dat$len[normwins], 2, mean)
}

# normalize scaffolds, so autosomal-like have value ~1
normMat <- t(t(depmat/dat$len) / denom)

# do test to compare mean betweent the two. it does seem to work!
## FOR NOW NOT USED
#pval <- apply(normMat,1,function(x) t.test(x~sex)$p.value)

# table to classify windows
outDF <- dat[sexwins,1:4]

# add mean depth and standard error of each window per sex
outDF$meanDepHomomorphic <- rowMeans(normMat[sexwins,sex=="homomorphic"])
outDF$seDepHomomorphic <- apply(normMat[sexwins,sex=="homomorphic"], 1, sd) / sqrt(sum(sex=="homomorphic"))
outDF$meanDepHeteromorphic <- rowMeans(normMat[sexwins,sex=="heteromorphic"])
outDF$seDepHeteromorphic <- apply(normMat[sexwins,sex=="heteromorphic"], 1, sd) / sqrt(sum(sex=="heteromorphic"))

outDF$diffDepthHomMinusHet <- outDF$meanDepHomomorphic - outDF$meanDepHeteromorphic

# NEEDS IMPROVING
# make thresholds as parameter or use t-test
outDF$winCategory <- "weird"
outDF$winCategory[outDF$diff < 0.15 & outDF$diff > -0.15] <- "autosomal" # needs extra condition, that both sexes are ~1
outDF$winCategory[outDF$diff < 0.65 & outDF$diff > 0.35] <- "sex-linked" # needs extra condition, heterogametic ~0.5 and homogametic ~1

# write output table with
outtsv <- paste0(outprefix, "_winClassification.tsv")
cat("Finished classification, will write results table to ", outtsv, " and plots to ", paste0(outprefix, "_plots/"), "\n")
write.table(x=outDF, file=outtsv, col.names=T, row.names=F, quote=F, sep="\t")


# make plot for each sex scaffold with window valeus and window classification
sexscaffs <- unique(dat$scaffold[sexwins])
for(scaff in sexscaffs){
outpng <- paste0(outprefix, "_plots/", scaff, ".png")
bitmap(outpng, w=8,h=4, res=300)
    par(oma=c(0,0,4,0))
    k <- outDF$scaffold == scaff
    midpoints <- apply(outDF[k,c("end", "start")], 1, mean)

    ylim <- c(min(min(outDF$meanDepHomomorphic[k] - outDF$seDepHomomorphic[k],
                     outDF$meanDepHeteromorphic[k] - outDF$seDepHeteromorphic[k], 
                     outDF$diffDepthHomMinusHet[k]), 0),
               max(max(outDF$meanDepHomomorphic[k] + outDF$seDepHomomorphic[k],
                     outDF$meanDepHeteromorphic[k] + outDF$seDepHeteromorphic[k], 
                     outDF$diffDepthHomMinusHet[k]),1))
    plot(x=midpoints, y=outDF$diffDepthHomMinusHet[k], type="l", lwd=2, main="", ylim=ylim,
            xlab=paste("Position in scaffold", scaff), ylab="Normalized depth", cex.lab=1.5)

    if(any(outDF$winCategory[k] == "autosomal")){
        rect(ybottom=ylim[1], ytop=ylim[2],
                    xleft=outDF$start[k][outDF$winCategory[k] == "autosomal"],
                    xright=outDF$end[k][outDF$winCategory[k] == "autosomal"], 
                    col="azure3", border=NA)
    }
    if(any(outDF$winCategory[k] == "weird")){
        rect(ybottom=ylim[1], ytop=ylim[2],
            xleft=outDF$start[k][outDF$winCategory[k] == "weird"],
            xright=outDF$end[k][outDF$winCategory[k] == "weird"], 
            col="antiquewhite", border=NA)
    }
    lines(x=midpoints, y=outDF$diffDepthHomMinusHet[k], lwd=2)

    lines(x=midpoints, y=outDF$meanDepHomomorphic[k], lwd=2, col = col12["homomorphic"])
    lines(x=midpoints, y=outDF$meanDepHomomorphic[k] + outDF$seDepHomomorphic[k], lwd=1, lty=2, col = col12["homomorphic"])
    lines(x=midpoints, y=outDF$meanDepHomomorphic[k] - outDF$seDepHomomorphic[k], lwd=1, lty=2, col = col12["homomorphic"])

    lines(x=midpoints, y=outDF$meanDepHeteromorphic[k], lwd=2, col = col12["heteromorphic"])
    lines(x=midpoints, y=outDF$meanDepHeteromorphic[k] + outDF$seDepHeteromorphic[k], lwd=1, lty=2, col = col12["heteromorphic"])
    lines(x=midpoints, y=outDF$meanDepHeteromorphic[k] - outDF$seDepHeteromorphic[k], lwd=1, lty=2, col = col12["heteromorphic"])


    legend(x=max(midpoints) * 0.25, y=ylim[2] * 1.3, bty="n",
            lty=1,lwd=2, col=c(col12, "black"), legend=c("Homomorphic mean", "Heteromorphic mean", "Difference Hom.-Het."), xpd=NA, cex=1.5)
    legend(x=max(midpoints) * 0.65, y=ylim[2] * 1.3, bty="n",
            fill=c("azure3", "antiquewhite"), legend=c("Autosomal-like (PAR)", "Weird"), xpd=NA, cex=1.5,
            title="Window classification")

dev.off()
}
