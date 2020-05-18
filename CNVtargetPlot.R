#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## script to make MLPA style plots based on WXS data ##
## uses GATK denoised cnratios and called segments as inputdata ##
## and a bed file that specifies the regions of interest ##
## requires BioConductor IRanges package ##


cnRatiosFile <- args[1]
cnRatiosFileNormal <- args[2]
calledSegFile <- args[3]
targetRegionFile <- args[4]

library(IRanges,quietly=TRUE)

cnRatios <- read.csv(cnRatiosFile,sep="\t",skip=1,header=F,stringsAsFactors = F)
colnames(cnRatios) <- c("Chr","Start","Stop","Feature","log2CNratio")
cnRatiosNormal <- read.csv(cnRatiosFileNormal,sep="\t",skip=1,header=F,stringsAsFactors = F)
colnames(cnRatiosNormal) <- c("Chr","Start","Stop","Feature","log2CNratio")

targetRegions <- read.csv(targetRegionFile,sep="\t",header=T,stringsAsFactors = F)
colnames(targetRegions) <- c("Chr","Start","Stop","Feature")


calledSeg <- read.csv(calledSegFile,sep="\t",header=T,stringsAsFactors = F)
diploidSegments <- calledSeg[calledSeg$Call == 0,]

tempList <- list()
for ( i in 1:nrow(diploidSegments)){
  tempList[[i]] <- cnRatios[cnRatios$Chr == diploidSegments$Chromosome[i] & cnRatios$Start >= diploidSegments$Start[i] & cnRatios$Stop <= diploidSegments$End[i],"log2CNratio"]
}
diploidSD <- sd(2**unlist(tempList))
diploidQuant99 <- quantile(x = 2**unlist(tempList),probs = c(0.005,0.995))
diploidQuant999 <- quantile(x = 2**unlist(tempList),probs = c(0.0005,0.9995))

targetRegions <- cbind(targetRegions,matrix(ncol=4,nrow=nrow(targetRegions),data = NA))
colnames(targetRegions)[(ncol(targetRegions)-3):ncol(targetRegions)] <- c("CNratioMean","CNratioSD","CNratioNormalMean","CNratioNormalSD")

for ( i in 1:nrow(targetRegions)){
  curChr <- cnRatios[cnRatios$Chr == targetRegions$Chr[i],]
  if (dim(curChr)[1] == 0){
    next
  }
  rangesA <- split(IRanges(curChr$Start, curChr$Stop), curChr$Chr)
  rangesB <- split(IRanges(targetRegions$Start[i], targetRegions$Stop[i]), targetRegions$Chr[i])
  ov <- countOverlaps(rangesA, rangesB, type="any")>0
  curRegion <- curChr[ov[[1]],]
  curChrNormal <- cnRatiosNormal[cnRatiosNormal$Chr == targetRegions$Chr[i],]
  rangesA <- split(IRanges(curChrNormal$Start, curChrNormal$Stop), curChrNormal$Chr)
  ov <- countOverlaps(rangesA, rangesB, type="any")>0
  curRegionNormal <- curChrNormal[ov[[1]],]
  
  #curRegionNormal <- curChrNormal[targetRegions$Start[i] > curChrNormal$Start & targetRegions$Start[i] < curChrNormal$Stop,]
  
  
  if (dim(curRegion)[1] == 1){
    targetRegions[i,"CNratioMean"] <- 2**curRegion$log2CNratio  
    targetRegions[i,"CNratioNormalMean"] <- 2**curRegionNormal$log2CNratio  
  }
  else if (dim(curRegion)[1] > 1){
    targetRegions[i,"CNratioMean"] <- mean(2**curRegion$log2CNratio)
    targetRegions[i,"CNratioSD"] <- sd(2**curRegion$log2CNratio)
    targetRegions[i,"CNratioNormalMean"] <- mean(2**curRegionNormal$log2CNratio)
    targetRegions[i,"CNratioNormalSD"] <- sd(2**curRegionNormal$log2CNratio)
  }
}

cols <- rep(c('darkgreen','magenta'),50)

targetchrs <- targetRegions$Chr
targetchrs <- gsub("chr","",targetchrs)
targetchrs <- gsub("X",23,targetchrs)
targetchrs <- as.numeric(gsub("Y",24,targetchrs))

targetRegions$col <- "lightgrey"
for ( i in 1:nrow(targetRegions)){
  if ( i == 1){
    curChr <- targetchrs[i]
    curCol <- 1
    targetRegions$col[i] <- cols[curCol]
  }else if ( targetchrs[i] == curChr ){
    targetRegions$col[i] <- cols[curCol]
  }else if (targetchrs[i] != curChr ){
    if (targetchrs[i] < curChr){
      break
    }
    curCol <- curCol+1
    targetRegions$col[i] <- cols[curCol]
    curChr <- targetchrs[i]
  }
}

sampleName <- strsplit(cnRatiosFile,"\\/")[[1]]
sampleName <- sampleName[length(sampleName)]
normalName <- strsplit(cnRatiosFileNormal,"\\/")[[1]]
normalName <- normalName[length(normalName)]
targetName <- strsplit(targetRegionFile,"\\/")[[1]]
targetName <- targetName[length(targetName)]

fileName <- paste(strsplit(sampleName,"\\.")[[1]][1],strsplit(targetName,"\\.")[[1]][1],"png",sep=".")
#targetRegions <- targetRegions[!is.na(targetRegions$CNratioMean),]

png(fileName, 12, 7, units="in", type="cairo", res=300, bg="white")
#png(filename = fileName,width = 1440,height = 960)
layout(mat=matrix(ncol=1,nrow=2,data=c(1,2)),height=c(1,3))
par(mar=c(0,5,2,2))
plotTitle <- paste("Tumor ",strsplit(sampleName,"\\_")[[1]][1]," - Normal ",strsplit(normalName,"\\_")[[1]][1]," ",strsplit(targetName,"\\.")[[1]][1],sep="")
plot(c(1,(nrow(targetRegions)+1)),targetRegions$CNratioMean[c(1,2)],ylim=c(0,2),main=plotTitle,axes=F,pch=20,xlab="",ylab="",cex=0,cex.lab=0.8)
for (i in 1:nrow(targetRegions)){
  text(x=i,y=0,labels = paste(targetRegions$Chr[i],":",targetRegions$Start[i],"-",targetRegions$Stop[i],sep=""),srt=90,adj=0,cex=0.6)
}
par(mar=c(6,5,0,2))
plot(c(1,(nrow(targetRegions)+1)),targetRegions$CNratioMean[c(1,2)],ylim=c(0,4),axes=F,pch=20,xlab="",ylab="Denoised Copy Ratio",cex=0,cex.lab=0.8)
for (i in 1:nrow(targetRegions)){
  polygon(x=c((i-0.5),(i+0.5),(i+0.5),(i-0.5)),y=c(0,0,4,4),col=adjustcolor(targetRegions$col[i],alpha.f = 0.15),border = F)
  #text(x=i,y=2.05,labels = paste(targetRegions$Chr[i],":",targetRegions$Start[i],"-",targetRegions$Stop[i],sep=""),srt=90,adj=0,cex=1.2)
}
points(c(1:nrow(targetRegions)),targetRegions$CNratioNormalMean,pch=20,col='grey',cex=1)
points(c(1:nrow(targetRegions))[targetRegions$CNratioNormalMean < diploidQuant99[1]],targetRegions$CNratioNormalMean[targetRegions$CNratioNormalMean < diploidQuant99[1]],pch=20,col='grey',cex=1.3)
points(c(1:nrow(targetRegions))[targetRegions$CNratioNormalMean > diploidQuant99[2]],targetRegions$CNratioNormalMean[targetRegions$CNratioNormalMean > diploidQuant99[2]],pch=20,col='grey',cex=1.3)
points(c(1:nrow(targetRegions)),targetRegions$CNratioMean,pch=20,cex=1)
#targetRegions$CNratioMean[which(targetRegions$CNratioMean > 4)] <- 4.001
xBar <- c(1:nrow(targetRegions))[!is.na(targetRegions$CNratioSD) & targetRegions$CNratioMean <= 4]
yBarMin <- targetRegions$CNratioMean[!is.na(targetRegions$CNratioSD) & targetRegions$CNratioMean <= 4] - targetRegions$CNratioSD[!is.na(targetRegions$CNratioSD) & targetRegions$CNratioMean <= 4]
yBarMax <- targetRegions$CNratioMean[!is.na(targetRegions$CNratioSD) & targetRegions$CNratioMean <= 4] + targetRegions$CNratioSD[!is.na(targetRegions$CNratioSD) & targetRegions$CNratioMean <= 4]


arrows(x0 = xBar,x1 = xBar,y0=yBarMin,y1=yBarMax,length = 0.02,angle=90,code = 3,lwd=1)
points(c(1:nrow(targetRegions))[targetRegions$CNratioMean < diploidQuant99[1]],targetRegions$CNratioMean[targetRegions$CNratioMean < diploidQuant99[1]],pch=20,col='red',cex=1.3)
points(c(1:nrow(targetRegions))[targetRegions$CNratioMean > diploidQuant99[2]],targetRegions$CNratioMean[targetRegions$CNratioMean > diploidQuant99[2]],pch=20,col='blue',cex=1.3)
#points(c(1:nrow(targetRegions))[targetRegions$CNratioMean < -2],rep(-2,sum(targetRegions$CNratioMean < -2,na.rm=T)),pch=20,col='red',cex=1.6)
points(c(1:nrow(targetRegions))[which(targetRegions$CNratioMean > 4)],rep(4,sum(targetRegions$CNratioMean > 4,na.rm=T)),pch=20,col='blue',cex=1.6)
if ( sum(targetRegions$CNratioMean > 4,na.rm=T) > 0){
	text(c(1:nrow(targetRegions))[which(targetRegions$CNratioMean > 4)],rep(3.85,sum(targetRegions$CNratioMean > 4,na.rm=T)),labels=paste(round(targetRegions$CNratioMean[which(targetRegions$CNratioMean > 4)],1),"±",round(targetRegions$CNratioSD[which(targetRegions$CNratioMean > 4)],1)),cex=0.8)
}
lines(c(0.5,(nrow(targetRegions)+0.5)),rep(diploidQuant99[1],2),lty=2,lwd=1)
lines(c(0.5,(nrow(targetRegions)+0.5)),rep(diploidQuant99[2],2),lty=2,lwd=1)
#lines(c(0.5,(nrow(targetRegions)+0.5)),rep(diploidQuant999[1],2),lty=2,lwd=2)
#lines(c(0.5,(nrow(targetRegions)+0.5)),rep(diploidQuant999[2],2),lty=2,lwd=2)
text(nrow(targetRegions)+1,diploidQuant99[1],labels ="99% CI",adj = c(0,0.5),cex=0.5)
text(nrow(targetRegions)+1,diploidQuant99[2],labels="99% CI",adj = c(0,0.5),cex=0.5)
axis(1,at=c(1:nrow(targetRegions)),labels = targetRegions$Feature,las=2,cex.axis=0.8)
axis(2,cex.axis=0.8)
legend("bottomright",legend=c("Normal","Tumor","Tumor - Deleted","Tumor - Amplified"),col=c("grey","black","red","blue"),cex=c(0.7),bty='n',pch=20)
dev.off()

