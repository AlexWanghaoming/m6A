library(Rsamtools)
library(GenomicFeatures)
library(data.table)
library(ChIPseeker)
# load("/home/wanghm/wanghm_R/methylation/methyFPKM.RData")
m4C_gff <- fread("/home/wanghm/whm/pacbio_data/PH-1/6mA/xulab/m4C.gff", skip = 1)
m6A_gff <- fread("/home/wanghm/whm/pacbio_data/PH-1/6mA/xulab/m6A.gff", skip = 1)
m4CGrange = GRanges(seqnames = m4C_gff$V1,ranges = IRanges(start = m4C_gff$V4, end = m4C_gff$V5),strand = m4C_gff$V7)
m6AGrange = GRanges(seqnames = m6A_gff$V1,ranges = IRanges(start = m6A_gff$V4, end = m6A_gff$V5),strand = m6A_gff$V7)
txdb <- makeTxDbFromGFF(file = "/home/wanghm/whm/1FG/PH1-fusion-version2_represent_transcript_rm_Mt.gff")
promoter <- promoters(genes(txdb), upstream = 500, downstream = 1000)
######## rm(list = c(-merged_df))
######## val:merged_df
# highfpkm_geneName <- merged_df[which(merged_df$V2.x>30),][,3]
# lowfpkm_geneName <- merged_df[which(merged_df$V2.x<1),][,3]
# ########get high/low fpkm genes promoter GRanges
# highfpkm_promoter <- promoter[which(names(promoter) %in% highfpkm_geneName),]
# lowfpkm_promoter <- promoter[which(names(promoter) %in% lowfpkm_geneName),]

############# add single base Methylation Levels to methylation Frequence
# methyList <- split(m4C_gff,f = m4C_gff$V1)
# fun1 <- function(x,strand){
#   methylationLevel <- as.numeric(sapply(strsplit(x[which(x$V7 == strand),]$V9,";|="), function(x) x[8]))
# }
# 
# positive_methylationLevel <- lapply(methyList, fun1, strand="+")
# negative_methylationLevel <- lapply(methyList, fun1, strand="-")
# peak.cov1 <- coverage(m4CGrange[strand(m4CGrange) == "+"],weight = unlist(positive_methylationLevel))
# peak.cov2 <- coverage(m4CGrange[strand(m4CGrange) == "-"],weight = unlist(negative_methylationLevel))

# load("/home/wanghm/wanghm_R/methylation/across_geneBody.RData")
m6AGrange2 <- m6AGrange
strand(m6AGrange2[strand(m6AGrange2)=="+"]) = "*"
strand(m6AGrange2[strand(m6AGrange2)=="-"]) = "+"
strand(m6AGrange2[strand(m6AGrange2)=="*"]) = "-"

m4CGrange2 <- m4CGrange
strand(m4CGrange2[strand(m4CGrange2)=="+"]) = "*"
strand(m4CGrange2[strand(m4CGrange2)=="-"]) = "+"
strand(m4CGrange2[strand(m4CGrange2)=="*"]) = "-"

forEachStrand <- function(y) {### operate Views list
  tagMatrixList <-lapply(y, function(x) t(viewApply(x, as.vector)))
  tagMatrix <- do.call("rbind", tagMatrixList)
  rownames(tagMatrix) <- 1:nrow(tagMatrix)
  tagMatrix <- tagMatrix[rowSums(tagMatrix) != 0,]   # remove no methylation genes
  return(tagMatrix)
}
##### Function getAroundTSSMatrix 
getAroundTssMatrix <- function(window,GRangeObj){
  peak.cov1 <- coverage(GRangeObj[strand(GRangeObj) == "+"]) # coverage return Rle Object
  peak.cov2 <- coverage(GRangeObj[strand(GRangeObj) == "-"]) # coverage return Rle Object
  cov.len1 <- elementNROWS(peak.cov1)
  cov.width1 <- GRanges(seqnames = names(cov.len1), IRanges(start = rep(1,length(cov.len1)), end = cov.len1))
  windows1 <- subsetByOverlaps(window, cov.width1, type = "within")  # reconfirm all promoter regions are in genome
  chr.idx1 <- intersect(names(peak.cov1), unique(as.character(seqnames(windows1))))
  peakView1 <- Views(peak.cov1[chr.idx1], as(windows1[strand(windows1)=="+"], "IntegerRangesList")[chr.idx1])## showMethods("coerce") can get all avaliable transforms
  
  
  # For peak.cov2
  cov.len2 <- elementNROWS(peak.cov2)
  cov.width2 <- GRanges(seqnames = names(cov.len2), IRanges(start = rep(1,length(cov.len2)), end = cov.len2))
  windows2 <- subsetByOverlaps(window, cov.width2, type = "within")  # reconfirm all promoter regions are in genome
  chr.idx2 <- intersect(names(peak.cov2), unique(as.character(seqnames(windows2))))
  peakView2 <- Views(peak.cov2[chr.idx2], as(windows2[strand(windows2)=="-"], "IntegerRangesList")[chr.idx2])## showMethods("coerce") can get all avaliable transforms

  diff_strandList <- list(peakView1,peakView2)
  diffStrandMatrixList <- lapply(diff_strandList, forEachStrand)
  return(diffStrandMatrixList)
}

aroundTss <- function(genomeFile, upstream = NULL,downstream = NULL, windowSize, methy) {
  k = 1
  ifelse(methy=="6mA", grlist <- list(m6AGrange,m6AGrange2), grlist <- list(m4CGrange,m4CGrange2))
  for(j in grlist) {
    geneMethylationCount <- countOverlaps(promoter,j)
    promoter_trimed <- promoter[names(promoter) %in% names(geneMethylationCount[geneMethylationCount!=0])]
    promoterSeq <- getSeq(FaFile(genomeFile), promoter_trimed)  # get sequences by a GRange
    tt <- getAroundTssMatrix(promoter, j)
    tt[[2]] <- t(apply(tt[[2]],1,rev))
    tt1 <- do.call(rbind,tt)
    ss <- colSums(tt1)
    if(k==1){
      ifelse(methy=="6mA",countBase <- consensusMatrix(promoterSeq)[1,], countBase <- consensusMatrix(promoterSeq)[2,]) 
    }else if(k==2){
      ifelse(methy=="6mA",countBase <- consensusMatrix(promoterSeq)[4,], countBase <- consensusMatrix(promoterSeq)[3,])
    }
    ss <- ss / countBase
    pos <- value <- NULL
    tagCount <- data.frame(pos = c(-499:1000), value = ss)
    winPos <- matrix(tagCount$pos, ncol = windowSize, byrow = TRUE)
    winValue <- matrix(tagCount$value, ncol = windowSize, byrow = TRUE)
    xx = rowMeans(winPos)
    yy = rowMeans(winValue)
    if (k == 1) {
      plot(xx,yy, type="l",col="#6b318a", xlab = "Distance from TSS(bp)",ylab = paste(methy,"methylation occupancy",sep = " "), ylim=c(0.005,0.012), lwd=2)
      legend("topright", c("Coding strand","Template strand"), lty = c(1,1,2), lwd = 2, 
             col = c("#6b318a", "#48bfd9"), bty = "n", cex = 0.8, xpd = TRUE)
    }else if(k == 2){
      lines(xx,yy,col="#48bfd9",lwd=2)
    }
    k=k+1
  }
}
aroundTss(genomeFile = "/home/wanghm/whm/1FG/xulab_ph1.fasta",windowSize = 50, methy = "6mA")





