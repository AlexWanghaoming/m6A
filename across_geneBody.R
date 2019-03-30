library(Rsamtools)
library(GenomicFeatures)
library(data.table)
txdb <- makeTxDbFromGFF(file = "/home/wanghm/whm/1FG/35genome/Fusarium_graminearum.RR1.35.chr.gff3")
ge <- genes(txdb)
start(ge) <- start(ge)-500; end(ge) <- end(ge) + 500

m6AGrange2 <- m6AGrange
strand(m6AGrange2[strand(m6AGrange2)=="+"]) = "*"
strand(m6AGrange2[strand(m6AGrange2)=="-"]) = "+"
strand(m6AGrange2[strand(m6AGrange2)=="*"]) = "-"

m4CGrange2 <- m4CGrange
strand(m4CGrange2[strand(m4CGrange2)=="+"]) = "*"
strand(m4CGrange2[strand(m4CGrange2)=="-"]) = "+"
strand(m4CGrange2[strand(m4CGrange2)=="*"]) = "-"

# getwd()
out_put_fasta <- function(){
  flag=0
  for(i in list(m4CGrange,m4CGrange2)){
    geneMethylationCount <- countOverlaps(ge, i)      # remove no methelation genes
    ge_trimed <- ge[names(ge) %in% names(geneMethylationCount[geneMethylationCount!=0])]
    ge_trimed <- ge_trimed[start(ge_trimed)>0]
    geneBodySeq <- getSeq(FaFile("/home/wanghm/whm/1FG/35genome/Fusarium_graminearum.RR1_with_mito.fa"), ge_trimed)
    writeXStringSet(geneBodySeq, paste0("/home/wanghm/wanghm_R/methylation/gebody",flag,".fasta"))
    # system("awk 'BEGIN{i=1}/>/{$0=\">\"i;i++;print$0;next}{print$0}' /home/wanghm/wanghm_R/methylation/gebody0.fasta > /home/wanghm/Desktop/gebody0.fasta")
    flag <- flag + 1
  }
}
out_put_fasta()

getGeneBodyMatrix <- function(window,GRangeObj){
  peak.cov1 <- coverage(GRangeObj[strand(GRangeObj) == "+"]) # coverage return Rle Object
  peak.cov2 <- coverage(GRangeObj[strand(GRangeObj) == "-"]) # for example with m6A methylation
  #### add methylationLevels info to peak.cov1 & peak.cov2

  # i=1;j=1
  # for (x in names(peak.cov1)) {
  #   vx <- as.vector(unlist(peak.cov1[x]))
  #   vx[vx != 0] <- vx[vx != 0]*unlist(positive_methylationLevel[x])
  #   peak.cov_1[[i]] <- Rle(vx)   ## when append elements to a list, must use num index, not characters index.
  #   i=i+1
  # }
  # for (y in names(peak.cov2)) {
  #   vy <- as.vector(unlist(peak.cov2[y]))
  #   vy[vy != 0] <- vy[vy != 0]*unlist(negative_methylationLevel[y])
  #   peak.cov_2[[j]] <- Rle(vy)
  #   j=j+1
  # }
  # names(peak.cov_1) <- names(peak.cov1)
  # names(peak.cov_2) <- names(peak.cov2)
  #
  # peak.cov1 <- peak.cov_1
  # peak.cov2 <- peak.cov_2
  # For peak.cov1
  cov.len1 <- elementNROWS(peak.cov1)
  cov.width1 <- GRanges(seqnames = names(cov.len1), IRanges(start = rep(1,length(cov.len1)), end = cov.len1))
  windows1 <- subsetByOverlaps(window, cov.width1, type = "within")  # reconfirm all promoter regions are in genome
  chr.idx1 <- intersect(names(peak.cov1), unique(as.character(seqnames(windows1))))

  # For peak.cov2
  cov.len2 <- elementNROWS(peak.cov2)
  cov.width2 <- GRanges(seqnames = names(cov.len2), IRanges(start = rep(1,length(cov.len2)), end = cov.len2))
  windows2 <- subsetByOverlaps(window, cov.width2, type = "within")  # reconfirm all promoter regions are in genome
  chr.idx2 <- intersect(names(peak.cov2), unique(as.character(seqnames(windows2))))

  peakView1 <- Views(peak.cov1[chr.idx1], as(windows1[strand(windows1)=="+"], "IntegerRangesList")[chr.idx1])## showMethods("coerce") can get all avaliable transforms
  peakView2 <- Views(peak.cov2[chr.idx2], as(windows2[strand(windows2)=="-"], "IntegerRangesList")[chr.idx2])## showMethods("coerce") can get all avaliable transforms
  diff_strandList <- list(peakView1,peakView2)
  diffStrandMatrixList = list()
  for(y in 1:2){

    forEachView <- function(z,binNum=70){   # for each gene: split gene length to 50 bins and calculate sum methylation counts
      as.vector(z)
      aa <- floor(length(z)/50)
      ifelse(y==1, split.factor <- c(rep(-10:-1, each=50), rep(1:50,each=floor((length(z)-1000)/50)), rep(50,each=(length(z)-1000)%%50), rep(51:60, each=50)),
                   split.factor <- c(rep(-10:-1, each=50), rep(1,each=(length(z)-1000)%%50), rep(1:50,each=floor((length(z)-1000)/50)), rep(51:60, each=50))
      )
      binCountSum <- sapply(split(z, split.factor), sum)
      return(binCountSum)
    }

    k=diff_strandList[[y]]
    library(BiocParallel)
    MulticoreParam = MulticoreParam(workers = 4)
    tagMatrixList <-bplapply(k, function(x) t(viewApply(x, forEachView)), BPPARAM = MulticoreParam)
    tagMatrix <- do.call("rbind", tagMatrixList)
    rownames(tagMatrix) <- 1:nrow(tagMatrix)
    tagMatrix <- tagMatrix[rowSums(tagMatrix) != 0,]   # remove no methylation genes
    diffStrandMatrixList[[y]] = tagMatrix
  }
  return(diffStrandMatrixList)
}
# load("/home/wanghm/wanghm_R/methylation/across_geneBody.RData")

m6A_matrix1 <- getGeneBodyMatrix(ge, m6AGrange)
m6A_matrix2 <- getGeneBodyMatrix(ge, m6AGrange2)
m4C_matrix1 <- getGeneBodyMatrix(ge, m4CGrange)
m4C_matrix2 <- getGeneBodyMatrix(ge, m4CGrange2)

acrossGeneBody <- function(genomeFile, upstream = NULL,downstream = NULL, methy) {
  for(j in 1:2) {
    if (j==1) {
      ifelse(methy=="m6A",tt <- m6A_matrix1, tt <- m4C_matrix1)
      ifelse(methy=="m6A",cc <- read.csv("/home/wanghm/wanghm_R/methylation/m6A_genebody0",header = F), 
                          cc <- read.csv("/home/wanghm/wanghm_R/methylation/m4C_genebody0",header = F))
    }else if(j==2){
      ifelse(methy=="m6A",tt <- m6A_matrix2, tt <- m4C_matrix2)
      ifelse(methy=="m6A",cc <- read.csv("/home/wanghm/wanghm_R/methylation/m6A_genebody1",header = F),
                          cc <- read.csv("/home/wanghm/wanghm_R/methylation/m4C_genebody1",header = F))

    }
    tt[[2]] <- t(apply(tt[[2]],1,rev))
    tt1 <- do.call(rbind,tt)
    ss <- colSums(tt1)  #  tt1 is a matrix with 70 cols, ss are 70 values, each value represent mean methy site counts in bins
    countBase <- cc$V1
    ss <- ss / countBase
    pos <- c(seq(-475,0,50), seq(20,2000,40), seq(2025,2500,50))
    if (j==1) {
      plot(pos,ss, type="l",col="red", xlab = "GeneBody",ylab = paste(methy,"Methylation Occupancy"," "),ylim = c(0.005,0.013),lwd=2,xaxt = 'n')
      axis(1,at = c(0,2000), labels = c("TSS", "TES"))
      legend("topright", c("modification in coding strand","modification in template strand"), lty = c(1,1,2), lwd = 2, 
             col = c("red", "green"), bty = "n", cex = 0.8, xpd = TRUE)
    }else if(j==2){
      lines(pos,ss,col="green",lwd=2)
    }
    # else{
    #   lines(xx,yy, col="black", lty=2)
    # }
  }
}
acrossGeneBody(genomeFile = "/home/wanghm/whm/1FG/35genome/Fusarium_graminearum.RR1_with_mito.fa",methy = "m4C")



















