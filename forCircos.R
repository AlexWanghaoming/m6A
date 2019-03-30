library(GenomicFeatures)
library(Rsamtools)
library(RColorBrewer)
library(circlize)
library(tidyverse)
m6A_gff <- data.table::fread("/Users/alexwang/data/pacbio/0.fg/genome/m6A.gff", skip = 1)
m6A_strand_specific <- GRanges(seqnames = m6A_gff$V1,ranges = IRanges(start = m6A_gff$V4, end = m6A_gff$V5), strand = m6A_gff$V7)
ph1_genome <- readDNAStringSet("/Users/alexwang/data/pacbio/0.fg/genome/Fusarium_graminearum.RR1_with_mito.fa",format = "fasta")
methy.level <- strsplit(m6A_gff$V9, ";|=") %>%
  sapply(function(x) x[8])
m6A_strand_specific$level <- as.numeric(methy.level)      ## 加入甲基化位点的levels信息

circos_strandSpecific6mA_pre <- function(windowSize,chrNum) {    # Circos图数据预处理
  chrLabels = unlist(strsplit(names(ph1_genome[chrNum])," "))[1]
  ran <- IRanges(start = seq(1, width(ph1_genome[chrNum]), by = windowSize),
            width = windowSize)
  end(ran[length(ran)]) <- width(ph1_genome[chrNum])
  ph1genome_bins <- GRanges(seqnames = rep(chrLabels, length(ran)), ranges = ran, strand = "*")   # 将全基因组划为窗口

  positive_strand_6mA <- m6A_strand_specific[which(strand(m6A_strand_specific) == "+")]
  m6A_counts1 <- countOverlaps(ph1genome_bins, positive_strand_6mA)
  
  negative_strand_6mA <- m6A_strand_specific[which(strand(m6A_strand_specific) == "-")]
  m6A_counts2 <- countOverlaps(ph1genome_bins, negative_strand_6mA)
  
  cov <- coverage(m6A_strand_specific, weight = "level")
  bin.view <- Views(cov[chrNum], as(ph1genome_bins, "IntegerRangesList"))## showMethods("coerce") can get all avaliable transforms
  bin_avg.levels <- bin.view[[1]]%>%
                      viewApply(sum) / (m6A_counts1+m6A_counts2)   # 计算窗口内的甲基化levels
  # binnedAverage(ph1genome_bins, cov[1:4], varname = "level")     # 也可以用binnedAverage（)计算
  bin_avg.levels[is.infinite(bin_avg.levels)] <- 0
  
  bins_seq <- getSeq(FaFile("/Users/alexwang/data/pacbio/0.fg/genome/Fusarium_graminearum.RR1_with_mito.fa"),
           ph1genome_bins)
  density <- countOverlaps(ph1genome_bins, m6A_strand_specific) / letterFrequency(bins_seq, "AT")  # 计算总甲基化频率
  Pos <- (m6A_counts1 / letterFrequency(bins_seq, "A"))
  Neg <- (m6A_counts2 / letterFrequency(bins_seq, "T"))
  GC_content <- (letterFrequency(bins_seq, "GC") / letterFrequency(bins_seq, "ATCG"))  # 计算GC含量

  mcols(ph1genome_bins) <- cbind(density,Pos,Neg, GC_content, bin_avg.levels)  ## 将以上信息作为 mcols 加入GRanges对象
  return(ph1genome_bins) 
}
get_plot_df <- function(win.size){
  total_info <- c(circos_strandSpecific6mA_pre(win.size,1),
                  circos_strandSpecific6mA_pre(win.size,2),
                  circos_strandSpecific6mA_pre(win.size,3),
                  circos_strandSpecific6mA_pre(win.size,4))
  bedDF <- as.data.frame(total_info)[,c(-4,-5)]  # remove 4th,5th cols
  bedDF[,1] <- as.character(bedDF[,1])  # factor -> charactors
  bedDF[,1] <- paste0("chr", bedDF[,1])
  return(bedDF)
}

# 绘制circos
chro <- c("chr1", "chr2", "chr3", "chr4")
starts <- c(0,0,0,0)
ends <- c(11760891, 8997558, 7792947, 9395062)
genoCir <- data.frame(chr=chro, start=starts, end=ends)

circos.clear()
circos.par(gap.degree = 8, start.degree = 30, track.height = 1.0, cell.padding = c(0,0,0,0))
circos.genomicInitialize(data = genoCir, plotType = "labels") # 不显示刻度，只显示labels
circos.genomicTrackPlotRegion(ylim = c(0,0.5), 
                              bg.col = c("#3c1240", "#6a205d", "#595959", "#39a2ae"),
                              bg.border = 1, track.height = 0.02)
## Track1: 链特异性甲基化频率
bedDF <- get_plot_df(win.size = 50000)
bedList <- list(bedDF[,c(1,2,3,5)],bedDF[,c(1,2,3,6)])
circos.genomicTrackPlotRegion(bedList,ylim = c(-0.02,0.02),track.height = 0.2,bg.border="white",
                              panel.fun = function(region,value,...){
                                i = getI(...)
                                if(i == 1){
                                  circos.genomicLines(region,value,lwd = 1.5,col = "#22b21a",...)
                                }else{
                                  circos.genomicLines(region,-value,lwd = 1.5, col = "#4281a4",...)
                                }
                                cell.xlim=get.cell.meta.data("cell.xlim")
                                circos.lines(cell.xlim,y = c(0,0))
                              })
## Track2: 总甲基化频率
f1 <- colorRamp2(breaks = c(0.001, 0.015),colors = c("white","#ffa630"))
circos.genomicTrackPlotRegion(bedDF[,c(1,2,3,4)],track.height = 0.1,
                              panel.fun = function(region,value,...){
                                circos.genomicRect(region,value,col = f1(value[[1]]),border = NA,...)
                                # cell.xlim=get.cell.meta.data("cell.xlim")
                                # cell.ylim=get.cell.meta.data("cell.ylim")
                                # circos.lines(cell.xlim,col = "white",y = c((cell.ylim[1]+cell.ylim[2])/2,
                                #                              (cell.ylim[1]+cell.ylim[2])/2))
                              })
## Track3: 甲基化水平
bedDF2 <- get_plot_df(win.size = 100000)  # 将窗口大小改为10万
m <- mean(bedDF2$bin_avg.levels)
bedDF2.low <- filter(bedDF2, bin_avg.levels<m) # dplyr 
bedDF2.high <- filter(bedDF2, bin_avg.levels>m)
list2 <- list(bedDF2.low, bedDF2.high)
circos.genomicTrackPlotRegion(list2,track.height = 0.1,ylim = c(0.2,0.8),bg.col="#84bc9c",
                              panel.fun = function(region, value,...){
                                i = getI(...)
                                if (i == 1) {
                                  circos.genomicPoints(region, 0.3, cex=0.1, col = "red")# 甲基化水平高的点用绿，低用红色表示
                                }else{
                                  circos.genomicPoints(region, 0.7, cex=0.1, col = "green")
                                }
                              })
## Track4: GC 含量
f2 <- colorRamp2(breaks = c(0.2,0.8), colors = c("white", "red"))
circos.genomicTrackPlotRegion(bedDF[,c(1,2,3,7)], track.height = 0.1,
                              panel.fun = function(region, value, ...){
                                circos.genomicRect(region, value, col = f2(value[[1]]), border = NA, ...)
                              })






