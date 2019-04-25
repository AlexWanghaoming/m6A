library(GenomicFeatures)
library(Rsamtools)
library(RColorBrewer)
library(rtracklayer)
txdb <-makeTxDbFromGFF(file = "/home/wanghm/whm/1FG/PH1-fusion-version2_represent_transcript_rm_Mt.gff")  # Here a warning 14387 transcripts totally, but only load 14336.
# calculating per gene cds methylation fraction
m6A_gff <-data.table::fread("/home/wanghm/whm/pacbio_data/PH-1/6mA/xulab/m6A.gff")
m6AGrange <- GRanges(seqnames = m6A_gff$V1,ranges = IRanges(start = m6A_gff$V4, end = m6A_gff$V5),strand = m6A_gff$V7)
  
m6AGrange_total <- GRanges(seqnames = m6A_gff$V1,ranges = IRanges(start = m6A_gff$V4, end = m6A_gff$V5),strand = "*")

# get all gene sequence
# features can be one of: gene, peak ....such GRanges object, methy can be m6A or m4C 
get_feature_methylation_ratio <- function(feature){
  allgeneSeq <- getSeq(FaFile("/home/wanghm/whm/1FG/xulab_ph1.fasta"), feature)
  ph1_gene <- feature # feature is GRanges Object
  
  # count each gene binding methylation sites
  allgeneCounts <- countOverlaps(ph1_gene, m6AGrange_total)
  
  # get each gene AT counts
  A_counts_perGenes <- letterFrequency(allgeneSeq, letters = "A")
  T_counts_perGenes <- letterFrequency(allgeneSeq, letters = "T")
  # total methylation ratio
  total_m6AFra <- allgeneCounts / (A_counts_perGenes + T_counts_perGenes)
  
  # calculate strand-specific methylation ratio
  sameDir_withGenes_6mACounts <- countOverlaps(ph1_gene, m6AGrange)
  diffDir_withGenes_6mACounts <- allgeneCounts - sameDir_withGenes_6mACounts
  samDir_6mAFra <- sameDir_withGenes_6mACounts / A_counts_perGenes
  diffDir_6mAFra <- diffDir_withGenes_6mACounts / T_counts_perGenes
  
  # make df merged
  m6A_fra_df <- data.frame(total=total_m6AFra, sameDir=samDir_6mAFra, diffDir=diffDir_6mAFra, row.names = names(ph1_gene))
  colnames(m6A_fra_df) <- c("total","sameDir","diffDir")
  return(m6A_fra_df)
}
m6A_fra_df <- get_feature_methylation_ratio(feature = genes(txdb))

# filter non-methylation sites genes
m6A_fra_df <- m6A_fra_df[m6A_fra_df$total!=0,]

m6A_fra_df$difference <- m6A_fra_df$sameDir - m6A_fra_df$diffDir
m6A_fra_df$foldchange <- log2(m6A_fra_df$sameDir/m6A_fra_df$diffDir)
squash::hist2(m6A_fra_df$difference, m6A_fra_df$foldchange)
## ********************************************************************* ##
## calculate each strands methylation across genome
ph1_genome <- readDNAStringSet("Fusarium_graminearum.RR1_with_mito.fa")

## plot strand_specific 6mA cross genome chromosomes
plot_strandSpecific6mA <- function(windowSize,chrNum) {
  chrLabels <- unlist(strsplit(names(ph1_genome[chrNum]), " "))[1]
  ran <-
    IRanges(start = seq(1, width(ph1_genome[chrNum]), by = windowSize),
            width = windowSize)
  end(ran[length(ran)]) <- width(ph1_genome[chrNum])
  ph1genome_bins <-
    GRanges(seqnames = rep(chrLabels, length(ran)),
            ranges = ran,
            strand = "*")
  positive_strand_6mA <-
    m6A_strand_specific[which(strand(m6AGrange) == "+")]
  m6A_counts1 <- countOverlaps(ph1genome_bins, positive_strand_6mA)
  negative_strand_6mA <-
    m6A_strand_specific[which(strand(m6AGrange) == "-")]
  m6A_counts2 <- countOverlaps(ph1genome_bins, negative_strand_6mA)
  bins_seq <-
    getSeq(FaFile("Fusarium_graminearum.RR1_with_mito.fa"),
           ph1genome_bins)
  Pos <- m6A_counts1 / letterFrequency(bins_seq, "A")-0.0081
  Neg <- m6A_counts2 / letterFrequency(bins_seq, "T")-0.0081
  rang <- seq(windowSize/2,width(ph1_genome[chrNum]),windowSize)[1:length(Pos)]
  plot(x=rang, y=Pos, col="#ee6b31", lwd=1, type="l",xlab=paste0("GenomicRange_chr",chrNum), 
       ylab="Strand Sepecifc 6mA Occupancy",
       ylim=c(-0.005,0.005))
  lines(x=rang, y=Neg, col="#73c2e7", lwd=1)
  abline(h=0,col="red")
  fast_slow <- read.table("/home/wanghm/Desktop/2speed-state.yellow-purple.txt")
  fast_slow2 <- split(fast_slow,fast_slow$V1)[[chrNum]]     # split return a list
  apply(fast_slow2,1,xx)
}
## add transparent rectangle to plot
xx <- function(x) {
  rect(xleft = x[2],xright = x[3], ybottom = -0.005,ytop = 0.005,
       col = ifelse(unlist(strsplit(x[4],"="))[2] == "yellow", 
                    rgb(1,1,0,alpha = 0.3),
                    rgb(0.63,0.125,0.94, alpha = 0.3)),border = NA)
}
par(mfrow = c(4,1))
plot_strandSpecific6mA(10000,1)
plot_strandSpecific6mA(10000,2)
plot_strandSpecific6mA(10000,3)
plot_strandSpecific6mA(10000,4)

## ********************************************************************* ##
# wilcox test indicate that 6mA enrich in one strand.
# load("/home/wanghm/wanghm_R/methylation/strand_specificAnalysis.RData")
group1 = c(df1[which(df1$SameDir > df1$DiffDir), ][, 2], df1[which(df1$SameDir <df1$DiffDir), ][, 3])
group2 = c(df1[which(df1$SameDir < df1$DiffDir), ][, 2], df1[which(df1$SameDir >df1$DiffDir), ][, 3])
x = data.frame(group1,group2)
library(ggstatsplot)
library(tidyverse)
xx <- gather(x);
xx$key <- as.factor(xx$key)
## use ggbetweenstats for stataitic analysis
ggbetweenstats(xx, x = key, y=value, notch = T) + xlab("groups") + ylab("Methylation Occupancy")
# boxplot(x, ylim=c(0,0.08),col=brewer.pal(2,"Greens"), ylab="6mA Occupancy", xlab="6mA always enrich in one strand")
# text(x=0.05,labels = "   pvalue < 2.2e-16",pos = 4)
# RColorBrewer::display.brewer.all()
### Non-dependence sample t-test and wilcox test for samDir_6mAFra & diffDir_6mAFra
iris
var.test(group1, group2) # before t-test judging whether var is homogeneous( p>0.05 )
t.test(df3$SameDir, df3$DiffDir, var.equal = F, paired = T)
wilcox.test(df3$SameDir, df3$DiffDir, paired = T)$p.value
wilcox.test(group1, group2, paired = T)
t.test(group1, group2, var.equal = T)

