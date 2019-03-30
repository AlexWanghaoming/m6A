library(GenomicFeatures)
library(Rsamtools)
library(RColorBrewer)
library(rtracklayer)
longestTranTxdb <-
  makeTxDbFromGFF(
    file = "/home/wanghm/whm/1FG/35genome/Fusarium_graminearum.RR1.35.chr.gff3",
    format = c("auto", "gff3", "gtf"),
    dataSource = "gff3 file for Fusarium_graminearum",
    organism = "Fusarium graminearum"
  )  # Here a warning 14387 transcripts totally, but only load 14336.
# calculating per gene cds methylation fraction
m6A_gff <-
  data.table::fread(
    "/home/wanghm/whm/pacbio_data/Pacbio/PH-1/6mA/official_methylationResult/m6A.gff",
    skip = 1
  )
read.table("/home/wanghm/whm/pacbio_data/Pacbio/PH-1/6mA/official_methylationResult/m6A.gff")
  m6AGrange <- GRanges(
    seqnames = m6A_gff$V1,
    ranges = IRanges(start = m6A_gff$V4, end = m6A_gff$V5),
    strand = m6A_gff$V7
  )
  m6AGrange_total <- GRanges(
    seqnames = m6A_gff$V1,
    ranges = IRanges(start = m6A_gff$V4, end = m6A_gff$V5),
    strand = "*"
  )
  allgeneSeq <-
    Rsamtools::getSeq(FaFile("Fusarium_graminearum.RR1_with_mito.fa"),
                      genes(longestTranTxdb))
ph1_gene <- genes(longestTranTxdb)
allgeneCounts <- countOverlaps(ph1_gene, m6AGrange_total)
A_counts_perGenes <- letterFrequency(allgeneSeq, letters = "A")
T_counts_perGenes <- letterFrequency(allgeneSeq, letters = "T")
total_m6AFra <-
  allgeneCounts / (A_counts_perGenes + T_counts_perGenes)
sameDir_withGenes_6mACounts <-
  countOverlaps(ph1_gene, m6AGrange)
diffDir_withGenes_6mACounts <-
  allgeneCounts - sameDir_withGenes_6mACounts

samDir_6mAFra <- sameDir_withGenes_6mACounts / A_counts_perGenes
diffDir_6mAFra <- diffDir_withGenes_6mACounts / T_counts_perGenes
m6A_fra_df <-
  data.frame(
    total = total_m6AFra,
    samDir = samDir_6mAFra,
    diffDir = diffDir_6mAFra,
    strand = strand(ph1_gene),row.names = names(ph1_gene)
  )

colnames(m6A_fra_df) <- c("total", "SameDir", "DiffDir", "strand")
m6A_fra_df <-
  m6A_fra_df[which(m6A_fra_df$total != 0), ]  # remove non-methylation gene

# df3 <- m6A_fra_df[which(order(m6A_fra_df$total)>length(m6A_fra_df$total)*0.75),] # top 1/4 methyFra genes
df1 <- m6A_fra_df[order(m6A_fra_df$total, decreasing = T), ]
df1[which(df1$SameDir == 0),]
list1 <- split(df1,rep(1:ceiling(nrow(df1) / 50), each=50, length.out=nrow(df1)))  # split dataframe into equal groups
groupTest <- function(x){
  samdir <- x[,2]
  diffdir <- x[,3]
  pVal <- wilcox.test(samdir,diffdir,paired = T)$p.value
  return(pVal)
}
sapply(list1, groupTest)   # calculate groups' pvalue

# library(squash)
# hist2(log2(m6A_fra_df$SameDir),log2(m6A_fra_df$DiffDir), xlim=c(-10,-3), ylim=c(-10,-3))

## ********************************************************************* ##
## calculate each strands methylation across genome
ph1_genome <-
  readDNAStringSet("Fusarium_graminearum.RR1_with_mito.fa")

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
  # sp1 = spline(x=rang, y=Pos,n = 2000)
  # sp2 = spline(x=rang, y=Neg,n = 2000)
  plot(x=rang, y=Pos, col="#ee6b31", lwd=1, type="l",xlab=paste0("GenomicRange_chr",chrNum), 
       ylab="Strand Sepecifc 6mA Occupancy",
       ylim=c(-0.005,0.005))
  lines(x=rang, y=Neg, col="#73c2e7", lwd=1)
  abline(h=0,col="red")
  fast_slow <- read.table("/home/wanghm/Desktop/2speed-state.yellow-purple.txt")
  fast_slow2 <- split(fast_slow,fast_slow$V1)[[chrNum]]     # split return a list
  apply(fast_slow2,1,xx)
  # legend(x = "topleft",
  #        legend = c("sense strand", "antisense strand"),
  #        fill = c("#ee6b31","#73c2e7"),cex=0.75)
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
group1 = c(df1[which(df1$SameDir > df1$DiffDir), ][, 2], df1[which(df1$SameDir <df1$DiffDir), ][, 3])
group2 = c(df1[which(df1$SameDir < df1$DiffDir), ][, 2], df1[which(df1$SameDir >df1$DiffDir), ][, 3])
boxplot(x = data.frame(group1,group2), ylim=c(0,0.08),col=brewer.pal(2,"Greens"), ylab="6mA Occupancy", xlab="6mA always enrich in one strand")
text(x=0.05,labels = "   pvalue < 2.2e-16",pos = 4)
# RColorBrewer::display.brewer.all()
### Non-dependence sample t-test and wilcox test for samDir_6mAFra & diffDir_6mAFra

var.test(group1, group2) # before t-test judging whether var is homogeneous( p>0.05 )
t.test(df3$SameDir, df3$DiffDir, var.equal = F, paired = T)
wilcox.test(df3$SameDir, df3$DiffDir, paired = T)$p.value
wilcox.test(group1, group2, paired = T)
t.test(group1, group2, var.equal = T)

