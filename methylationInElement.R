## 
## these script for m6A sites annotation and calculate m6A ratio on each genomic elements
## Wang Haoming

library(tidyverse)
library(GenomicFeatures)
library(Rsamtools)
library(RColorBrewer)
library(rtracklayer)
txdb <- makeTxDbFromGFF(file = "/home/wanghm/whm/1FG/PH1-fusion-version2_represent_transcript_rm_Mt.gff") 
m6A_gff <- data.table::fread("/home/wanghm/whm/pacbio_data/PH-1/6mA/xulab/m6A.gff")
m6AGrange <- GRanges(seqnames = m6A_gff$V1,ranges = IRanges(start = m6A_gff$V4, end = m6A_gff$V5),strand = m6A_gff$V7)
m6AGrange_total <- GRanges(seqnames = m6A_gff$V1,ranges = IRanges(start = m6A_gff$V4, end = m6A_gff$V5),strand = "*") # m6A_total do not distinct pos/neg strands
### ************** count 6mA in Genes
ge <- genes(txdb)
siteCountsOngenes <- countOverlaps(ge,m6AGrange_total)
table(siteCountsOngenes)
# concentrate on methylated genes
target_gene.list <- names(siteCountsOngenes[siteCountsOngenes!=0])  # 15598 genes totally 
## plot barplot for 6mA-marked genes counts, 
  ggpubr::ggbarplot(reshape2::melt(data.frame(1425, NROW(siteCountsOngenes)-1425)), x="variable", y= "value", label = T, fill = "lightblue", width = 0.6)+
  scale_y_continuous(expand = c(0,0),limits = c(0,17000)) +
    scale_x_discrete(breaks = c("X1425", "NROW.siteCountsOngenes....1425"), labels = c("non-6mA genes","6mA_marked genes")) + 
    xlab("") + ylab("Number of genes") + theme(axis.text.x = element_text(vjust = -1,size = 10))
  
## ********************************************************************* ##
library(ChIPseeker)
library(tidyverse)
library(tidyselect)
m6A_anno.gr <- annotatePeak(m6AGrange_total, tssRegion = c(-500,300), TxDb = txdb, 
                         genomicAnnotationPriority=c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Intergenic"), 
                         addFlankGeneInfo = T) %>% as.GRanges()
pm <- promoters(txdb, upstream = 500, downstream = 300)

m6A_anno <- m6A_anno.gr$annotation

# vars selection
a1 <- length(starts_with("Exon", vars = m6A_anno))
a2 <- length(starts_with("Intron", vars = m6A_anno))
a3 <- length(starts_with("Promoter", vars = m6A_anno))
a4 <- length(starts_with("3' UTR", vars = m6A_anno))
a5 <- length(starts_with("5' UTR", vars = m6A_anno))
a6 <- length(starts_with("Downstream", vars = m6A_anno))
a7 <- length(starts_with("Distal Intergenic",vars = m6A_anno))

## pie plot
pie.data <- data.frame(group = c("Promoter", "Exon", "Intron", "5'UTR", "3'UTR", "Intergenic"), 
                       counts=c(a3,a1,a2,a5,a4,a6+a7))
aa <- pie.data$counts/sum(pie.data$counts)
labs <- paste0(round(aa,3)*100,"%")
p <- ggpubr::ggpie(pie.data, x = "counts", label = labs[c(1,3,6,2,4,5)], fill = "group", lab.pos = "in",
                   color="white", palette = RColorBrewer::brewer.pal(6, "Set3"), legend.title="")
p1 <- ggpar(p, legend = "right", tickslab = F)

## *******************************  
calculate_feature_methyFra <- function(x) {
    allcdsSeq <- getSeq(FaFile("/home/wanghm/whm/1FG/xulab_ph1.fasta"), x)
    methylationFraction <- sum(countOverlaps(x, m6AGrange_total)) / sum(letterFrequency(allcdsSeq, letters = c("A", "T")))
    return(methylationFraction)
}

get_feature_methylation_ratio <- function(feature){
  # feature <- feature[names(feature) %in% target_gene.list]
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
  # m6A_fra_df <- m6A_fra_df[m6A_fra_df[,1]!=0,]
  # m6A_fra_df <- na.omit(m6A_fra_df)
  return(m6A_fra_df)
}   # *** strand_specific methylation ratio

# make features' GRange obj
pmt.gr <- promoters(txdb, upstream = 500, downstream = 300)
cds.gr <- cds(txdb)
exon.gr <- exons(txdb)
fiveUTR.gr <- unlist(fiveUTRsByTranscript(txdb),use.names = F)
intron.gr <- unlist(intronsByTranscript(txdb), use.names = F)
threeUTR.gr <- unlist(threeUTRsByTranscript(txdb), use.names = F)
## extract genic and intergenic region as GRanges object 
genic.gr <- reduce(ge, ignore.strand=T)
intergenic <- gaps(genic.gr)
intergenic.gr <- intergenic[strand(intergenic)=="*"]


b1 <- calculate_feature_methyFra(pmt.gr)
b2 <- calculate_feature_methyFra(cds.gr)
b3 <- calculate_feature_methyFra(exon.gr)
b4 <- calculate_feature_methyFra(fiveUTR.gr)
b5 <- calculate_feature_methyFra(threeUTR.gr)
b6 <- calculate_feature_methyFra(intron.gr)
b7 <- calculate_feature_methyFra(genic.gr)
b8 <- calculate_feature_methyFra(intergenic.gr)

methyRatio <- data.frame(type=c("Promoter", "CDS", "Exon", "5UTR", "3UTR", "Intron", "Genic", "Intergenic"), 
                         ratio=c(b1,b2,b3,b4,b5,b6,b7,b8))
## plot feature's methylation ratio
p1 <- ggpubr::ggdotchart(data = methyRatio, x = "type", y = "ratio", color = "type", 
                   add = "segments", rotate = T, label = round(methyRatio$ratio, 4), dot.size = 4)
ggpar(p1, legend = "right", tickslab = F, ylab = "6mA ratio", xlab = "")
  



