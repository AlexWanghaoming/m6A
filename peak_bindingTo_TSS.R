library("ChIPseeker")
library("ggplot2")
library("GenomicFeatures")
txdb <-
  makeTxDbFromGFF(
    file = "/home/wanghm/whm/1FG/fg.gtf",
    format = c("auto", "gff3", "gtf"),
    dataSource = "gtf file for Fusarium_graminearum",
    organism = "Fusarium graminearum"
  )
rtracklayer::export(txdb, "sss.gff", "gff")  # export gff file
rtracklayer::readGFF("/home/wanghm/whm/1FG/FG-fix_ORF-fix_nofusion5-rm_repeat.gtf")
txdb_noCDS <-
  makeTxDbFromGFF(
    file = "/home/wanghm/whm/1FG/FG_noCDS.gtf",
    format = c("auto", "gff3", "gtf"),
    dataSource = "gtf file for Fusarium_graminearum",
    organism = "Fusarium graminearum"
  )
txdb_noCDS
# any(transcripts(txdb)$tx_name == "FGRRES_00036.10")
# transcripts(txdb_noCDS)
#
# genes(txdb)$gene_id

###########txdb object To GRange object
# genes(txdb)
# exon_txdb <- exons(txdb)
# as.data.frame(transcripts(txdb))
# txdb
# transcriptsBy(txdb,by = "gene")
############## operate on GRange objext
# seqnames(exon_txd)
# transcript_txdb<- transcripts(txdb_noCDS)
# mcols(transcript_txdb)
# Rle("chr1",10)
# as.vector(Rle("chr3",19))
# ranges(exon_txdb)
# strand(exon_txdb)
getwd()

promoter <- getPromoters(TxDb = txdb,
                         upstream = 1000,
                         downstream = 1000)
promoter
# bio_region1 <- getBioRegion(TxDb = txdb_noCDS,upstream = 300,downstream = 50, by = "e")

# genes(txdb_noCDS)
# promoters(genes(txdb_noCDS))
################ Start Plot
gRange_peaks <-
  GenomicRanges::GRangesList(
    B = readPeakFile("/home/wanghm/wanghm_R/B_BinQuasiPeaks.bed"),
    S = readPeakFile("/home/wanghm/wanghm_R/S_BinQuasiPeaks.bed")
  )
############ cut bed file 1st col "chr1" to "1"
system(
  "cut -b 4- /home/wanghm/wanghm_R/B_BinQuasiPeaks.bed > /home/wanghm/wanghm_R/B_BinQuasiPeaks_rmChr.bed"
)
system(
  "cut -b 4- /home/wanghm/wanghm_R/S_BinQuasiPeaks.bed > /home/wanghm/wanghm_R/S_BinQuasiPeaks_rmChr.bed"
)
tagMatrixList <-
  lapply(
    c(B_peaks = "/home/wanghm/wanghm_R/B_BinQuasiPeaks_rmChr.bed",
      S_peaks = "/home/wanghm/wanghm_R/S_BinQuasiPeaks_rmChr.bed"),
    getTagMatrix,
    windows = promoter
  )

############### High dpi picture
pdf(
  file = "B&S_500&500_profile.pdf",
  width = 18,
  height = 15,
  bg = "white"
)
# plotAnnoPie(B_Annotation)
# upsetplot(B_Annotation)

dev.off()
###################################
###################################  Peaks Annotations
B_Annotation <-
  annotatePeak(peak = "/home/wanghm/wanghm_R/B_BinQuasiPeaks_rmChr.bed",
               TxDb = txdb_noCDS,
               tssRegion = c(-1000, 1000))
S_Annotation <-
  annotatePeak(peak = "/home/wanghm/wanghm_R/S_BinQuasiPeaks_rmChr.bed",
               TxDb = txdb_noCDS,
               tssRegion = c(-2000, 2000))

peakAnnoList <- lapply(
  c(B_peaks = "/home/wanghm/wanghm_R/B_BinQuasiPeaks_rmChr.bed",
    S_peaks = "/home/wanghm/wanghm_R/S_BinQuasiPeaks_rmChr.bed"),
  annotatePeak,
  TxDb = txdb_noCDS,
  tssRegion = c(-2000, 2000),
  verbose = F
)
# rm(geneBodyRegion)
geneInfo <- genes(txdb_noCDS)
median(width(ranges(geneInfo)))
??chipseeker
install.packages("clusterProfiler")
biocLite("clusterProfiler")

