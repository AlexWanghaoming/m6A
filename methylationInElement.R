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

## ********************************************************************* ##

pergene_cds <-
  function(x) {
    # calculate methylation Rre in each cds of genes in FG
    allcdsSeq <-
      Rsamtools::getSeq(FaFile("Fusarium_graminearum.RR1_with_mito.fa"), x)
    methylationFraction <-
      sum(countOverlaps(x, m4CGrange)) / sum(letterFrequency(allcdsSeq, letters = c("A", "T")))
    return(methylationFraction)
  }
aa <- cdsBy(longestTranTxdb, "gene")
result <- sapply(aa, pergene_cds)
sort(result, decreasing = T)[1:10]   # get top 10 genes whose CDS regions has most methyFra
length(result[which(result > 0.05)])

## ********************************************************************* ##
## calculate 6mA Fra in each features: exon, promoter, intron, UTR
methyFrac_element <-
  function(longestTranTxdb,
           genomeFile = "Fusarium_graminearum.RR1_with_mito.fa",
           gffFile = "/home/wanghm/whm/1FG/35genome/Fusarium_graminearum.RR1.35.chr.gff3",
           methylationGff = "/home/wanghm/whm/pacbio_data/Pacbio/PH-1/6mA/official_methylationResult/m6A.gff",
           element) {
    library(Rsamtools)
    library(ChIPseeker)
    library(Biostrings)
    library(GenomicFeatures)
    m4C_gff <- data.table::fread(methylationGff, skip = 1)
    m4CGrange <- GRanges(
      seqnames = m4C_gff$V1,
      ranges = IRanges(start = m4C_gff$V4, end = m4C_gff$V5),
      strand = "*"
    )
    element <-
      match.arg(element,
                c("promoter", "exon", "intron", "5utr", "3utr", "cds"))
    if (element == "promoter") {
      ph1Pmt <- promoters(longestTranTxdb,
                          upstream = 200,
                          downstream = 0)
      ph1Pmt <- ph1Pmt[start(ph1Pmt) > 0]
      pmtSeq <- Rsamtools::getSeq(FaFile(genomeFile), ph1Pmt)
      C_Inpromoters <-
        letterFrequency(pmtSeq, letters = c("A", "T")) # For 14336 Genes, calculate pmt
      # methylation Frequency per Gene.
      methylationFraction <-
        sum(countOverlaps(ph1Pmt, m4CGrange)) / sum(C_Inpromoters)
    }
    if (element == "exon") {
      all_exon <- exons(longestTranTxdb)
      allExonSeq <- Rsamtools::getSeq(FaFile(genomeFile), all_exon)
      methylationFraction <-
        sum(countOverlaps(all_exon, m4CGrange)) / sum(letterFrequency(allExonSeq, letters = c("A", "T")))
    }
    if (element == "cds") {
      all_cds <- cds(longestTranTxdb)
      allCdsSeq <- Rsamtools::getSeq(FaFile(genomeFile), all_cds)
      methylationFraction <-
        sum(countOverlaps(all_cds, m4CGrange)) / sum(letterFrequency(allCdsSeq, letters = c("A", "T")))
    }
    if (element == "intron") {
      perTxIntron <- intronsByTranscript(longestTranTxdb)
      perTxIntronDF <-
        as.data.frame(perTxIntron)   # transform list -> dataframe -> GRange Obj
      totalIntron <-
        GRanges(
          seqnames = perTxIntronDF$seqnames,
          ranges = IRanges(start = perTxIntronDF$start, end = perTxIntronDF$end),
          strand = perTxIntronDF$strand
        )
      allIntronSeq <- Rsamtools::getSeq(FaFile(genomeFile), totalIntron)
      methylationFraction <-
        sum(countOverlaps(totalIntron, m4CGrange)) / sum(letterFrequency(allIntronSeq, letters = c("A", "T")))
    }
    if (element == "5utr") {
      # perTx5Utr <- fiveUTRsByTranscript(longestTranTxdb)
      # perTx5UtrDF <- as.data.frame(perTx5Utr)
      # total5Utr <- GRanges(seqnames = perTx5UtrDF$seqnames,ranges = IRanges(start = perTx5UtrDF$start,end = perTx5UtrDF$end),strand = perTx5UtrDF$strand)
      a <- cdsBy(longestTranTxdb, "tx")
      b <- sapply(a, function(x)
        x[1])
      d <- Reduce(c, b)
      total5Utr <- flank(d,
                         both = F,
                         width = 200,
                         start = T)
      all5Utr <- Rsamtools::getSeq(FaFile(genomeFile), total5Utr)
      methylationFraction <-
        sum(countOverlaps(total5Utr, m4CGrange)) / sum(letterFrequency(all5Utr, letters = c("A", "T")))
    }
    if (element == "3utr") {
      # perTx3Utr<- threeUTRsByTranscript(longestTranTxdb)
      # perTx3UtrDF <- as.data.frame(perTx3Utr)
      # total3Utr <- GRanges(seqnames = perTx3UtrDF$seqnames,ranges = IRanges(start = perTx3UtrDF$start,end = perTx3UtrDF$end),strand = perTx3UtrDF$strand)
      a = cdsBy(longestTranTxdb, "tx")
      b = sapply(a, function(x)
        tail(x, 1))
      d = Reduce(c, b)
      total3Utr = flank(d,
                        both = F,
                        width = 200,
                        start = F)
      all3Utr <-
        Rsamtools::getSeq(FaFile("Fusarium_graminearum.RR1_with_mito.fa"),
                          total3Utr)
      methylationFraction <-
        sum(countOverlaps(total3Utr, m4CGrange)) / sum(letterFrequency(all3Utr, letters = c("A", "T")))
    }
    return(methylationFraction)
  }
methyFrac_element(longestTranTxdb = longestTranTxdb, element = "5utr")

## calculating exon perGene
# exonByGeneList <- exonsBy(longestTranTxdb,by = "tx")
# func1 <- function(x){
#   methy_m4C_perGene <- countOverlaps(x[1],m4CGrange)
#   return(methy_m4C_perGene)
# }
# exon_m4C_perGene <- sapply(exonByGeneList, func1)
# length(exon_m4C_perGene)

## ********************************************************************* ##

m4CFreInelement <-
  data.frame(
    MethyFrequency = c(
      total_methylationFreInPmt,
      total_methylationFreInExon,
      total_methylationFreInIntron,
      total_methylationFreIn5Utr,
      total_methylationFreIn3Utr
    ),
    row.names = c("Promoter", "Exon", "Intron", "5Utr", "3Utr")
  )
m4CFre <- t(m4CFreInelement)
barplot(m4CFre, ylab = "m4C Frequency")
# axis(2,at=pretty(m4CFre),labels = pretty(m4CFre))
barplot(m4CFreInelement[1, ])

# ____________________________________________________________________

## High methylation Fraction
library(data.table)

m4CHighFraction_gff <-
  fread(
    "/home/wanghm/whm/pacbio_data/Pacbio/PH-1/6mA/official_methylationResult/m4C_withHignmethyFra.gff",
    skip = 1
  )
m4CHighFraction_GRange = GRanges(
  seqnames = m4CHighFraction_gff$V1,
  ranges = IRanges(start = m4CHighFraction_gff$V4, end = m4CHighFraction_gff$V5),
  strand = "*"
)
## get Promoters and calculate m4C Frequency in region
ph1Pmt <- promoters(longestTranTxdb, upstream = 500, downstream = 50)
ph1Pmt <- ph1Pmt[start(ph1Pmt) > 0]
library(Rsamtools)
pmtSeq <-
  Rsamtools::getSeq(FaFile("Fusarium_graminearum.RR1_with_mito.fa"), ph1Pmt)
C_Inpromoters <-
  letterFrequency(pmtSeq, letters = "C") # For 14336 Genes, calculate pmt
# methylation Frequency per Gene.
total_methylationFreInPmt <-
  sum(countOverlaps(ph1Pmt, m4CHighFraction_GRange)) / sum(countOverlaps(ph1Pmt, m4CGrange))
## get Exons
all_exon <- exons(longestTranTxdb)
allExonSeq <-
  Rsamtools::getSeq(FaFile("Fusarium_graminearum.RR1_with_mito.fa"), all_exon)
total_methylationFreInExon <-
  sum(countOverlaps(all_exon, m4CHighFraction_GRange)) / sum(countOverlaps(all_exon, m4CGrange))
