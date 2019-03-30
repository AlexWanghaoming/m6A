library(Biostrings)
library(seqLogo)
# seqlogo plots can not be combination
## use biostrings & seqLogo
load("/home/wanghm/wanghm_R/methylation/methy1.RData")
comparedSeqLogo <- function(methylationFile,genomeFile){
  library(seqLogo)
  library(ggseqlogo)
  ph1_m4C <- readDNAStringSet(filepath = methylationFile, format = "fasta")
  matrixForSeqLogo1 <- consensusMatrix(DNAStringSet(ph1_m4C,start=19,end=25),as.prob = T)[1:4,]
  ph1_genome <- readDNAStringSet(filepath = genomeFile,format = "fasta")
  merged_genome <- Reduce(c,ph1_genome)
  samp <- sample(matchPattern("NNCNNNN",fixed = F,subject = merged_genome),length(ph1_m4C)) # fixed = False can match N
  matrixForSeqLogo2 <- consensusMatrix(samp, as.prob = T)[1:4,]
  p1 <- ggseqlogo(matrixForSeqLogo1)
  p2 <- ggseqlogo(matrixForSeqLogo2)
  gridExtra::grid.arrange(p1, p2)
  # library(seqLogo)
  # ss <- seqLogo(matrixForSeqLogo2,ic.scale = T)
}
comparedSeqLogo(methylationFile = "/home/wanghm/whm/pacbio_data/Pacbio/PH-1/6mA/official_methylationResult/m4C_Seq.fasta",genomeFile = "/home/wanghm/whm/1FG/35genome/Fusarium_graminearum.RR1_with_mito.fa")


# ________________________ calculate 6mA frequency!
# FG_genome <-
#   readDNAStringSet(filepath = "/home/wanghm/whm/pacbio_data/Pacbio/PH-1/genome_analysis/polished_genome/fg_polished.fasta", format = "fasta")
#
# m6AFrenquency <- 86508 / sum(letterFrequency(FG_genome, "A"))
#
# m5CFrenquenct <- 118261 / sum(letterFrequency(FG_genome, "C"))

ydj <-
  readDNAStringSet("/home/wanghm/whm/pacbio_data/Pacbio/YDJ-3/6mA/ydj.filted.fasta",
                   'fasta')
scaffold_names <- substr(names(ydj), 1, 11)
library(data.table)
ydj_m6A_gff <-
  fread("/home/wanghm/whm/pacbio_data/Pacbio/YDJ-3/6mA/m6A.gff",
        skip = 1)
ydj_m4C_gff <-
  fread("/home/wanghm/whm/pacbio_data/Pacbio/YDJ-3/6mA/m4C.gff",
        skip = 1)
# nrow(m6A_gff[which(m6A_gff$V1 == 1],) # filter dataframe
nrow(ydj_m6A_gff)
m6A_count = ydj_m6A_gff[, .N, by = .(V1)] # using data.table package for counting
m4C_count = ydj_m4C_gff[, .N, by = .(V1)]

# cbind(m6A_count, scaffold_names) # add a column to a dataframe
colnames(m6A_count) <- c("scaffold_names", "N")
aligned_m6A_count = merge(as.data.frame(scaffold_names),
                          m6A_count,
                          by = "scaffold_names",
                          all.x = TRUE)
colnames(m4C_count) <- c("scaffold_names", "N")
aligned_m4C_count <-
  merge(as.data.frame(scaffold_names),
        m4C_count,
        by = "scaffold_names",
        all.x = TRUE)

A_count = letterFrequency(ydj, "A")
C_count = letterFrequency(ydj, "C")
aFreq <- aligned_m6A_count$N / A_count
cFreq <- aligned_m4C_count$N / C_count
freqDF <- data.frame(aFreq, cFreq, row.names = scaffold_names)
colnames(freqDF) <- c("6mA", "m4C")
freqDF2 = t(freqDF) # inverse matrix

# plot m6A frequency in numbers of scaffolds
freqDF[is.na(freqDF)] <- 0 # sub NA in dataFrame
library(ggplot2)
p <- ggplot(freqDF) # discrete Var must be factor()
p + geom_line(aes(
  x = as.factor(rownames(freqDF)),
  y = freqDF$`6mA`,
  group = 1,
  colour = "red"
)) + geom_line(aes(
  x = as.factor(rownames(freqDF)),
  y = freqDF$m4C,
  group = 1,
  color = "yellow"
)) + xlab("68Scaffolds") + ylab("MethylationFrequency") + 
scale_colour_hue("MethylationType",labels=c("6mA", "m4C")) + theme(axis.text.x = element_blank())


# plot m6A frequency in each 4 chromosomers in FG
barplot(freqDF2,
        beside = TRUE,
        col = rainbow(2))
# legend = )
legend(x = "topleft",
       legend = rownames(freqDF2),
       fill = rainbow(2))  # if we want set lengend details, do not use "legend = " option in barplot()
title(main = "6mA and m4C Frequency in each chromosomes", font.main = 4)

########################################################################################################
## Source: biostar
## usage: plot methylation around TSS employing RRBS data
if (FALSE)
{
  #if we do not have reduced representation bisulphite sequencing (RRBS) data, we can not use method below to plot around TSS
  math = read.table(
    "/home/wanghm/wanghm_R/methylation/forR_dataframe",
    sep = " ",
    header = FALSE
  )
  m_dens <- math[which(math$V3 < 500 & math$V3 > -500), 2:3]
  # order() return a position index that the element in the vector
  m_dens <- m_dens[order(m_dens[, 2]),]
  results_dist <-
    apply(m_dens, 2, function(x)
      ave(x, m_dens[, 2], FUN = mean))
  results_dist2 <- results_dist[!(duplicated(results_dist[, 2])),]
  
  plot(
    x = results_dist2[, 2],
    y = results_dist2[, 1],
    type = "n",
    ylim = c(0, 100),
    xlab = "Dist to TSS",
    ylab = "Average Methylation"
  )
  fit <-
    lm(results_dist2[, 1] ~ poly(results_dist2[, 2], 6, raw = T))
  lines(results_dist2[, 2],
        predict(fit, data.frame(x = results_dist2[, 1])),
        lwd = 5,
        col = "black")
}
save.image() # save current workplace
##### ____________________________________________________













