library(Biostrings)
library(seqLogo)
# seqlogo plots can not be combination
## use biostrings & seqLogo
m6A_context <- readDNAStringSet(filepath = "/home/wanghm/whm/pacbio_data/PH-1/6mA/xulab/m6A_context.fasta", format = "fasta")
matrixForSeqLogo1 <- consensusMatrix(DNAStringSet(m6A_context, start = 2, end = 10), as.prob = T)[1:4,]
seqLogo(matrixForSeqLogo1, ic.scale = T)
comparedSeqLogo <- function(methylationFile){
  library(seqLogo)
  library(ggseqlogo)
  library(Biostrings)
  ph1_m6A <- readDNAStringSet(filepath = methylationFile, format = "fasta")
  matrixForSeqLogo1 <- consensusMatrix(DNAStringSet(ph1_m6A,start = 2, end = 10),as.prob = T)[1:4,]   # cut -2 +4 motifs 
  ph1_genome <- readDNAStringSet(filepath = "/home/wanghm/whm/1FG/xulab_ph1.fasta",format = "fasta")
  merged_genome <- Reduce(c,ph1_genome)
  samp <- sample(matchPattern("NNNNNANNNNN",fixed = F,subject = merged_genome),length(ph1_m6A)) # fixed = False can match N
  matrixForSeqLogo2 <- consensusMatrix(samp, as.prob = T)[1:4,]
  p1 <- ggseqlogo(matrixForSeqLogo1)
  p2 <- ggseqlogo(matrixForSeqLogo2)
  gridExtra::grid.arrange(p1, p2)
  # library(seqLogo)
  # ss <- seqLogo(matrixForSeqLogo2,ic.scale = T)
}
comparedSeqLogo(methylationFile = "/home/wanghm/whm/pacbio_data/PH-1/6mA/xulab/m6A_context.fasta")

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

