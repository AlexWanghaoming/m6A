library(GenomicFeatures)
library(tidyverse)
txdb <- makeTxDbFromGFF(file = "/home/wanghm/whm/1FG/PH1-fusion-version2_represent_transcript_rm_Mt.gff") 
m6A_gff <- data.table::fread("/home/wanghm/whm/pacbio_data/PH-1/6mA/xulab/m6A.gff")
m6AGrange <- GRanges(seqnames = m6A_gff$V1,ranges = IRanges(start = m6A_gff$V4, end = m6A_gff$V5),strand = m6A_gff$V7)
m6AGrange_total <- GRanges(seqnames = m6A_gff$V1,ranges = IRanges(start = m6A_gff$V4, end = m6A_gff$V5),strand = "*") 

# functions
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
grouped_calculation <- function(ungrouped.gr, label){
  library(dplyr)
  nam <- as.factor(names(ungrouped.gr))
  names(ungrouped.gr) <- NULL
  dd1 <- get_feature_methylation_ratio(ungrouped.gr)
  dd2 <- add_column(dd1, nam)
  df <- summarise_at(group_by(dd2, dd2$nam), c(1,2,3),mean)
  colnames(df) <- c("keeped_geneID", paste(label,"total", sep = "_"),
                    paste(label, "sameDir", sep = "_"),
                    paste(label, "diffDir", sep = "_"))
  return(df)
} # return grouped methylation ratio dataframe


######
####
### calculating promoter, exon, intron, 5UTR, 3UTR, feature methylaton ratio Separately
pmt.gr <- promoters(txdb, upstream = 500, downstream = 300)
pmt.gr <- pmt.gr[sort(names(pmt.gr))]
pmt.df <- get_feature_methylation_ratio(pmt.gr) %>% dplyr::rename(., promoter_total=total, 
                                                    promoter_sameDir=sameDir, 
                                                    promoter_diffDir=diffDir)
pmt.df <- pmt.df[rownames(pmt.df) %in% paste0(target_gene.list, "T0"),] %>% rownames_to_column(., var = "keeped_geneID") # 15598


exon.gr <- unlist(exonsBy(txdb, "gene"), use.names = T) %>% unique()
exon.df <- grouped_calculation(exon.gr, label = "exon") %>% filter(., keeped_geneID %in% target_gene.list) # 15598
exon.df$keeped_geneID <- paste0(exon.df$keeped_geneID,"T0")

intron.gr <- unlist(intronsByTranscript(txdb, use.names = T), use.names = T) %>% unique()
intron.df <- grouped_calculation(intron.gr, label = "intron") %>% filter(., keeped_geneID %in% paste0(target_gene.list,"T0")) # 15598

fiveUTR.gr <- unlist(fiveUTRsByTranscript(txdb, use.names =T),use.names = T) %>% unique()
ss <- fiveUTR.gr[!(duplicated(names(fiveUTR.gr)) | duplicated(names(fiveUTR.gr), fromLast=T))]
fiveUTR.gr <- ss[sort(names(ss))]
fiveUTR.df <- get_feature_methylation_ratio(fiveUTR.gr) %>% dplyr::rename(., fiveUTR_total=total, 
                                                                          fiveUTR_sameDir=sameDir, 
                                                                          fiveUTR_diffDir=diffDir)
fiveUTR.df <- fiveUTR.df[rownames(fiveUTR.df) %in% paste0(target_gene.list,"T0"),] %>% rownames_to_column(., var = "keeped_geneID")


threeUTR.gr <- unlist(threeUTRsByTranscript(txdb, use.names = T), use.names = T) %>% unique()
tt <- threeUTR.gr[!(duplicated(names(threeUTR.gr)) | duplicated(names(threeUTR.gr), fromLast=T))]
threeUTR.gr <- tt[sort(names(tt))]
threeUTR.df <- get_feature_methylation_ratio(threeUTR.gr) %>% dplyr::rename(., threeUTR_total=total, 
                                                                            threeUTR_sameDir=sameDir, 
                                                                            threeUTR_diffDir=diffDir)
threeUTR.df <- threeUTR.df[rownames(threeUTR.df) %in% paste0(target_gene.list,"T0"),] %>% rownames_to_column(., var = "keeped_geneID")

#### merge all features(promoters, introns, exons, CDS, 5UTR, 3UTR) strand_specific 
genic.df <- get_feature_methylation_ratio(genes(txdb)) %>% 
  rownames_to_column(., var = "keeped_geneID") %>%
  filter(., keeped_geneID %in% target_gene.list)
genic.df$keeped_geneID <- paste0(genic.df$keeped_geneID, "T0")

## merge all 
res.df <- merge(genic.df,exon.df, by="keeped_geneID", all=T)%>%
            merge(.,pmt.df, by="keeped_geneID",all=T) %>% 
            merge(., intron.df, by="keeped_geneID",all=T) %>% 
            merge(., fiveUTR.df, by="keeped_geneID",all=T) %>%
            merge(., threeUTR.df, by="keeped_geneID",all=T)
res.df$keeped_geneID <- gsub("T0","", res.df$keeped_geneID)
res.df[,2:ncol(res.df)] <- round(res.df[,2:ncol(res.df)],4)
# output table
write.table(res.df, file = "1.res/m6A_ratio.txt" ,sep = "\t", quote = F, row.names = F)


g1 <- d1[d1$key=="exon_sameDir",][,2]
g2 <- d1[d1$key=="exon_diffDir",][,2]
wilcox.test(g1,g2)
### plot boxplot 
d1 <- gather(res.df[,c(6,7)]) %>% add_column(feature = "exon")
d2 <- gather(res.df[,c(9,10)]) %>% add_column(feature = "promoter")
d3 <- gather(res.df[,c(12,13)]) %>% add_column(feature = "intron")
d4 <- gather(res.df[,c(15,16)]) %>% add_column(feature = "fiveUTR")
d5 <- gather(res.df[,c(18,19)]) %>% add_column(feature = "threeUTR")

### make a multi-variable dataframe
feature_methyFra.df <- rbind(d1,d2,d3,d4,d5)
feature_methyFra.df$key <- gsub(".*_","",feature_methyFra.df$key)
## rm NA and NaN in dataframe
boxplot.df <- feature_methyFra.df[!is.na(feature_methyFra.df$value),]
boxplot.df <- feature_methyFra.df[!is.nan(feature_methyFra.df$value),]
boxplot.df <- boxplot.df[is.finite(boxplot.df$value),]

library(ggstatsplot)
ylim(c(0,0.05))
#### it is difficult for base ggplot to add significant labels
# p <- ggplot(data = boxplot.df) +
#   geom_boxplot(aes(x=feature, y=value, fill=key))+
#   stat_summary(aes(x=feature, y=value, group=key),
#                fun.y = mean, geom="point", shape=20, size=3, color="red", fill="red", 
#                position = position_dodge(width = 0.75), na.rm = T)

p1 <- ggboxplot(boxplot.df, x="feature", y="value", color = "key", palette = "jco", ylab = "6mA ratio") + 
  stat_compare_means(aes(group=key), label = "p.signif") + ylim(c(0,0.04))+
  stat_summary(aes(x=feature, y=value, group=key),fun.y = mean, geom="point", shape=20, size=3, color="red", fill="red", 
                                                  position = position_dodge(width = 0.75), na.rm = T)
ggpar(p1, legend = "right", legend.title = "")

