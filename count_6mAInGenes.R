library(GenomicFeatures)
txdb <- makeTxDbFromGFF(file = "/Users/alexwang/data/pacbio/0.fg/genome/Fusarium_graminearum.RR1.35.chr.gff3") 
m6A_gff <- data.table::fread("/Users/alexwang/data/pacbio/0.fg/genome/m6A.gff", skip = 1)
## m6A_gff 需要去除mt 和 contigs的信息
m6A_gff <- m6A_gff[which(m6A_gff$V1 %in% c(1,2,3,4)),]
m6AGrange <- GRanges(seqnames = m6A_gff$V1, ranges = IRanges(start = m6A_gff$V4, end = m6A_gff$V5), strand = "*") # 将gff文件读取为GRanges格式
gene <- genes(txdb)
counts <- countOverlaps(gene, m6AGrange)

a <- counts[counts==1] %>%length()
b <- counts[counts==2] %>%length()
c <- counts[counts==3] %>%length()
d <- counts[counts>3] %>%length()
counts.df <- tibble::tibble(type=c("1","2","3",">3"), count=c(a,b,c,d))
library(ggpubr)
ggbarplot(counts.df, x = "type", y = 'count', combine = T, fill = "#8E2C31",
          xlab = "Number of Methylation Sites", ylab = "Number of Genes", label = T)
