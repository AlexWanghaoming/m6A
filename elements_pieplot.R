library(GenomicFeatures)
txdb <- makeTxDbFromGFF(file = "/Users/alexwang/data/pacbio/0.fg/genome/Fusarium_graminearum.RR1.35.chr.gff3") 
m6A_gff <- data.table::fread("/Users/alexwang/data/pacbio/0.fg/genome/m6A.gff", skip = 1)
## m6A_gff 需要去除mt 和 contigs的信息
m6A_gff <- m6A_gff[which(m6A_gff$V1 %in% c(1,2,3,4)),]
m6AGrange <- GRanges(seqnames = m6A_gff$V1, ranges = IRanges(start = m6A_gff$V4, end = m6A_gff$V5), strand = "*") # 将gff文件读取为GRanges格式
gene <- genes(txdb)
pmt <- promoters(txdb, upstream = 500, downstream = 500) 
intron <- intronicParts(txdb,linked.to.single.gene.only = F)  # 提取内含子和外显子位置，若有重叠基因包含同一外显子，计算多次
exon <- exonicParts(txdb, linked.to.single.gene.only = F)

# 获取6mA甲基化位点在各个elements上的个数，此处不可用sum(countOverlaps(m6AGrange, pmt))
methylsites.pmt <- subsetByOverlaps(m6AGrange, pmt)
## exon
methylsites.exon <- subsetByOverlaps(m6AGrange, exon)
exon_pmt.intersect <- intersect(methylsites.exon, methylsites.pmt)
methyl.onlyExon <- setdiff(methylsites.exon, exon_pmt.intersect)  # 从在外显子上的6mA位点中过滤掉与启动子区overlap的部分
## intron
methylsites.intron <- subsetByOverlaps(m6AGrange, intron)
intron_pmt.intersect <- intersect(methylsites.intron, methylsites.pmt)
methyl.onlyIntron <- setdiff(methylsites.intron, intron_pmt.intersect)
## intergenic
methyl.intergenic <- setdiff(m6AGrange, subsetByOverlaps(m6AGrange, gene))
intergenic_pmt.intersect <- intersect(methyl.intergenic, methylsites.pmt)
methyl.onlyIntergenic <- setdiff(methyl.intergenic, intergenic_pmt.intersect)
## 画饼图
pie.data <- data.frame(group = c("Promoter", "Exon", "Intron", "Intergenic"), counts = c(length(methylsites.pmt), 
                                                                                 length(methyl.onlyExon), 
                                                                                 length(methyl.onlyIntron), 
                                                                                 length(methyl.onlyIntergenic)))
library(ggpubr)
# ?ggpar
aa <- pie.data$counts/sum(pie.data$counts)  # 计算百分比
labs <- paste0(round(aa,2)*100,"%")
p <- ggpubr::ggpie(data = pie.data, x = "counts", label = labs[c(1,3,4,2)], fill = "group", lab.pos = "in", color = "white", 
              palette = RColorBrewer::brewer.pal(4,"Set3"), legend.title="")
p1 <- ggpar(p,legend = "right",tickslab = F)
# - - - - - - - - - - - - - - - - - - - - - 
# 为GRanges对象加一列标签，用于后续判断其分类
methyl.intergenic$type = "intergenic"
methylsites.pmt$type = "promoter"
methylsites.exon$type = "exon"
methylsites.intron$type = "intron"
# methyl.intergenic$type = NULL
# methylsites.pmt$type = NULL
# methylsites.exon$type = NULL
# methylsites.intron$type = NULL
library(tidyverse)
aa <- as_tibble(methyl.intergenic)
bb <- as_tibble(methylsites.pmt)
cc <- as_tibble(methylsites.exon)
dd <- as_tibble(methylsites.intron)
f1 <- function(x,y){
  merge(x,y, by.x=c("seqnames","start"), by.y=c("seqnames","start"), all=TRUE)
}
df <- Reduce(f1, list(aa,bb,cc,dd))   ## 将四个tibble数据库昂按seqnames和start 合并
df <- df[,c(6,10,14,18)]
# 创建一个elements.list空列表，迭代df的每一行，利用列表索引将其加入空列表中
elements.list=list()
for (i in 1:nrow(df)) {
  xx <- as.character(df[i,])
  elements.list[[i]] = xx[xx!=0]
}

tidy_elements <- tibble(sitename=as.character(1:length(elements.list)), genres=elements.list)
library(ggupset)
p2 <- tidy_elements %>%    # 画upsetplot
  ggplot(aes(x=genres)) +
  geom_bar(fill= "#085F5F", color="black") + ylab("Methylation Site Counts")+
  theme(panel.grid.major = element_blank())+
  scale_x_upset(order_by = "degree", name="") +
  theme_combmatrix(combmatrix.panel.point.color.fill = "#E15505",
                   combmatrix.panel.line.size = 0,
                   combmatrix.label.make_space = FALSE)
library(cowplot)
res <- plot_grid(p1, p2, labels = c("a", "b"), align = "v",ncol = 2)
save_plot("pie&upset_plot.pdf",plot = res,base_height = 10)
# RColorBrewer::display.brewer.all()
