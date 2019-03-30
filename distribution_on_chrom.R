library("Biostrings")
fo_genome <- readDNAStringSet("/Users/alexwang/data/pacbio/fo/Fot1_assembly.fasta",'fasta')
scaffold_names <- names(fo_genome)
library(data.table)
fo_m6A_gff <- fread("/Users/alexwang/data/pacbio/fo/fo_m6A.gff")
fo_m4C_gff <- fread("/Users/alexwang/data/pacbio/fo/fo_m4C.gff")
m6A_count = fo_m6A_gff[, .N, by = .(V1)] # data.table 分组统计
m4C_count = fo_m4C_gff[, .N, by = .(V1)]
colnames(m6A_count) <- c("scaffold_names", "N")
colnames(m4C_count) <- c("scaffold_names", "N")
df <- merge(m6A_count,m4C_count,by="scaffold_names", all=TRUE)# 合并两个dataframe 没有甲基化位点的染色体用NA补上
 # 将缺失值赋0
AT_count = letterFrequency(fo_genome, "AT");rownames(AT_count) <- scaffold_names
CG_count = letterFrequency(fo_genome, "CG");rownames(CG_count) <- scaffold_names
df_sort <- df[rownames(AT_count),]
df_sort <- as.data.frame(df_sort)
df_sort[is.na(df_sort)] <- 0  # 将scaffolds name 排序
at.ratio <- df_sort$N.x / AT_count # 计算6mA， m4C的比率
cg.ratio <- df_sort$N.y / CG_count
ratio.df <- data.frame(at.ratio, cg.ratio, row.names = scaffold_names)
colnames(ratio.df) <- c("m6A", "m4C")
ratio.df.t <- t(ratio.df)
library(prettyB)   # prettyB包美化base plot
bar <- prettyB::barplot.default(ratio.df.t,beside = T,col =c("orange","green"),xaxt = "n", legend=rownames(ratio.df.t))
# legend(x = "topleft",
#        legend = rownames(freqDF2),
#        fill = rainbow(2))  # if we want set lengend details, do not use "legend = " option in barplot()
## 自定义x轴标签
labs <- scaffold_names
text(cex=1, x=colMeans(bar)-0.25, y=-0.0017, labels = labs, xpd=TRUE, srt=65)  # 当barplot reside=TRUE时，x=colMeans(bar)
## 也可以用ggplot2实现
library(ggpubr)
library(reshape2)
ratio.df$name <- rownames(ratio.df)
ratio.df.ggplot <- melt(ratio.df)  # 宽表变长表
## ggbarplot更优秀，x.text.angle参数可以替代 theme(axis.text.x = element_text(angle = 60, hjust = 0.5, vjust = 0.5))
ggbarplot(ratio.df.ggplot, x = "name", y="value",fill = "variable",color="white",position = position_dodge2(), 
          palette = "jco",x.text.angle=60)


