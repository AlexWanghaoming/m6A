library(ggseqlogo)
library(Biostrings)
Args <- commandArgs(trailingOnly = TRUE)
gff <- Args[1]
shell.cmd1 <- paste0("awk '$3~/m6A/' ",as.character(gff), "| awk -F\";|=\" 'BEGIN{i=1}{print \">\"i;print$4;i++}' > m6A_context.fasta")
shell.cmd2 <- paste0("awk '$3~/m4C/' ",gff, "| awk -F\";|=\" 'BEGIN{i=1}{print \">\"i;print$4;i++}' > m4C_context.fasta")
print(shell.cmd2)
system(command = shell.cmd1)
system(command = shell.cmd2)

m6A.context <- readDNAStringSet("m6A_context.fasta", format = "fasta")
m4C.context <- readDNAStringSet("m4C_context.fasta", format = "fasta")

m6A_seqLogo.mat <- consensusMatrix(m6A.context, as.prob = T)[1:4,16:26]
p1 <- ggseqlogo(data = m6A_seqLogo.mat, method="probability",col_scheme = "clustalx")
m4C_seqLogo.mat <- consensusMatrix(m4C.context, as.prob = T)[1:4,16:26]
p2 <- ggseqlogo(data = m4C_seqLogo.mat, method="probability",col_scheme = "clustalx")
library(cowplot)
res <- plot_grid(p1, p2, labels = c("A", "B"), align = "h",nrow = 2)
save_plot("seqlogo.pdf",plot = res)
