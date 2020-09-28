NumMCLallRNASeq <- iterations14[match(rownames(allRNASeq), iterations14$name), ]$X__mclCluster
NumMCLallRNASeq <- data.frame(allRNASeq,
                              NumMCL = NumMCLallRNASeq)
NumMCLallRNASeq <- na.omit(NumMCLallRNASeq)
write.table(NumMCLallRNASeq, file = "~/Nakano_RNAseq/network_analysis/base/NumMCLallRNASeq.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
