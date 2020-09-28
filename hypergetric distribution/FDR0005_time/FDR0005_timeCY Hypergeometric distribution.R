NumMCL <- na.omit(read.table("~/Nakano_RNAseq/network_analysis/cytoscape/NumMCL.txt", sep = "\t", header = T))
rownames(NumMCL) <- NumMCL$name
####CY15_1h enrichment####
CY15_1h_DEGs_FDR0005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY15_1h_FDR0005.txt", sep = "\t", row.names = 1, header = T)
#全クラスターのCY15_1hのDEGs数
NumCY15_1h_DEGs_FDR0005 <- sum(!is.na(match(rownames(CY15_1h_DEGs_FDR0005), NumMCL$name)))
#全クラスターのCY15_1h以外のDEGs数
Numothergenes <- c(length(NumMCL$name)-NumCY15_1h_DEGs_FDR0005)
maxNumMCL <- max(NumMCL$X__mclCluster)
n <- 1
CY15_1h_enrichment_pvalue <- c()
check <- c()
CY15_1h_check_MCL <- list()

for(n in n:maxNumMCL){
  test <- NumMCL[NumMCL$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(rownames(CY15_1h_DEGs_FDR0005), rownames(test)))
  Numsubcluster <- nrow(test)
  CY15_1h_enrichment_pvalue <- c(CY15_1h_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumCY15_1h_DEGs_FDR0005, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(CY15_1h_enrichment_pvalue)[n] <- n
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(rownames(CY15_1h_DEGs_FDR0005), rownames(test))
    names(T_INT) <- n
  }
  CY15_1h_check_MCL <- c(CY15_1h_check_MCL, list(T_INT))
  print(maxNumMCL - n)
  n <- n + 1
}
CY15_1h_enrichment_qvalue <- p.adjust(CY15_1h_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
CY15_1h_enriched_Node <- CY15_1h_check_MCL[CY15_1h_enrichment_qvalue < 0.05]
NumCY15_1h_subcluster <- c()
i <- 1
for(i in i:length(CY15_1h_enriched_Node)){
  NumCY15_1h_subcluster <- c(NumCY15_1h_subcluster, names(CY15_1h_enriched_Node[[i]])[1])
  i <- i+1
}
NumCY15_1h_subcluster <- as.numeric(NumCY15_1h_subcluster)


#CY15_3h enrichment#
CY15_3h_DEGs_FDR0005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY15_3h_FDR0005.txt", sep = "\t", row.names = 1, header = T)
#全クラスターのCY15_3hのDEGs数
NumCY15_3h_DEGs_FDR0005 <- sum(!is.na(match(rownames(CY15_3h_DEGs_FDR0005), NumMCL$name)))
#全クラスターのCY15_3h以外のDEGs数
Numothergenes <- c(length(NumMCL$name)-NumCY15_3h_DEGs_FDR0005)
maxNumMCL <- max(NumMCL$X__mclCluster)
n <- 1
CY15_3h_enrichment_pvalue <- c()
check <- c()
CY15_3h_check_MCL <- list()

for(n in n:maxNumMCL){
  test <- NumMCL[NumMCL$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(rownames(CY15_3h_DEGs_FDR0005), rownames(test)))
  Numsubcluster <- nrow(test)
  CY15_3h_enrichment_pvalue <- c(CY15_3h_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumCY15_3h_DEGs_FDR0005, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(CY15_3h_enrichment_pvalue)[n] <- n
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(rownames(CY15_3h_DEGs_FDR0005), rownames(test))
    names(T_INT) <- n
  }
  CY15_3h_check_MCL <- c(CY15_3h_check_MCL, list(T_INT))
  print(maxNumMCL - n)
  n <- n + 1
}
CY15_3h_enrichment_qvalue <- p.adjust(CY15_3h_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
CY15_3h_enriched_Node <- CY15_3h_check_MCL[CY15_3h_enrichment_qvalue < 0.05]
NumCY15_3h_subcluster <- c()
i <- 1
for(i in i:length(CY15_3h_enriched_Node)){
  NumCY15_3h_subcluster <- c(NumCY15_3h_subcluster, names(CY15_3h_enriched_Node[[i]])[1])
  i <- i+1
}
NumCY15_3h_subcluster <- as.numeric(NumCY15_3h_subcluster)



#CY15_12h enrichment#
CY15_12h_DEGs_FDR0005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY15_12h_FDR0005.txt", sep = "\t", row.names = 1, header = T)
#全クラスターのCY15_12hのDEGs数
NumCY15_12h_DEGs_FDR0005 <- sum(!is.na(match(rownames(CY15_12h_DEGs_FDR0005), NumMCL$name)))
#全クラスターのCY15_12h以外のDEGs数
Numothergenes <- c(length(NumMCL$name)-NumCY15_12h_DEGs_FDR0005)
maxNumMCL <- max(NumMCL$X__mclCluster)
n <- 1
CY15_12h_enrichment_pvalue <- c()
check <- c()
CY15_12h_check_MCL <- list()

for(n in n:maxNumMCL){
  test <- NumMCL[NumMCL$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(rownames(CY15_12h_DEGs_FDR0005), rownames(test)))
  Numsubcluster <- nrow(test)
  CY15_12h_enrichment_pvalue <- c(CY15_12h_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumCY15_12h_DEGs_FDR0005, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(CY15_12h_enrichment_pvalue)[n] <- n
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(rownames(CY15_12h_DEGs_FDR0005), rownames(test))
    names(T_INT) <- n
  }
  CY15_12h_check_MCL <- c(CY15_12h_check_MCL, list(T_INT))
  print(maxNumMCL - n)
  n <- n + 1
}
CY15_12h_enrichment_qvalue <- p.adjust(CY15_12h_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
CY15_12h_enriched_Node <- CY15_12h_check_MCL[CY15_12h_enrichment_qvalue < 0.05]
NumCY15_12h_subcluster <- c()
i <- 1
for(i in i:length(CY15_12h_enriched_Node)){
  NumCY15_12h_subcluster <- c(NumCY15_12h_subcluster, names(CY15_12h_enriched_Node[[i]])[1])
  i <- i+1
}
NumCY15_12h_subcluster <- as.numeric(NumCY15_12h_subcluster)



#CY15_24h enrichment#
CY15_24h_DEGs_FDR0005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY15_24h_FDR0005.txt", sep = "\t", row.names = 1, header = T)
#全クラスターのCY15_24hのDEGs数
NumCY15_24h_DEGs_FDR0005 <- sum(!is.na(match(rownames(CY15_24h_DEGs_FDR0005), NumMCL$name)))
#全クラスターのCY15_24h以外のDEGs数
Numothergenes <- c(length(NumMCL$name)-NumCY15_24h_DEGs_FDR0005)
maxNumMCL <- max(NumMCL$X__mclCluster)
n <- 1
CY15_24h_enrichment_pvalue <- c()
check <- c()
CY15_24h_check_MCL <- list()

for(n in n:maxNumMCL){
  test <- NumMCL[NumMCL$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(rownames(CY15_24h_DEGs_FDR0005), rownames(test)))
  Numsubcluster <- nrow(test)
  CY15_24h_enrichment_pvalue <- c(CY15_24h_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumCY15_24h_DEGs_FDR0005, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(CY15_24h_enrichment_pvalue)[n] <- n
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(rownames(CY15_24h_DEGs_FDR0005), rownames(test))
    names(T_INT) <- n
  }
  CY15_24h_check_MCL <- c(CY15_24h_check_MCL, list(T_INT))
  print(maxNumMCL - n)
  n <- n + 1
}
CY15_24h_enrichment_qvalue <- p.adjust(CY15_24h_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
CY15_24h_enriched_Node <- CY15_24h_check_MCL[CY15_24h_enrichment_qvalue < 0.05]
NumCY15_24h_subcluster <- c()
i <- 1
for(i in i:length(CY15_24h_enriched_Node)){
  NumCY15_24h_subcluster <- c(NumCY15_24h_subcluster, names(CY15_24h_enriched_Node[[i]])[1])
  i <- i+1
}
NumCY15_24h_subcluster <- as.numeric(NumCY15_24h_subcluster)



#CY15_48h enrichment#
CY15_48h_DEGs_FDR0005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY15_48h_FDR0005.txt", sep = "\t", row.names = 1, header = T)
#全クラスターのCY15_48hのDEGs数
NumCY15_48h_DEGs_FDR0005 <- sum(!is.na(match(rownames(CY15_48h_DEGs_FDR0005), NumMCL$name)))
#全クラスターのCY15_48h以外のDEGs数
Numothergenes <- c(length(NumMCL$name)-NumCY15_48h_DEGs_FDR0005)
maxNumMCL <- max(NumMCL$X__mclCluster)
n <- 1
CY15_48h_enrichment_pvalue <- c()
check <- c()
CY15_48h_check_MCL <- list()

for(n in n:maxNumMCL){
  test <- NumMCL[NumMCL$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(rownames(CY15_48h_DEGs_FDR0005), rownames(test)))
  Numsubcluster <- nrow(test)
  CY15_48h_enrichment_pvalue <- c(CY15_48h_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumCY15_48h_DEGs_FDR0005, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(CY15_48h_enrichment_pvalue)[n] <- n
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(rownames(CY15_48h_DEGs_FDR0005), rownames(test))
    names(T_INT) <- n
  }
  CY15_48h_check_MCL <- c(CY15_48h_check_MCL, list(T_INT))
  print(maxNumMCL - n)
  n <- n + 1
}
CY15_48h_enrichment_qvalue <- p.adjust(CY15_48h_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
CY15_48h_enriched_Node <- CY15_48h_check_MCL[CY15_48h_enrichment_qvalue < 0.05]
NumCY15_48h_subcluster <- c()
i <- 1
for(i in i:length(CY15_48h_enriched_Node)){
  NumCY15_48h_subcluster <- c(NumCY15_48h_subcluster, names(CY15_48h_enriched_Node[[i]])[1])
  i <- i+1
}
NumCY15_48h_subcluster <- as.numeric(NumCY15_48h_subcluster)


#names(CY20_enriched_Node[[11]])[1]でクラスター番号とってこれる。
table_CY15_1h_enriched_Node <- data.frame(unlist(CY15_1h_enriched_Node))
write.table(table_CY15_1h_enriched_Node, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY15FDR0005_time/table_CY15_1h_enriched_Node.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
table_CY15_3h_enriched_Node <- data.frame(unlist(CY15_3h_enriched_Node))
write.table(table_CY15_3h_enriched_Node, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY15FDR0005_time/table_CY15_3h_enriched_Node.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
table_CY15_12h_enriched_Node <- data.frame(unlist(CY15_12h_enriched_Node))
write.table(table_CY15_12h_enriched_Node, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY15FDR0005_time/table_CY15_12h_enriched_Node.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
table_CY15_24h_enriched_Node <- data.frame(unlist(CY15_24h_enriched_Node))
write.table(table_CY15_24h_enriched_Node, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY15FDR0005_time/table_CY15_24h_enriched_Node.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
table_CY15_48h_enriched_Node <- data.frame(unlist(CY15_48h_enriched_Node))
write.table(table_CY15_48h_enriched_Node, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY15FDR0005_time/table_CY15_48h_enriched_Node.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)

NumMCL <- na.omit(read.table("~/Nakano_RNAseq/network_analysis/cytoscape/NumMCL.txt", sep = "\t", header = T))
rownames(NumMCL) <- NumMCL$name
####CY16_1h enrichment####
CY16_1h_DEGs_FDR0005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY16_1h_FDR0005.txt", sep = "\t", row.names = 1, header = T)
#全クラスターのCY16_1hのDEGs数
NumCY16_1h_DEGs_FDR0005 <- sum(!is.na(match(rownames(CY16_1h_DEGs_FDR0005), NumMCL$name)))
#全クラスターのCY16_1h以外のDEGs数
Numothergenes <- c(length(NumMCL$name)-NumCY16_1h_DEGs_FDR0005)
maxNumMCL <- max(NumMCL$X__mclCluster)
n <- 1
CY16_1h_enrichment_pvalue <- c()
check <- c()
CY16_1h_check_MCL <- list()

for(n in n:maxNumMCL){
  test <- NumMCL[NumMCL$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(rownames(CY16_1h_DEGs_FDR0005), rownames(test)))
  Numsubcluster <- nrow(test)
  CY16_1h_enrichment_pvalue <- c(CY16_1h_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumCY16_1h_DEGs_FDR0005, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(CY16_1h_enrichment_pvalue)[n] <- n
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(rownames(CY16_1h_DEGs_FDR0005), rownames(test))
    names(T_INT) <- n
  }
  CY16_1h_check_MCL <- c(CY16_1h_check_MCL, list(T_INT))
  print(maxNumMCL - n)
  n <- n + 1
}
CY16_1h_enrichment_qvalue <- p.adjust(CY16_1h_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
CY16_1h_enriched_Node <- CY16_1h_check_MCL[CY16_1h_enrichment_qvalue < 0.05]
NumCY16_1h_subcluster <- c()
i <- 1
for(i in i:length(CY16_1h_enriched_Node)){
  NumCY16_1h_subcluster <- c(NumCY16_1h_subcluster, names(CY16_1h_enriched_Node[[i]])[1])
  i <- i+1
}
NumCY16_1h_subcluster <- as.numeric(NumCY16_1h_subcluster)


####CY16_3h enrichment####
CY16_3h_DEGs_FDR0005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY16_3h_FDR0005.txt", sep = "\t", row.names = 1, header = T)
#全クラスターのCY16_3hのDEGs数
NumCY16_3h_DEGs_FDR0005 <- sum(!is.na(match(rownames(CY16_3h_DEGs_FDR0005), NumMCL$name)))
#全クラスターのCY16_3h以外のDEGs数
Numothergenes <- c(length(NumMCL$name)-NumCY16_3h_DEGs_FDR0005)
maxNumMCL <- max(NumMCL$X__mclCluster)
n <- 1
CY16_3h_enrichment_pvalue <- c()
check <- c()
CY16_3h_check_MCL <- list()

for(n in n:maxNumMCL){
  test <- NumMCL[NumMCL$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(rownames(CY16_3h_DEGs_FDR0005), rownames(test)))
  Numsubcluster <- nrow(test)
  CY16_3h_enrichment_pvalue <- c(CY16_3h_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumCY16_3h_DEGs_FDR0005, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(CY16_3h_enrichment_pvalue)[n] <- n
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(rownames(CY16_3h_DEGs_FDR0005), rownames(test))
    names(T_INT) <- n
  }
  CY16_3h_check_MCL <- c(CY16_3h_check_MCL, list(T_INT))
  print(maxNumMCL - n)
  n <- n + 1
}
CY16_3h_enrichment_qvalue <- p.adjust(CY16_3h_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
CY16_3h_enriched_Node <- CY16_3h_check_MCL[CY16_3h_enrichment_qvalue < 0.05]
NumCY16_3h_subcluster <- c()
i <- 1
for(i in i:length(CY16_3h_enriched_Node)){
  NumCY16_3h_subcluster <- c(NumCY16_3h_subcluster, names(CY16_3h_enriched_Node[[i]])[1])
  i <- i+1
}
NumCY16_3h_subcluster <- as.numeric(NumCY16_3h_subcluster)



####CY16_12h enrichment####
CY16_12h_DEGs_FDR0005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY16_12h_FDR0005.txt", sep = "\t", row.names = 1, header = T)
#全クラスターのCY16_12hのDEGs数
NumCY16_12h_DEGs_FDR0005 <- sum(!is.na(match(rownames(CY16_12h_DEGs_FDR0005), NumMCL$name)))
#全クラスターのCY16_12h以外のDEGs数
Numothergenes <- c(length(NumMCL$name)-NumCY16_12h_DEGs_FDR0005)
maxNumMCL <- max(NumMCL$X__mclCluster)
n <- 1
CY16_12h_enrichment_pvalue <- c()
check <- c()
CY16_12h_check_MCL <- list()

for(n in n:maxNumMCL){
  test <- NumMCL[NumMCL$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(rownames(CY16_12h_DEGs_FDR0005), rownames(test)))
  Numsubcluster <- nrow(test)
  CY16_12h_enrichment_pvalue <- c(CY16_12h_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumCY16_12h_DEGs_FDR0005, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(CY16_12h_enrichment_pvalue)[n] <- n
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(rownames(CY16_12h_DEGs_FDR0005), rownames(test))
    names(T_INT) <- n
  }
  CY16_12h_check_MCL <- c(CY16_12h_check_MCL, list(T_INT))
  print(maxNumMCL - n)
  n <- n + 1
}
CY16_12h_enrichment_qvalue <- p.adjust(CY16_12h_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
CY16_12h_enriched_Node <- CY16_12h_check_MCL[CY16_12h_enrichment_qvalue < 0.05]
NumCY16_12h_subcluster <- c()
i <- 1
for(i in i:length(CY16_12h_enriched_Node)){
  NumCY16_12h_subcluster <- c(NumCY16_12h_subcluster, names(CY16_12h_enriched_Node[[i]])[1])
  i <- i+1
}
NumCY16_12h_subcluster <- as.numeric(NumCY16_12h_subcluster)



####CY16_24h enrichment####
CY16_24h_DEGs_FDR0005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY16_24h_FDR0005.txt", sep = "\t", row.names = 1, header = T)
#全クラスターのCY16_24hのDEGs数
NumCY16_24h_DEGs_FDR0005 <- sum(!is.na(match(rownames(CY16_24h_DEGs_FDR0005), NumMCL$name)))
#全クラスターのCY16_24h以外のDEGs数
Numothergenes <- c(length(NumMCL$name)-NumCY16_24h_DEGs_FDR0005)
maxNumMCL <- max(NumMCL$X__mclCluster)
n <- 1
CY16_24h_enrichment_pvalue <- c()
check <- c()
CY16_24h_check_MCL <- list()

for(n in n:maxNumMCL){
  test <- NumMCL[NumMCL$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(rownames(CY16_24h_DEGs_FDR0005), rownames(test)))
  Numsubcluster <- nrow(test)
  CY16_24h_enrichment_pvalue <- c(CY16_24h_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumCY16_24h_DEGs_FDR0005, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(CY16_24h_enrichment_pvalue)[n] <- n
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(rownames(CY16_24h_DEGs_FDR0005), rownames(test))
    names(T_INT) <- n
  }
  CY16_24h_check_MCL <- c(CY16_24h_check_MCL, list(T_INT))
  print(maxNumMCL - n)
  n <- n + 1
}
CY16_24h_enrichment_qvalue <- p.adjust(CY16_24h_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
CY16_24h_enriched_Node <- CY16_24h_check_MCL[CY16_24h_enrichment_qvalue < 0.05]
NumCY16_24h_subcluster <- c()
i <- 1
for(i in i:length(CY16_24h_enriched_Node)){
  NumCY16_24h_subcluster <- c(NumCY16_24h_subcluster, names(CY16_24h_enriched_Node[[i]])[1])
  i <- i+1
}
NumCY16_24h_subcluster <- as.numeric(NumCY16_24h_subcluster)



####CY16_48h enrichment####
CY16_48h_DEGs_FDR0005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY16_48h_FDR0005.txt", sep = "\t", row.names = 1, header = T)
#全クラスターのCY16_48hのDEGs数
NumCY16_48h_DEGs_FDR0005 <- sum(!is.na(match(rownames(CY16_48h_DEGs_FDR0005), NumMCL$name)))
#全クラスターのCY16_48h以外のDEGs数
Numothergenes <- c(length(NumMCL$name)-NumCY16_48h_DEGs_FDR0005)
maxNumMCL <- max(NumMCL$X__mclCluster)
n <- 1
CY16_48h_enrichment_pvalue <- c()
check <- c()
CY16_48h_check_MCL <- list()

for(n in n:maxNumMCL){
  test <- NumMCL[NumMCL$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(rownames(CY16_48h_DEGs_FDR0005), rownames(test)))
  Numsubcluster <- nrow(test)
  CY16_48h_enrichment_pvalue <- c(CY16_48h_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumCY16_48h_DEGs_FDR0005, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(CY16_48h_enrichment_pvalue)[n] <- n
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(rownames(CY16_48h_DEGs_FDR0005), rownames(test))
    names(T_INT) <- n
  }
  CY16_48h_check_MCL <- c(CY16_48h_check_MCL, list(T_INT))
  print(maxNumMCL - n)
  n <- n + 1
}
CY16_48h_enrichment_qvalue <- p.adjust(CY16_48h_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
CY16_48h_enriched_Node <- CY16_48h_check_MCL[CY16_48h_enrichment_qvalue < 0.05]
NumCY16_48h_subcluster <- c()
i <- 1
for(i in i:length(CY16_48h_enriched_Node)){
  NumCY16_48h_subcluster <- c(NumCY16_48h_subcluster, names(CY16_48h_enriched_Node[[i]])[1])
  i <- i+1
}
NumCY16_48h_subcluster <- as.numeric(NumCY16_48h_subcluster)


#names(CY20_enriched_Node[[11]])[1]でクラスター番号とってこれる。
table_CY16_1h_enriched_Node <- data.frame(unlist(CY16_1h_enriched_Node))
write.table(table_CY16_1h_enriched_Node, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY16FDR0005_time/table_CY16_1h_enriched_Node.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
table_CY16_3h_enriched_Node <- data.frame(unlist(CY16_3h_enriched_Node))
write.table(table_CY16_3h_enriched_Node, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY16FDR0005_time/table_CY16_3h_enriched_Node.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
table_CY16_12h_enriched_Node <- data.frame(unlist(CY16_12h_enriched_Node))
write.table(table_CY16_12h_enriched_Node, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY16FDR0005_time/table_CY16_12h_enriched_Node.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
table_CY16_24h_enriched_Node <- data.frame(unlist(CY16_24h_enriched_Node))
write.table(table_CY16_24h_enriched_Node, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY16FDR0005_time/table_CY16_24h_enriched_Node.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
table_CY16_48h_enriched_Node <- data.frame(unlist(CY16_48h_enriched_Node))
write.table(table_CY16_48h_enriched_Node, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY16FDR0005_time/table_CY16_48h_enriched_Node.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)


NumMCL <- na.omit(read.table("~/Nakano_RNAseq/network_analysis/cytoscape/NumMCL.txt", sep = "\t", header = T))
rownames(NumMCL) <- NumMCL$name
####CY20_1h enrichment####
CY20_1h_DEGs_FDR0005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY20_1h_FDR0005.txt", sep = "\t", row.names = 1, header = T)
#全クラスターのCY20_1hのDEGs数
NumCY20_1h_DEGs_FDR0005 <- sum(!is.na(match(rownames(CY20_1h_DEGs_FDR0005), NumMCL$name)))
#全クラスターのCY20_1h以外のDEGs数
Numothergenes <- c(length(NumMCL$name)-NumCY20_1h_DEGs_FDR0005)
maxNumMCL <- max(NumMCL$X__mclCluster)
n <- 1
CY20_1h_enrichment_pvalue <- c()
check <- c()
CY20_1h_check_MCL <- list()

for(n in n:maxNumMCL){
  test <- NumMCL[NumMCL$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(rownames(CY20_1h_DEGs_FDR0005), rownames(test)))
  Numsubcluster <- nrow(test)
  CY20_1h_enrichment_pvalue <- c(CY20_1h_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumCY20_1h_DEGs_FDR0005, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(CY20_1h_enrichment_pvalue)[n] <- n
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(rownames(CY20_1h_DEGs_FDR0005), rownames(test))
    names(T_INT) <- n
  }
  CY20_1h_check_MCL <- c(CY20_1h_check_MCL, list(T_INT))
  print(maxNumMCL - n)
  n <- n + 1
}
CY20_1h_enrichment_qvalue <- p.adjust(CY20_1h_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
CY20_1h_enriched_Node <- CY20_1h_check_MCL[CY20_1h_enrichment_qvalue < 0.05]
NumCY20_1h_subcluster <- c()
i <- 1
for(i in i:length(CY20_1h_enriched_Node)){
  NumCY20_1h_subcluster <- c(NumCY20_1h_subcluster, names(CY20_1h_enriched_Node[[i]])[1])
  i <- i+1
}
NumCY20_1h_subcluster <- as.numeric(NumCY20_1h_subcluster)


####CY20_3h enrichment####
CY20_3h_DEGs_FDR0005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY20_3h_FDR0005.txt", sep = "\t", row.names = 1, header = T)
#全クラスターのCY20_3hのDEGs数
NumCY20_3h_DEGs_FDR0005 <- sum(!is.na(match(rownames(CY20_3h_DEGs_FDR0005), NumMCL$name)))
#全クラスターのCY20_3h以外のDEGs数
Numothergenes <- c(length(NumMCL$name)-NumCY20_3h_DEGs_FDR0005)
maxNumMCL <- max(NumMCL$X__mclCluster)
n <- 1
CY20_3h_enrichment_pvalue <- c()
check <- c()
CY20_3h_check_MCL <- list()

for(n in n:maxNumMCL){
  test <- NumMCL[NumMCL$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(rownames(CY20_3h_DEGs_FDR0005), rownames(test)))
  Numsubcluster <- nrow(test)
  CY20_3h_enrichment_pvalue <- c(CY20_3h_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumCY20_3h_DEGs_FDR0005, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(CY20_3h_enrichment_pvalue)[n] <- n
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(rownames(CY20_3h_DEGs_FDR0005), rownames(test))
    names(T_INT) <- n
  }
  CY20_3h_check_MCL <- c(CY20_3h_check_MCL, list(T_INT))
  print(maxNumMCL - n)
  n <- n + 1
}
CY20_3h_enrichment_qvalue <- p.adjust(CY20_3h_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
CY20_3h_enriched_Node <- CY20_3h_check_MCL[CY20_3h_enrichment_qvalue < 0.05]
NumCY20_3h_subcluster <- c()
i <- 1
for(i in i:length(CY20_3h_enriched_Node)){
  NumCY20_3h_subcluster <- c(NumCY20_3h_subcluster, names(CY20_3h_enriched_Node[[i]])[1])
  i <- i+1
}
NumCY20_3h_subcluster <- as.numeric(NumCY20_3h_subcluster)



####CY20_12h enrichment####
CY20_12h_DEGs_FDR0005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY20_12h_FDR0005.txt", sep = "\t", row.names = 1, header = T)
#全クラスターのCY20_12hのDEGs数
NumCY20_12h_DEGs_FDR0005 <- sum(!is.na(match(rownames(CY20_12h_DEGs_FDR0005), NumMCL$name)))
#全クラスターのCY20_12h以外のDEGs数
Numothergenes <- c(length(NumMCL$name)-NumCY20_12h_DEGs_FDR0005)
maxNumMCL <- max(NumMCL$X__mclCluster)
n <- 1
CY20_12h_enrichment_pvalue <- c()
check <- c()
CY20_12h_check_MCL <- list()

for(n in n:maxNumMCL){
  test <- NumMCL[NumMCL$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(rownames(CY20_12h_DEGs_FDR0005), rownames(test)))
  Numsubcluster <- nrow(test)
  CY20_12h_enrichment_pvalue <- c(CY20_12h_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumCY20_12h_DEGs_FDR0005, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(CY20_12h_enrichment_pvalue)[n] <- n
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(rownames(CY20_12h_DEGs_FDR0005), rownames(test))
    names(T_INT) <- n
  }
  CY20_12h_check_MCL <- c(CY20_12h_check_MCL, list(T_INT))
  print(maxNumMCL - n)
  n <- n + 1
}
CY20_12h_enrichment_qvalue <- p.adjust(CY20_12h_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
CY20_12h_enriched_Node <- CY20_12h_check_MCL[CY20_12h_enrichment_qvalue < 0.05]
NumCY20_12h_subcluster <- c()
i <- 1
for(i in i:length(CY20_12h_enriched_Node)){
  NumCY20_12h_subcluster <- c(NumCY20_12h_subcluster, names(CY20_12h_enriched_Node[[i]])[1])
  i <- i+1
}
NumCY20_12h_subcluster <- as.numeric(NumCY20_12h_subcluster)



####CY20_24h enrichment####
CY20_24h_DEGs_FDR0005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY20_24h_FDR0005.txt", sep = "\t", row.names = 1, header = T)
#全クラスターのCY20_24hのDEGs数
NumCY20_24h_DEGs_FDR0005 <- sum(!is.na(match(rownames(CY20_24h_DEGs_FDR0005), NumMCL$name)))
#全クラスターのCY20_24h以外のDEGs数
Numothergenes <- c(length(NumMCL$name)-NumCY20_24h_DEGs_FDR0005)
maxNumMCL <- max(NumMCL$X__mclCluster)
n <- 1
CY20_24h_enrichment_pvalue <- c()
check <- c()
CY20_24h_check_MCL <- list()

for(n in n:maxNumMCL){
  test <- NumMCL[NumMCL$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(rownames(CY20_24h_DEGs_FDR0005), rownames(test)))
  Numsubcluster <- nrow(test)
  CY20_24h_enrichment_pvalue <- c(CY20_24h_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumCY20_24h_DEGs_FDR0005, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(CY20_24h_enrichment_pvalue)[n] <- n
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(rownames(CY20_24h_DEGs_FDR0005), rownames(test))
    names(T_INT) <- n
  }
  CY20_24h_check_MCL <- c(CY20_24h_check_MCL, list(T_INT))
  print(maxNumMCL - n)
  n <- n + 1
}
CY20_24h_enrichment_qvalue <- p.adjust(CY20_24h_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
CY20_24h_enriched_Node <- CY20_24h_check_MCL[CY20_24h_enrichment_qvalue < 0.05]
NumCY20_24h_subcluster <- c()
i <- 1
for(i in i:length(CY20_24h_enriched_Node)){
  NumCY20_24h_subcluster <- c(NumCY20_24h_subcluster, names(CY20_24h_enriched_Node[[i]])[1])
  i <- i+1
}
NumCY20_24h_subcluster <- as.numeric(NumCY20_24h_subcluster)



####CY20_48h enrichment####
CY20_48h_DEGs_FDR0005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY20_48h_FDR0005.txt", sep = "\t", row.names = 1, header = T)
#全クラスターのCY20_48hのDEGs数
NumCY20_48h_DEGs_FDR0005 <- sum(!is.na(match(rownames(CY20_48h_DEGs_FDR0005), NumMCL$name)))
#全クラスターのCY20_48h以外のDEGs数
Numothergenes <- c(length(NumMCL$name)-NumCY20_48h_DEGs_FDR0005)
maxNumMCL <- max(NumMCL$X__mclCluster)
n <- 1
CY20_48h_enrichment_pvalue <- c()
check <- c()
CY20_48h_check_MCL <- list()

for(n in n:maxNumMCL){
  test <- NumMCL[NumMCL$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(rownames(CY20_48h_DEGs_FDR0005), rownames(test)))
  Numsubcluster <- nrow(test)
  CY20_48h_enrichment_pvalue <- c(CY20_48h_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumCY20_48h_DEGs_FDR0005, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(CY20_48h_enrichment_pvalue)[n] <- n
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(rownames(CY20_48h_DEGs_FDR0005), rownames(test))
    names(T_INT) <- n
  }
  CY20_48h_check_MCL <- c(CY20_48h_check_MCL, list(T_INT))
  print(maxNumMCL - n)
  n <- n + 1
}
CY20_48h_enrichment_qvalue <- p.adjust(CY20_48h_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
CY20_48h_enriched_Node <- CY20_48h_check_MCL[CY20_48h_enrichment_qvalue < 0.05]
NumCY20_48h_subcluster <- c()
i <- 1
for(i in i:length(CY20_48h_enriched_Node)){
  NumCY20_48h_subcluster <- c(NumCY20_48h_subcluster, names(CY20_48h_enriched_Node[[i]])[1])
  i <- i+1
}
NumCY20_48h_subcluster <- as.numeric(NumCY20_48h_subcluster)


#names(CY20_enriched_Node[[11]])[1]でクラスター番号とってこれる。
table_CY20_1h_enriched_Node <- data.frame(unlist(CY20_1h_enriched_Node))
write.table(table_CY20_1h_enriched_Node, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY20FDR0005_time/table_CY20_1h_enriched_Node.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
table_CY20_3h_enriched_Node <- data.frame(unlist(CY20_3h_enriched_Node))
write.table(table_CY20_3h_enriched_Node, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY20FDR0005_time/table_CY20_3h_enriched_Node.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
table_CY20_12h_enriched_Node <- data.frame(unlist(CY20_12h_enriched_Node))
write.table(table_CY20_12h_enriched_Node, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY20FDR0005_time/table_CY20_12h_enriched_Node.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
table_CY20_24h_enriched_Node <- data.frame(unlist(CY20_24h_enriched_Node))
write.table(table_CY20_24h_enriched_Node, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY20FDR0005_time/table_CY20_24h_enriched_Node.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
table_CY20_48h_enriched_Node <- data.frame(unlist(CY20_48h_enriched_Node))
write.table(table_CY20_48h_enriched_Node, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY20FDR0005_time/table_CY20_48h_enriched_Node.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
