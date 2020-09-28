iterations14 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/iteration14_node.txt", sep = "\t", header = T)
rownames(iterations14) <- iterations14$name
####CY15 enrichment####
CY15_DEGs_FDR0005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/union/FDR0005/CY15_FDR0005_time.txt", sep = "\t", row.names = 1, header = T)
#全クラスターのCY15のDEGs数
NumCY15_DEGs_FDR0005 <- sum(!is.na(match(rownames(CY15_DEGs_FDR0005), iterations14$name)))
#全クラスターのCY15以外のDEGs数
Numothergenes <- c(length(iterations14$name)-NumCY15_DEGs_FDR0005)
maxiterations14 <- max(iterations14$X__mclCluster)
n <- 1
CY15_enrichment_pvalue <- c()
check <- c()
CY15_check_MCL <- list()
subnodes <- list()
for(n in n:maxiterations14){
  test <- iterations14[iterations14$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(rownames(CY15_DEGs_FDR0005), rownames(test)))
  Numsubcluster <- nrow(test)
  CY15_enrichment_pvalue <- c(CY15_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumCY15_DEGs_FDR0005, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(CY15_enrichment_pvalue)[n] <- n
  subnodes <- c(subnodes, list(Numsubcluster))
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(rownames(CY15_DEGs_FDR0005), rownames(test))
    names(T_INT) <- rep(n, times = length(T_INT))
    CY15_check_MCL <- c(CY15_check_MCL, list(T_INT))
  }
  print(maxiterations14 - n)
  n <- n + 1
}
#各サブクラスターにおいてp-value < 0.05であるAGIを引っ張れている。
CY15_enriched_Node <- CY15_check_MCL[CY15_enrichment_pvalue < 0.05]
NumCY15_subcluster <- c()
i <- 1
for(i in i:length(CY15_enriched_Node)){
  NumCY15_subcluster <- c(NumCY15_subcluster, names(CY15_enriched_Node[[i]])[1])
  i <- i+1
}
NumCY15_subcluster <- as.numeric(NumCY15_subcluster)
####CY16 enrichment####
CY16_DEGs_FDR0005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/union/FDR0005/CY16_FDR0005_time.txt", sep = "\t", row.names = 1, header = T)
#全クラスターのCY16のDEGs数
NumCY16_DEGs_FDR0005 <- sum(!is.na(match(rownames(CY16_DEGs_FDR0005), iterations14$name)))
#全クラスターのCY16以外のDEGs数
Numothergenes <- c(length(iterations14$name)-NumCY16_DEGs_FDR0005)
maxiterations14 <- max(iterations14$X__mclCluster)
n <- 1
CY16_enrichment_pvalue <- c()
check <- c()
CY16_check_MCL <- list()
subnodes <- list()
for(n in n:maxiterations14){
  test <- iterations14[iterations14$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(rownames(CY16_DEGs_FDR0005), rownames(test)))
  Numsubcluster <- nrow(test)
  CY16_enrichment_pvalue <- c(CY16_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumCY16_DEGs_FDR0005, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(CY16_enrichment_pvalue)[n] <- n
  subnodes <- c(subnodes, list(Numsubcluster))
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(rownames(CY16_DEGs_FDR0005), rownames(test))
    names(T_INT) <- rep(n, times = length(T_INT))
    CY16_check_MCL <- c(CY16_check_MCL, list(T_INT))
  }
  print(maxiterations14 - n)
  n <- n + 1
}
#各サブクラスターにおいてp-value < 0.05であるAGIを引っ張れている。
CY16_enriched_Node <- CY16_check_MCL[CY16_enrichment_pvalue < 0.05]
NumCY16_subcluster <- c()
i <- 1
for(i in i:length(CY16_enriched_Node)){
  NumCY16_subcluster <- c(NumCY16_subcluster, names(CY16_enriched_Node[[i]])[1])
  i <- i+1
}
NumCY16_subcluster <- as.numeric(NumCY16_subcluster)
####CY20 enrichment####
CY20_DEGs_FDR0005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/union/FDR0005/CY20_FDR0005_time.txt", sep = "\t", row.names = 1, header = T)
#全クラスターのCY20のDEGs数
NumCY20_DEGs_FDR0005 <- sum(!is.na(match(rownames(CY20_DEGs_FDR0005), iterations14$name)))
#全クラスターのCY20以外のDEGs数
Numothergenes <- c(length(iterations14$name)-NumCY20_DEGs_FDR0005)
maxiterations14 <- max(iterations14$X__mclCluster)
n <- 1
CY20_enrichment_pvalue <- c()
check <- c()
CY20_check_MCL <- list()
subnodes <- list()
for(n in n:maxiterations14){
  test <- iterations14[iterations14$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(rownames(CY20_DEGs_FDR0005), rownames(test)))
  Numsubcluster <- nrow(test)
  CY20_enrichment_pvalue <- c(CY20_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumCY20_DEGs_FDR0005, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(CY20_enrichment_pvalue)[n] <- n
  subnodes <- c(subnodes, list(Numsubcluster))
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(rownames(CY20_DEGs_FDR0005), rownames(test))
    names(T_INT) <- rep(n, times = length(T_INT))
    CY20_check_MCL <- c(CY20_check_MCL, list(T_INT))
  }
  print(maxiterations14 - n)
  n <- n + 1
}
#各サブクラスターにおいてp-value < 0.05であるAGIを引っ張れている。
CY20_enriched_Node <- CY20_check_MCL[CY20_enrichment_pvalue < 0.05]
NumCY20_subcluster <- c()
i <- 1
for(i in i:length(CY20_enriched_Node)){
  NumCY20_subcluster <- c(NumCY20_subcluster, names(CY20_enriched_Node[[i]])[1])
  i <- i+1
}
NumCY20_subcluster <- as.numeric(NumCY20_subcluster)
####output####
#names(CY20_enriched_Node[[11]])[1]でクラスター番号とってこれる。
CY20_enriched_Node <- data.frame(unlist(CY20_enriched_Node))
CY16_enriched_Node <- data.frame(unlist(CY16_enriched_Node))
CY15_enriched_Node <- data.frame(unlist(CY15_enriched_Node))
write.table(CY20_enriched_Node, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/iterations14/CY20_enriched_Node.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
write.table(CY16_enriched_Node, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/iterations14/CY16_enriched_Node.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
write.table(CY15_enriched_Node, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/iterations14/CY15_enriched_Node.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)

write.table(CY20_enrichment_pvalue, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/iterations14/CY_subcluster/CY20_enrichment_pvalue.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
write.table(CY16_enrichment_pvalue, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/iterations14/CY_subcluster/CY16_enrichment_pvalue.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
write.table(CY15_enrichment_pvalue, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/iterations14/CY_subcluster/CY15_enrichment_pvalue.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)