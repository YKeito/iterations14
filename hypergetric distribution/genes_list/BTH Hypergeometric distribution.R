iterations14 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/iteration14_node.txt", sep = "\t", header = T)
rownames(iterations14) <- iterations14$name
####BTH enrichment####
BTH_DEGs <- read.table("~/Nakano_RNAseq/network_analysis/base/BTH_affected_AGI.txt", sep = "\t", header = T)
rownames(BTH_DEGs) <- BTH_DEGs$AGI
#全クラスターのBTHのDEGs数
NumBTH_DEGs <- sum(!is.na(match(BTH_DEGs$AGI, iterations14$name)))
#全クラスターのBTH以外のDEGs数
Numothergenes <- c(length(iterations14$name)-NumBTH_DEGs)
maxiterations14 <- max(iterations14$X__mclCluster)
n <- 1
BTH_enrichment_pvalue <- c()
check <- c()
BTH_check_MCL <- list()
for(n in n:maxiterations14){
  test <- iterations14[iterations14$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(rownames(BTH_DEGs), rownames(test)))
  Numsubcluster <- nrow(test)
  BTH_enrichment_pvalue <- c(BTH_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumBTH_DEGs, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(BTH_enrichment_pvalue)[n] <- n
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(rownames(BTH_DEGs), rownames(test))
    names(T_INT) <- rep(n, times = length(T_INT))
    BTH_check_MCL <- c(BTH_check_MCL, list(T_INT))
  }
  print(maxiterations14 - n)
  n <- n + 1
}
BTH_enrichment_qvalue <- p.adjust(BTH_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
BTH_enriched_Node <- BTH_check_MCL[BTH_enrichment_qvalue < 0.05]
NumBTH_subcluster <- c()
i <- 1
for(i in i:length(BTH_enriched_Node)){
  NumBTH_subcluster <- c(NumBTH_subcluster, names(BTH_enriched_Node[[i]])[1])
  i <- i+1
}
NumBTH_subcluster <- as.numeric(NumBTH_subcluster)
#names(BTH_enriched_Node[[11]])[1]でクラスター番号とってこれる。
BTH_enriched_Node <- data.frame(unlist(BTH_enriched_Node))
BTH_statistics <- data.frame(NumBRT_DEGs = check, 
                             NumSubcluster = unlist(subnodes), 
                             p_value = BTH_enrichment_pvalue, 
                             q_value = BTH_enrichment_qvalue)

write.table(BTH_enriched_Node, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/iterations14/genes_set/BTH_enriched_Node.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
save(BTH_enrichment_qvalue, file = "~/Nakano_RNAseq/network_analysis/.RData/iterations14/BTH_enrichment_qvalue.RData")
write.table(BTH_statistics, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/iterations14/genes_set/BTH_statistics.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
