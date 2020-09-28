iterations14 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/iteration14_node.txt", sep = "\t", header = T)
rownames(iterations14) <- iterations14$name
####TF_family enrichment####
TF_family <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/arabidopsis_TF_family.txt", sep = "\t", header = T)
maxiterations14 <- max(iterations14$X__mclCluster)
maxTF <- length(unique(TF_family$TF))
i <- 1
temp <- c()
all_TF <- as.character(unique(TF_family$TF))
allTF_enrichment_pvalue <- c()
allTF_check_MCL <- list()
for(i in i:maxTF){
  temp <- TF_family[grep(unique(TF_family$TF)[i], TF_family$TF), ]$AGI
  #母集団成功数
  NumTF <- length(intersect(temp, iterations14$name))
  #母集団成功数-成功数
  Numothergenes <- c(length(iterations14$name)-NumTF)
  n <- 1
  TF_enrichment_pvalue <- c()
  check <- c()
  TF_check_MCL <- c()
  subnodes <- list()
  for(n in n:maxiterations14){
    test <- iterations14[iterations14$X__mclCluster == n, ]
    #標本数
    Numsubcluster <- nrow(test)
    #
    NumTF_DEGs <- length(intersect(temp, test$name))
    TF_enrichment_pvalue <- c(TF_enrichment_pvalue, phyper(c(NumTF_DEGs - 1), NumTF, Numothergenes, Numsubcluster, lower.tail = FALSE))
    names(TF_enrichment_pvalue)[n] <- n
    subnodes <- c(subnodes, list(Numsubcluster))
    check <- c(check, NumTF_DEGs)
    T_INT <- intersect(rownames(test), iterations14$name)
    names(T_INT) <- rep(n, times = length(T_INT))
    T_INT <- paste(T_INT, collapse = " | ")
    T_INT[which(T_INT %in% "")] <- "NA"
    TF_check_MCL <- c(TF_check_MCL, T_INT)
    n <- n + 1
  }
  allTF_enrichment_pvalue <- cbind(allTF_enrichment_pvalue, TF_enrichment_pvalue)
  TF_enrichment_pvalue <- c()
  allTF_check_MCL <- c(allTF_check_MCL, list(TF_check_MCL))
  TF_check_MCL <- c()
  print(maxTF - i)
  i <- i+1
}
SumTF_subcluster <- data.frame(allTF_enrichment_pvalue)
colnames(SumTF_subcluster) <- unique(TF_family$TF)
SumTF_subcluster <- as.data.frame(t(SumTF_subcluster))
allTF_DEGs <- data.frame(DEGs = unlist(allTF_check_MCL[1]))
####output####
write.table(SumTF_subcluster, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/iterations14/allsubcluster/SumTF_subcluster.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
write.table(allTF_DEGs, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/iterations14/allsubcluster/allTF_DEGs.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
