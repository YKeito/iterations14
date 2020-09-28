iterations14 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/iteration14_node.txt", sep = "\t", header = T)
rownames(iterations14) <- iterations14$name
####MAPK enrichment####
MAPK <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/MAPK_signaling_pathway.txt", sep = "\t", header = T)
MAPK <- strtrim(MAPK$AGI, 9)
#全クラスターのMAPKのDEGs数
NumMAPK <- sum(!is.na(match(MAPK, iterations14$name)))
#全クラスターのMAPK以外のDEGs数
Numothergenes <- c(length(iterations14$name)-NumMAPK)
maxiterations14 <- max(iterations14$X__mclCluster)
n <- 1
MAPK_enrichment_pvalue <- c()
check <- c()
MAPK_check_MCL <- list()
subnodes <- list()
for(n in n:maxiterations14){
  test <- iterations14[iterations14$X__mclCluster == n, ]
  NumMAPK_DEGs <- length(intersect(MAPK, rownames(test)))
  Numsubcluster <- nrow(test)
  MAPK_enrichment_pvalue <- c(MAPK_enrichment_pvalue, phyper(c(NumMAPK_DEGs - 1), NumMAPK, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(MAPK_enrichment_pvalue)[n] <- n
  subnodes <- c(subnodes, list(Numsubcluster))
  check <- c(check, NumMAPK_DEGs)
  T_INT <- intersect(MAPK, rownames(test))
  names(T_INT) <- rep(n, times = length(T_INT))
  MAPK_check_MCL <- c(MAPK_check_MCL, list(T_INT))
  print(maxiterations14 - n)
  n <- n + 1
}
#各サブクラスターにおいてp-value < 0.05であるAGIを引っ張れている。
data <- c()
AGI <- list()
i <- 1
for(i in i:length(MAPK_check_MCL)){
  data <- unlist(MAPK_check_MCL[i])
  data <- paste(data, collapse = " | ")
  data[which(data %in% "")] <- "NA"
  AGI <- c(AGI, list(data))
  data <- c()
  m <- m+1
  i <- i+1
}

NumSubCluster <- c()
i <- 1
for(i in i:maxiterations14){
  NumSubCluster <- c(NumSubCluster, paste0("SubCluster", i))
  i <- i+1
}
MAPK_statistics <- data.frame(NumSubcluster = NumSubCluster,
                             NumMAPK_DEGs = check, 
                             MAPK_DEGs_AGI = unlist(AGI),
                             NumSubcluster = unlist(subnodes), 
                             p_value = MAPK_enrichment_pvalue, 
                             enrichment_score = -log2(MAPK_enrichment_pvalue)
                             )

####output####
write.table(MAPK_statistics, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/iterations14/genes_set/MAPK_statistics.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
