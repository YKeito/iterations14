iterations14 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/iteration14_node.txt", sep = "\t", header = T)
rownames(iterations14) <- iterations14$name
####TF_target enrichment####
TF_target <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/arabidopsis_TF_target.txt", sep = "\t", header = T)
TF_target$TargetLocus <- toupper(TF_target$TargetLocus)
TF_target$TFLocus <- toupper(TF_target$TFLocus)
population <- iterations14$name
allTF <- unique(TF_target$TFLocus)
####
allTF_target_pvalue <- c()
total_i <- length(unique(TF_target$TFLocus))
i <- 1
for(i in i:total_i){
  target <- allTF[i]
  population_target <- intersect(target, population)
  n <- 1
  total_n <- max(iterations14$X__mclCluster)
  TF_target_pvalue <- c()
  for(n in n:total_n){
    sample <- iterations14[iterations14$X__mclCluster == n, ]$name
    sample_target <- intersect(population_target, sample)
    a <- length(population)
    b <- length(population_target)
    c <- length(sample)
    d <- length(sample_target)
    TF_target_pvalue <- c(TF_target_pvalue, phyper(c(d-1), b, c(a-b), c, lower.tail = F))
    n <- n+1
  }
  allTF_target_pvalue <- cbind(allTF_target_pvalue, TF_target_pvalue)
  print(i)
  i < i+1
}
allTF_target_enrichment <- t(data.frame(allTF_target_pvalue))
rownames(allTF_target_enrichment) <- allTF
colnames(allTF_target_enrichment) <- 1:442
allTF_target_enrichment