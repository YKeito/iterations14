iterations14 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/iteration14_node.txt", sep = "\t", header = T, stringsAsFactors = F)
rownames(iterations14) <- iterations14$name
####TF_family enrichment####
TF_family <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/arabidopsis_TF_family.txt", sep = "\t", header = T, stringsAsFactors = F)
N <- nrow(allRNASeq)
TF_AGI <- TF_family$AGI
TF_RNASeq <- intersect(rownames(allRNASeq), TF_AGI)
M <- length(TF_RNASeq)
enrichment_score <- c()
i <- 1
for(i in i:442){
  subcluster_node <- iterations14$name[iterations14$X__mclCluster == i]
  n <- length(subcluster_node)
  x <- length(intersect(subcluster_node, TF_RNASeq))
  enrichment_score <- c(enrichment_score, phyper(x-1, M, N-M, n, lower.tail = F))
  print(i)
  i <- i+1
}