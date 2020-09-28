#
test <- intersect(CY_hub, rownames(allRNASeq))
allRNASeq_test <- allRNASeq[match(test, rownames(allRNASeq)), ] #sum(is.na(allRNASeq_test))
allRNASeq_test <- as.matrix(allRNASeq_test)
n <- 1
m <- 2
allRNASeq_PCC <- c()
allRNASeq_PCC_all <- list()
for(n in n:c(nrow(allRNASeq_test)-1)){
  for(m in m:nrow(allRNASeq_test)){
    allRNASeq_PCC <- c(allRNASeq_PCC, cor(allRNASeq_test[n, ], allRNASeq_test[m, ], method = "pearson"))
    m <- m+1
  }
  allRNASeq_PCC_all <- c(allRNASeq_PCC_all, list(allRNASeq_PCC))
  allRNASeq_PCC <- c()
  print(n)
  n <- n+1
  m <- n+1
}
#####Cytoscape_format####
#source genes
rowtotal <- nrow(allRNASeq_test)
m <- 1
source_genes <- c()
for(m in m:c(rowtotal-1)){
  source_genes <- c(source_genes, rep(rownames(allRNASeq_test)[m], times = c(rowtotal-m)))
  print(m)
  m <- m+1
}
#target genes
n <- 2
target_genes <- c()
for(n in n:rowtotal){
  target_genes <- c(target_genes, rownames(allRNASeq_test)[n:rowtotal])
  print(n)
  n <- n+1
}

allRNASeq_cytoscape <- data.frame(source_genes = unlist(source_genes),
                              interaction_value = unlist(allRNASeq_PCC_all),
                              target_genes = unlist(target_genes)
)


save(allRNASeq_cytoscape, file = "~/Nakano_RNAseq/network_analysis/base/Subcluster_Information/allRNASeq_hub_cytoscape.RData")