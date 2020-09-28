#ATTED <- read.table("~/bigdata/yasue/ATTEDII/GeneExp_v3", fill = T, stringsAsFactors = F)
rownames(ATTED) <- toupper(rownames(ATTED))
rownames(allRNASeq) <- toupper(rownames(allRNASeq))
test <- intersect(CY_hub, rownames(ATTED))
ATTED_test <- ATTED[match(test, rownames(ATTED)), ] #sum(is.na(ATTED_test))
ATTED_test <- as.matrix(ATTED_test)
n <- 1
m <- 2
ATTED_PCC <- c()
ATTED_PCC_all <- list()
for(n in n:c(nrow(ATTED_test)-1)){
  for(m in m:nrow(ATTED_test)){
    ATTED_PCC <- c(ATTED_PCC, cor(ATTED_test[n, ], ATTED_test[m, ], method = "pearson"))
    m <- m+1
  }
  ATTED_PCC_all <- c(ATTED_PCC_all, list(ATTED_PCC))
  ATTED_PCC <- c()
  print(n)
  n <- n+1
  m <- n+1
}
#####Cytoscape_format####
#source genes
rowtotal <- nrow(ATTED_test)
m <- 1
source_genes <- c()
for(m in m:c(rowtotal-1)){
  source_genes <- c(source_genes, rep(rownames(ATTED_test)[m], times = c(rowtotal-m)))
  print(m)
  m <- m+1
}
#target genes
n <- 2
target_genes <- c()
for(n in n:rowtotal){
  target_genes <- c(target_genes, rownames(ATTED_test)[n:rowtotal])
  print(n)
  n <- n+1
}

ATTED_cytoscape <- data.frame(source_genes = unlist(source_genes),
                              interaction_value = unlist(ATTED_PCC_all),
                              target_genes = unlist(target_genes)
                              )
save(ATTED_cytoscape, file = "~/bigdata/yasue/ATTEDII/ATTED_hub_cytoscape_th.RData")