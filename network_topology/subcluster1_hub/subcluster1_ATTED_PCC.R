#ATTED <- read.table("~/bigdata/yasue/ATTEDII/GeneExp_v3", fill = T, stringsAsFactors = F)
rownames(ATTED) <- toupper(rownames(ATTED))
rownames(allRNASeq) <- toupper(rownames(allRNASeq))
test <- intersect(rownames(allRNASeq), rownames(ATTED))
ATTED_test <- ATTED[match(test, rownames(ATTED)), ] #sum(is.na(ATTED_test))
ATTED_test <- as.matrix(ATTED_test)
n <- 1
m <- 2
base <- c()
PCC <- c()
PCC_all <- list()
PCC_pvalue <- c()
PCC_pvalue_all <- list()
n <- 1
m <- 2
base <- c()
PCC <- c()
PCC_all <- list()
PCC_pvalue <- c()
PCC_pvalue_all <- list()
for(n in n:c(nrow(ATTED_test)-1)){
  for(m in m:nrow(ATTED_test)){
    base <- cor.test(ATTED_test[n, ], ATTED_test[m, ], method = "pearson")
    PCC <- c(PCC, base$estimate)
    PCC_pvalue <- c(PCC_pvalue, base$p.value)
    m <- m+1
  }
  PCC_all <- c(PCC_all, list(PCC))
  PCC_pvalue_all <- c(PCC_pvalue_all, list(PCC_pvalue))
  PCC <- c()
  PCC_pvalue <- c()
  print(n)
  n <- n+1
  m <- n+1
}
####PCC q-value####
PCC_qvalue_all <- p.adjust(unlist(PCC_pvalue_all), method = "BH")
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
                              interaction_value = unlist(PCC_all),
                              target_genes = unlist(target_genes),
                              p_value = unlist(PCC_pvalue_all),
                              q_value = PCC_qvalue_all
)
ATTED_cytoscape_th <- ATTED_cytoscape[ATTED_cytoscape$q_value < 0.05, ]
save(ATTED_cytoscape_th, file = "~/bigdata/yasue/ATTEDII/ATTED_cytoscape_th.RData")