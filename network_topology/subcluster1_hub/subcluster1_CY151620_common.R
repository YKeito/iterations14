CY151620_FDR005 <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/CYall_FDR005.txt", sep = "\t", header = T)
CY151620_common <- CY151620_FDR005[CY151620_FDR005$CYall_Venn == "CY15:CY16:CY20", ]$AGI
subcluster1_node <- iterations14[iterations14$X__mclCluster == 1, ]$name
subcluster1_CY151620_common <- intersect(CY151620_common, subcluster1_node)
subcluster1_CY151620_common <- allRNASeq[match(subcluster1_CY151620_common, rownames(allRNASeq)), 1:15]
n <- 1
m <- 2
PCC <- c()
PCC_all <- list()
for(n in n:c(nrow(subcluster1_CY151620_common)-1)){
  for(m in m:nrow(subcluster1_CY151620_common)){
    PCC <- c(PCC, cor(as.numeric(subcluster1_CY151620_common[n, ]), as.numeric(subcluster1_CY151620_common[m, ]), method = "pearson"))
    m <- m+1
  }
  PCC_all <- c(PCC_all, list(PCC))
  PCC <- c()
  print(n)
  n <- n+1
  m <- n+1
}
####PCC p-value####
n <- 1
m <- 2
PCC_pvalue <- c()
PCC_pvalue_all <- list()
for(n in n:c(nrow(subcluster1_CY151620_common)-1)){
  for(m in m:nrow(subcluster1_CY151620_common)){
    PCC_pvalue <- c(PCC_pvalue, cor.test(as.numeric(subcluster1_CY151620_common[n, ]), as.numeric(subcluster1_CY151620_common[m, ]), method = "pearson")$p.value)
    m <- m+1
  }
  PCC_pvalue_all <- c(PCC_pvalue_all, list(PCC_pvalue))
  PCC_pvalue <- c()
  print(n)
  n <- n+1
  m <- n+1
}
####PCC q-value####
PCC_qvalue_all <- p.adjust(unlist(PCC_pvalue_all), method = "BH")
#####Cytoscape_format####
#source genes
rowtotal <- nrow(subcluster1_CY151620_common)
m <- 1
source_genes <- c()
for(m in m:c(rowtotal-1)){
  source_genes <- c(source_genes, rep(rownames(subcluster1_CY151620_common)[m], times = c(rowtotal-m)))
  print(m)
  m <- m+1
}
#target genes
n <- 2
target_genes <- c()
for(n in n:rowtotal){
  target_genes <- c(target_genes, rownames(subcluster1_CY151620_common)[n:rowtotal])
  print(n)
  n <- n+1
}
subcluster1_CY151620_common_cytoscape <- data.frame(source_genes = unlist(source_genes), 
                                             interaction_value = 1/abs(unlist(PCC_all)), 
                                             target_genes = unlist(target_genes),
                                             p_value = unlist(PCC_pvalue_all),
                                             q_value = PCC_qvalue_all
)
subcluster1_CY151620_common_cytoscape_th <- subcluster1_CY151620_common_cytoscape[subcluster1_CY151620_common_cytoscape$q_value < 0.05, ]
write.table(subcluster1_CY151620_common_cytoscape_th, "~/Nakano_RNAseq/network_analysis/base/subcluster1_CY151620_common_cytoscape_th.txt", append=F, quote = F, sep = "\t", row.names = F)
