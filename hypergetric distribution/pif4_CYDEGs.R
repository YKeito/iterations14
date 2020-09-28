CY15_DEGs_FDR005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/union/FDR005/CY15_FDR005_time.txt", sep = "\t", row.names = 1, header = T, stringsAsFactors = F)
CY16_DEGs_FDR005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/union/FDR005/CY16_FDR005_time.txt", sep = "\t", row.names = 1, header = T, stringsAsFactors = F)
CY20_DEGs_FDR005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/union/FDR005/CY20_FDR005_time.txt", sep = "\t", row.names = 1, header = T, stringsAsFactors = F)
iterations14 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/iteration14_node.txt", sep = "\t", header = T, stringsAsFactors = F)
rownames(iterations14) <- iterations14$name
####pif4 enrichment####
pif4_DEGs <- read.table("~/pif4_DEGs.txt", sep = "\t", header = T, stringsAsFactors = F)
rownames(pif4_DEGs) <- pif4_DEGs$pif4_DEGs
#CY15
N <- 33566
M <- nrow(pif4_DEGs)
n <- nrow(CY15_DEGs_FDR005)
x <- length(intersect(rownames(CY15_DEGs_FDR005), rownames(pif4_DEGs)))
phyper(c(x-1), M, c(N - M), n, lower.tail = F)
#CY16
N <- 33566
M <- nrow(pif4_DEGs)
n <- nrow(CY16_DEGs_FDR005)
x <- length(intersect(rownames(CY16_DEGs_FDR005), rownames(pif4_DEGs)))
phyper(c(x-1), M, c(N - M), n, lower.tail = F)
#CY20
N <- 33566
M <- nrow(pif4_DEGs)
n <- nrow(CY20_DEGs_FDR005)
x <- length(intersect(rownames(CY20_DEGs_FDR005), rownames(pif4_DEGs)))
phyper(c(x-1), M, c(N - M), n, lower.tail = F)

#subcluster1
N <- nrow(allRNASeq)
M <- nrow(pif4_DEGs)
n <- length(iterations14$name[iterations14$X__mclCluster == 1])
x <- length(intersect(iterations14$name[iterations14$X__mclCluster == 1], rownames(pif4_DEGs)))
phyper(c(x-1), M, c(N - M), n, lower.tail = F)
