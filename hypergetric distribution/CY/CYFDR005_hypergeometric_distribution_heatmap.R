NumSubCluster <- c()
i <- 1
for(i in i:total){
  NumSubCluster <- c(NumSubCluster, paste0("Subcluster", i))
  i <- i+1
}
iterations14 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/iteration14_node.txt", sep = "\t", header = T, stringsAsFactors = F)
rownames(iterations14) <- iterations14$name
####CY15 enrichment####
CY15_DEGs_FDR005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/union/FDR005/CY15_FDR005_time.txt", sep = "\t", row.names = 1, header = T, stringsAsFactors = F)
#全クラスターのCY15のDEGs数
N <- nrow(allRNASeq)
#全クラスターのCY15以外のDEGs数
M <- length(intersect(rownames(allRNASeq), rownames(CY15_DEGs_FDR005)))
i <- 1
CY15_enrichment_pvalue <- c()
check <- c()
CY_DEGs_AGI <- list()
subnodes <- list()
total <- 10
for(i in i:total){
  subcluster_node <- iterations14$name[iterations14$X__mclCluster == i]
  n <- length(subcluster_node)
  overlap_AGI <- intersect(subcluster_node, rownames(CY15_DEGs_FDR005))
  x <- length(overlap_AGI)
  CY15_enrichment_pvalue <- c(CY15_enrichment_pvalue, phyper(x-1, M, N-M, n, lower.tail = F))
  names(CY15_enrichment_pvalue)[i] <- i
  subnodes <- c(subnodes, list(n))
  CY_DEGs_AGI <- c(CY_DEGs_AGI, list(overlap_AGI))
  check <- c(check, x)
  print(total - i)
  i <- i + 1
}
data <- c()
AGI <- list()
i <- 1
for(i in i:length(CY_DEGs_AGI)){
  data <- unlist(CY_DEGs_AGI[i])
  data <- paste(data, collapse = " | ")
  data[which(data %in% "")] <- "NA"
  AGI <- c(AGI, list(data))
  data <- c()
  m <- m+1
  i <- i+1
}


CY15_statistics <- data.frame(NumSubcluster = NumSubCluster,
                              NumCY_DEGs = check, 
                              CY_DEGs_AGI = unlist(AGI),
                              NumSubcluster = unlist(subnodes), 
                              p_value = CY15_enrichment_pvalue, 
                              enrichment_score = -log10(CY15_enrichment_pvalue),
                              stringsAsFactors = F)
####CY16 enrichment####
####CY16 enrichment####
CY16_DEGs_FDR005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/union/FDR005/CY16_FDR005_time.txt", sep = "\t", row.names = 1, header = T, stringsAsFactors = F)
#全クラスターのCY16のDEGs数
N <- nrow(allRNASeq)
#全クラスターのCY16以外のDEGs数
M <- length(intersect(rownames(allRNASeq), rownames(CY16_DEGs_FDR005)))
i <- 1
CY16_enrichment_pvalue <- c()
check <- c()
CY_DEGs_AGI <- list()
subnodes <- list()
total <- 10
for(i in i:total){
  subcluster_node <- iterations14$name[iterations14$X__mclCluster == i]
  n <- length(subcluster_node)
  overlap_AGI <- intersect(subcluster_node, rownames(CY16_DEGs_FDR005))
  x <- length(overlap_AGI)
  CY16_enrichment_pvalue <- c(CY16_enrichment_pvalue, phyper(x-1, M, N-M, n, lower.tail = F))
  names(CY16_enrichment_pvalue)[i] <- i
  subnodes <- c(subnodes, list(n))
  CY_DEGs_AGI <- c(CY_DEGs_AGI, list(overlap_AGI))
  check <- c(check, x)
  print(total - i)
  i <- i + 1
}
data <- c()
AGI <- list()
i <- 1
for(i in i:length(CY_DEGs_AGI)){
  data <- unlist(CY_DEGs_AGI[i])
  data <- paste(data, collapse = " | ")
  data[which(data %in% "")] <- "NA"
  AGI <- c(AGI, list(data))
  data <- c()
  m <- m+1
  i <- i+1
}


CY16_statistics <- data.frame(NumSubcluster = NumSubCluster,
                              NumCY_DEGs = check, 
                              CY_DEGs_AGI = unlist(AGI),
                              NumSubcluster = unlist(subnodes), 
                              p_value = CY16_enrichment_pvalue, 
                              enrichment_score = -log10(CY16_enrichment_pvalue),
                              stringsAsFactors = F)
####CY20 enrichment####
CY20_DEGs_FDR005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/union/FDR005/CY20_FDR005_time.txt", sep = "\t", row.names = 1, header = T, stringsAsFactors = F)
#全クラスターのCY20のDEGs数
N <- nrow(allRNASeq)
#全クラスターのCY20以外のDEGs数
M <- length(intersect(rownames(allRNASeq), rownames(CY20_DEGs_FDR005)))
i <- 1
CY20_enrichment_pvalue <- c()
check <- c()
CY_DEGs_AGI <- list()
subnodes <- list()
total <- 10
for(i in i:total){
  subcluster_node <- iterations14$name[iterations14$X__mclCluster == i]
  n <- length(subcluster_node)
  overlap_AGI <- intersect(subcluster_node, rownames(CY20_DEGs_FDR005))
  x <- length(overlap_AGI)
  CY20_enrichment_pvalue <- c(CY20_enrichment_pvalue, phyper(x-1, M, N-M, n, lower.tail = F))
  names(CY20_enrichment_pvalue)[i] <- i
  subnodes <- c(subnodes, list(n))
  CY_DEGs_AGI <- c(CY_DEGs_AGI, list(overlap_AGI))
  check <- c(check, x)
  print(total - i)
  i <- i + 1
}
data <- c()
AGI <- list()
i <- 1
for(i in i:length(CY_DEGs_AGI)){
  data <- unlist(CY_DEGs_AGI[i])
  data <- paste(data, collapse = " | ")
  data[which(data %in% "")] <- "NA"
  AGI <- c(AGI, list(data))
  data <- c()
  m <- m+1
  i <- i+1
}


CY20_statistics <- data.frame(NumSubcluster = NumSubCluster,
                              NumCY_DEGs = check, 
                              CY_DEGs_AGI = unlist(AGI),
                              NumSubcluster = unlist(subnodes), 
                              p_value = CY20_enrichment_pvalue, 
                              enrichment_score = -log10(CY20_enrichment_pvalue),
                              stringsAsFactors = F)
####BTH enrichment####
BTH <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/BTH_affected_AGI.txt", sep = "\t", header = T, stringsAsFactors = F)
rownames(BTH) <- BTH$AGI
N <- nrow(allRNASeq)
#全クラスターのBTH以外のDEGs数
M <- length(intersect(rownames(allRNASeq), rownames(BTH)))
i <- 1
BTH_enrichment_pvalue <- c()
check <- c()
BTH_DEGs_AGI <- list()
subnodes <- list()
total <- 10
for(i in i:total){
  subcluster_node <- iterations14$name[iterations14$X__mclCluster == i]
  n <- length(subcluster_node)
  overlap_AGI <- intersect(subcluster_node, rownames(BTH))
  x <- length(overlap_AGI)
  BTH_enrichment_pvalue <- c(BTH_enrichment_pvalue, phyper(x-1, M, N-M, n, lower.tail = F))
  names(BTH_enrichment_pvalue)[i] <- i
  subnodes <- c(subnodes, list(n))
  BTH_DEGs_AGI <- c(BTH_DEGs_AGI, list(overlap_AGI))
  check <- c(check, x)
  print(total - i)
  i <- i + 1
}
data <- c()
AGI <- list()
i <- 1
for(i in i:length(BTH_DEGs_AGI)){
  data <- unlist(BTH_DEGs_AGI[i])
  data <- paste(data, collapse = " | ")
  data[which(data %in% "")] <- "NA"
  AGI <- c(AGI, list(data))
  data <- c()
  m <- m+1
  i <- i+1
}


BTH_statistics <- data.frame(NumSubcluster = NumSubCluster,
                             NumBTH_DEGs = check, 
                             BTH_DEGs_AGI = unlist(AGI),
                             NumSubcluster = unlist(subnodes), 
                             p_value = BTH_enrichment_pvalue, 
                             enrichment_score = -log10(BTH_enrichment_pvalue),
                             stringsAsFactors = F)

####heat map####
n <- 1
total <- 10
T_all <- c()
for(n in n:total){
  temp <- "0000"
  temp <- paste0(substr(temp, 1, nchar(temp)-nchar(n)), n)
  T_all <- c(T_all, temp)
}

library(ggplot2)
plantactivator_statics <- data.frame(SubCluster = rep(T_all, times = 4),
                                     plantactivators = rep(c("BTH", "CY15", "CY16", "CY20"), each = 4*nrow(CY15_statistics)),
                                     p_value = c(BTH_statistics$p_value, CY15_statistics$p_value, CY16_statistics$p_value, CY20_statistics$p_value)
)
library(reshape2)
temp <- dcast(plantactivator_statics, SubCluster ~ plantactivators)
temp$BTH <- BTH_statistics$p_value
temp$CY15 <- CY15_statistics$p_value
temp$CY16 <- CY16_statistics$p_value
temp$CY20 <- CY20_statistics$p_value

test <- temp[, 2:5] < 5e-2
temp <- temp[apply(test, 1, sum) != 0, ]
plantactivator_statics <- melt(temp)

colnames(plantactivator_statics) <- c("SubCluster", "plantactivators", "p_value")
plantactivator_statics <- data.frame(plantactivator_statics, 
                                     log2 = -log2(plantactivator_statics$p_value),
                                     log10 = -log10(plantactivator_statics$p_value)
)

test <- plantactivator_statics
test$p_value[test$p_value >= 5e-2] <- 1.0
test$p_value[test$p_value < 5e-2 & test$p_value >= 5e-5] <- 0.8
test$p_value[test$p_value < 5e-5 & test$p_value >= 5e-10] <- 0.6
test$p_value[test$p_value < 5e-10 & test$p_value >= 5e-20] <- 0.4
test$p_value[test$p_value < 5e-20 & test$p_value >= 5e-40] <- 0.15
test$p_value[test$p_value < 5e-40] <- 0
ghm <- ggplot(test, aes(plantactivators, SubCluster))+
  geom_tile(aes(fill = p_value), color = "white") +
  scale_fill_gradient(low="red",high="white") +
  coord_flip()
ghm <- ghm + labs(x = "",y = "") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +  theme(legend.position = "none",axis.ticks = element_blank(), axis.text.x = element_text(size = 60, angle = 280, hjust = 0, colour = "grey50"))
ghm <- ghm + theme(axis.text=element_text(size=60), axis.title=element_text(size=60))
plot(ghm)
ggsave(file = "~/Nakano_RNAseq/network_analysis/geneslist_heatmap/iterations14/plantactivators_FDR005.png", plot = ghm, dpi = 100, width = 40, height = 16.2)

####output####
write.table(CY20_statistics, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/iterations14/CY_subcluster/CY20_statistics.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
write.table(CY16_statistics, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/iterations14/CY_subcluster/CY16_statistics.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
write.table(CY15_statistics, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/iterations14/CY_subcluster/CY15_statistics.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
write.table(BTH_statistics, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/iterations14/genes_set/BTH_statics.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)