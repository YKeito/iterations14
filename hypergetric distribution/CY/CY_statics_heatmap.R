iterations14 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/iteration14_node.txt", sep = "\t", header = T)
rownames(iterations14) <- iterations14$name
####CY15 enrichment####
CY15_DEGs_FDR0005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/union/FDR0005/CY15_FDR0005_time.txt", sep = "\t", row.names = 1, header = T)
#全クラスターのCY15のDEGs数
NumCY15_DEGs_FDR0005 <- sum(!is.na(match(rownames(CY15_DEGs_FDR0005), iterations14$name)))
#全クラスターのCY15以外のDEGs数
Numothergenes <- c(length(iterations14$name)-NumCY15_DEGs_FDR0005)
maxiterations14 <- max(iterations14$X__mclCluster)
n <- 1
CY15_enrichment_pvalue <- c()
check <- c()
CY15_check_MCL <- list()
subnodes <- list()
for(n in n:maxiterations14){
  test <- iterations14[iterations14$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(rownames(CY15_DEGs_FDR0005), rownames(test)))
  Numsubcluster <- nrow(test)
  CY15_enrichment_pvalue <- c(CY15_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumCY15_DEGs_FDR0005, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(CY15_enrichment_pvalue)[n] <- n
  subnodes <- c(subnodes, list(Numsubcluster))
  check <- c(check, NumCY_DEGs)
  T_INT <- intersect(rownames(CY15_DEGs_FDR0005), rownames(test))
  names(T_INT) <- rep(n, times = length(T_INT))
  CY15_check_MCL <- c(CY15_check_MCL, list(T_INT))
  print(maxiterations14 - n)
  n <- n + 1
}
#各サブクラスターにおいてp-value < 0.05であるAGIを引っ張れている。
data <- c()
AGI <- list()
i <- 1
for(i in i:length(CY15_check_MCL)){
  data <- unlist(CY15_check_MCL[i])
  data <- paste(data, collapse = " | ")
  data[which(data %in% "")] <- "NA"
  AGI <- c(AGI, list(data))
  data <- c()
  m <- m+1
  i <- i+1
}
NumCY15_subcluster <- as.numeric(NumCY15_subcluster)

NumSubCluster <- c()
i <- 1
for(i in i:maxiterations14){
  NumSubCluster <- c(NumSubCluster, paste0("SubCluster", i))
  i <- i+1
}
CY15_statistics <- data.frame(NumSubcluster = NumSubCluster,
                              NumCY_DEGs = check, 
                              CY_DEGs_AGI = unlist(AGI),
                              NumSubcluster = unlist(subnodes), 
                              p_value = CY15_enrichment_pvalue, 
                              enrichment_score = -log10(CY15_enrichment_pvalue)
                              )
####CY16 enrichment####
####CY16 enrichment####
CY16_DEGs_FDR0005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/union/FDR0005/CY16_FDR0005_time.txt", sep = "\t", row.names = 1, header = T)
#全クラスターのCY16のDEGs数
NumCY16_DEGs_FDR0005 <- sum(!is.na(match(rownames(CY16_DEGs_FDR0005), iterations14$name)))
#全クラスターのCY16以外のDEGs数
Numothergenes <- c(length(iterations14$name)-NumCY16_DEGs_FDR0005)
maxiterations14 <- max(iterations14$X__mclCluster)
n <- 1
CY16_enrichment_pvalue <- c()
check <- c()
CY16_check_MCL <- list()
subnodes <- list()
for(n in n:maxiterations14){
  test <- iterations14[iterations14$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(rownames(CY16_DEGs_FDR0005), rownames(test)))
  Numsubcluster <- nrow(test)
  CY16_enrichment_pvalue <- c(CY16_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumCY16_DEGs_FDR0005, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(CY16_enrichment_pvalue)[n] <- n
  subnodes <- c(subnodes, list(Numsubcluster))
  check <- c(check, NumCY_DEGs)
  T_INT <- intersect(rownames(CY16_DEGs_FDR0005), rownames(test))
  names(T_INT) <- rep(n, times = length(T_INT))
  CY16_check_MCL <- c(CY16_check_MCL, list(T_INT))
  print(maxiterations14 - n)
  n <- n + 1
}
#各サブクラスターにおいてp-value < 0.05であるAGIを引っ張れている。
data <- c()
AGI <- list()
i <- 1
for(i in i:length(CY16_check_MCL)){
  data <- unlist(CY16_check_MCL[i])
  data <- paste(data, collapse = " | ")
  data[which(data %in% "")] <- "NA"
  AGI <- c(AGI, list(data))
  data <- c()
  m <- m+1
  i <- i+1
}
NumCY16_subcluster <- as.numeric(NumCY16_subcluster)

CY16_statistics <- data.frame(NumSubcluster = NumSubCluster,
                              NumCY_DEGs = check, 
                              CY_DEGs_AGI = unlist(AGI),
                              NumSubcluster = unlist(subnodes), 
                              p_value = CY16_enrichment_pvalue, 
                              enrichment_score = -log10(CY16_enrichment_pvalue)
                              )
####CY20 enrichment####
CY20_DEGs_FDR0005 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/union/FDR0005/CY20_FDR0005_time.txt", sep = "\t", row.names = 1, header = T)
#全クラスターのCY20のDEGs数
NumCY20_DEGs_FDR0005 <- sum(!is.na(match(rownames(CY20_DEGs_FDR0005), iterations14$name)))
#全クラスターのCY20以外のDEGs数
Numothergenes <- c(length(iterations14$name)-NumCY20_DEGs_FDR0005)
maxiterations14 <- max(iterations14$X__mclCluster)
n <- 1
CY20_enrichment_pvalue <- c()
check <- c()
CY20_check_MCL <- list()
subnodes <- list()
for(n in n:maxiterations14){
  test <- iterations14[iterations14$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(rownames(CY20_DEGs_FDR0005), rownames(test)))
  Numsubcluster <- nrow(test)
  CY20_enrichment_pvalue <- c(CY20_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumCY20_DEGs_FDR0005, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(CY20_enrichment_pvalue)[n] <- n
  subnodes <- c(subnodes, list(Numsubcluster))
  check <- c(check, NumCY_DEGs)
  T_INT <- intersect(rownames(CY20_DEGs_FDR0005), rownames(test))
  names(T_INT) <- rep(n, times = length(T_INT))
  CY20_check_MCL <- c(CY20_check_MCL, list(T_INT))
  print(maxiterations14 - n)
  n <- n + 1
}
#各サブクラスターにおいてp-value < 0.05であるAGIを引っ張れている。
data <- c()
AGI <- list()
i <- 1
for(i in i:length(CY20_check_MCL)){
  data <- unlist(CY20_check_MCL[i])
  data <- paste(data, collapse = " | ")
  data[which(data %in% "")] <- "NA"
  AGI <- c(AGI, list(data))
  data <- c()
  m <- m+1
  i <- i+1
}
NumCY20_subcluster <- as.numeric(NumCY20_subcluster)

CY20_statistics <- data.frame(NumSubcluster = NumSubCluster, 
                              NumCY_DEGs = check, 
                              CY_DEGs_AGI = unlist(AGI),
                              NumSubcluster = unlist(subnodes), 
                              p_value = CY20_enrichment_pvalue, 
                              enrichment_score = -log10(CY20_enrichment_pvalue)
                              )
####BTH enrichment####
BTH <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/BTH_affected_AGI.txt", sep = "\t", row.names = 1, header = T)
#全クラスターのBTHのDEGs数
NumBTH <- sum(!is.na(match(rownames(BTH), iterations14$name)))
#全クラスターのBTH以外のDEGs数
Numothergenes <- c(length(iterations14$name)-NumBTH)
maxiterations14 <- max(iterations14$X__mclCluster)
n <- 1
BTH_enrichment_pvalue <- c()
check <- c()
BTH_check_MCL <- list()
subnodes <- list()
for(n in n:maxiterations14){
  test <- iterations14[iterations14$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(rownames(BTH), rownames(test)))
  Numsubcluster <- nrow(test)
  BTH_enrichment_pvalue <- c(BTH_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumBTH, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(BTH_enrichment_pvalue)[n] <- n
  subnodes <- c(subnodes, list(Numsubcluster))
  check <- c(check, NumCY_DEGs)
  T_INT <- intersect(rownames(BTH), rownames(test))
  names(T_INT) <- rep(n, times = length(T_INT))
  BTH_check_MCL <- c(BTH_check_MCL, list(T_INT))
  print(maxiterations14 - n)
  n <- n + 1
}
#各サブクラスターにおいてp-value < 0.05であるAGIを引っ張れている。
data <- c()
AGI <- list()
i <- 1
for(i in i:length(BTH_check_MCL)){
  data <- unlist(BTH_check_MCL[i])
  data <- paste(data, collapse = " | ")
  data[which(data %in% "")] <- "NA"
  AGI <- c(AGI, list(data))
  data <- c()
  m <- m+1
  i <- i+1
}
NumBTH_subcluster <- as.numeric(NumBTH_subcluster)

BTH_statistics <- data.frame(NumSubcluster = NumSubCluster, 
                             NumCY_DEGs = check, 
                             CY_DEGs_AGI = unlist(AGI),
                             NumSubcluster = unlist(subnodes), 
                             p_value = BTH_enrichment_pvalue, 
                             enrichment_score = -log10(BTH_enrichment_pvalue)
                             )


####heat map####
n <- 1
total <- 442
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

test <- temp[, 2:5] < 0.05
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
test$p_value[test$p_value < 5e-10 & test$p_value >= 5e-20] <- 0.3
test$p_value[test$p_value < 5e-20 & test$p_value >= 5e-40] <- 0.2
test$p_value[test$p_value < 5e-40] <- 0.1

ghm <- ggplot(test, aes(plantactivators, SubCluster))+
  geom_tile(aes(fill = p_value), color = "white") +
  scale_fill_gradient(low="red",high="white") +
  coord_flip()
ghm <- ghm + labs(x = "",y = "") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +  theme(legend.position = "none",axis.ticks = element_blank(), axis.text.x = element_text(size = 60, angle = 280, hjust = 0, colour = "grey50"))
ghm <- ghm + theme(axis.text=element_text(size=60), axis.title=element_text(size=60))
plot(ghm)
ggsave(file = "~/Nakano_RNAseq/network_analysis/geneslist_heatmap/iterations14/plantactivators.png", plot = ghm, dpi = 100, width = 40, height = 16.2)

####output####
write.table(CY20_statistics, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/iterations14/CY_subcluster/CY20_statistics.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
write.table(CY16_statistics, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/iterations14/CY_subcluster/CY16_statistics.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
write.table(CY15_statistics, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/iterations14/CY_subcluster/CY15_statistics.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)
write.table(BTH_statistics, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/iterations14/genes_set/BTH_statics.txt",sep="\t", col.names=T, row.names=T, append=F, quote=F)