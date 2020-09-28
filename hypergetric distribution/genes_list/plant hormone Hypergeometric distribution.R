####JA list####
JA <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/jasmonic_acid.txt", sep = "\t", header = T)
JA <- as.character(unlist(JA))
JA <- strtrim(JA, 9)
####SA list####
SA <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/salicylic_acid.txt", sep = "\t", header = T)
SA <- as.character(unlist(SA))
SA <- strtrim(SA, 9)
####ethylene####
ET <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/ethlene.txt", sep = "\t", header = T)
ET <- as.character(unlist(ET))
ET <- strtrim(ET, 9)
####abscisic acid list####
ABA <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/abscisic_acid.txt", sep = "\t", header = T)
ABA <- as.character(unlist(ABA))
ABA <- strtrim(ABA, 9)
####auzin list####
AUX <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/auxin.txt", sep = "\t", header = T)
AUX <- as.character(unlist(AUX))
AUX <- strtrim(AUX, 9)
####brassinosteroid list####
BR <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/brassinosteroid.txt", sep = "\t", header = T)
BR <- as.character(unlist(BR))
BR <- strtrim(BR, 9)
####cytokine list####
Cytokinine <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/cytokinine.txt", sep = "\t", header = T)
Cytokinine <- as.character(unlist(Cytokinine))
Cytokinine <- strtrim(Cytokinine, 9)
####gibberellin list####
GA <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/gibberellin.txt", sep = "\t", header = T)
GA <- as.character(unlist(GA))
GA <- strtrim(GA, 9)

iterations14 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/iteration14_node.txt", sep = "\t", header = T)

####plant hormone####
#SA#
#全クラスターのSAのDEGs数
NumSA <- sum(!is.na(match(SA, iterations14$name)))
#全クラスターのSA以外のDEGs数
Numothergenes <- c(length(iterations14$name)-NumSA)
maxiterations14 <- max(iterations14$X__mclCluster)
n <- 1
SA_enrichment_pvalue <- c()
check <- c()
SA_check_MCL <- list()
for(n in n:maxiterations14){
  test <- iterations14[iterations14$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(SA, rownames(test)))
  Numsubcluster <- nrow(test)
  SA_enrichment_pvalue <- c(SA_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumSA, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(SA_enrichment_pvalue)[n] <- n
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(SA, rownames(test))
    names(T_INT) <- n
  }
  SA_check_MCL <- c(SA_check_MCL, list(T_INT))
  print(maxiterations14 - n)
  n <- n + 1
}
SA_enrichment_qvalue <- p.adjust(SA_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
SA_enriched_Node <- SA_check_MCL[SA_enrichment_qvalue < 0.05]
NumSA_subcluster <- c()
i <- 1
for(i in i:length(SA_enriched_Node)){
  NumSA_subcluster <- c(NumSA_subcluster, names(SA_enriched_Node[[i]])[1])
  i <- i+1
}
NumSA_subcluster <- as.numeric(NumSA_subcluster)
#names(SA_enriched_Node[[11]])[1]でクラスター番号とってこれる。
#JA#
#全クラスターのJAのDEGs数
NumJA <- sum(!is.na(match(JA, iterations14$name)))
#全クラスターのJA以外のDEGs数
Numothergenes <- c(length(iterations14$name)-NumJA)
maxiterations14 <- max(iterations14$X__mclCluster)
n <- 1
JA_enrichment_pvalue <- c()
check <- c()
JA_check_MCL <- list()
for(n in n:maxiterations14){
  test <- iterations14[iterations14$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(JA, rownames(test)))
  Numsubcluster <- nrow(test)
  JA_enrichment_pvalue <- c(JA_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumJA, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(JA_enrichment_pvalue)[n] <- n
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(JA, rownames(test))
    names(T_INT) <- n
  }
  JA_check_MCL <- c(JA_check_MCL, list(T_INT))
  print(maxiterations14 - n)
  n <- n + 1
}
JA_enrichment_qvalue <- p.adjust(JA_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
JA_enriched_Node <- JA_check_MCL[JA_enrichment_qvalue < 0.05]
NumJA_subcluster <- c()
i <- 1
for(i in i:length(JA_enriched_Node)){
  NumJA_subcluster <- c(NumJA_subcluster, names(JA_enriched_Node[[i]])[1])
  i <- i+1
}
NumJA_subcluster <- as.numeric(NumJA_subcluster)
#names(JA_enriched_Node[[11]])[1]でクラスター番号とってこれる。
#ET#
#全クラスターのETのDEGs数
NumET <- sum(!is.na(match(ET, iterations14$name)))
#全クラスターのET以外のDEGs数
Numothergenes <- c(length(iterations14$name)-NumET)
maxiterations14 <- max(iterations14$X__mclCluster)
n <- 1
ET_enrichment_pvalue <- c()
check <- c()
ET_check_MCL <- list()
for(n in n:maxiterations14){
  test <- iterations14[iterations14$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(ET, rownames(test)))
  Numsubcluster <- nrow(test)
  ET_enrichment_pvalue <- c(ET_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumET, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(ET_enrichment_pvalue)[n] <- n
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(ET, rownames(test))
    names(T_INT) <- n
  }
  ET_check_MCL <- c(ET_check_MCL, list(T_INT))
  print(maxiterations14 - n)
  n <- n + 1
}
ET_enrichment_qvalue <- p.adjust(ET_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
ET_enriched_Node <- ET_check_MCL[ET_enrichment_qvalue < 0.05]
NumET_subcluster <- c()
i <- 1
for(i in i:length(ET_enriched_Node)){
  NumET_subcluster <- c(NumET_subcluster, names(ET_enriched_Node[[i]])[1])
  i <- i+1
}
NumET_subcluster <- as.numeric(NumET_subcluster)
#names(ET_enriched_Node[[11]])[1]でクラスター番号とってこれる。
#ABA#
#全クラスターのABAのDEGs数
NumABA <- sum(!is.na(match(ABA, iterations14$name)))
#全クラスターのABA以外のDEGs数
Numothergenes <- c(length(iterations14$name)-NumABA)
maxiterations14 <- max(iterations14$X__mclCluster)
n <- 1
ABA_enrichment_pvalue <- c()
check <- c()
ABA_check_MCL <- list()
for(n in n:maxiterations14){
  test <- iterations14[iterations14$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(ABA, rownames(test)))
  Numsubcluster <- nrow(test)
  ABA_enrichment_pvalue <- c(ABA_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumABA, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(ABA_enrichment_pvalue)[n] <- n
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(ABA, rownames(test))
    names(T_INT) <- n
  }
  ABA_check_MCL <- c(ABA_check_MCL, list(T_INT))
  print(maxiterations14 - n)
  n <- n + 1
}
ABA_enrichment_qvalue <- p.adjust(ABA_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
ABA_enriched_Node <- ABA_check_MCL[ABA_enrichment_qvalue < 0.05]
NumABA_subcluster <- c()
i <- 1
for(i in i:length(ABA_enriched_Node)){
  NumABA_subcluster <- c(NumABA_subcluster, names(ABA_enriched_Node[[i]])[1])
  i <- i+1
}
NumABA_subcluster <- as.numeric(NumABA_subcluster)
#names(ABA_enriched_Node[[11]])[1]でクラスター番号とってこれる。
#AUX#
#全クラスターのAUXのDEGs数
NumAUX <- sum(!is.na(match(AUX, iterations14$name)))
#全クラスターのAUX以外のDEGs数
Numothergenes <- c(length(iterations14$name)-NumAUX)
maxiterations14 <- max(iterations14$X__mclCluster)
n <- 1
AUX_enrichment_pvalue <- c()
check <- c()
AUX_check_MCL <- list()
for(n in n:maxiterations14){
  test <- iterations14[iterations14$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(AUX, rownames(test)))
  Numsubcluster <- nrow(test)
  AUX_enrichment_pvalue <- c(AUX_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumAUX, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(AUX_enrichment_pvalue)[n] <- n
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(AUX, rownames(test))
    names(T_INT) <- n
  }
  AUX_check_MCL <- c(AUX_check_MCL, list(T_INT))
  print(maxiterations14 - n)
  n <- n + 1
}
AUX_enrichment_qvalue <- p.adjust(AUX_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
AUX_enriched_Node <- AUX_check_MCL[AUX_enrichment_qvalue < 0.05]
NumAUX_subcluster <- c()
i <- 1
for(i in i:length(AUX_enriched_Node)){
  NumAUX_subcluster <- c(NumAUX_subcluster, names(AUX_enriched_Node[[i]])[1])
  i <- i+1
}
NumAUX_subcluster <- as.numeric(NumAUX_subcluster)
#names(AUX_enriched_Node[[11]])[1]でクラスター番号とってこれる。
#BR#
#全クラスターのBRのDEGs数
NumBR <- sum(!is.na(match(BR, iterations14$name)))
#全クラスターのBR以外のDEGs数
Numothergenes <- c(length(iterations14$name)-NumBR)
maxiterations14 <- max(iterations14$X__mclCluster)
n <- 1
BR_enrichment_pvalue <- c()
check <- c()
BR_check_MCL <- list()
for(n in n:maxiterations14){
  test <- iterations14[iterations14$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(BR, rownames(test)))
  Numsubcluster <- nrow(test)
  BR_enrichment_pvalue <- c(BR_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumBR, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(BR_enrichment_pvalue)[n] <- n
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(BR, rownames(test))
    names(T_INT) <- n
  }
  BR_check_MCL <- c(BR_check_MCL, list(T_INT))
  print(maxiterations14 - n)
  n <- n + 1
}
BR_enrichment_qvalue <- p.adjust(BR_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
BR_enriched_Node <- BR_check_MCL[BR_enrichment_qvalue < 0.05]
NumBR_subcluster <- c()
i <- 1
for(i in i:length(BR_enriched_Node)){
  NumBR_subcluster <- c(NumBR_subcluster, names(BR_enriched_Node[[i]])[1])
  i <- i+1
}
NumBR_subcluster <- as.numeric(NumBR_subcluster)
#names(BR_enriched_Node[[11]])[1]でクラスター番号とってこれる。
#GA#
#全クラスターのGAのDEGs数
NumGA <- sum(!is.na(match(GA, iterations14$name)))
#全クラスターのGA以外のDEGs数
Numothergenes <- c(length(iterations14$name)-NumGA)
maxiterations14 <- max(iterations14$X__mclCluster)
n <- 1
GA_enrichment_pvalue <- c()
check <- c()
GA_check_MCL <- list()
for(n in n:maxiterations14){
  test <- iterations14[iterations14$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(GA, rownames(test)))
  Numsubcluster <- nrow(test)
  GA_enrichment_pvalue <- c(GA_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumGA, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(GA_enrichment_pvalue)[n] <- n
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(GA, rownames(test))
    names(T_INT) <- n
  }
  GA_check_MCL <- c(GA_check_MCL, list(T_INT))
  print(maxiterations14 - n)
  n <- n + 1
}
GA_enrichment_qvalue <- p.adjust(GA_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
GA_enriched_Node <- GA_check_MCL[GA_enrichment_qvalue < 0.05]
NumGA_subcluster <- c()
i <- 1
for(i in i:length(GA_enriched_Node)){
  NumGA_subcluster <- c(NumGA_subcluster, names(GA_enriched_Node[[i]])[1])
  i <- i+1
}
NumGA_subcluster <- as.numeric(NumGA_subcluster)
#names(GA_enriched_Node[[11]])[1]でクラスター番号とってこれる。
#Cytokinine#
#全クラスターのCytokinineのDEGs数
NumCytokinine <- sum(!is.na(match(Cytokinine, iterations14$name)))
#全クラスターのCytokinine以外のDEGs数
Numothergenes <- c(length(iterations14$name)-NumCytokinine)
maxiterations14 <- max(iterations14$X__mclCluster)
n <- 1
Cytokinine_enrichment_pvalue <- c()
check <- c()
Cytokinine_check_MCL <- list()
for(n in n:maxiterations14){
  test <- iterations14[iterations14$X__mclCluster == n, ]
  NumCY_DEGs <- length(intersect(Cytokinine, rownames(test)))
  Numsubcluster <- nrow(test)
  Cytokinine_enrichment_pvalue <- c(Cytokinine_enrichment_pvalue, phyper(c(NumCY_DEGs - 1), NumCytokinine, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(Cytokinine_enrichment_pvalue)[n] <- n
  check <- c(check, NumCY_DEGs)
  if(NumCY_DEGs != 0){
    T_INT <- intersect(Cytokinine, rownames(test))
    names(T_INT) <- n
  }
  Cytokinine_check_MCL <- c(Cytokinine_check_MCL, list(T_INT))
  print(maxiterations14 - n)
  n <- n + 1
}
Cytokinine_enrichment_qvalue <- p.adjust(Cytokinine_enrichment_pvalue, method = "BH")
#各サブクラスターにおいてq-value < 0.05であるAGIを引っ張れている。
Cytokinine_enriched_Node <- Cytokinine_check_MCL[Cytokinine_enrichment_qvalue < 0.05]
NumCytokinine_subcluster <- c()
i <- 1
for(i in i:length(Cytokinine_enriched_Node)){
  NumCytokinine_subcluster <- c(NumCytokinine_subcluster, names(Cytokinine_enriched_Node[[i]])[1])
  i <- i+1
}
NumCytokinine_subcluster <- as.numeric(NumCytokinine_subcluster)
#names(Cytokinine_enriched_Node[[11]])[1]でクラスター番号とってこれる。
