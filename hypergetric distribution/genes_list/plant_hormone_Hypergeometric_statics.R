#write(JA, "~/JA.txt", sep = "\t")
#write(SA, "~/SA.txt", sep = "\t")
#write(ET, "~/ET.txt", sep = "\t")
#write(ABA, "~/ABA.txt", sep = "\t")
#write(AUX, "~/AUX.txt", sep = "\t")
#write(BR, "~/BR.txt", sep = "\t")
#write(GA, "~/GA.txt", sep = "\t")
#write(CK, "~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/new/CK.txt", sep = "\t")
####JA list####
#JA <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/jasmonic_acid.txt", sep = "\t", header = T)
#JA <- as.character(unlist(JA))
#JA <- strtrim(JA, 9)
JA <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/new/JA.txt", sep = "\t", header = T)
JA <- JA$AGI

####SA list####
#SA <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/salicylic_acid.txt", sep = "\t", header = T)
#SA <- as.character(unlist(SA))
#SA <- strtrim(SA, 9)
SA <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/new/SA.txt", sep = "\t", header = T)
SA <- SA$AGI

####ethylene####
#ET <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/ethlene.txt", sep = "\t", header = T)
#ET <- as.character(unlist(ET))
#ET <- strtrim(ET, 9)
ET <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/new/ET.txt", sep = "\t", header = T)
ET <- ET$AGI

####abscisic acid list####
#ABA <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/abscisic_acid.txt", sep = "\t", header = T)
#ABA <- as.character(unlist(ABA))
#ABA <- strtrim(ABA, 9)
ABA <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/new/ABA.txt", sep = "\t", header = T)
ABA <- ABA$AGI


####auxin list####
#AUX <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/auxin.txt", sep = "\t", header = T)
#AUX <- as.character(unlist(AUX))
#AUX <- strtrim(AUX, 9)
AUX <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/new/Auxin.txt", sep = "\t", header = T)
AUX <- AUX$AGI

####brassinosteroid list####
#BR <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/brassinosteroid.txt", sep = "\t", header = T)
#BR <- as.character(unlist(BR))
#BR <- strtrim(BR, 9)
BR <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/new/BR.txt", sep = "\t", header = T)
BR <- BR$AGI

####cytokine list####
#CK <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/cytokinine.txt", sep = "\t", header = T)
#CK <- as.character(unlist(CK))
#CK <- strtrim(CK, 9)
CK <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/new/CK.txt", sep = "\t", header = T)
CK <- CK$AGI

####gibberellin list####
#GA <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/gibberellin.txt", sep = "\t", header = T)
#GA <- as.character(unlist(GA))
#GA <- strtrim(GA, 9)
GA <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/plant_hormone/new/GA.txt", sep = "\t", header = T)
GA <- GA$AGI

####plant hormone####
iterations14 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/iteration14_node.txt", sep = "\t", header = T)

####SA####
####SA enrichment####
#全クラスターのSAのDEGs数
NumSA <- sum(!is.na(match(SA, iterations14$name)))
#全クラスターのSA以外のDEGs数
Numothergenes <- c(length(iterations14$name)-NumSA)
maxiterations14 <- max(iterations14$X__mclCluster)
n <- 1
SA_enrichment_pvalue <- c()
check <- c()
SA_check_MCL <- list()
subnodes <- list()
for(n in n:maxiterations14){
  test <- iterations14[iterations14$X__mclCluster == n, ]
  NumSA_DEGs <- length(intersect(SA, test$name))
  Numsubcluster <- nrow(test)
  SA_enrichment_pvalue <- c(SA_enrichment_pvalue, phyper(c(NumSA_DEGs - 1), NumSA, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(SA_enrichment_pvalue)[n] <- n
  subnodes <- c(subnodes, list(Numsubcluster))
  check <- c(check, NumSA_DEGs)
  T_INT <- intersect(SA, test$name)
  names(T_INT) <- rep(n, times = length(T_INT))
  SA_check_MCL <- c(SA_check_MCL, list(T_INT))
  print(maxiterations14 - n)
  n <- n + 1
}
#各サブクラスターにおいてp-value < 0.05であるAGIを引っ張れている。
data <- c()
AGI <- list()
i <- 1
for(i in i:length(SA_check_MCL)){
  data <- unlist(SA_check_MCL[i])
  data <- paste(data, collapse = " | ")
  data[which(data %in% "")] <- "NA"
  AGI <- c(AGI, list(data))
  data <- c()
  m <- m+1
  i <- i+1
}

SA_statistics <- data.frame(NumSubcluster = NumSubCluster, 
                            NumSA_DEGs = check,
                            SA_DEGs_AGI = unlist(AGI),
                            NumNodes = unlist(subnodes),
                            p_value = SA_enrichment_pvalue, 
                            enrichment_score = -log2(SA_enrichment_pvalue)
)
####JA####
####JA enrichment####
#全クラスターのJAのDEGs数
NumJA <- sum(!is.na(match(JA, iterations14$name)))
#全クラスターのJA以外のDEGs数
Numothergenes <- c(length(iterations14$name)-NumJA)
maxiterations14 <- max(iterations14$X__mclCluster)
n <- 1
JA_enrichment_pvalue <- c()
check <- c()
JA_check_MCL <- list()
subnodes <- list()
for(n in n:maxiterations14){
  test <- iterations14[iterations14$X__mclCluster == n, ]
  NumJA_DEGs <- length(intersect(JA, test$name))
  Numsubcluster <- nrow(test)
  JA_enrichment_pvalue <- c(JA_enrichment_pvalue, phyper(c(NumJA_DEGs - 1), NumJA, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(JA_enrichment_pvalue)[n] <- n
  subnodes <- c(subnodes, list(Numsubcluster))
  check <- c(check, NumJA_DEGs)
  T_INT <- intersect(JA, test$name)
  names(T_INT) <- rep(n, times = length(T_INT))
  JA_check_MCL <- c(JA_check_MCL, list(T_INT))
  print(maxiterations14 - n)
  n <- n + 1
}
#各サブクラスターにおいてp-value < 0.05であるAGIを引っ張れている。
data <- c()
AGI <- list()
i <- 1
for(i in i:length(JA_check_MCL)){
  data <- unlist(JA_check_MCL[i])
  data <- paste(data, collapse = " | ")
  data[which(data %in% "")] <- "NA"
  AGI <- c(AGI, list(data))
  data <- c()
  m <- m+1
  i <- i+1
}

JA_statistics <- data.frame(NumSubcluster = NumSubCluster, 
                            NumJA_DEGs = check,
                            JA_DEGs_AGI = unlist(AGI),
                            NumNodes = unlist(subnodes),
                            p_value = JA_enrichment_pvalue, 
                            enrichment_score = -log2(JA_enrichment_pvalue)
)
####ET####
####ET enrichment####
#全クラスターのETのDEGs数
NumET <- sum(!is.na(match(ET, iterations14$name)))
#全クラスターのET以外のDEGs数
Numothergenes <- c(length(iterations14$name)-NumET)
maxiterations14 <- max(iterations14$X__mclCluster)
n <- 1
ET_enrichment_pvalue <- c()
check <- c()
ET_check_MCL <- list()
subnodes <- list()
for(n in n:maxiterations14){
  test <- iterations14[iterations14$X__mclCluster == n, ]
  NumET_DEGs <- length(intersect(ET, test$name))
  Numsubcluster <- nrow(test)
  ET_enrichment_pvalue <- c(ET_enrichment_pvalue, phyper(c(NumET_DEGs - 1), NumET, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(ET_enrichment_pvalue)[n] <- n
  subnodes <- c(subnodes, list(Numsubcluster))
  check <- c(check, NumET_DEGs)
  T_INT <- intersect(ET, test$name)
  names(T_INT) <- rep(n, times = length(T_INT))
  ET_check_MCL <- c(ET_check_MCL, list(T_INT))
  print(maxiterations14 - n)
  n <- n + 1
}
#各サブクラスターにおいてp-value < 0.05であるAGIを引っ張れている。
data <- c()
AGI <- list()
i <- 1
for(i in i:length(ET_check_MCL)){
  data <- unlist(ET_check_MCL[i])
  data <- paste(data, collapse = " | ")
  data[which(data %in% "")] <- "NA"
  AGI <- c(AGI, list(data))
  data <- c()
  m <- m+1
  i <- i+1
}

ET_statistics <- data.frame(NumSubcluster = NumSubCluster, 
                            NumET_DEGs = check,
                            ET_DEGs_AGI = unlist(AGI),
                            NumNodes = unlist(subnodes),
                            p_value = ET_enrichment_pvalue, 
                            enrichment_score = -log2(ET_enrichment_pvalue)
)
####ABA####
####ABA enrichment####
#全クラスターのABAのDEGs数
NumABA <- sum(!is.na(match(ABA, iterations14$name)))
#全クラスターのABA以外のDEGs数
Numothergenes <- c(length(iterations14$name)-NumABA)
maxiterations14 <- max(iterations14$X__mclCluster)
n <- 1
ABA_enrichment_pvalue <- c()
check <- c()
ABA_check_MCL <- list()
subnodes <- list()
for(n in n:maxiterations14){
  test <- iterations14[iterations14$X__mclCluster == n, ]
  NumABA_DEGs <- length(intersect(ABA, test$name))
  Numsubcluster <- nrow(test)
  ABA_enrichment_pvalue <- c(ABA_enrichment_pvalue, phyper(c(NumABA_DEGs - 1), NumABA, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(ABA_enrichment_pvalue)[n] <- n
  subnodes <- c(subnodes, list(Numsubcluster))
  check <- c(check, NumABA_DEGs)
  T_INT <- intersect(ABA, test$name)
  names(T_INT) <- rep(n, times = length(T_INT))
  ABA_check_MCL <- c(ABA_check_MCL, list(T_INT))
  print(maxiterations14 - n)
  n <- n + 1
}
#各サブクラスターにおいてp-value < 0.05であるAGIを引っ張れている。
data <- c()
AGI <- list()
i <- 1
for(i in i:length(ABA_check_MCL)){
  data <- unlist(ABA_check_MCL[i])
  data <- paste(data, collapse = " | ")
  data[which(data %in% "")] <- "NA"
  AGI <- c(AGI, list(data))
  data <- c()
  m <- m+1
  i <- i+1
}

ABA_statistics <- data.frame(NumSubcluster = NumSubCluster, 
                             NumABA_DEGs = check,
                             ABA_DEGs_AGI = unlist(AGI),
                             NumNodes = unlist(subnodes),
                             p_value = ABA_enrichment_pvalue, 
                             enrichment_score = -log2(ABA_enrichment_pvalue)
)
####AUX####
####AUX enrichment####
#全クラスターのAUXのDEGs数
NumAUX <- sum(!is.na(match(AUX, iterations14$name)))
#全クラスターのAUX以外のDEGs数
Numothergenes <- c(length(iterations14$name)-NumAUX)
maxiterations14 <- max(iterations14$X__mclCluster)
n <- 1
AUX_enrichment_pvalue <- c()
check <- c()
AUX_check_MCL <- list()
subnodes <- list()
for(n in n:maxiterations14){
  test <- iterations14[iterations14$X__mclCluster == n, ]
  NumAUX_DEGs <- length(intersect(AUX, test$name))
  Numsubcluster <- nrow(test)
  AUX_enrichment_pvalue <- c(AUX_enrichment_pvalue, phyper(c(NumAUX_DEGs - 1), NumAUX, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(AUX_enrichment_pvalue)[n] <- n
  subnodes <- c(subnodes, list(Numsubcluster))
  check <- c(check, NumAUX_DEGs)
  T_INT <- intersect(AUX, test$name)
  names(T_INT) <- rep(n, times = length(T_INT))
  AUX_check_MCL <- c(AUX_check_MCL, list(T_INT))
  print(maxiterations14 - n)
  n <- n + 1
}
#各サブクラスターにおいてp-value < 0.05であるAGIを引っ張れている。
data <- c()
AGI <- list()
i <- 1
for(i in i:length(AUX_check_MCL)){
  data <- unlist(AUX_check_MCL[i])
  data <- paste(data, collapse = " | ")
  data[which(data %in% "")] <- "NA"
  AGI <- c(AGI, list(data))
  data <- c()
  m <- m+1
  i <- i+1
}

AUX_statistics <- data.frame(NumSubcluster = NumSubCluster, 
                             NumAUX_DEGs = check,
                             AUX_DEGs_AGI = unlist(AGI),
                             NumNodes = unlist(subnodes),
                             p_value = AUX_enrichment_pvalue, 
                             enrichment_score = -log2(AUX_enrichment_pvalue)
)
####BR####
####BR enrichment####
#全クラスターのBRのDEGs数
NumBR <- sum(!is.na(match(BR, iterations14$name)))
#全クラスターのBR以外のDEGs数
Numothergenes <- c(length(iterations14$name)-NumBR)
maxiterations14 <- max(iterations14$X__mclCluster)
n <- 1
BR_enrichment_pvalue <- c()
check <- c()
BR_check_MCL <- list()
subnodes <- list()
for(n in n:maxiterations14){
  test <- iterations14[iterations14$X__mclCluster == n, ]
  NumBR_DEGs <- length(intersect(BR, test$name))
  Numsubcluster <- nrow(test)
  BR_enrichment_pvalue <- c(BR_enrichment_pvalue, phyper(c(NumBR_DEGs - 1), NumBR, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(BR_enrichment_pvalue)[n] <- n
  subnodes <- c(subnodes, list(Numsubcluster))
  check <- c(check, NumBR_DEGs)
  T_INT <- intersect(BR, test$name)
  names(T_INT) <- rep(n, times = length(T_INT))
  BR_check_MCL <- c(BR_check_MCL, list(T_INT))
  print(maxiterations14 - n)
  n <- n + 1
}
#各サブクラスターにおいてp-value < 0.05であるAGIを引っ張れている。
data <- c()
AGI <- list()
i <- 1
for(i in i:length(BR_check_MCL)){
  data <- unlist(BR_check_MCL[i])
  data <- paste(data, collapse = " | ")
  data[which(data %in% "")] <- "NA"
  AGI <- c(AGI, list(data))
  data <- c()
  m <- m+1
  i <- i+1
}

BR_statistics <- data.frame(NumSubcluster = NumSubCluster, 
                            NumBR_DEGs = check,
                            BR_DEGs_AGI = unlist(AGI),
                            NumNodes = unlist(subnodes),
                            p_value = BR_enrichment_pvalue, 
                            enrichment_score = -log2(BR_enrichment_pvalue)
)
####GA####
####GA enrichment####
#全クラスターのGAのDEGs数
NumGA <- sum(!is.na(match(GA, iterations14$name)))
#全クラスターのGA以外のDEGs数
Numothergenes <- c(length(iterations14$name)-NumGA)
maxiterations14 <- max(iterations14$X__mclCluster)
n <- 1
GA_enrichment_pvalue <- c()
check <- c()
GA_check_MCL <- list()
subnodes <- list()
for(n in n:maxiterations14){
  test <- iterations14[iterations14$X__mclCluster == n, ]
  NumGA_DEGs <- length(intersect(GA, test$name))
  Numsubcluster <- nrow(test)
  GA_enrichment_pvalue <- c(GA_enrichment_pvalue, phyper(c(NumGA_DEGs - 1), NumGA, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(GA_enrichment_pvalue)[n] <- n
  subnodes <- c(subnodes, list(Numsubcluster))
  check <- c(check, NumGA_DEGs)
  T_INT <- intersect(GA, test$name)
  names(T_INT) <- rep(n, times = length(T_INT))
  GA_check_MCL <- c(GA_check_MCL, list(T_INT))
  print(maxiterations14 - n)
  n <- n + 1
}
#各サブクラスターにおいてp-value < 0.05であるAGIを引っ張れている。
data <- c()
AGI <- list()
i <- 1
for(i in i:length(GA_check_MCL)){
  data <- unlist(GA_check_MCL[i])
  data <- paste(data, collapse = " | ")
  data[which(data %in% "")] <- "NA"
  AGI <- c(AGI, list(data))
  data <- c()
  m <- m+1
  i <- i+1
}

GA_statistics <- data.frame(NumSubcluster = NumSubCluster, 
                            NumGA_DEGs = check,
                            GA_DEGs_AGI = unlist(AGI),
                            NumNodes = unlist(subnodes),
                            p_value = GA_enrichment_pvalue, 
                            enrichment_score = -log2(GA_enrichment_pvalue)
)
####CK####
####CK enrichment####
#全クラスターのCKのDEGs数
NumCK <- sum(!is.na(match(CK, iterations14$name)))
#全クラスターのCK以外のDEGs数
Numothergenes <- c(length(iterations14$name)-NumCK)
maxiterations14 <- max(iterations14$X__mclCluster)
n <- 1
CK_enrichment_pvalue <- c()
check <- c()
CK_check_MCL <- list()
subnodes <- list()
for(n in n:maxiterations14){
  test <- iterations14[iterations14$X__mclCluster == n, ]
  NumCK_DEGs <- length(intersect(CK, test$name))
  Numsubcluster <- nrow(test)
  CK_enrichment_pvalue <- c(CK_enrichment_pvalue, phyper(c(NumCK_DEGs - 1), NumCK, Numothergenes, Numsubcluster, lower.tail = FALSE))
  names(CK_enrichment_pvalue)[n] <- n
  subnodes <- c(subnodes, list(Numsubcluster))
  check <- c(check, NumCK_DEGs)
  T_INT <- intersect(CK, test$name)
  names(T_INT) <- rep(n, times = length(T_INT))
  CK_check_MCL <- c(CK_check_MCL, list(T_INT))
  print(maxiterations14 - n)
  n <- n + 1
}
#各サブクラスターにおいてp-value < 0.05であるAGIを引っ張れている。
data <- c()
AGI <- list()
i <- 1
for(i in i:length(CK_check_MCL)){
  data <- unlist(CK_check_MCL[i])
  data <- paste(data, collapse = " | ")
  data[which(data %in% "")] <- "NA"
  AGI <- c(AGI, list(data))
  data <- c()
  m <- m+1
  i <- i+1
}

CK_statistics <- data.frame(NumSubcluster = NumSubCluster, 
                            NumCK_DEGs = check,
                            CK_DEGs_AGI = unlist(AGI),
                            NumNodes = unlist(subnodes),
                            p_value = CK_enrichment_pvalue, 
                            enrichment_score = -log2(CK_enrichment_pvalue)
                            )
###########################################################################
planthormone_statics <- data.frame(SubCluster = rep(T_all, times = 8),
                                   planthormone = rep(c("SA", "JA", "ET", "ABA", "AUX", "BR", "GA", "CK"), each = 8*length(T_all)),
                                   p_value = c(SA_statistics$p_value,
                                               JA_statistics$p_value,
                                               ET_statistics$p_value,
                                               ABA_statistics$p_value,
                                               AUX_statistics$p_value,
                                               BR_statistics$p_value,
                                               GA_statistics$p_value,
                                               CK_statistics$p_value
                                               )
                                   )
library(reshape2)
temp <- dcast(planthormone_statics, SubCluster ~ planthormone)
temp$SA <- SA_statistics$p_value
temp$JA <- JA_statistics$p_value
temp$ET <- ET_statistics$p_value
temp$ABA <- ABA_statistics$p_value
temp$AUX <- AUX_statistics$p_value
temp$BR <- BR_statistics$p_value
temp$GA <- GA_statistics$p_value
temp$CK <- CK_statistics$p_value
test <- temp[, 2:9] < 0.05
temp <- temp[apply(test, 1, sum) != 0, ]
planthormone_statics <- melt(temp)
colnames(planthormone_statics) <- c("SubCluster", "planthormone", "p_value")
planthormone_statics <- data.frame(planthormone_statics, 
                                     log2 = -log2(planthormone_statics$p_value),
                                     log10 = -log10(planthormone_statics$p_value)
)

test <- planthormone_statics
test$p_value[test$p_value >= 0.05] <- 1.0
test$p_value[test$p_value < 0.05 & test$p_value >= 0.01] <- 0.8
test$p_value[test$p_value < 1e-2 & test$p_value >= 1e-4] <- 0.4
test$p_value[test$p_value < 1e-4] <- 0
g <- ggplot(test, aes(planthormone, SubCluster))+ 
  geom_tile(aes(fill = p_value), color = "white") + 
  scale_fill_gradient(low="red",high="white") +
  coord_flip()
g <- g + labs(x = "",y = "") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +  theme(legend.position = "none",axis.ticks = element_blank(), axis.text.x = element_text(size = 40, angle = 280, hjust = 0, colour = "grey50"))
g <- g + theme(axis.text=element_text(size=40), axis.title=element_text(size=40))
g <- g + theme(panel.grid.minor = element_line(colour = "brack"))
plot(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/geneslist_heatmap/iterations14/planthormone.png", plot = g, dpi = 100, width = 20.8, height = 7.2)
