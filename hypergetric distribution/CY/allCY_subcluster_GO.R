load("~/bigdata/yasue/GO_analysis/GO_List.RData")
####CY15####
i <- 1
total <- length(CY15_enriched_Node)
T_GO_Full <- unlist(GO_List)
EnrichedGO_CY15 <- list()
EnrichedGO_pvalue_CY15 <- list()

CY15_Node <- list()
t <- proc.time()
for (i in 1:total) {
  temp_all_GO <- c()
  temp_all_pV <- c()
  T_AGI <- unlist(CY15_enriched_Node[[i]])
  T_GO <- unlist(GO_List[match(T_AGI, names(GO_List))])
  k <- length(T_GO)
  T_GO_Uni <- unique(T_GO)
  
  j <- 1
  total_j <- length(T_GO_Uni)
  for (j in j:total_j) {
    AGI_name <- c()
    q <- sum(T_GO == T_GO_Uni[j])
    m <- sum(T_GO_Full == T_GO_Uni[j])
    n <- length(T_GO_Full) - m
    
    temp <- phyper(q-1, m, n, k, lower.tail=FALSE)
    temp_all_pV <- c(temp_all_pV, temp)
    temp_all_GO <- c(temp_all_GO, T_GO_Uni[j])
    j <- j+1
  }
  temp_CY15 <- GO_List[match(T_AGI, names(GO_List))]
  AGI_name <- c(AGI_name, intersect(T_AGI, names(temp_CY15)))
  names(AGI_name) <- rep(NumCY15_subcluster[i], times = length(AGI_name))
  CY15_Node <- c(CY15_Node, list(AGI_name))
  AGI_name <- c()
  EnrichedGO_pvalue_CY15 <- c(EnrichedGO_pvalue_CY15, list(temp_all_pV))
  names(EnrichedGO_pvalue_CY15)[i] <- NumCY15_subcluster[i]
  EnrichedGO_CY15 <- c(EnrichedGO_CY15, list(temp_all_GO))
  names(EnrichedGO_CY15)[i] <- NumCY15_subcluster[i]
  print(total - i)
  i <- i+1
}
proc.time()-t
save(CY15_Node, file = "~/Nakano_RNAseq/network_analysis/.RData/iterations14/CY15_Node.RData")
####CY16####
i <- 1
total <- length(CY16_enriched_Node)
T_GO_Full <- unlist(GO_List)
EnrichedGO_CY16 <- list()
EnrichedGO_pvalue_CY16 <- list()

CY16_Node <- list()
t <- proc.time()

#各サブクラスターにおいてp-value < 0.05であるAGIを引っ張れている。
CY16_enriched_Node <- CY16_check_MCL[CY16_enrichment_pvalue < 0.05]
NumCY16_subcluster <- c()
i <- 1
for(i in i:length(CY16_enriched_Node)){
  NumCY16_subcluster <- c(NumCY16_subcluster, names(CY16_enriched_Node[[i]])[1])
  i <- i+1
}
NumCY16_subcluster <- as.numeric(NumCY16_subcluster)
for (i in 1:total) {
  temp_all_GO <- c()
  temp_all_pV <- c()
  T_AGI <- unlist(CY16_enriched_Node[[i]])
  T_GO <- unlist(GO_List[match(T_AGI, names(GO_List))])
  k <- length(T_GO)
  T_GO_Uni <- unique(T_GO)
  
  j <- 1
  total_j <- length(T_GO_Uni)
  for (j in j:total_j) {
    AGI_name <- c()
    q <- sum(T_GO == T_GO_Uni[j])
    m <- sum(T_GO_Full == T_GO_Uni[j])
    n <- length(T_GO_Full) - m
    
    temp <- phyper(q-1, m, n, k, lower.tail=FALSE)
    temp_all_pV <- c(temp_all_pV, temp)
    temp_all_GO <- c(temp_all_GO, T_GO_Uni[j])
    j <- j+1
  }
  temp_CY16 <- GO_List[match(T_AGI, names(GO_List))]
  AGI_name <- c(AGI_name, intersect(T_AGI, names(temp_CY16)))
  names(AGI_name) <- rep(NumCY16_subcluster[i], times = length(AGI_name))
  CY16_Node <- c(CY16_Node, list(AGI_name))
  AGI_name <- c()
  EnrichedGO_pvalue_CY16 <- c(EnrichedGO_pvalue_CY16, list(temp_all_pV))
  names(EnrichedGO_pvalue_CY16)[i] <- NumCY16_subcluster[i]
  
  EnrichedGO_CY16 <- c(EnrichedGO_CY16, list(temp_all_GO))
  names(EnrichedGO_CY16)[i] <- NumCY16_subcluster[i]
  
  print(total - i)
  i <- i+1
}
proc.time()-t
save(CY16_Node, file = "~/Nakano_RNAseq/network_analysis/.RData/iterations14/CY16_Node.RData")
####CY20####
i <- 1
total <- length(CY20_enriched_Node)
T_GO_Full <- unlist(GO_List)
EnrichedGO_CY20 <- list()
EnrichedGO_pvalue_CY20 <- list()
CY20_Node <- list()
t <- proc.time()

#各サブクラスターにおいてp-value < 0.05であるAGIを引っ張れている。
CY20_enriched_Node <- CY20_check_MCL[CY20_enrichment_pvalue < 0.05]
NumCY20_subcluster <- c()
i <- 1
for(i in i:length(CY20_enriched_Node)){
  NumCY20_subcluster <- c(NumCY20_subcluster, names(CY20_enriched_Node[[i]])[1])
  i <- i+1
}
NumCY20_subcluster <- as.numeric(NumCY20_subcluster)


for (i in 1:total) {
  temp_all_GO <- c()
  temp_all_pV <- c()
  T_AGI <- unlist(CY20_enriched_Node[[i]])
  T_GO <- unlist(GO_List[match(T_AGI, names(GO_List))])
  k <- length(T_GO)
  T_GO_Uni <- unique(T_GO)
  
  j <- 1
  total_j <- length(T_GO_Uni)
  for (j in j:total_j) {
    AGI_name <- c()
    q <- sum(T_GO == T_GO_Uni[j])
    m <- sum(T_GO_Full == T_GO_Uni[j])
    n <- length(T_GO_Full) - m
    
    temp <- phyper(q-1, m, n, k, lower.tail=FALSE)
    temp_all_pV <- c(temp_all_pV, temp)
    temp_all_GO <- c(temp_all_GO, T_GO_Uni[j])
    j <- j+1
  }
  temp_CY20 <- GO_List[match(T_AGI, names(GO_List))]
  AGI_name <- c(AGI_name, intersect(T_AGI, names(temp_CY20)))
  names(AGI_name) <- rep(NumCY20_subcluster[i], times = length(AGI_name))
  CY20_Node <- c(CY20_Node, list(AGI_name))
  AGI_name <- c()
  EnrichedGO_pvalue_CY20 <- c(EnrichedGO_pvalue_CY20, list(temp_all_pV))
  names(EnrichedGO_pvalue_CY20)[i] <- NumCY20_subcluster[i]
  
  EnrichedGO_CY20 <- c(EnrichedGO_CY20, list(temp_all_GO))
  names(EnrichedGO_CY20)[i] <- NumCY20_subcluster[i]
  
  print(total - i)
  i <- i+1
}
proc.time()-t
save(CY20_Node, file = "~/Nakano_RNAseq/network_analysis/.RData/iterations14/CY20_Node.RData")