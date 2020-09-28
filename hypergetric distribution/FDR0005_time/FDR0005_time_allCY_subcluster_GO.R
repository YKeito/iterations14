####CY15_1h####
i <- 1
total <- length(CY15_1h_enriched_Node)
T_GO_Full <- unlist(GO_List)
EnrichedGO_CY15_1h <- list()
EnrichedGO_pvalue_CY15_1h <- list()

CY15_1h_Node <- list()
t <- proc.time()
for (i in 1:total) {
  temp_all_GO <- c()
  temp_all_pV <- c()
  T_AGI <- unlist(CY15_1h_enriched_Node[[i]])
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
  temp_CY15_1h <- GO_List[match(T_AGI, names(GO_List))]
  AGI_name <- c(AGI_name, intersect(T_AGI, names(temp_CY15_1h)))
  names(AGI_name) <- NumCY15_1h_subcluster[i]
  CY15_1h_Node <- c(CY15_1h_Node, list(AGI_name))
  AGI_name <- c()
  EnrichedGO_pvalue_CY15_1h <- c(EnrichedGO_pvalue_CY15_1h, list(temp_all_pV))
  names(EnrichedGO_pvalue_CY15_1h)[i] <- NumCY15_1h_subcluster[i]
  
  EnrichedGO_CY15_1h <- c(EnrichedGO_CY15_1h, list(temp_all_GO))
  names(EnrichedGO_CY15_1h)[i] <- NumCY15_1h_subcluster[i]
  
  print(total - i)
  i <- i+1
}
proc.time()-t
save(CY15_1h_Node, file = "~/Nakano_RNAseq/network_analysis/.RData/CY15_1h_Node.RData")

####CY15_3h####
i <- 1
total <- length(CY15_3h_enriched_Node)
T_GO_Full <- unlist(GO_List)
EnrichedGO_CY15_3h <- list()
EnrichedGO_pvalue_CY15_3h <- list()

CY15_3h_Node <- list()
t <- proc.time()
for (i in 1:total) {
  temp_all_GO <- c()
  temp_all_pV <- c()
  T_AGI <- unlist(CY15_3h_enriched_Node[[i]])
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
  temp_CY15_3h <- GO_List[match(T_AGI, names(GO_List))]
  AGI_name <- c(AGI_name, intersect(T_AGI, names(temp_CY15_3h)))
  names(AGI_name) <- NumCY15_3h_subcluster[i]
  CY15_3h_Node <- c(CY15_3h_Node, list(AGI_name))
  AGI_name <- c()
  EnrichedGO_pvalue_CY15_3h <- c(EnrichedGO_pvalue_CY15_3h, list(temp_all_pV))
  names(EnrichedGO_pvalue_CY15_3h)[i] <- NumCY15_3h_subcluster[i]
  
  EnrichedGO_CY15_3h <- c(EnrichedGO_CY15_3h, list(temp_all_GO))
  names(EnrichedGO_CY15_3h)[i] <- NumCY15_3h_subcluster[i]
  
  print(total - i)
  i <- i+1
}
proc.time()-t
save(CY15_3h_Node, file = "~/Nakano_RNAseq/network_analysis/.RData/CY15_3h_Node.RData")

####CY15_12h####
i <- 1
total <- length(CY15_12h_enriched_Node)
T_GO_Full <- unlist(GO_List)
EnrichedGO_CY15_12h <- list()
EnrichedGO_pvalue_CY15_12h <- list()

CY15_12h_Node <- list()
t <- proc.time()
for (i in 1:total) {
  temp_all_GO <- c()
  temp_all_pV <- c()
  T_AGI <- unlist(CY15_12h_enriched_Node[[i]])
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
  temp_CY15_12h <- GO_List[match(T_AGI, names(GO_List))]
  AGI_name <- c(AGI_name, intersect(T_AGI, names(temp_CY15_12h)))
  names(AGI_name) <- NumCY15_12h_subcluster[i]
  CY15_12h_Node <- c(CY15_12h_Node, list(AGI_name))
  AGI_name <- c()
  EnrichedGO_pvalue_CY15_12h <- c(EnrichedGO_pvalue_CY15_12h, list(temp_all_pV))
  names(EnrichedGO_pvalue_CY15_12h)[i] <- NumCY15_12h_subcluster[i]
  
  EnrichedGO_CY15_12h <- c(EnrichedGO_CY15_12h, list(temp_all_GO))
  names(EnrichedGO_CY15_12h)[i] <- NumCY15_12h_subcluster[i]
  
  print(total - i)
  i <- i+1
}
proc.time()-t
save(CY15_12h_Node, file = "~/Nakano_RNAseq/network_analysis/.RData/CY15_12h_Node.RData")

####CY15_24h####
i <- 1
total <- length(CY15_24h_enriched_Node)
T_GO_Full <- unlist(GO_List)
EnrichedGO_CY15_24h <- list()
EnrichedGO_pvalue_CY15_24h <- list()

CY15_24h_Node <- list()
t <- proc.time()
for (i in 1:total) {
  temp_all_GO <- c()
  temp_all_pV <- c()
  T_AGI <- unlist(CY15_24h_enriched_Node[[i]])
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
  temp_CY15_24h <- GO_List[match(T_AGI, names(GO_List))]
  AGI_name <- c(AGI_name, intersect(T_AGI, names(temp_CY15_24h)))
  names(AGI_name) <- NumCY15_24h_subcluster[i]
  CY15_24h_Node <- c(CY15_24h_Node, list(AGI_name))
  AGI_name <- c()
  EnrichedGO_pvalue_CY15_24h <- c(EnrichedGO_pvalue_CY15_24h, list(temp_all_pV))
  names(EnrichedGO_pvalue_CY15_24h)[i] <- NumCY15_24h_subcluster[i]
  
  EnrichedGO_CY15_24h <- c(EnrichedGO_CY15_24h, list(temp_all_GO))
  names(EnrichedGO_CY15_24h)[i] <- NumCY15_24h_subcluster[i]
  
  print(total - i)
  i <- i+1
}
proc.time()-t
save(CY15_24h_Node, file = "~/Nakano_RNAseq/network_analysis/.RData/CY15_24h_Node.RData")

####CY15_48h####
i <- 1
total <- length(CY15_48h_enriched_Node)
T_GO_Full <- unlist(GO_List)
EnrichedGO_CY15_48h <- list()
EnrichedGO_pvalue_CY15_48h <- list()

CY15_48h_Node <- list()
t <- proc.time()
for (i in 1:total) {
  temp_all_GO <- c()
  temp_all_pV <- c()
  T_AGI <- unlist(CY15_48h_enriched_Node[[i]])
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
  temp_CY15_48h <- GO_List[match(T_AGI, names(GO_List))]
  AGI_name <- c(AGI_name, intersect(T_AGI, names(temp_CY15_48h)))
  names(AGI_name) <- NumCY15_48h_subcluster[i]
  CY15_48h_Node <- c(CY15_48h_Node, list(AGI_name))
  AGI_name <- c()
  EnrichedGO_pvalue_CY15_48h <- c(EnrichedGO_pvalue_CY15_48h, list(temp_all_pV))
  names(EnrichedGO_pvalue_CY15_48h)[i] <- NumCY15_48h_subcluster[i]
  
  EnrichedGO_CY15_48h <- c(EnrichedGO_CY15_48h, list(temp_all_GO))
  names(EnrichedGO_CY15_48h)[i] <- NumCY15_48h_subcluster[i]
  
  print(total - i)
  i <- i+1
}
proc.time()-t
save(CY15_48h_Node, file = "~/Nakano_RNAseq/network_analysis/.RData/CY15_48h_Node.RData")

####CY16_1h####
i <- 1
total <- length(CY16_1h_enriched_Node)
T_GO_Full <- unlist(GO_List)
EnrichedGO_CY16_1h <- list()
EnrichedGO_pvalue_CY16_1h <- list()

CY16_1h_Node <- list()
t <- proc.time()
for (i in 1:total) {
  temp_all_GO <- c()
  temp_all_pV <- c()
  T_AGI <- unlist(CY16_1h_enriched_Node[[i]])
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
  temp_CY16_1h <- GO_List[match(T_AGI, names(GO_List))]
  AGI_name <- c(AGI_name, intersect(T_AGI, names(temp_CY16_1h)))
  names(AGI_name) <- NumCY16_1h_subcluster[i]
  CY16_1h_Node <- c(CY16_1h_Node, list(AGI_name))
  AGI_name <- c()
  EnrichedGO_pvalue_CY16_1h <- c(EnrichedGO_pvalue_CY16_1h, list(temp_all_pV))
  names(EnrichedGO_pvalue_CY16_1h)[i] <- NumCY16_1h_subcluster[i]
  
  EnrichedGO_CY16_1h <- c(EnrichedGO_CY16_1h, list(temp_all_GO))
  names(EnrichedGO_CY16_1h)[i] <- NumCY16_1h_subcluster[i]
  
  print(total - i)
  i <- i+1
}
proc.time()-t
save(CY16_1h_Node, file = "~/Nakano_RNAseq/network_analysis/.RData/CY16_1h_Node.RData")

####CY16_3h####
i <- 1
total <- length(CY16_3h_enriched_Node)
T_GO_Full <- unlist(GO_List)
EnrichedGO_CY16_3h <- list()
EnrichedGO_pvalue_CY16_3h <- list()

CY16_3h_Node <- list()
t <- proc.time()
for (i in 1:total) {
  temp_all_GO <- c()
  temp_all_pV <- c()
  T_AGI <- unlist(CY16_3h_enriched_Node[[i]])
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
  temp_CY16_3h <- GO_List[match(T_AGI, names(GO_List))]
  AGI_name <- c(AGI_name, intersect(T_AGI, names(temp_CY16_3h)))
  names(AGI_name) <- NumCY16_3h_subcluster[i]
  CY16_3h_Node <- c(CY16_3h_Node, list(AGI_name))
  AGI_name <- c()
  EnrichedGO_pvalue_CY16_3h <- c(EnrichedGO_pvalue_CY16_3h, list(temp_all_pV))
  names(EnrichedGO_pvalue_CY16_3h)[i] <- NumCY16_3h_subcluster[i]
  
  EnrichedGO_CY16_3h <- c(EnrichedGO_CY16_3h, list(temp_all_GO))
  names(EnrichedGO_CY16_3h)[i] <- NumCY16_3h_subcluster[i]
  
  print(total - i)
  i <- i+1
}
proc.time()-t
save(CY16_3h_Node, file = "~/Nakano_RNAseq/network_analysis/.RData/CY16_3h_Node.RData")

####CY16_12h####
i <- 1
total <- length(CY16_12h_enriched_Node)
T_GO_Full <- unlist(GO_List)
EnrichedGO_CY16_12h <- list()
EnrichedGO_pvalue_CY16_12h <- list()

CY16_12h_Node <- list()
t <- proc.time()
for (i in 1:total) {
  temp_all_GO <- c()
  temp_all_pV <- c()
  T_AGI <- unlist(CY16_12h_enriched_Node[[i]])
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
  temp_CY16_12h <- GO_List[match(T_AGI, names(GO_List))]
  AGI_name <- c(AGI_name, intersect(T_AGI, names(temp_CY16_12h)))
  names(AGI_name) <- NumCY16_12h_subcluster[i]
  CY16_12h_Node <- c(CY16_12h_Node, list(AGI_name))
  AGI_name <- c()
  EnrichedGO_pvalue_CY16_12h <- c(EnrichedGO_pvalue_CY16_12h, list(temp_all_pV))
  names(EnrichedGO_pvalue_CY16_12h)[i] <- NumCY16_12h_subcluster[i]
  
  EnrichedGO_CY16_12h <- c(EnrichedGO_CY16_12h, list(temp_all_GO))
  names(EnrichedGO_CY16_12h)[i] <- NumCY16_12h_subcluster[i]
  
  print(total - i)
  i <- i+1
}
proc.time()-t
save(CY16_12h_Node, file = "~/Nakano_RNAseq/network_analysis/.RData/CY16_12h_Node.RData")

####CY16_24h####
i <- 1
total <- length(CY16_24h_enriched_Node)
T_GO_Full <- unlist(GO_List)
EnrichedGO_CY16_24h <- list()
EnrichedGO_pvalue_CY16_24h <- list()

CY16_24h_Node <- list()
t <- proc.time()
for (i in 1:total) {
  temp_all_GO <- c()
  temp_all_pV <- c()
  T_AGI <- unlist(CY16_24h_enriched_Node[[i]])
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
  temp_CY16_24h <- GO_List[match(T_AGI, names(GO_List))]
  AGI_name <- c(AGI_name, intersect(T_AGI, names(temp_CY16_24h)))
  names(AGI_name) <- NumCY16_24h_subcluster[i]
  CY16_24h_Node <- c(CY16_24h_Node, list(AGI_name))
  AGI_name <- c()
  EnrichedGO_pvalue_CY16_24h <- c(EnrichedGO_pvalue_CY16_24h, list(temp_all_pV))
  names(EnrichedGO_pvalue_CY16_24h)[i] <- NumCY16_24h_subcluster[i]
  
  EnrichedGO_CY16_24h <- c(EnrichedGO_CY16_24h, list(temp_all_GO))
  names(EnrichedGO_CY16_24h)[i] <- NumCY16_24h_subcluster[i]
  
  print(total - i)
  i <- i+1
}
proc.time()-t
save(CY16_24h_Node, file = "~/Nakano_RNAseq/network_analysis/.RData/CY16_24h_Node.RData")

####CY16_48h####
i <- 1
total <- length(CY16_48h_enriched_Node)
T_GO_Full <- unlist(GO_List)
EnrichedGO_CY16_48h <- list()
EnrichedGO_pvalue_CY16_48h <- list()

CY16_48h_Node <- list()
t <- proc.time()
for (i in 1:total) {
  temp_all_GO <- c()
  temp_all_pV <- c()
  T_AGI <- unlist(CY16_48h_enriched_Node[[i]])
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
  temp_CY16_48h <- GO_List[match(T_AGI, names(GO_List))]
  AGI_name <- c(AGI_name, intersect(T_AGI, names(temp_CY16_48h)))
  names(AGI_name) <- NumCY16_48h_subcluster[i]
  CY16_48h_Node <- c(CY16_48h_Node, list(AGI_name))
  AGI_name <- c()
  EnrichedGO_pvalue_CY16_48h <- c(EnrichedGO_pvalue_CY16_48h, list(temp_all_pV))
  names(EnrichedGO_pvalue_CY16_48h)[i] <- NumCY16_48h_subcluster[i]
  
  EnrichedGO_CY16_48h <- c(EnrichedGO_CY16_48h, list(temp_all_GO))
  names(EnrichedGO_CY16_48h)[i] <- NumCY16_48h_subcluster[i]
  
  print(total - i)
  i <- i+1
}
proc.time()-t
save(CY16_48h_Node, file = "~/Nakano_RNAseq/network_analysis/.RData/CY16_48h_Node.RData")


####CY20_1h####
i <- 1
total <- length(CY20_1h_enriched_Node)
T_GO_Full <- unlist(GO_List)
EnrichedGO_CY20_1h <- list()
EnrichedGO_pvalue_CY20_1h <- list()

CY20_1h_Node <- list()
t <- proc.time()
for (i in 1:total) {
  temp_all_GO <- c()
  temp_all_pV <- c()
  T_AGI <- unlist(CY20_1h_enriched_Node[[i]])
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
  temp_CY20_1h <- GO_List[match(T_AGI, names(GO_List))]
  AGI_name <- c(AGI_name, intersect(T_AGI, names(temp_CY20_1h)))
  names(AGI_name) <- NumCY20_1h_subcluster[i]
  CY20_1h_Node <- c(CY20_1h_Node, list(AGI_name))
  AGI_name <- c()
  EnrichedGO_pvalue_CY20_1h <- c(EnrichedGO_pvalue_CY20_1h, list(temp_all_pV))
  names(EnrichedGO_pvalue_CY20_1h)[i] <- NumCY20_1h_subcluster[i]
  
  EnrichedGO_CY20_1h <- c(EnrichedGO_CY20_1h, list(temp_all_GO))
  names(EnrichedGO_CY20_1h)[i] <- NumCY20_1h_subcluster[i]
  
  print(total - i)
  i <- i+1
}
proc.time()-t
save(CY20_1h_Node, file = "~/Nakano_RNAseq/network_analysis/.RData/CY20_1h_Node.RData")

####CY20_3h####
i <- 1
total <- length(CY20_3h_enriched_Node)
T_GO_Full <- unlist(GO_List)
EnrichedGO_CY20_3h <- list()
EnrichedGO_pvalue_CY20_3h <- list()

CY20_3h_Node <- list()
t <- proc.time()
for (i in 1:total) {
  temp_all_GO <- c()
  temp_all_pV <- c()
  T_AGI <- unlist(CY20_3h_enriched_Node[[i]])
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
  temp_CY20_3h <- GO_List[match(T_AGI, names(GO_List))]
  AGI_name <- c(AGI_name, intersect(T_AGI, names(temp_CY20_3h)))
  names(AGI_name) <- NumCY20_3h_subcluster[i]
  CY20_3h_Node <- c(CY20_3h_Node, list(AGI_name))
  AGI_name <- c()
  EnrichedGO_pvalue_CY20_3h <- c(EnrichedGO_pvalue_CY20_3h, list(temp_all_pV))
  names(EnrichedGO_pvalue_CY20_3h)[i] <- NumCY20_3h_subcluster[i]
  
  EnrichedGO_CY20_3h <- c(EnrichedGO_CY20_3h, list(temp_all_GO))
  names(EnrichedGO_CY20_3h)[i] <- NumCY20_3h_subcluster[i]
  
  print(total - i)
  i <- i+1
}
proc.time()-t
save(CY20_3h_Node, file = "~/Nakano_RNAseq/network_analysis/.RData/CY20_3h_Node.RData")

####CY20_12h####
i <- 1
total <- length(CY20_12h_enriched_Node)
T_GO_Full <- unlist(GO_List)
EnrichedGO_CY20_12h <- list()
EnrichedGO_pvalue_CY20_12h <- list()

CY20_12h_Node <- list()
t <- proc.time()
for (i in 1:total) {
  temp_all_GO <- c()
  temp_all_pV <- c()
  T_AGI <- unlist(CY20_12h_enriched_Node[[i]])
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
  temp_CY20_12h <- GO_List[match(T_AGI, names(GO_List))]
  AGI_name <- c(AGI_name, intersect(T_AGI, names(temp_CY20_12h)))
  names(AGI_name) <- NumCY20_12h_subcluster[i]
  CY20_12h_Node <- c(CY20_12h_Node, list(AGI_name))
  AGI_name <- c()
  EnrichedGO_pvalue_CY20_12h <- c(EnrichedGO_pvalue_CY20_12h, list(temp_all_pV))
  names(EnrichedGO_pvalue_CY20_12h)[i] <- NumCY20_12h_subcluster[i]
  
  EnrichedGO_CY20_12h <- c(EnrichedGO_CY20_12h, list(temp_all_GO))
  names(EnrichedGO_CY20_12h)[i] <- NumCY20_12h_subcluster[i]
  
  print(total - i)
  i <- i+1
}
proc.time()-t
save(CY20_12h_Node, file = "~/Nakano_RNAseq/network_analysis/.RData/CY20_12h_Node.RData")

####CY20_24h####
i <- 1
total <- length(CY20_24h_enriched_Node)
T_GO_Full <- unlist(GO_List)
EnrichedGO_CY20_24h <- list()
EnrichedGO_pvalue_CY20_24h <- list()

CY20_24h_Node <- list()
t <- proc.time()
for (i in 1:total) {
  temp_all_GO <- c()
  temp_all_pV <- c()
  T_AGI <- unlist(CY20_24h_enriched_Node[[i]])
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
  temp_CY20_24h <- GO_List[match(T_AGI, names(GO_List))]
  AGI_name <- c(AGI_name, intersect(T_AGI, names(temp_CY20_24h)))
  names(AGI_name) <- NumCY20_24h_subcluster[i]
  CY20_24h_Node <- c(CY20_24h_Node, list(AGI_name))
  AGI_name <- c()
  EnrichedGO_pvalue_CY20_24h <- c(EnrichedGO_pvalue_CY20_24h, list(temp_all_pV))
  names(EnrichedGO_pvalue_CY20_24h)[i] <- NumCY20_24h_subcluster[i]
  
  EnrichedGO_CY20_24h <- c(EnrichedGO_CY20_24h, list(temp_all_GO))
  names(EnrichedGO_CY20_24h)[i] <- NumCY20_24h_subcluster[i]
  
  print(total - i)
  i <- i+1
}
proc.time()-t
save(CY20_24h_Node, file = "~/Nakano_RNAseq/network_analysis/.RData/CY20_24h_Node.RData")

####CY20_48h####
i <- 1
total <- length(CY20_48h_enriched_Node)
T_GO_Full <- unlist(GO_List)
EnrichedGO_CY20_48h <- list()
EnrichedGO_pvalue_CY20_48h <- list()

CY20_48h_Node <- list()
t <- proc.time()
for (i in 1:total) {
  temp_all_GO <- c()
  temp_all_pV <- c()
  T_AGI <- unlist(CY20_48h_enriched_Node[[i]])
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
  temp_CY20_48h <- GO_List[match(T_AGI, names(GO_List))]
  AGI_name <- c(AGI_name, intersect(T_AGI, names(temp_CY20_48h)))
  names(AGI_name) <- NumCY20_48h_subcluster[i]
  CY20_48h_Node <- c(CY20_48h_Node, list(AGI_name))
  AGI_name <- c()
  EnrichedGO_pvalue_CY20_48h <- c(EnrichedGO_pvalue_CY20_48h, list(temp_all_pV))
  names(EnrichedGO_pvalue_CY20_48h)[i] <- NumCY20_48h_subcluster[i]
  
  EnrichedGO_CY20_48h <- c(EnrichedGO_CY20_48h, list(temp_all_GO))
  names(EnrichedGO_CY20_48h)[i] <- NumCY20_48h_subcluster[i]
  
  print(total - i)
  i <- i+1
}
proc.time()-t
save(CY20_48h_Node, file = "~/Nakano_RNAseq/network_analysis/.RData/CY20_48h_Node.RData")
