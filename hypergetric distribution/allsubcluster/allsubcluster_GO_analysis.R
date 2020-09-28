load("~/bigdata/yasue/GO_analysis/GO_List.RData")
i <- 1
total <- max(iterations14$X__mclCluster)
T_GO_Full <- unlist(GO_List)
EnrichedGO_MCL <- list()
EnrichedGO_pvalue_MCL <- list()

check_AGI <- list()
t <- proc.time()
for (i in 1:total) {
  temp_all_GO <- c()
  temp_all_pV <- c()
  T_AGI <- iterations14$name[iterations14$X__mclCluster == i]
  T_GO <- unlist(GO_List[match(T_AGI, names(GO_List))])
  k <- length(T_GO)
  T_GO_Uni <- unique(T_GO)
  AGI_name <- c()
  j <- 1
  total_j <- length(T_GO_Uni)
  for (j in j:total_j) {
    q <- sum(T_GO == T_GO_Uni[j])
    m <- sum(T_GO_Full == T_GO_Uni[j])
    n <- length(T_GO_Full) - m
    temp <- phyper(q-1, m, n, k, lower.tail=FALSE)
    temp_all_pV <- c(temp_all_pV, temp)
    temp_all_GO <- c(temp_all_GO, T_GO_Uni[j])
    j <- j+1
  }
  temp50 <- GO_List[match(T_AGI, names(GO_List))]
  AGI_name <- c(AGI_name, names(temp50))
  names(AGI_name) <- i
  check_AGI <- c(check_AGI, list(AGI_name))
  AGI_name <- c()

  EnrichedGO_pvalue_MCL <- c(EnrichedGO_pvalue_MCL, list(temp_all_pV))
  names(EnrichedGO_pvalue_MCL)[i] <- i
  
  EnrichedGO_MCL <- c(EnrichedGO_MCL, list(temp_all_GO))
  names(EnrichedGO_MCL)[i] <- i
  
  print(total - i)
  i <- i+1
}
proc.time()-t