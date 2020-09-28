load("~/bigdata/yasue/GO_analysis/ATH_GO_GOSLIM.RData")
####CY15####
GO_SLIM <- data.frame(AGI=ATH_GO_GOSLIM$`locus name`, ID=ATH_GO_GOSLIM$`GO ID`, Term=ATH_GO_GOSLIM$`GO term`)

total <- length(EnrichedGO_CY15)
T_AGI <- c()
temp <- c()
temp1 <- c()
temp2 <- list()
temp3 <- c()
temp4 <- c()
temp5 <- c()
test2 <- c()
T_12 <- rep("2", times=length(unlist(EnrichedGO_pvalue_CY15)))
aaaa <- p.adjust(unlist(EnrichedGO_pvalue_CY15), method = "BH")
T_12[aaaa < 0.05] <- "1"
n <- 1
for (n in 1:total) {
  T_GO <- EnrichedGO_CY15[[n]]
  T_Term <- GO_SLIM$Term[match(T_GO, GO_SLIM$ID)]
  m <- 1
  for(m in m:length(T_GO)){
    test <- GO_SLIM$AGI[T_GO[[m]] == GO_SLIM$ID]
    data <- c(data, intersect(CY15_Node[[n]], test))
    data <- paste(data, collapse = " | ")
    data[which(data %in% "")] <- "NA"
    test2 <- c(test2, list(data))
    data <- c()
    m <- m+1
  }
  significance <- switch(T_12[n],
                         "1"="Y",
                         "2"="N")
  temp <- c(temp, rep(names(EnrichedGO_CY15)[n], times=length(T_GO)))
  temp1 <- c(temp1, T_GO)
  test2 <- c(test2, as.character(unlist(temp2)))
  temp3 <- c(temp3, as.character(unlist(T_Term)))
  temp4 <- c(temp4, as.character(unlist(EnrichedGO_pvalue_CY15[[n]])))
  temp5 <- c(temp5, rep(significance, times=length(T_GO)))
  print(total - n)
  n <- n+1
}
temp6 <- -log2(as.numeric(temp4))

EnrichedGO_Summ_CY15 <- data.frame(CluNum = temp,
                                  GO_ID = temp1,
                                  Node = unlist(test2),
                                  GO_term = temp3,
                                  qvalue = as.numeric(temp4),
                                  enrichment_score = temp6
                                  )
####CY16####
total <- length(EnrichedGO_CY16)
T_AGI <- c()
temp <- c()
temp1 <- c()
temp2 <- list()
temp3 <- c()
temp4 <- c()
temp5 <- c()
test2 <- c()
T_12 <- rep("2", times=length(unlist(EnrichedGO_pvalue_CY16)))
aaaa <- p.adjust(unlist(EnrichedGO_pvalue_CY16), method = "BH")
T_12[aaaa < 0.05] <- "1"
n <- 1
for (n in 1:total) {
  T_GO <- EnrichedGO_CY16[[n]]
  T_Term <- GO_SLIM$Term[match(T_GO, GO_SLIM$ID)]
  m <- 1
  for(m in m:length(T_GO)){
    test <- GO_SLIM$AGI[T_GO[[m]] == GO_SLIM$ID]
    data <- c(data, intersect(CY16_Node[[n]], test))
    data <- paste(data, collapse = " | ")
    data[which(data %in% "")] <- "NA"
    test2 <- c(test2, list(data))
    data <- c()
    m <- m+1
  }
  significance <- switch(T_12[n],
                         "1"="Y",
                         "2"="N")
  temp <- c(temp, rep(names(EnrichedGO_CY16)[n], times=length(T_GO)))
  temp1 <- c(temp1, T_GO)
  test2 <- c(test2, as.character(unlist(temp2)))
  temp3 <- c(temp3, as.character(unlist(T_Term)))
  temp4 <- c(temp4, as.character(unlist(EnrichedGO_pvalue_CY16[[n]])))
  temp5 <- c(temp5, rep(significance, times=length(T_GO)))
  print(total - n)
  n <- n+1
}
temp6 <- -log2(as.numeric(temp4))

EnrichedGO_Summ_CY16 <- data.frame(CluNum = temp,
                                   GO_ID = temp1,
                                   Node = unlist(test2),
                                   GO_term = temp3,
                                   qvalue = as.numeric(temp4),
                                   enrichment_score = temp6
                                   )

####CY20####
total <- length(EnrichedGO_CY20)
T_AGI <- c()
temp <- c()
temp1 <- c()
temp2 <- list()
temp3 <- c()
temp4 <- c()
temp5 <- c()
test2 <- c()
T_12 <- rep("2", times=length(unlist(EnrichedGO_pvalue_CY20)))
aaaa <- p.adjust(unlist(EnrichedGO_pvalue_CY20), method = "BH")
T_12[aaaa < 0.05] <- "1"
n <- 1
for (n in 1:total) {
  T_GO <- EnrichedGO_CY20[[n]]
  T_Term <- GO_SLIM$Term[match(T_GO, GO_SLIM$ID)]
  m <- 1
  for(m in m:length(T_GO)){
    test <- GO_SLIM$AGI[T_GO[[m]] == GO_SLIM$ID]
    data <- c(data, intersect(CY20_Node[[n]], test))
    data <- paste(data, collapse = " | ")
    data[which(data %in% "")] <- "NA"
    test2 <- c(test2, list(data))
    data <- c()
    m <- m+1
  }
  significance <- switch(T_12[n],
                         "1"="Y",
                         "2"="N")
  temp <- c(temp, rep(names(EnrichedGO_CY20)[n], times=length(T_GO)))
  temp1 <- c(temp1, T_GO)
  test2 <- c(test2, as.character(unlist(temp2)))
  temp3 <- c(temp3, as.character(unlist(T_Term)))
  temp4 <- c(temp4, as.character(unlist(EnrichedGO_pvalue_CY20[[n]])))
  temp5 <- c(temp5, rep(significance, times=length(T_GO)))
  print(total - n)
  n <- n+1
}
temp6 <- -log2(as.numeric(temp4))

EnrichedGO_Summ_CY20 <- data.frame(CluNum = temp,
                                   GO_ID = temp1,
                                   Node = unlist(test2),
                                   GO_term = temp3,
                                   qvalue = as.numeric(temp4),
                                   enrichment_score = temp6
                                   )

EnrichedGO_Summ_CY15FDR005 <- EnrichedGO_Summ_CY15[EnrichedGO_Summ_CY15$qvalue < 0.05, ]
EnrichedGO_Summ_CY16FDR005 <- EnrichedGO_Summ_CY15[EnrichedGO_Summ_CY16$qvalue < 0.05, ]
EnrichedGO_Summ_CY20FDR005 <- EnrichedGO_Summ_CY15[EnrichedGO_Summ_CY20$qvalue < 0.05, ]

save(EnrichedGO_Summ_CY20FDR005, file = "~/Nakano_RNAseq/network_analysis/.RData/iterations14/EnrichedGO_Summ_CY20.RData")
write.table(EnrichedGO_Summ_CY20FDR005, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/iterations14/CY_subcluster/GO/EnrichedGO_Summ_CY20FDR005.txt",sep="\t", col.names=T, row.names=F, append=F, quote=F)
save(EnrichedGO_Summ_CY16FDR005, file = "~/Nakano_RNAseq/network_analysis/.RData/iterations14/EnrichedGO_Summ_CY16.RData")
write.table(EnrichedGO_Summ_CY16FDR005, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/iterations14/CY_subcluster/GO/EnrichedGO_Summ_CY16FDR005.txt",sep="\t", col.names=T, row.names=F, append=F, quote=F)
save(EnrichedGO_Summ_CY15FDR005, file = "~/Nakano_RNAseq/network_analysis/.RData/iterations14/EnrichedGO_Summ_CY15.RData")
write.table(EnrichedGO_Summ_CY15FDR005, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/iterations14/CY_subcluster/GO/EnrichedGO_Summ_CY15FDR005.txt",sep="\t", col.names=T, row.names=F, append=F, quote=F)
