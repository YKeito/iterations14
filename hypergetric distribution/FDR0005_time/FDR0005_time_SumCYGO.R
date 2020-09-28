####CY15_1h####
GO_SLIM <- data.frame(AGI=ATH_GO_GOSLIM$`locus name`, ID=ATH_GO_GOSLIM$`GO ID`, Term=ATH_GO_GOSLIM$`GO term`)

total <- length(EnrichedGO_CY15_1h)
T_AGI <- c()
temp <- c()
temp1 <- c()
temp2 <- list()
temp3 <- c()
temp4 <- c()
temp5 <- c()
test2 <- c()
T_12 <- rep("2", times=length(unlist(EnrichedGO_pvalue_CY15_1h)))
aaaa <- p.adjust(unlist(EnrichedGO_pvalue_CY15_1h), method = "BH")
T_12[aaaa < 0.05] <- "1"
n <- 1
for (n in 1:total) {
  T_GO <- EnrichedGO_CY15_1h[[n]]
  T_Term <- GO_SLIM$Term[match(T_GO, GO_SLIM$ID)]
  m <- 1
  for(m in m:length(T_GO)){
    test <- GO_SLIM$AGI[T_GO[[m]] == GO_SLIM$ID]
    data <- c(data, intersect(CY15_1h_Node[[n]], test))
    data <- paste(data, collapse = " | ")
    data[which(data %in% "")] <- "NA"
    test2 <- c(test2, list(data))
    data <- c()
    m <- m+1
  }
  significance <- switch(T_12[n],
                         "1"="Y",
                         "2"="N")
  temp <- c(temp, rep(names(EnrichedGO_CY15_1h)[n], times=length(T_GO)))
  temp1 <- c(temp1, T_GO)
  test2 <- c(test2, as.character(unlist(temp2)))
  temp3 <- c(temp3, as.character(unlist(T_Term)))
  temp4 <- c(temp4, as.character(unlist(EnrichedGO_pvalue_CY15_1h[[n]])))
  temp5 <- c(temp5, rep(significance, times=length(T_GO)))
  print(total - n)
  n <- n+1
}
temp6 <- -log2(as.numeric(temp4))

EnrichedGO_Summ_CY15_1h <- data.frame(CluNum = temp,
                                   GO_ID = temp1,
                                   Node = unlist(test2),
                                   GO_term = temp3,
                                   qvalue = as.numeric(temp4),
                                   enrichment_score = temp6
                                   )

#CY15_3h#
GO_SLIM <- data.frame(AGI=ATH_GO_GOSLIM$`locus name`, ID=ATH_GO_GOSLIM$`GO ID`, Term=ATH_GO_GOSLIM$`GO term`)

total <- length(EnrichedGO_CY15_3h)
T_AGI <- c()
temp <- c()
temp1 <- c()
temp2 <- list()
temp3 <- c()
temp4 <- c()
temp5 <- c()
test2 <- c()
T_12 <- rep("2", times=length(unlist(EnrichedGO_pvalue_CY15_3h)))
aaaa <- p.adjust(unlist(EnrichedGO_pvalue_CY15_3h), method = "BH")
T_12[aaaa < 0.05] <- "1"
n <- 1
for (n in 1:total) {
  T_GO <- EnrichedGO_CY15_3h[[n]]
  T_Term <- GO_SLIM$Term[match(T_GO, GO_SLIM$ID)]
  m <- 1
  for(m in m:length(T_GO)){
    test <- GO_SLIM$AGI[T_GO[[m]] == GO_SLIM$ID]
    data <- c(data, intersect(CY15_3h_Node[[n]], test))
    data <- paste(data, collapse = " | ")
    data[which(data %in% "")] <- "NA"
    test2 <- c(test2, list(data))
    data <- c()
    m <- m+1
  }
  significance <- switch(T_12[n],
                         "1"="Y",
                         "2"="N")
  temp <- c(temp, rep(names(EnrichedGO_CY15_3h)[n], times=length(T_GO)))
  temp1 <- c(temp1, T_GO)
  test2 <- c(test2, as.character(unlist(temp2)))
  temp3 <- c(temp3, as.character(unlist(T_Term)))
  temp4 <- c(temp4, as.character(unlist(EnrichedGO_pvalue_CY15_3h[[n]])))
  temp5 <- c(temp5, rep(significance, times=length(T_GO)))
  print(total - n)
  n <- n+1
}
temp6 <- -log2(as.numeric(temp4))

EnrichedGO_Summ_CY15_3h <- data.frame(CluNum = temp,
                                      GO_ID = temp1,
                                      Node = unlist(test2),
                                      GO_term = temp3,
                                      qvalue = as.numeric(temp4),
                                      enrichment_score = temp6
                                      )

#CY15_12h#
GO_SLIM <- data.frame(AGI=ATH_GO_GOSLIM$`locus name`, ID=ATH_GO_GOSLIM$`GO ID`, Term=ATH_GO_GOSLIM$`GO term`)

total <- length(EnrichedGO_CY15_12h)
T_AGI <- c()
temp <- c()
temp1 <- c()
temp2 <- list()
temp3 <- c()
temp4 <- c()
temp5 <- c()
test2 <- c()
T_12 <- rep("2", times=length(unlist(EnrichedGO_pvalue_CY15_12h)))
aaaa <- p.adjust(unlist(EnrichedGO_pvalue_CY15_12h), method = "BH")
T_12[aaaa < 0.05] <- "1"
n <- 1
for (n in 1:total) {
  T_GO <- EnrichedGO_CY15_12h[[n]]
  T_Term <- GO_SLIM$Term[match(T_GO, GO_SLIM$ID)]
  m <- 1
  for(m in m:length(T_GO)){
    test <- GO_SLIM$AGI[T_GO[[m]] == GO_SLIM$ID]
    data <- c(data, intersect(CY15_12h_Node[[n]], test))
    data <- paste(data, collapse = " | ")
    data[which(data %in% "")] <- "NA"
    test2 <- c(test2, list(data))
    data <- c()
    m <- m+1
  }
  significance <- switch(T_12[n],
                         "1"="Y",
                         "2"="N")
  temp <- c(temp, rep(names(EnrichedGO_CY15_12h)[n], times=length(T_GO)))
  temp1 <- c(temp1, T_GO)
  test2 <- c(test2, as.character(unlist(temp2)))
  temp3 <- c(temp3, as.character(unlist(T_Term)))
  temp4 <- c(temp4, as.character(unlist(EnrichedGO_pvalue_CY15_12h[[n]])))
  temp5 <- c(temp5, rep(significance, times=length(T_GO)))
  print(total - n)
  n <- n+1
}
temp6 <- -log2(as.numeric(temp4))

EnrichedGO_Summ_CY15_12h <- data.frame(CluNum = temp,
                                      GO_ID = temp1,
                                      Node = unlist(test2),
                                      GO_term = temp3,
                                      qvalue = as.numeric(temp4),
                                      enrichment_score = temp6
                                      )

#CY15_24h#
GO_SLIM <- data.frame(AGI=ATH_GO_GOSLIM$`locus name`, ID=ATH_GO_GOSLIM$`GO ID`, Term=ATH_GO_GOSLIM$`GO term`)

total <- length(EnrichedGO_CY15_24h)
T_AGI <- c()
temp <- c()
temp1 <- c()
temp2 <- list()
temp3 <- c()
temp4 <- c()
temp5 <- c()
test2 <- c()
T_12 <- rep("2", times=length(unlist(EnrichedGO_pvalue_CY15_24h)))
aaaa <- p.adjust(unlist(EnrichedGO_pvalue_CY15_24h), method = "BH")
T_12[aaaa < 0.05] <- "1"
n <- 1
for (n in 1:total) {
  T_GO <- EnrichedGO_CY15_24h[[n]]
  T_Term <- GO_SLIM$Term[match(T_GO, GO_SLIM$ID)]
  m <- 1
  for(m in m:length(T_GO)){
    test <- GO_SLIM$AGI[T_GO[[m]] == GO_SLIM$ID]
    data <- c(data, intersect(CY15_24h_Node[[n]], test))
    data <- paste(data, collapse = " | ")
    data[which(data %in% "")] <- "NA"
    test2 <- c(test2, list(data))
    data <- c()
    m <- m+1
  }
  significance <- switch(T_12[n],
                         "1"="Y",
                         "2"="N")
  temp <- c(temp, rep(names(EnrichedGO_CY15_24h)[n], times=length(T_GO)))
  temp1 <- c(temp1, T_GO)
  test2 <- c(test2, as.character(unlist(temp2)))
  temp3 <- c(temp3, as.character(unlist(T_Term)))
  temp4 <- c(temp4, as.character(unlist(EnrichedGO_pvalue_CY15_24h[[n]])))
  temp5 <- c(temp5, rep(significance, times=length(T_GO)))
  print(total - n)
  n <- n+1
}
temp6 <- -log2(as.numeric(temp4))

EnrichedGO_Summ_CY15_24h <- data.frame(CluNum = temp,
                                      GO_ID = temp1,
                                      Node = unlist(test2),
                                      GO_term = temp3,
                                      qvalue = as.numeric(temp4),
                                      enrichment_score = temp6
                                      )

#CY15_48h#
GO_SLIM <- data.frame(AGI=ATH_GO_GOSLIM$`locus name`, ID=ATH_GO_GOSLIM$`GO ID`, Term=ATH_GO_GOSLIM$`GO term`)

total <- length(EnrichedGO_CY15_48h)
T_AGI <- c()
temp <- c()
temp1 <- c()
temp2 <- list()
temp3 <- c()
temp4 <- c()
temp5 <- c()
test2 <- c()
T_12 <- rep("2", times=length(unlist(EnrichedGO_pvalue_CY15_48h)))
aaaa <- p.adjust(unlist(EnrichedGO_pvalue_CY15_48h), method = "BH")
T_12[aaaa < 0.05] <- "1"
n <- 1
for (n in 1:total) {
  T_GO <- EnrichedGO_CY15_48h[[n]]
  T_Term <- GO_SLIM$Term[match(T_GO, GO_SLIM$ID)]
  m <- 1
  for(m in m:length(T_GO)){
    test <- GO_SLIM$AGI[T_GO[[m]] == GO_SLIM$ID]
    data <- c(data, intersect(CY15_48h_Node[[n]], test))
    data <- paste(data, collapse = " | ")
    data[which(data %in% "")] <- "NA"
    test2 <- c(test2, list(data))
    data <- c()
    m <- m+1
  }
  significance <- switch(T_12[n],
                         "1"="Y",
                         "2"="N")
  temp <- c(temp, rep(names(EnrichedGO_CY15_48h)[n], times=length(T_GO)))
  temp1 <- c(temp1, T_GO)
  test2 <- c(test2, as.character(unlist(temp2)))
  temp3 <- c(temp3, as.character(unlist(T_Term)))
  temp4 <- c(temp4, as.character(unlist(EnrichedGO_pvalue_CY15_48h[[n]])))
  temp5 <- c(temp5, rep(significance, times=length(T_GO)))
  print(total - n)
  n <- n+1
}
temp6 <- -log2(as.numeric(temp4))

EnrichedGO_Summ_CY15_48h <- data.frame(CluNum = temp,
                                      GO_ID = temp1,
                                      Node = unlist(test2),
                                      GO_term = temp3,
                                      qvalue = as.numeric(temp4),
                                      enrichment_score = temp6
                                      )




###################################################################################################################################
EnrichedGO_Summ_CY15_1h_FDR005 <- EnrichedGO_Summ_CY15_1h[EnrichedGO_Summ_CY15_1h$qvalue < 0.05, ]
EnrichedGO_Summ_CY15_3h_FDR005 <- EnrichedGO_Summ_CY15_3h[EnrichedGO_Summ_CY15_3h$qvalue < 0.05, ]
EnrichedGO_Summ_CY15_12h_FDR005 <- EnrichedGO_Summ_CY15_12h[EnrichedGO_Summ_CY15_12h$qvalue < 0.05, ]
EnrichedGO_Summ_CY15_24h_FDR005 <- EnrichedGO_Summ_CY15_24h[EnrichedGO_Summ_CY15_24h$qvalue < 0.05, ]
EnrichedGO_Summ_CY15_48h_FDR005 <- EnrichedGO_Summ_CY15_48h[EnrichedGO_Summ_CY15_48h$qvalue < 0.05, ]
#降順に並び替え#
EnrichedGO_Summ_CY15_1h_FDR005 <- EnrichedGO_Summ_CY15_1h_FDR005[order(EnrichedGO_Summ_CY15_1h_FDR005$enrichment_score, decreasing = T), ]
EnrichedGO_Summ_CY15_3h_FDR005 <- EnrichedGO_Summ_CY15_3h_FDR005[order(EnrichedGO_Summ_CY15_3h_FDR005$enrichment_score, decreasing = T), ]
EnrichedGO_Summ_CY15_12h_FDR005 <- EnrichedGO_Summ_CY15_12h_FDR005[order(EnrichedGO_Summ_CY15_12h_FDR005$enrichment_score, decreasing = T), ]
EnrichedGO_Summ_CY15_24h_FDR005 <- EnrichedGO_Summ_CY15_24h_FDR005[order(EnrichedGO_Summ_CY15_24h_FDR005$enrichment_score, decreasing = T), ]
EnrichedGO_Summ_CY15_48h_FDR005 <- EnrichedGO_Summ_CY15_48h_FDR005[order(EnrichedGO_Summ_CY15_48h_FDR005$enrichment_score, decreasing = T), ]

save(EnrichedGO_Summ_CY15_1h_FDR005, file = "~/Nakano_RNAseq/network_analysis/.RData/EnrichedGO_Summ_CY15_1h_FDR005.RData")
write.table(EnrichedGO_Summ_CY15_1h_FDR005, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY15FDR0005_time/EnrichedGO_Summ_CY15_1h_FDR005.txt",sep="\t", col.names=T, row.names=F, append=F, quote=F)
save(EnrichedGO_Summ_CY15_3h_FDR005, file = "~/Nakano_RNAseq/network_analysis/.RData/EnrichedGO_Summ_CY15_3h.RData")
write.table(EnrichedGO_Summ_CY15_3h_FDR005, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY15FDR0005_time/EnrichedGO_Summ_CY15_3h_FDR005.txt",sep="\t", col.names=T, row.names=F, append=F, quote=F)
save(EnrichedGO_Summ_CY15_12h_FDR005, file = "~/Nakano_RNAseq/network_analysis/.RData/EnrichedGO_Summ_CY15_12h.RData")
write.table(EnrichedGO_Summ_CY15_12h_FDR005, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY15FDR0005_time/EnrichedGO_Summ_CY15_12h_FDR005.txt",sep="\t", col.names=T, row.names=F, append=F, quote=F)
save(EnrichedGO_Summ_CY15_24h_FDR005, file = "~/Nakano_RNAseq/network_analysis/.RData/EnrichedGO_Summ_CY15_24h.RData")
write.table(EnrichedGO_Summ_CY15_24h_FDR005, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY15FDR0005_time/EnrichedGO_Summ_CY15_24h_FDR005.txt",sep="\t", col.names=T, row.names=F, append=F, quote=F)
save(EnrichedGO_Summ_CY15_48h_FDR005, file = "~/Nakano_RNAseq/network_analysis/.RData/EnrichedGO_Summ_CY15_48h.RData")
write.table(EnrichedGO_Summ_CY15_48h_FDR005, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY15FDR0005_time/EnrichedGO_Summ_CY15_48h_FDR005.txt",sep="\t", col.names=T, row.names=F, append=F, quote=F)
####CY16_1h####
GO_SLIM <- data.frame(AGI=ATH_GO_GOSLIM$`locus name`, ID=ATH_GO_GOSLIM$`GO ID`, Term=ATH_GO_GOSLIM$`GO term`)

total <- length(EnrichedGO_CY16_1h)
T_AGI <- c()
temp <- c()
temp1 <- c()
temp2 <- list()
temp3 <- c()
temp4 <- c()
temp5 <- c()
test2 <- c()
T_12 <- rep("2", times=length(unlist(EnrichedGO_pvalue_CY16_1h)))
aaaa <- p.adjust(unlist(EnrichedGO_pvalue_CY16_1h), method = "BH")
T_12[aaaa < 0.05] <- "1"
n <- 1
for (n in 1:total) {
  T_GO <- EnrichedGO_CY16_1h[[n]]
  T_Term <- GO_SLIM$Term[match(T_GO, GO_SLIM$ID)]
  m <- 1
  for(m in m:length(T_GO)){
    test <- GO_SLIM$AGI[T_GO[[m]] == GO_SLIM$ID]
    data <- c(data, intersect(CY16_1h_Node[[n]], test))
    data <- paste(data, collapse = " | ")
    data[which(data %in% "")] <- "NA"
    test2 <- c(test2, list(data))
    data <- c()
    m <- m+1
  }
  significance <- switch(T_12[n],
                         "1"="Y",
                         "2"="N")
  temp <- c(temp, rep(names(EnrichedGO_CY16_1h)[n], times=length(T_GO)))
  temp1 <- c(temp1, T_GO)
  test2 <- c(test2, as.character(unlist(temp2)))
  temp3 <- c(temp3, as.character(unlist(T_Term)))
  temp4 <- c(temp4, as.character(unlist(EnrichedGO_pvalue_CY16_1h[[n]])))
  temp5 <- c(temp5, rep(significance, times=length(T_GO)))
  print(total - n)
  n <- n+1
}
temp6 <- -log2(as.numeric(temp4))

EnrichedGO_Summ_CY16_1h <- data.frame(CluNum = temp,
                                      GO_ID = temp1,
                                      Node = unlist(test2),
                                      GO_term = temp3,
                                      qvalue = as.numeric(temp4),
                                      enrichment_score = temp6
)

####CY16_3h####
GO_SLIM <- data.frame(AGI=ATH_GO_GOSLIM$`locus name`, ID=ATH_GO_GOSLIM$`GO ID`, Term=ATH_GO_GOSLIM$`GO term`)

total <- length(EnrichedGO_CY16_3h)
T_AGI <- c()
temp <- c()
temp1 <- c()
temp2 <- list()
temp3 <- c()
temp4 <- c()
temp5 <- c()
test2 <- c()
T_12 <- rep("2", times=length(unlist(EnrichedGO_pvalue_CY16_3h)))
aaaa <- p.adjust(unlist(EnrichedGO_pvalue_CY16_3h), method = "BH")
T_12[aaaa < 0.05] <- "1"
n <- 1
for (n in 1:total) {
  T_GO <- EnrichedGO_CY16_3h[[n]]
  T_Term <- GO_SLIM$Term[match(T_GO, GO_SLIM$ID)]
  m <- 1
  for(m in m:length(T_GO)){
    test <- GO_SLIM$AGI[T_GO[[m]] == GO_SLIM$ID]
    data <- c(data, intersect(CY16_3h_Node[[n]], test))
    data <- paste(data, collapse = " | ")
    data[which(data %in% "")] <- "NA"
    test2 <- c(test2, list(data))
    data <- c()
    m <- m+1
  }
  significance <- switch(T_12[n],
                         "1"="Y",
                         "2"="N")
  temp <- c(temp, rep(names(EnrichedGO_CY16_3h)[n], times=length(T_GO)))
  temp1 <- c(temp1, T_GO)
  test2 <- c(test2, as.character(unlist(temp2)))
  temp3 <- c(temp3, as.character(unlist(T_Term)))
  temp4 <- c(temp4, as.character(unlist(EnrichedGO_pvalue_CY16_3h[[n]])))
  temp5 <- c(temp5, rep(significance, times=length(T_GO)))
  print(total - n)
  n <- n+1
}
temp6 <- -log2(as.numeric(temp4))

EnrichedGO_Summ_CY16_3h <- data.frame(CluNum = temp,
                                      GO_ID = temp1,
                                      Node = unlist(test2),
                                      GO_term = temp3,
                                      qvalue = as.numeric(temp4),
                                      enrichment_score = temp6
)

####CY16_12h####
GO_SLIM <- data.frame(AGI=ATH_GO_GOSLIM$`locus name`, ID=ATH_GO_GOSLIM$`GO ID`, Term=ATH_GO_GOSLIM$`GO term`)

total <- length(EnrichedGO_CY16_12h)
T_AGI <- c()
temp <- c()
temp1 <- c()
temp2 <- list()
temp3 <- c()
temp4 <- c()
temp5 <- c()
test2 <- c()
T_12 <- rep("2", times=length(unlist(EnrichedGO_pvalue_CY16_12h)))
aaaa <- p.adjust(unlist(EnrichedGO_pvalue_CY16_12h), method = "BH")
T_12[aaaa < 0.05] <- "1"
n <- 1
for (n in 1:total) {
  T_GO <- EnrichedGO_CY16_12h[[n]]
  T_Term <- GO_SLIM$Term[match(T_GO, GO_SLIM$ID)]
  m <- 1
  for(m in m:length(T_GO)){
    test <- GO_SLIM$AGI[T_GO[[m]] == GO_SLIM$ID]
    data <- c(data, intersect(CY16_12h_Node[[n]], test))
    data <- paste(data, collapse = " | ")
    data[which(data %in% "")] <- "NA"
    test2 <- c(test2, list(data))
    data <- c()
    m <- m+1
  }
  significance <- switch(T_12[n],
                         "1"="Y",
                         "2"="N")
  temp <- c(temp, rep(names(EnrichedGO_CY16_12h)[n], times=length(T_GO)))
  temp1 <- c(temp1, T_GO)
  test2 <- c(test2, as.character(unlist(temp2)))
  temp3 <- c(temp3, as.character(unlist(T_Term)))
  temp4 <- c(temp4, as.character(unlist(EnrichedGO_pvalue_CY16_12h[[n]])))
  temp5 <- c(temp5, rep(significance, times=length(T_GO)))
  print(total - n)
  n <- n+1
}
temp6 <- -log2(as.numeric(temp4))

EnrichedGO_Summ_CY16_12h <- data.frame(CluNum = temp,
                                       GO_ID = temp1,
                                       Node = unlist(test2),
                                       GO_term = temp3,
                                       qvalue = as.numeric(temp4),
                                       enrichment_score = temp6
)

####CY16_24h####
GO_SLIM <- data.frame(AGI=ATH_GO_GOSLIM$`locus name`, ID=ATH_GO_GOSLIM$`GO ID`, Term=ATH_GO_GOSLIM$`GO term`)

total <- length(EnrichedGO_CY16_24h)
T_AGI <- c()
temp <- c()
temp1 <- c()
temp2 <- list()
temp3 <- c()
temp4 <- c()
temp5 <- c()
test2 <- c()
T_12 <- rep("2", times=length(unlist(EnrichedGO_pvalue_CY16_24h)))
aaaa <- p.adjust(unlist(EnrichedGO_pvalue_CY16_24h), method = "BH")
T_12[aaaa < 0.05] <- "1"
n <- 1
for (n in 1:total) {
  T_GO <- EnrichedGO_CY16_24h[[n]]
  T_Term <- GO_SLIM$Term[match(T_GO, GO_SLIM$ID)]
  m <- 1
  for(m in m:length(T_GO)){
    test <- GO_SLIM$AGI[T_GO[[m]] == GO_SLIM$ID]
    data <- c(data, intersect(CY16_24h_Node[[n]], test))
    data <- paste(data, collapse = " | ")
    data[which(data %in% "")] <- "NA"
    test2 <- c(test2, list(data))
    data <- c()
    m <- m+1
  }
  significance <- switch(T_12[n],
                         "1"="Y",
                         "2"="N")
  temp <- c(temp, rep(names(EnrichedGO_CY16_24h)[n], times=length(T_GO)))
  temp1 <- c(temp1, T_GO)
  test2 <- c(test2, as.character(unlist(temp2)))
  temp3 <- c(temp3, as.character(unlist(T_Term)))
  temp4 <- c(temp4, as.character(unlist(EnrichedGO_pvalue_CY16_24h[[n]])))
  temp5 <- c(temp5, rep(significance, times=length(T_GO)))
  print(total - n)
  n <- n+1
}
temp6 <- -log2(as.numeric(temp4))

EnrichedGO_Summ_CY16_24h <- data.frame(CluNum = temp,
                                       GO_ID = temp1,
                                       Node = unlist(test2),
                                       GO_term = temp3,
                                       qvalue = as.numeric(temp4),
                                       enrichment_score = temp6
)

####CY16_48h####
GO_SLIM <- data.frame(AGI=ATH_GO_GOSLIM$`locus name`, ID=ATH_GO_GOSLIM$`GO ID`, Term=ATH_GO_GOSLIM$`GO term`)

total <- length(EnrichedGO_CY16_48h)
T_AGI <- c()
temp <- c()
temp1 <- c()
temp2 <- list()
temp3 <- c()
temp4 <- c()
temp5 <- c()
test2 <- c()
T_12 <- rep("2", times=length(unlist(EnrichedGO_pvalue_CY16_48h)))
aaaa <- p.adjust(unlist(EnrichedGO_pvalue_CY16_48h), method = "BH")
T_12[aaaa < 0.05] <- "1"
n <- 1
for (n in 1:total) {
  T_GO <- EnrichedGO_CY16_48h[[n]]
  T_Term <- GO_SLIM$Term[match(T_GO, GO_SLIM$ID)]
  m <- 1
  for(m in m:length(T_GO)){
    test <- GO_SLIM$AGI[T_GO[[m]] == GO_SLIM$ID]
    data <- c(data, intersect(CY16_48h_Node[[n]], test))
    data <- paste(data, collapse = " | ")
    data[which(data %in% "")] <- "NA"
    test2 <- c(test2, list(data))
    data <- c()
    m <- m+1
  }
  significance <- switch(T_12[n],
                         "1"="Y",
                         "2"="N")
  temp <- c(temp, rep(names(EnrichedGO_CY16_48h)[n], times=length(T_GO)))
  temp1 <- c(temp1, T_GO)
  test2 <- c(test2, as.character(unlist(temp2)))
  temp3 <- c(temp3, as.character(unlist(T_Term)))
  temp4 <- c(temp4, as.character(unlist(EnrichedGO_pvalue_CY16_48h[[n]])))
  temp5 <- c(temp5, rep(significance, times=length(T_GO)))
  print(total - n)
  n <- n+1
}
temp6 <- -log2(as.numeric(temp4))

EnrichedGO_Summ_CY16_48h <- data.frame(CluNum = temp,
                                       GO_ID = temp1,
                                       Node = unlist(test2),
                                       GO_term = temp3,
                                       qvalue = as.numeric(temp4),
                                       enrichment_score = temp6
)




###################################################################################################################################
EnrichedGO_Summ_CY16_1h_FDR005 <- EnrichedGO_Summ_CY16_1h[EnrichedGO_Summ_CY16_1h$qvalue < 0.05, ]
EnrichedGO_Summ_CY16_3h_FDR005 <- EnrichedGO_Summ_CY16_3h[EnrichedGO_Summ_CY16_3h$qvalue < 0.05, ]
EnrichedGO_Summ_CY16_12h_FDR005 <- EnrichedGO_Summ_CY16_12h[EnrichedGO_Summ_CY16_12h$qvalue < 0.05, ]
EnrichedGO_Summ_CY16_24h_FDR005 <- EnrichedGO_Summ_CY16_24h[EnrichedGO_Summ_CY16_24h$qvalue < 0.05, ]
EnrichedGO_Summ_CY16_48h_FDR005 <- EnrichedGO_Summ_CY16_48h[EnrichedGO_Summ_CY16_48h$qvalue < 0.05, ]
#降順に並び替え#
EnrichedGO_Summ_CY16_1h_FDR005 <- EnrichedGO_Summ_CY16_1h_FDR005[order(EnrichedGO_Summ_CY16_1h_FDR005$enrichment_score, decreasing = T), ]
EnrichedGO_Summ_CY16_3h_FDR005 <- EnrichedGO_Summ_CY16_3h_FDR005[order(EnrichedGO_Summ_CY16_3h_FDR005$enrichment_score, decreasing = T), ]
EnrichedGO_Summ_CY16_12h_FDR005 <- EnrichedGO_Summ_CY16_12h_FDR005[order(EnrichedGO_Summ_CY16_12h_FDR005$enrichment_score, decreasing = T), ]
EnrichedGO_Summ_CY16_24h_FDR005 <- EnrichedGO_Summ_CY16_24h_FDR005[order(EnrichedGO_Summ_CY16_24h_FDR005$enrichment_score, decreasing = T), ]
EnrichedGO_Summ_CY16_48h_FDR005 <- EnrichedGO_Summ_CY16_48h_FDR005[order(EnrichedGO_Summ_CY16_48h_FDR005$enrichment_score, decreasing = T), ]

save(EnrichedGO_Summ_CY16_1h_FDR005, file = "~/Nakano_RNAseq/network_analysis/.RData/EnrichedGO_Summ_CY16_1h_FDR005.RData")
write.table(EnrichedGO_Summ_CY16_1h_FDR005, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY16FDR0005_time/EnrichedGO_Summ_CY16_1h_FDR005.txt",sep="\t", col.names=T, row.names=F, append=F, quote=F)
save(EnrichedGO_Summ_CY16_3h_FDR005, file = "~/Nakano_RNAseq/network_analysis/.RData/EnrichedGO_Summ_CY16_3h.RData")
write.table(EnrichedGO_Summ_CY16_3h_FDR005, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY16FDR0005_time/EnrichedGO_Summ_CY16_3h_FDR005.txt",sep="\t", col.names=T, row.names=F, append=F, quote=F)
save(EnrichedGO_Summ_CY16_12h_FDR005, file = "~/Nakano_RNAseq/network_analysis/.RData/EnrichedGO_Summ_CY16_12h.RData")
write.table(EnrichedGO_Summ_CY16_12h_FDR005, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY16FDR0005_time/EnrichedGO_Summ_CY16_12h_FDR005.txt",sep="\t", col.names=T, row.names=F, append=F, quote=F)
save(EnrichedGO_Summ_CY16_24h_FDR005, file = "~/Nakano_RNAseq/network_analysis/.RData/EnrichedGO_Summ_CY16_24h.RData")
write.table(EnrichedGO_Summ_CY16_24h_FDR005, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY16FDR0005_time/EnrichedGO_Summ_CY16_24h_FDR005.txt",sep="\t", col.names=T, row.names=F, append=F, quote=F)
save(EnrichedGO_Summ_CY16_48h_FDR005, file = "~/Nakano_RNAseq/network_analysis/.RData/EnrichedGO_Summ_CY16_48h.RData")
write.table(EnrichedGO_Summ_CY16_48h_FDR005, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY16FDR0005_time/EnrichedGO_Summ_CY16_48h_FDR005.txt",sep="\t", col.names=T, row.names=F, append=F, quote=F)

####CY20_1h####
GO_SLIM <- data.frame(AGI=ATH_GO_GOSLIM$`locus name`, ID=ATH_GO_GOSLIM$`GO ID`, Term=ATH_GO_GOSLIM$`GO term`)

total <- length(EnrichedGO_CY20_1h)
T_AGI <- c()
temp <- c()
temp1 <- c()
temp2 <- list()
temp3 <- c()
temp4 <- c()
temp5 <- c()
test2 <- c()
T_12 <- rep("2", times=length(unlist(EnrichedGO_pvalue_CY20_1h)))
aaaa <- p.adjust(unlist(EnrichedGO_pvalue_CY20_1h), method = "BH")
T_12[aaaa < 0.05] <- "1"
n <- 1
for (n in 1:total) {
  T_GO <- EnrichedGO_CY20_1h[[n]]
  T_Term <- GO_SLIM$Term[match(T_GO, GO_SLIM$ID)]
  m <- 1
  for(m in m:length(T_GO)){
    test <- GO_SLIM$AGI[T_GO[[m]] == GO_SLIM$ID]
    data <- c(data, intersect(CY20_1h_Node[[n]], test))
    data <- paste(data, collapse = " | ")
    data[which(data %in% "")] <- "NA"
    test2 <- c(test2, list(data))
    data <- c()
    m <- m+1
  }
  significance <- switch(T_12[n],
                         "1"="Y",
                         "2"="N")
  temp <- c(temp, rep(names(EnrichedGO_CY20_1h)[n], times=length(T_GO)))
  temp1 <- c(temp1, T_GO)
  test2 <- c(test2, as.character(unlist(temp2)))
  temp3 <- c(temp3, as.character(unlist(T_Term)))
  temp4 <- c(temp4, as.character(unlist(EnrichedGO_pvalue_CY20_1h[[n]])))
  temp5 <- c(temp5, rep(significance, times=length(T_GO)))
  print(total - n)
  n <- n+1
}
temp6 <- -log2(as.numeric(temp4))

EnrichedGO_Summ_CY20_1h <- data.frame(CluNum = temp,
                                      GO_ID = temp1,
                                      Node = unlist(test2),
                                      GO_term = temp3,
                                      qvalue = as.numeric(temp4),
                                      enrichment_score = temp6
)

####CY20_3h####
GO_SLIM <- data.frame(AGI=ATH_GO_GOSLIM$`locus name`, ID=ATH_GO_GOSLIM$`GO ID`, Term=ATH_GO_GOSLIM$`GO term`)

total <- length(EnrichedGO_CY20_3h)
T_AGI <- c()
temp <- c()
temp1 <- c()
temp2 <- list()
temp3 <- c()
temp4 <- c()
temp5 <- c()
test2 <- c()
T_12 <- rep("2", times=length(unlist(EnrichedGO_pvalue_CY20_3h)))
aaaa <- p.adjust(unlist(EnrichedGO_pvalue_CY20_3h), method = "BH")
T_12[aaaa < 0.05] <- "1"
n <- 1
for (n in 1:total) {
  T_GO <- EnrichedGO_CY20_3h[[n]]
  T_Term <- GO_SLIM$Term[match(T_GO, GO_SLIM$ID)]
  m <- 1
  for(m in m:length(T_GO)){
    test <- GO_SLIM$AGI[T_GO[[m]] == GO_SLIM$ID]
    data <- c(data, intersect(CY20_3h_Node[[n]], test))
    data <- paste(data, collapse = " | ")
    data[which(data %in% "")] <- "NA"
    test2 <- c(test2, list(data))
    data <- c()
    m <- m+1
  }
  significance <- switch(T_12[n],
                         "1"="Y",
                         "2"="N")
  temp <- c(temp, rep(names(EnrichedGO_CY20_3h)[n], times=length(T_GO)))
  temp1 <- c(temp1, T_GO)
  test2 <- c(test2, as.character(unlist(temp2)))
  temp3 <- c(temp3, as.character(unlist(T_Term)))
  temp4 <- c(temp4, as.character(unlist(EnrichedGO_pvalue_CY20_3h[[n]])))
  temp5 <- c(temp5, rep(significance, times=length(T_GO)))
  print(total - n)
  n <- n+1
}
temp6 <- -log2(as.numeric(temp4))

EnrichedGO_Summ_CY20_3h <- data.frame(CluNum = temp,
                                      GO_ID = temp1,
                                      Node = unlist(test2),
                                      GO_term = temp3,
                                      qvalue = as.numeric(temp4),
                                      enrichment_score = temp6
)

####CY20_12h####
GO_SLIM <- data.frame(AGI=ATH_GO_GOSLIM$`locus name`, ID=ATH_GO_GOSLIM$`GO ID`, Term=ATH_GO_GOSLIM$`GO term`)

total <- length(EnrichedGO_CY20_12h)
T_AGI <- c()
temp <- c()
temp1 <- c()
temp2 <- list()
temp3 <- c()
temp4 <- c()
temp5 <- c()
test2 <- c()
T_12 <- rep("2", times=length(unlist(EnrichedGO_pvalue_CY20_12h)))
aaaa <- p.adjust(unlist(EnrichedGO_pvalue_CY20_12h), method = "BH")
T_12[aaaa < 0.05] <- "1"
n <- 1
for (n in 1:total) {
  T_GO <- EnrichedGO_CY20_12h[[n]]
  T_Term <- GO_SLIM$Term[match(T_GO, GO_SLIM$ID)]
  m <- 1
  for(m in m:length(T_GO)){
    test <- GO_SLIM$AGI[T_GO[[m]] == GO_SLIM$ID]
    data <- c(data, intersect(CY20_12h_Node[[n]], test))
    data <- paste(data, collapse = " | ")
    data[which(data %in% "")] <- "NA"
    test2 <- c(test2, list(data))
    data <- c()
    m <- m+1
  }
  significance <- switch(T_12[n],
                         "1"="Y",
                         "2"="N")
  temp <- c(temp, rep(names(EnrichedGO_CY20_12h)[n], times=length(T_GO)))
  temp1 <- c(temp1, T_GO)
  test2 <- c(test2, as.character(unlist(temp2)))
  temp3 <- c(temp3, as.character(unlist(T_Term)))
  temp4 <- c(temp4, as.character(unlist(EnrichedGO_pvalue_CY20_12h[[n]])))
  temp5 <- c(temp5, rep(significance, times=length(T_GO)))
  print(total - n)
  n <- n+1
}
temp6 <- -log2(as.numeric(temp4))

EnrichedGO_Summ_CY20_12h <- data.frame(CluNum = temp,
                                       GO_ID = temp1,
                                       Node = unlist(test2),
                                       GO_term = temp3,
                                       qvalue = as.numeric(temp4),
                                       enrichment_score = temp6
)

#CY20_24h#
GO_SLIM <- data.frame(AGI=ATH_GO_GOSLIM$`locus name`, ID=ATH_GO_GOSLIM$`GO ID`, Term=ATH_GO_GOSLIM$`GO term`)

total <- length(EnrichedGO_CY20_24h)
T_AGI <- c()
temp <- c()
temp1 <- c()
temp2 <- list()
temp3 <- c()
temp4 <- c()
temp5 <- c()
test2 <- c()
T_12 <- rep("2", times=length(unlist(EnrichedGO_pvalue_CY20_24h)))
aaaa <- p.adjust(unlist(EnrichedGO_pvalue_CY20_24h), method = "BH")
T_12[aaaa < 0.05] <- "1"
n <- 1
for (n in 1:total) {
  T_GO <- EnrichedGO_CY20_24h[[n]]
  T_Term <- GO_SLIM$Term[match(T_GO, GO_SLIM$ID)]
  m <- 1
  for(m in m:length(T_GO)){
    test <- GO_SLIM$AGI[T_GO[[m]] == GO_SLIM$ID]
    data <- c(data, intersect(CY20_24h_Node[[n]], test))
    data <- paste(data, collapse = " | ")
    data[which(data %in% "")] <- "NA"
    test2 <- c(test2, list(data))
    data <- c()
    m <- m+1
  }
  significance <- switch(T_12[n],
                         "1"="Y",
                         "2"="N")
  temp <- c(temp, rep(names(EnrichedGO_CY20_24h)[n], times=length(T_GO)))
  temp1 <- c(temp1, T_GO)
  test2 <- c(test2, as.character(unlist(temp2)))
  temp3 <- c(temp3, as.character(unlist(T_Term)))
  temp4 <- c(temp4, as.character(unlist(EnrichedGO_pvalue_CY20_24h[[n]])))
  temp5 <- c(temp5, rep(significance, times=length(T_GO)))
  print(total - n)
  n <- n+1
}
temp6 <- -log2(as.numeric(temp4))

EnrichedGO_Summ_CY20_24h <- data.frame(CluNum = temp,
                                       GO_ID = temp1,
                                       Node = unlist(test2),
                                       GO_term = temp3,
                                       qvalue = as.numeric(temp4),
                                       enrichment_score = temp6
)

#CY20_48h#
GO_SLIM <- data.frame(AGI=ATH_GO_GOSLIM$`locus name`, ID=ATH_GO_GOSLIM$`GO ID`, Term=ATH_GO_GOSLIM$`GO term`)

total <- length(EnrichedGO_CY20_48h)
T_AGI <- c()
temp <- c()
temp1 <- c()
temp2 <- list()
temp3 <- c()
temp4 <- c()
temp5 <- c()
test2 <- c()
T_12 <- rep("2", times=length(unlist(EnrichedGO_pvalue_CY20_48h)))
aaaa <- p.adjust(unlist(EnrichedGO_pvalue_CY20_48h), method = "BH")
T_12[aaaa < 0.05] <- "1"
n <- 1
for (n in 1:total) {
  T_GO <- EnrichedGO_CY20_48h[[n]]
  T_Term <- GO_SLIM$Term[match(T_GO, GO_SLIM$ID)]
  m <- 1
  for(m in m:length(T_GO)){
    test <- GO_SLIM$AGI[T_GO[[m]] == GO_SLIM$ID]
    data <- c(data, intersect(CY20_48h_Node[[n]], test))
    data <- paste(data, collapse = " | ")
    data[which(data %in% "")] <- "NA"
    test2 <- c(test2, list(data))
    data <- c()
    m <- m+1
  }
  significance <- switch(T_12[n],
                         "1"="Y",
                         "2"="N")
  temp <- c(temp, rep(names(EnrichedGO_CY20_48h)[n], times=length(T_GO)))
  temp1 <- c(temp1, T_GO)
  test2 <- c(test2, as.character(unlist(temp2)))
  temp3 <- c(temp3, as.character(unlist(T_Term)))
  temp4 <- c(temp4, as.character(unlist(EnrichedGO_pvalue_CY20_48h[[n]])))
  temp5 <- c(temp5, rep(significance, times=length(T_GO)))
  print(total - n)
  n <- n+1
}
temp6 <- -log2(as.numeric(temp4))

EnrichedGO_Summ_CY20_48h <- data.frame(CluNum = temp,
                                       GO_ID = temp1,
                                       Node = unlist(test2),
                                       GO_term = temp3,
                                       qvalue = as.numeric(temp4),
                                       enrichment_score = temp6
)




###################################################################################################################################
EnrichedGO_Summ_CY20_1h_FDR005 <- EnrichedGO_Summ_CY20_1h[EnrichedGO_Summ_CY20_1h$qvalue < 0.05, ]
EnrichedGO_Summ_CY20_3h_FDR005 <- EnrichedGO_Summ_CY20_3h[EnrichedGO_Summ_CY20_3h$qvalue < 0.05, ]
EnrichedGO_Summ_CY20_12h_FDR005 <- EnrichedGO_Summ_CY20_12h[EnrichedGO_Summ_CY20_12h$qvalue < 0.05, ]
EnrichedGO_Summ_CY20_24h_FDR005 <- EnrichedGO_Summ_CY20_24h[EnrichedGO_Summ_CY20_24h$qvalue < 0.05, ]
EnrichedGO_Summ_CY20_48h_FDR005 <- EnrichedGO_Summ_CY20_48h[EnrichedGO_Summ_CY20_48h$qvalue < 0.05, ]
#降順に並び替え#
EnrichedGO_Summ_CY20_1h_FDR005 <- EnrichedGO_Summ_CY20_1h_FDR005[order(EnrichedGO_Summ_CY20_1h_FDR005$enrichment_score, decreasing = T), ]
EnrichedGO_Summ_CY20_3h_FDR005 <- EnrichedGO_Summ_CY20_3h_FDR005[order(EnrichedGO_Summ_CY20_3h_FDR005$enrichment_score, decreasing = T), ]
EnrichedGO_Summ_CY20_12h_FDR005 <- EnrichedGO_Summ_CY20_12h_FDR005[order(EnrichedGO_Summ_CY20_12h_FDR005$enrichment_score, decreasing = T), ]
EnrichedGO_Summ_CY20_24h_FDR005 <- EnrichedGO_Summ_CY20_24h_FDR005[order(EnrichedGO_Summ_CY20_24h_FDR005$enrichment_score, decreasing = T), ]
EnrichedGO_Summ_CY20_48h_FDR005 <- EnrichedGO_Summ_CY20_48h_FDR005[order(EnrichedGO_Summ_CY20_48h_FDR005$enrichment_score, decreasing = T), ]
####out put####
save(EnrichedGO_Summ_CY20_1h_FDR005, file = "~/Nakano_RNAseq/network_analysis/.RData/EnrichedGO_Summ_CY20_1h_FDR005.RData")
write.table(EnrichedGO_Summ_CY20_1h_FDR005, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY20FDR0005_time/EnrichedGO_Summ_CY20_1h_FDR005.txt",sep="\t", col.names=T, row.names=F, append=F, quote=F)
save(EnrichedGO_Summ_CY20_3h_FDR005, file = "~/Nakano_RNAseq/network_analysis/.RData/EnrichedGO_Summ_CY20_3h.RData")
write.table(EnrichedGO_Summ_CY20_3h_FDR005, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY20FDR0005_time/EnrichedGO_Summ_CY20_3h_FDR005.txt",sep="\t", col.names=T, row.names=F, append=F, quote=F)
save(EnrichedGO_Summ_CY20_12h_FDR005, file = "~/Nakano_RNAseq/network_analysis/.RData/EnrichedGO_Summ_CY20_12h.RData")
write.table(EnrichedGO_Summ_CY20_12h_FDR005, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY20FDR0005_time/EnrichedGO_Summ_CY20_12h_FDR005.txt",sep="\t", col.names=T, row.names=F, append=F, quote=F)
save(EnrichedGO_Summ_CY20_24h_FDR005, file = "~/Nakano_RNAseq/network_analysis/.RData/EnrichedGO_Summ_CY20_24h.RData")
write.table(EnrichedGO_Summ_CY20_24h_FDR005, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY20FDR0005_time/EnrichedGO_Summ_CY20_24h_FDR005.txt",sep="\t", col.names=T, row.names=F, append=F, quote=F)
save(EnrichedGO_Summ_CY20_48h_FDR005, file = "~/Nakano_RNAseq/network_analysis/.RData/EnrichedGO_Summ_CY20_48h.RData")
write.table(EnrichedGO_Summ_CY20_48h_FDR005, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY20FDR0005_time/EnrichedGO_Summ_CY20_48h_FDR005.txt",sep="\t", col.names=T, row.names=F, append=F, quote=F)
