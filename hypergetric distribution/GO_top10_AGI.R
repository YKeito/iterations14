#subcluster1
test <- c()
test2 <- list()
Number <- c()
i <- 1
for(i in i:10){
  a <- unlist(strsplit(as.character(subcluster1_top10$Node[i]), " | "))
  a <- a[seq(1, length(a), 2)]
  expression <- allRNASeq[match(a, rownames(allRNASeq)), ]
  m <- 1
  for(m in m:nrow(expression)){
    temp <- expression[m, ]
    temp <- temp[temp$CY15_1h >= 3 | temp$CY15_3h >= 3 | temp$CY15_12h >= 3 | temp$CY15_24h >=3 | temp$CY15_48h >=3 | 
                 temp$CY16_1h >= 3 | temp$CY16_3h >= 3 | temp$CY16_12h >= 3 | temp$CY16_24h >=3 | temp$CY16_48h >=3 |
                 temp$CY20_1h >= 3 | temp$CY20_3h >= 3 | temp$CY20_12h >= 3 | temp$CY20_24h >=3 | temp$CY20_48h >=3, ]
    test <- c(test, rownames(temp))
  }
  test2 <- c(test2, list(test))
  test <- c()
  Number <- c(Number, rep(i, times = length(unlist(test2[i]))))
}
names(test2) <- c(1:10)

subcluster1 <- data.frame(top10 = Number,
                          AGI = unlist(test2))
write.table(subcluster1, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/iterations14/allsubcluster/subcluster1_expression.txt",sep="\t", col.names=T, row.names=F, append=F, quote=F)
#subcluster2
test <- c()
test2 <- list()
Number <- c()
i <- 1
for(i in i:10){
  a <- unlist(strsplit(as.character(subcluster2_top10$Node[i]), " | "))
  a <- a[seq(1, length(a), 2)]
  expression <- allRNASeq[match(a, rownames(allRNASeq)), ]
  m <- 1
  for(m in m:nrow(expression)){
    temp <- expression[m, ]
    temp <- temp[temp$CY15_1h >= 3 | temp$CY15_3h >= 3 | temp$CY15_12h >= 3 | temp$CY15_24h >=3 | temp$CY15_48h >=3 | 
                   temp$CY16_1h >= 3 | temp$CY16_3h >= 3 | temp$CY16_12h >= 3 | temp$CY16_24h >=3 | temp$CY16_48h >=3 |
                   temp$CY20_1h >= 3 | temp$CY20_3h >= 3 | temp$CY20_12h >= 3 | temp$CY20_24h >=3 | temp$CY20_48h >=3, ]
    test <- c(test, rownames(temp))
  }
  test2 <- c(test2, list(test))
  test <- c()
  Number <- c(Number, rep(i, times = length(unlist(test2[i]))))
}
names(test2) <- c(1:10)

subcluster2 <- data.frame(top10 = Number,
                          AGI = unlist(test2))
#subcluster8
test <- c()
test2 <- list()
Number <- c()
i <- 1
for(i in i:10){
  a <- unlist(strsplit(as.character(subcluster8_top10$Node[i]), " | "))
  a <- a[seq(1, length(a), 2)]
  expression <- allRNASeq[match(a, rownames(allRNASeq)), ]
  m <- 1
  for(m in m:nrow(expression)){
    temp <- expression[m, ]
    temp <- temp[temp$CY15_1h >= 3 | temp$CY15_3h >= 3 | temp$CY15_12h >= 3 | temp$CY15_24h >=3 | temp$CY15_48h >=3 | 
                   temp$CY16_1h >= 3 | temp$CY16_3h >= 3 | temp$CY16_12h >= 3 | temp$CY16_24h >=3 | temp$CY16_48h >=3 |
                   temp$CY20_1h >= 3 | temp$CY20_3h >= 3 | temp$CY20_12h >= 3 | temp$CY20_24h >=3 | temp$CY20_48h >=3, ]
    test <- c(test, rownames(temp))
  }
  test2 <- c(test2, list(test))
  test <- c()
  Number <- c(Number, rep(i, times = length(unlist(test2[i]))))
}
names(test2) <- c(1:10)

subcluster8 <- data.frame(top10 = Number,
                          AGI = unlist(test2))
#subcluster40
test <- c()
test2 <- list()
Number <- c()
i <- 1
for(i in i:10){
  a <- unlist(strsplit(as.character(subcluster40_top10$Node[i]), " | "))
  a <- a[seq(1, length(a), 2)]
  expression <- allRNASeq[match(a, rownames(allRNASeq)), ]
  m <- 1
  for(m in m:nrow(expression)){
    temp <- expression[m, ]
    temp <- temp[temp$CY15_1h >= 3 | temp$CY15_3h >= 3 | temp$CY15_12h >= 3 | temp$CY15_24h >=3 | temp$CY15_48h >=3 | 
                   temp$CY16_1h >= 3 | temp$CY16_3h >= 3 | temp$CY16_12h >= 3 | temp$CY16_24h >=3 | temp$CY16_48h >=3 |
                   temp$CY20_1h >= 3 | temp$CY20_3h >= 3 | temp$CY20_12h >= 3 | temp$CY20_24h >=3 | temp$CY20_48h >=3, ]
    test <- c(test, rownames(temp))
  }
  test2 <- c(test2, list(test))
  test <- c()
  Number <- c(Number, rep(i, times = length(unlist(test2[i]))))
}
names(test2) <- c(1:10)

subcluster40 <- data.frame(top10 = Number,
                           AGI = unlist(test2))
