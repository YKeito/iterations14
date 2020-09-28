####CY_hub combn####
CY_hub <-  c(subcluster1_CY151620_hub, subcluster1_CY15_hub, subcluster1_CY1516_hub, subcluster1_CY1520_hub, subcluster1_CY16_hub, subcluster1_CY1620_hub, subcluster1_CY20_hub)
CY_hub <- CY_hub[order(CY_hub)]
test <- combn(CY_hub, 2)
m <- 1
CYpair1 <- c()
for(m in m:ncol(test)){
  CYpair1 <- c(CYpair1, paste0(test[1, m], test[2, m]))
  m <- m+1
}
#test <- combn(CY_hub, 2)
#m <- 1
#CYpair2 <- c()
#for(m in m:ncol(test)){
#  CYpair2 <- c(CYpair2, paste0(test[2, m], test[1, m]))
#  m <- m+1
#}
####allRNASeq_pairall####
n <- 1
m <- 1
allRNASeqpair <- c()
allRNASeqpairall <- c()
for(n in n:length(CY_hub)){
  for(m in m:length(CY_hub)){
    test <- CY_hub[n] == allRNASeq_cytoscape_th$source_genes & CY_hub[m] == allRNASeq_cytoscape_th$target_genes
    allRNASeqpair <- rbind(allRNASeqpair, allRNASeq_cytoscape_th[test, ])
    m <- m+1
  }
  allRNASeqpairall <- rbind(allRNASeqpairall, allRNASeqpair)
  allRNASeqpair <- c()
  print(n)
  n <- n+1
  m <- 1
}

n <- 1
allRNASeq_PCC <- c()
for(n in n:nrow(allRNASeqpairall)){
  allRNASeq_PCC <- c(allRNASeq_PCC, paste0(allRNASeqpairall$source_genes[n], allRNASeqpairall$target_genes[n]))
  print(n)
  n <- 1
}
#####CY_hubPCC####
temp4 <- match(CYpair1, allRNASeq_PCC)
temp4 <- temp4[!is.na(temp4)]

CY_hubPCC <- allRNASeqpairall[temp4, ]

####ATTED_pairall####
n <- 1
m <- 1
ATTED_pair <- c()
ATTED_pairall <- c()
for(n in n:length(CY_hub)){
  for(m in m:length(CY_hub)){
    test <- CY_hub[n] == ATTED_cytoscape_th$source_genes & CY_hub[m] == ATTED_cytoscape_th$target_genes
    ATTED_pair <- rbind(ATTED_pair, ATTED_cytoscape_th[test, ])
    m <- m+1
  }
  ATTED_pairall <- rbind(ATTED_pairall, ATTED_pair)
  ATTED_pair <- c()
  print(n)
  n <- n+1
  m <- 1
}
n <- 1
ATTED_PCC <- c()
for(n in n:nrow(ATTED_pairall)){
  ATTED_PCC <- c(ATTED_PCC, paste0(ATTED_pairall$source_genes[n], ATTED_pairall$target_genes[n]))
  print(n)
  n <- 1
}
####ATTED_PCC####
temp <- match(CYpair1, ATTED_PCC)
temp <- temp[!is.na(temp)]
ATTED_hubPCC <- ATTED_pairall[temp, ]
####CY_heatmap####
CY_heatmap <- data.frame(AGI = CYpair1[order(CYpair1)],
                         PCC = rep(0, times = length(CYpair1)),
                         sample = rep(0, times = length(CYpair1))
)
#ハブ遺伝子の総当たりCYpair1
#実際のエッジ組み合わせallRNASeq_PCC
#ATTEDの組み合わせATTED_PCC
CY_unique_pair <- setdiff(allRNASeq_PCC, ATTED_PCC)
ATTED_unique_pair <- setdiff(ATTED_PCC, allRNASeq_PCC)
CYATTED_common <- intersect(ATTED_PCC, allRNASeq_PCC)
not_edges <- setdiff(CYpair1, allRNASeq_PCC)

#sample
a <- match(CY_unique_pair, CY_heatmap$AGI)
temp <- a[!is.na(a)]
b <- match(CYATTED_common, CY_heatmap$AGI)
temp2 <- b[!is.na(b)]
CY_heatmap[temp, ]$sample <- "CY"
CY_heatmap[temp2, ]$sample <- "CYATTED"
#PCC
CY_heatmap[temp, ]$PCC <- allRNASeqpairall[match(CY_unique_pair, allRNASeq_PCC), ]$interaction_value
CY_heatmap[temp2, ]$PCC <- allRNASeqpairall[match(CYATTED_common, allRNASeq_PCC), ]$interaction_value
#length(CY_heatmap[temp, ]$PCC) == length(allRNASeqpairall[match(CY_unique_pair, allRNASeq_PCC), ]$interaction_value)
#length(CY_heatmap[temp2, ]$PCC) == length(allRNASeqpairall[match(CYATTED_common, allRNASeq_PCC), ]$interaction_value)
####ATTED_heatmap####
ATTED_heatmap <- data.frame(AGI = CYpair1[order(CYpair1)],
                            PCC = rep(0, times = length(CYpair1)),
                            sample = rep(0, times = length(CYpair1))
)
#sample
c <- match(ATTED_unique_pair, ATTED_heatmap$AGI)
temp3 <- c[!is.na(c)]
d <- match(CYATTED_common, ATTED_heatmap$AGI)
temp4 <- d[!is.na(d)]
ATTED_heatmap[temp3, ]$sample <- "ATTED"
ATTED_heatmap[temp4, ]$sample <- "CYATTED"
ATTED_heatmap[temp3, ]$PCC <- ATTED_pairall[match(ATTED_unique_pair, ATTED_PCC), ]$interaction_value
ATTED_heatmap[temp4, ]$PCC <- ATTED_pairall[match(CYATTED_common, ATTED_PCC), ]$interaction_value

#check length(ATTED_heatmap[temp3, ]$PCC) == length(ATTED_pairall[match(ATTED_unique_pair, ATTED_PCC), ]$interaction_value)
####unique_edges####
CY_uniqueedges <- na.omit(CY_heatmap[CY_heatmap$sample == "CY", ])#204行
CY_source <- substring(CY_uniqueedges$AGI, 1, 9)
CY_target <- substring(CY_uniqueedges$AGI, 10, 18)
CY_AGIunion <- union(CY_source, CY_target)
m <- 1
CY_test <- c()
for(m in m:length(CY_AGIunion)){
  CY_test <- c(CY_test, length(grep(CY_AGIunion[m], CY_uniqueedges$AGI)))
  names(CY_test)[m] <- CY_AGIunion[m]
}
ATTED_uniqueedges <- na.omit(ATTED_heatmap[ATTED_heatmap$sample == "ATTED", ])#176行
ATTED_source <- substring(ATTED_uniqueedges$AGI, 1, 9)
ATTED_target <- substring(ATTED_uniqueedges$AGI, 10, 18)
ATTED_AGIunion <- union(ATTED_source, ATTED_target)
m <- 1
ATTED_test <- c()
for(m in m:length(ATTED_AGIunion)){
  ATTED_test <- c(ATTED_test, length(grep(ATTED_AGIunion[m], ATTED_uniqueedges$AGI)))
  names(ATTED_test)[m] <- ATTED_AGIunion[m]
}
test <- c(ATTED_test, rep(0, times = c(length(CY_test)-length(ATTED_test))))
names(test)[31:39] <- "NA"
i <- 1
CYATTED_defference <- c()
CYATTEDunion <- union(names(CY_test), names(ATTED_test))
aa <- c()
bb <- c()
for(i in i:length(CYATTEDunion)){
  aa <- c(aa, CY_test[names(CY_test) == CYATTEDunion[i]])
  bb <- c(bb, test[names(test) == CYATTEDunion[i]])
  i <- i+1
}
cc <- rep(0, times = 9)
names(cc) <- setdiff(CYATTEDunion, names(test))
dd <- c(bb, cc)
e <- 1
ww <- c()
for(e in e:39){
 ww <- c(ww, CY_test[order(names(CY_test))][e] - dd[order(names(dd))][e]) 
}
CY_ATTED_edges <- data.frame(edges =ww, 
                             abs_edges = abs(ww))
write.table(CY_ATTED_edges, "~/Nakano_RNAseq/network_analysis/subcluster1_results/CYATTED_edges.txt", sep = "\t")