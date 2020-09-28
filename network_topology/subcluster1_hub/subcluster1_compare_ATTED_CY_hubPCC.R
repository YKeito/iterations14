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
####sum####
CY_heatmap[CY_heatmap$PCC == 0, ]$PCC <- NA
CY_heatmap[CY_heatmap$sample == 0, ]$sample <- NA

source_genes <-  as.character(unlist(c(substring(CY_heatmap$AGI, 1, 9), substring(CY_heatmap$AGI, 10, 18))))
PCC <-  c(rep(CY_heatmap$PCC, times = 2))
target_genes <- as.character(unlist(c(substring(CY_heatmap$AGI, 10, 18), substring(CY_heatmap$AGI, 1, 9))))
sample <- as.character(unlist(c(rep(CY_heatmap$sample, times = 2))))
CY_heatmap <- data.frame(source_genes = c(source_genes, "ZZ", "Z"),
                         PCC = c(PCC, -1, 1),
                         target_genes = c(target_genes, "Z", "ZZ"),
                         sample = c(sample, "AA", "AA")
)

ATTED_heatmap[ATTED_heatmap$PCC == 0, ]$PCC <- NA
ATTED_heatmap[ATTED_heatmap$sample == 0, ]$sample <- NA
source_genes <- as.character(unlist(c(substring(ATTED_heatmap$AGI, 1, 9), substring(ATTED_heatmap$AGI, 10, 18))))
PCC <- c(rep(ATTED_heatmap$PCC, times = 2))
target_genes <- as.character(unlist(c(substring(ATTED_heatmap$AGI, 10, 18), substring(ATTED_heatmap$AGI, 1, 9))))
sample <- as.character(unlist(c(rep(ATTED_heatmap$sample, times = 2))))

ATTED_heatmap <- data.frame(source_genes = c(source_genes, "ZZ", "Z"),
                            PCC = c(PCC, -1, 1),
                            target_genes = c(target_genes, "Z", "ZZ"),
                            sample = c(sample, "AA", "AA")
)

####plot heatmap ggplot2####
#library(RColorBrewer)
g <- ggplot(CY_heatmap, aes(x = target_genes, y = source_genes, file = PCC), na.rm = TRUE)
g <- g + geom_tile(aes(fill = PCC), color = "white")
g <- g + scale_fill_gradientn("value", colours = rev(brewer.pal(9, "RdBu")), na.value = "white")
g <- g + labs(x = "",y = "")
g <- g + scale_x_discrete(expand = c(0, 0))
g <- g + scale_y_discrete(expand = c(0, 0))
g <- g + theme(legend.position = "none",axis.ticks = element_blank(), axis.text.x = element_text(size = 7, angle = 280, hjust = 0, colour = "grey50"))
plot(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/results_Fig/ATTED_hubPCC_Fig.png", plot = g, dpi = 100)

ATTED_heatmap <- ATTED_heatmap[order(ATTED_heatmap$source_genes), ]
g <- ggplot(ATTED_heatmap, aes(x = target_genes, source_genes, file = PCC), na.rm = TRUE)
g <- g + geom_tile(aes(fill = PCC), color = "white")
g <- g + scale_fill_gradientn("value", colours = rev(brewer.pal(9, "RdBu")), na.value = "white")
g <- g + labs(x = "",y = "") + scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +  theme(legend.position = "none",axis.ticks = element_blank(), axis.text.x = element_text(size = 7, angle = 280, hjust = 0, colour = "grey50"))
plot(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/results_Fig/CY_hubPCC_Fig.png", plot = g, dpi = 100)

####PCC of PCC####
ggplotCY_PCC <- CY_heatmap$PCC
ggplotATTED_PCC <- ATTED_heatmap$PCC

ggplotCYPCC_ATTED_PCC <- data.frame(CY_PCC = ggplotCY_PCC,
                                    ATTED_PCC = ggplotATTED_PCC
)

#ggplotCYPCC_ATTED_PCC <- na.omit(ggplotCYPCC_ATTED_PCC)
i <- 1
yasue <- c()
for(i in i:nrow(ggplotCYPCC_ATTED_PCC)){
  yasue <- c(yasue, is.na(ggplotCYPCC_ATTED_PCC$CY_PCC[i]) == is.na(ggplotCYPCC_ATTED_PCC$ATTED_PCC[i]))
  i <- i+1
}
keito <- ggplotCYPCC_ATTED_PCC[!yasue, ]

keito2 <- data.frame(rbind(ggplotCYPCC_ATTED_PCC, keito))
keito2$CY_PCC[is.na(keito2$CY_PCC)] <- 0
keito2$ATTED_PCC[is.na(keito2$ATTED_PCC)] <- 0

i <- 1
yasue2 <- c()
for(i in i:nrow(ggplotCYPCC_ATTED_PCC)){
  yasue2 <- c(yasue2, keito2$CY_PCC[i] == keito2$ATTED_PCC[i])
  i <- i+1
}


library(ggplot2)
keito2 <- keito2[keito2$CY_PCC != 1, ]
keito2 <- keito2[keito2$CY_PCC != -1, ]
keito2 <- keito2[keito2$ATTED_PCC != 1, ]
keito2 <- keito2[keito2$ATTED_PCC != -1, ]




p <- ggplot(keito2, aes(x=CY_PCC, y=ATTED_PCC))
p <- p + ylim(-1, 1)
p <- p + geom_point()
plot(p)

q <- p + stat_smooth(method = "lm", se = FALSE, colour = "red", size = 1)
plot(q)
a <- lm(formula = CY_PCC~ATTED_PCC, data = keito2)