#/home/yasue/Nakano_RNAseq/network_analysis/script
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
####PCC of PCC data.frame作成####
ATTED_heatmap <- data.frame(AGI = CYpair1[order(CYpair1)],
                            PCC = rep(0, times = length(CYpair1))
                            )

ATTED_pair <- paste0(ATTED_cytoscape$source_genes, ATTED_cytoscape$target_genes)
CY_pair <- paste0(allRNASeq_cytoscape$source_genes, allRNASeq_cytoscape$target_genes)
CYATTED_commonedges <- intersect(CY_pair, ATTED_pair)


ATTED_temp <- data.frame(AGI = ATTED_pair, 
                     PCC = ATTED_cytoscape$interaction_value
                     )

ATTED_heatmap[match(CYATTED_commonedges, ATTED_heatmap$AGI), ]$PCC <- ATTED_temp$PCC
ATTED_heatmap[ATTED_heatmap$PCC == 0, ]$PCC <- NA
ATTED_temp <- na.omit(ATTED_heatmap)


CY_heatmap <- data.frame(AGI = CYpair1[order(CYpair1)],
                         PCC = rep(0, times = length(CYpair1)),
                         )

CY_temp <- data.frame(AGI = CY_pair,
                      PCC = allRNASeq_cytoscape$interaction_value
                      )
CY_temp <- CY_temp[match(CYATTED_commonedges, CY_heatmap$AGI), ]

CY_heatmap[match(CYATTED_commonedges, CY_heatmap$AGI), ]$PCC <- CY_temp$PCC
CY_heatmap[CY_heatmap$PCC == 0, ]$PCC <- NA
CY_temp <- na.omit(CY_heatmap)

#WRKY50なし
#grep("AT5G26170", CY_temp$AGI)#WRKY50
#grep("AT5G26170", CY_pair)

grep("AT1G06160", CY_temp$AGI)#ORA59
grep("AT2G47520", CY_temp$AGI)#ERF71
grep("AT4G37750", CY_temp$AGI)#ANT

sample <- rep("black", times = nrow(keito3))
sample[grep("AT1G06160", CY_temp$AGI)] <- "red" #ORA59
sample[grep("AT2G47520", CY_temp$AGI)] <- "green"#ERF71
sample[grep("AT4G37750", CY_temp$AGI)] <- "blue"#ANT


keito3 <- data.frame(ATTED_PCC = ATTED_temp$PCC,
                     CY_PCC = CY_temp$PCC,
                     sample = sample)

####PCC of PCC plot####
p <- ggplot(keito3, aes(x=CY_PCC, y=ATTED_PCC, color = sample))
p <- p
p <- p + geom_point()
p <- p +  scale_colour_identity(guide = "legend")
p <- p
plot(p)

q <- p + stat_smooth(method = "lm", se = FALSE, colour = "red", size = 1)
plot(q)
a <- lm(formula = CY_PCC~ATTED_PCC, data = keito3)
summary(a)