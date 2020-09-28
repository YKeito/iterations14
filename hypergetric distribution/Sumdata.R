#significance FDR < 0.05
EnrichedGO_Summ_MCLFDR005 <- EnrichedGO_Summ_MCL[EnrichedGO_Summ_MCL$qvalue < 0.05, ]
EnrichedGO_Summ_CY15FDR005 <- EnrichedGO_Summ_CY15[EnrichedGO_Summ_CY15$qvalue < 0.05, ]
EnrichedGO_Summ_CY16FDR005 <- EnrichedGO_Summ_CY16[EnrichedGO_Summ_CY16$qvalue < 0.05, ]
EnrichedGO_Summ_CY20FDR005 <- EnrichedGO_Summ_CY20[EnrichedGO_Summ_CY20$qvalue < 0.05, ]
write.table(EnrichedGO_Summ_MCLFDR005, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/EnrichedGO_Summ_MCLFDR005.txt",sep="\t", col.names=T, row.names=F, append=F, quote=F)
write.table(EnrichedGO_Summ_CY15FDR005, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/EnrichedGO_Summ_CY15FDR005.txt",sep="\t", col.names=T, row.names=F, append=F, quote=F)
write.table(EnrichedGO_Summ_CY16FDR005, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/EnrichedGO_Summ_CY16FDR005.txt",sep="\t", col.names=T, row.names=F, append=F, quote=F)
write.table(EnrichedGO_Summ_CY20FDR005, file = "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/EnrichedGO_Summ_CY20FDR005.txt",sep="\t", col.names=T, row.names=F, append=F, quote=F)


#
#allRNASeq <- read.table("~/Nakano_RNAseq/network_analysis/base/allRNASeq_union.txt", sep = "\t", row.names = 1, header = T)
#allRNASeq <- allRNASeq[, 1:36]
attribute <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/attribute/allattribute.csv", sep = ",", row.names = 1, header = T)
rownames(attribute) <- attribute$name
test <- allRNASeq[match(attribute$name, rownames(allRNASeq)), ]
test2 <- data.frame(test, attribute)
write.table(test2, "~/Nakano_RNAseq/network_analysis/test2.txt", append=F, quote = F, sep = "\t", row.names = T)
#####Venn Diagram CY151620_Numsubcluster####
library(gplots)
#CY15_enriched_Node <- CY15_check_MCL[CY15_enrichment_qvalue < 0.05]
#CY16_enriched_Node <- CY16_check_MCL[CY16_enrichment_qvalue < 0.05]
#CY20_enriched_Node <- CY20_check_MCL[CY20_enrichment_qvalue < 0.05]
CYNode <- list(CY15_SubCluster1_Node = CY15_enriched_Node[[1]], CY16_SubCluster1_Node = CY16_enriched_Node[[1]], CY20_SubCluster1_Node = CY20_enriched_Node[[1]])
CY_Venn <- venn(CYNode)
CY_SubCluster1_Node <- data.frame(unlist(attr(CY_Venn, "intersections")))
write.table(CY_SubCluster1_Node, "~/Nakano_RNAseq/network_analysis/Hypergeometric distribution/CY_SubCluster1_Node.txt", append=F, quote = F, sep = "\t", row.names = T)
