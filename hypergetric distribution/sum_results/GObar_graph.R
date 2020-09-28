####nottimecource####
#unique(EnrichedGO_Summ_CY15$CluNum)
#subcluster1
subcluster1_top10 <- EnrichedGO_Summ_MCLFDR005[EnrichedGO_Summ_MCLFDR005$CluNum == 1, ]
subcluster1_top10 <- subcluster1_top10[order(subcluster1_top10$enrichment_score, decreasing = T), ]
subcluster1_top10 <- subcluster1_top10[1:10, ]
library(ggplot2)
g <- ggplot(
  subcluster1_top10,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(FDR)")
g <- g + xlab("")
g <- g + labs(title = "subcluster1 GO")
g <- g + theme(axis.text=element_text(size=60), axis.title=element_text(size=60,face="bold"))
g <- g + theme_bw(base_size = 60)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/iterations14/subcluster/subcluster1_top10_GO.png", plot = g, dpi = 100, width = 21.6, height = 12)
#subcluster2
subcluster2_top10 <- EnrichedGO_Summ_MCLFDR005[EnrichedGO_Summ_MCLFDR005$CluNum == 2, ]
subcluster2_top10 <- subcluster2_top10[order(subcluster2_top10$enrichment_score, decreasing = T), ]
subcluster2_top10 <- subcluster2_top10[1:10, ]
library(ggplot2)
g <- ggplot(
  subcluster2_top10,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(FDR)")
g <- g + xlab("")
g <- g + labs(title = "subcluster2 GO")
g <- g + theme(axis.text=element_text(size=60), axis.title=element_text(size=60,face="bold"))
g <- g + theme_bw(base_size = 60)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/iterations14/subcluster/subcluster2_top10_GO.png", plot = g, dpi = 100, width = 28.8, height = 20)
#subcluster3
subcluster3_top10 <- EnrichedGO_Summ_MCLFDR005[EnrichedGO_Summ_MCLFDR005$CluNum == 3, ]
subcluster3_top10 <- subcluster3_top10[order(subcluster3_top10$enrichment_score, decreasing = T), ]
subcluster3_top10 <- subcluster3_top10[1:10, ]
library(ggplot2)
subcluster3_top10$GO_term <- c("nitrile biosynthetic process", 
                               "cytoplasmic stress granule", 
                               "response to molecule of fungal origin", 
                               "defense response signaling pathway", 
                               "glucosinolate catabolic process", 
                               "monovalent inorganic cation transport", 
                               "regulation of pH", 
                               "monovalent cation:proton antiporter activity", 
                               "solute:proton antiporter activity", 
                               "response to ozone")
g <- ggplot(
  subcluster3_top10,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(FDR)")
g <- g + xlab("")
g <- g + labs(title = "subcluster3 GO")
g <- g + theme(axis.text=element_text(size=60), axis.title=element_text(size=60,face="bold"))
g <- g + theme_bw(base_size = 60)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/iterations14/subcluster/subcluster3_top10_GO.png", plot = g, dpi = 100, width = 25, height = 15)
#subcluster7
subcluster7_top10 <- EnrichedGO_Summ_MCLFDR005[EnrichedGO_Summ_MCLFDR005$CluNum == 7, ]
subcluster7_top10 <- subcluster7_top10[order(subcluster7_top10$enrichment_score, decreasing = T), ]
subcluster7_top10 <- subcluster7_top10[1:10, ]
library(ggplot2)
g <- ggplot(
  subcluster7_top10,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(FDR)")
g <- g + xlab("")
g <- g + labs(title = "subcluster7 GO")
g <- g + theme(axis.text=element_text(size=60), axis.title=element_text(size=60,face="bold"))
g <- g + theme_bw(base_size = 60)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/iterations14/subcluster/subcluster7_top10_GO.png", plot = g, dpi = 100, width = 25, height = 18)
