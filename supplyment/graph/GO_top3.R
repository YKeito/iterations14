####nottimecource####
#unique(EnrichedGO_Summ_CY15$CluNum)
#subcluster1
subcluster1_top3 <- EnrichedGO_Summ_MCLFDR005[EnrichedGO_Summ_MCLFDR005$CluNum == 1, ]
subcluster1_top3 <- subcluster1_top3[order(subcluster1_top3$enrichment_score, decreasing = T), ]
subcluster1_top3 <- subcluster1_top3[1:3, ]
library(ggplot2)
g <- ggplot(
  subcluster1_top3,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(FDR)")
g <- g + xlab("")
g <- g + theme(axis.text=element_text(size=60, face="bold"), axis.title=element_text(size=60, face="bold"))
g <- g + theme_bw(base_size = 60)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/iterations14/subcluster/subcluster1_top3_GO.png", plot = g, dpi = 100, width = 21.6, height = 10)
#subcluster2
subcluster2_top3 <- EnrichedGO_Summ_MCLFDR005[EnrichedGO_Summ_MCLFDR005$CluNum == 2, ]
subcluster2_top3 <- subcluster2_top3[order(subcluster2_top3$enrichment_score, decreasing = T), ]
subcluster2_top3 <- subcluster2_top3[1:3, ]
library(ggplot2)
g <- ggplot(
  subcluster2_top3,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(FDR)")
g <- g + xlab("")
g <- g + theme(axis.text=element_text(size=60, face="bold"), axis.title=element_text(size=60, face="bold"))
g <- g + theme_bw(base_size = 60)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/iterations14/subcluster/subcluster2_top3_GO.png", plot = g, dpi = 100, width = 18, height = 10)
#subcluster3
subcluster3_top3 <- EnrichedGO_Summ_MCLFDR005[EnrichedGO_Summ_MCLFDR005$CluNum == 3, ]
subcluster3_top3 <- subcluster3_top3[order(subcluster3_top3$enrichment_score, decreasing = T), ]
subcluster3_top3 <- subcluster3_top3[1:3, ]
library(ggplot2)
g <- ggplot(
  subcluster3_top3,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(FDR)")
g <- g + xlab("")
g <- g + theme(axis.text=element_text(size=60, face="bold"), axis.title=element_text(size=60, face="bold"))
g <- g + theme_bw(base_size = 60)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/iterations14/subcluster/subcluster3_top3_GO.png", plot = g, dpi = 100, width = 22.5, height = 10)
#subcluster4
subcluster4_top3 <- EnrichedGO_Summ_MCLFDR005[EnrichedGO_Summ_MCLFDR005$CluNum == 4, ]
subcluster4_top3 <- subcluster4_top3[order(subcluster4_top3$enrichment_score, decreasing = T), ]
subcluster4_top3 <- subcluster4_top3[1:3, ]
library(ggplot2)
g <- ggplot(
  subcluster4_top3,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(FDR)")
g <- g + xlab("")
g <- g + theme(axis.text=element_text(size=60, face="bold"), axis.title=element_text(size=60, face="bold"))
g <- g + theme_bw(base_size = 60)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/iterations14/subcluster/subcluster4_top3_GO.png", plot = g, dpi = 100, width = 25, height = 10)

#subcluster5
subcluster5_top3 <- EnrichedGO_Summ_MCLFDR005[EnrichedGO_Summ_MCLFDR005$CluNum == 5, ]
subcluster5_top3 <- subcluster5_top3[order(subcluster5_top3$enrichment_score, decreasing = T), ]
subcluster5_top3 <- subcluster5_top3[1:3, ]
library(ggplot2)
g <- ggplot(
  subcluster5_top3,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(FDR)")
g <- g + xlab("")
g <- g + theme(axis.text=element_text(size=60, face="bold"), axis.title=element_text(size=60, face="bold"))
g <- g + theme_bw(base_size = 60)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/iterations14/subcluster/subcluster5_top3_GO.png", plot = g, dpi = 100, width = 25, height = 10)
#subcluster6
subcluster6_top3 <- EnrichedGO_Summ_MCLFDR005[EnrichedGO_Summ_MCLFDR005$CluNum == 6, ]
subcluster6_top3 <- subcluster6_top3[order(subcluster6_top3$enrichment_score, decreasing = T), ]
subcluster6_top3 <- subcluster6_top3[1:3, ]
library(ggplot2)
g <- ggplot(
  subcluster6_top3,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(FDR)")
g <- g + xlab("")
g <- g + theme(axis.text=element_text(size=60, face="bold"), axis.title=element_text(size=60, face="bold"))
g <- g + theme_bw(base_size = 60)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/iterations14/subcluster/subcluster6_top3_GO.png", plot = g, dpi = 100, width = 25, height = 10)

#subcluster7
subcluster7_top3 <- EnrichedGO_Summ_MCLFDR005[EnrichedGO_Summ_MCLFDR005$CluNum == 7, ]
subcluster7_top3 <- subcluster7_top3[order(subcluster7_top3$enrichment_score, decreasing = T), ]
subcluster7_top3 <- subcluster7_top3[1:3, ]
library(ggplot2)
g <- ggplot(
  subcluster7_top3,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(FDR)")
g <- g + xlab("")
g <- g + theme(axis.text=element_text(size=60, face="bold"), axis.title=element_text(size=60, face="bold"))
g <- g + theme_bw(base_size = 60)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/iterations14/subcluster/subcluster7_top3_GO.png", plot = g, dpi = 100, width = 25, height = 10)
#subcluster9
subcluster9_top3 <- EnrichedGO_Summ_MCLFDR005[EnrichedGO_Summ_MCLFDR005$CluNum == 9, ]
subcluster9_top3 <- subcluster9_top3[order(subcluster9_top3$enrichment_score, decreasing = T), ]
subcluster9_top3 <- subcluster9_top3[1:3, ]
library(ggplot2)
g <- ggplot(
  subcluster9_top3,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(FDR)")
g <- g + xlab("")
g <- g + theme(axis.text=element_text(size=60, face="bold"), axis.title=element_text(size=60, face="bold"))
g <- g + theme_bw(base_size = 60)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/iterations14/subcluster/subcluster9_top3_GO.png", plot = g, dpi = 100, width = 25, height = 10)
