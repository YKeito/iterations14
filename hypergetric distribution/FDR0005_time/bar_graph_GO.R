#####bar graph####
####CY15####
#subcluster1_CY15_1h
library(ggplot2)
#NumCY15_1h_subcluster
subcluster1_CY15_1h <- EnrichedGO_Summ_CY15_1h_FDR005[EnrichedGO_Summ_CY15_1h_FDR005$CluNum == 1, ]
subcluster1_CY15_1h <- subcluster1_CY15_1h[1:10, ]
g <- ggplot(
  subcluster1_CY15_1h,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(q-value)")
g <- g + xlab("GO term")
g <- g + labs(title = "subcluster1 CY15 1h GO")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/subcluster1_CY15_1h_GO.png", plot = g)
#subcluster86_CY15_1h
library(ggplot2)
#NumCY15_1h_subcluster
subcluster86_CY15_1h <- EnrichedGO_Summ_CY15_1h_FDR005[EnrichedGO_Summ_CY15_1h_FDR005$CluNum == 86, ]
subcluster86_CY15_1h <- subcluster86_CY15_1h[1:10, ]
g <- ggplot(
  subcluster86_CY15_1h,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(q-value)")
g <- g + xlab("GO term")
g <- g + labs(title = "subcluster86 CY15 1h GO")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/subcluster86_CY15_1h_GO.png", plot = g)

#subcluster1_CY15_3h
#NumCY15_3h_subcluster
library(ggplot2)
subcluster1_CY15_3h <- EnrichedGO_Summ_CY15_3h_FDR005[EnrichedGO_Summ_CY15_3h_FDR005$CluNum == 1, ]
subcluster1_CY15_3h <- subcluster1_CY15_3h[1:10, ]
g <- ggplot(
  subcluster1_CY15_3h,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(q-value)")
g <- g + xlab("GO term")
g <- g + labs(title = "subcluster1 CY15 3h GO")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/subcluster1_CY15_3h_GO.png", plot = g)
#subcluster62_CY15_3h
#NumCY15_3h_subcluster
library(ggplot2)
subcluster62_CY15_3h <- EnrichedGO_Summ_CY15_3h_FDR005[EnrichedGO_Summ_CY15_3h_FDR005$CluNum == 62, ]
subcluster62_CY15_3h <- subcluster62_CY15_3h[1:10, ]
g <- ggplot(
  subcluster62_CY15_3h,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(q-value)")
g <- g + xlab("GO term")
g <- g + labs(title = "subcluster62 CY15 3h GO")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/subcluster62_CY15_3h_GO.png", plot = g)


#subcluster1_CY15_12h
#NumCY15_12h_subcluster
library(ggplot2)
subcluster1_CY15_12h <- EnrichedGO_Summ_CY15_12h_FDR005[EnrichedGO_Summ_CY15_12h_FDR005$CluNum == 1, ]
subcluster1_CY15_12h <- subcluster1_CY15_12h[1:10, ]
g <- ggplot(
  subcluster1_CY15_12h,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(q-value)")
g <- g + xlab("GO term")
g <- g + labs(title = "subcluster1 CY15 12h GO")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/subcluster1_CY15_12h_GO.png", plot = g)
#subcluster12_CY15_12h
#NumCY15_12h_subcluster
library(ggplot2)
subcluster12_CY15_12h <- EnrichedGO_Summ_CY15_12h_FDR005[EnrichedGO_Summ_CY15_12h_FDR005$CluNum == 12, ]
subcluster12_CY15_12h <- subcluster12_CY15_12h[1:10, ]
g <- ggplot(
  subcluster12_CY15_12h,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(q-value)")
g <- g + xlab("GO term")
g <- g + labs(title = "subcluster12 CY15 12h GO")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/subcluster12_CY15_12h_GO.png", plot = g)


#subcluster1_CY15_24h
#NumCY15_24h_subcluster
library(ggplot2)
subcluster1_CY15_24h <- EnrichedGO_Summ_CY15_24h_FDR005[EnrichedGO_Summ_CY15_24h_FDR005$CluNum == 1, ]
subcluster1_CY15_24h <- subcluster1_CY15_24h[1:10, ]
g <- ggplot(
  subcluster1_CY15_24h,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(q-value)")
g <- g + xlab("GO term")
g <- g + labs(title = "subcluster1 CY15 24h GO")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/subcluster1_CY15_24h_GO.png", plot = g)

####CY16####
#subcluster2_CY16_1h
library(ggplot2)
#NumCY16_1h_subcluster
subcluster2_CY16_1h <- EnrichedGO_Summ_CY16_1h_FDR005[EnrichedGO_Summ_CY16_1h_FDR005$CluNum == 2, ]
subcluster2_CY16_1h <- subcluster2_CY16_1h[1:10, ]
g <- ggplot(
  subcluster2_CY16_1h,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(q-value)")
g <- g + xlab("GO term")
g <- g + labs(title = "subcluster2 CY16 1h GO")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/subcluster2_CY16_1h_GO.png", plot = g)
#subcluster1_CY16_3h
#NumCY16_3h_subcluster
library(ggplot2)
subcluster1_CY16_3h <- EnrichedGO_Summ_CY16_3h_FDR005[EnrichedGO_Summ_CY16_3h_FDR005$CluNum == 1, ]
subcluster1_CY16_3h <- subcluster1_CY16_3h[1:10, ]
g <- ggplot(
  subcluster1_CY16_3h,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(q-value)")
g <- g + xlab("GO term")
g <- g + labs(title = "subcluster1 CY16 3h GO")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/subcluster1_CY16_3h_GO.png", plot = g)
#subcluster1_CY16_12h
#NumCY16_12h_subcluster
library(ggplot2)
subcluster1_CY16_12h <- EnrichedGO_Summ_CY16_12h_FDR005[EnrichedGO_Summ_CY16_12h_FDR005$CluNum == 1, ]
subcluster1_CY16_12h <- subcluster1_CY16_12h[1:10, ]
g <- ggplot(
  subcluster1_CY16_12h,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(q-value)")
g <- g + xlab("GO term")
g <- g + labs(title = "subcluster1 CY16 12h GO")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/subcluster1_CY16_12h_GO.png", plot = g)
#subcluster1_CY16_24h
#NumCY16_24h_subcluster
library(ggplot2)
subcluster1_CY16_24h <- EnrichedGO_Summ_CY16_24h_FDR005[EnrichedGO_Summ_CY16_24h_FDR005$CluNum == 1, ]
subcluster1_CY16_24h <- subcluster1_CY16_24h[1:10, ]
g <- ggplot(
  subcluster1_CY16_24h,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(q-value)")
g <- g + xlab("GO term")
g <- g + labs(title = "subcluster1 CY16 24h GO")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/subcluster1_CY16_24h_GO.png", plot = g)
#subcluster14_CY16_24h
#NumCY16_24h_subcluster
library(ggplot2)
subcluster14_CY16_24h <- EnrichedGO_Summ_CY16_24h_FDR005[EnrichedGO_Summ_CY16_24h_FDR005$CluNum == 14, ]
subcluster14_CY16_24h <- subcluster14_CY16_24h[1:10, ]
g <- ggplot(
  subcluster14_CY16_24h,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(q-value)")
g <- g + xlab("GO term")
g <- g + labs(title = "subcluster14 CY16 24h GO")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/subcluster14_CY16_24h_GO.png", plot = g)



####CY20####
#subcluster1_CY20_24h
#NumCY20_24h_subcluster
library(ggplot2)
subcluster1_CY20_24h <- EnrichedGO_Summ_CY20_24h_FDR005[EnrichedGO_Summ_CY20_24h_FDR005$CluNum == 1, ]
subcluster1_CY20_24h <- subcluster1_CY20_24h[1:10, ]
g <- ggplot(
  subcluster1_CY20_24h,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(q-value)")
g <- g + xlab("GO term")
g <- g + labs(title = "subcluster1 CY20 24h GO")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/subcluster1_CY20_24h_GO.png", plot = g)
####CY20_48h####
#subcluster3_CY20_48h
#NumCY20_48h_subcluster
library(ggplot2)
subcluster3_CY20_48h <- EnrichedGO_Summ_CY20_48h_FDR005[EnrichedGO_Summ_CY20_48h_FDR005$CluNum == 3, ]
subcluster3_CY20_48h <- subcluster3_CY20_48h[1:10, ]
g <- ggplot(
  subcluster3_CY20_48h,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(q-value)")
g <- g + xlab("GO term")
g <- g + labs(title = "subcluster3 CY20 48h GO")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/subcluster3_CY20_48h_GO.png", plot = g)

#subcluster7_CY20_48h
#NumCY20_48h_subcluster
library(ggplot2)
subcluster7_CY20_48h <- EnrichedGO_Summ_CY20_48h_FDR005[EnrichedGO_Summ_CY20_48h_FDR005$CluNum == 3, ]
subcluster7_CY20_48h <- subcluster7_CY20_48h[1:10, ]
g <- ggplot(
  subcluster7_CY20_48h,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(q-value)")
g <- g + xlab("GO term")
g <- g + labs(title = "subcluster7 CY20 48h GO")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/subcluster7_CY20_48h_GO.png", plot = g)

#subcluster17_CY20_48h
#NumCY20_48h_subcluster
library(ggplot2)
subcluster17_CY20_48h <- EnrichedGO_Summ_CY20_48h_FDR005[EnrichedGO_Summ_CY20_48h_FDR005$CluNum == 17, ]
subcluster17_CY20_48h <- subcluster17_CY20_48h[1:10, ]
g <- ggplot(
  subcluster17_CY20_48h,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(q-value)")
g <- g + xlab("GO term")
g <- g + labs(title = "subcluster17 CY20 48h GO")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/subcluster17_CY20_48h_GO.png", plot = g)

#subcluster35_CY20_48h
#NumCY20_48h_subcluster
library(ggplot2)
subcluster35_CY20_48h <- EnrichedGO_Summ_CY20_48h_FDR005[EnrichedGO_Summ_CY20_48h_FDR005$CluNum == 35, ]
subcluster35_CY20_48h <- subcluster35_CY20_48h[1:10, ]
g <- ggplot(
  subcluster35_CY20_48h,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(q-value)")
g <- g + xlab("GO term")
g <- g + labs(title = "subcluster35 CY20 48h GO")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/subcluster35_CY20_48h_GO.png", plot = g)

#subcluster42_CY20_48h
#NumCY20_48h_subcluster
library(ggplot2)
subcluster42_CY20_48h <- EnrichedGO_Summ_CY20_48h_FDR005[EnrichedGO_Summ_CY20_48h_FDR005$CluNum == 35, ]
subcluster42_CY20_48h <- subcluster42_CY20_48h[1:10, ]
g <- ggplot(
  subcluster42_CY20_48h,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(q-value)")
g <- g + xlab("GO term")
g <- g + labs(title = "subcluster42 CY20 48h GO")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/subcluster42_CY20_48h_GO.png", plot = g)

#subcluster51_CY20_48h
#NumCY20_48h_subcluster
library(ggplot2)
subcluster51_CY20_48h <- EnrichedGO_Summ_CY20_48h_FDR005[EnrichedGO_Summ_CY20_48h_FDR005$CluNum == 51, ]
subcluster51_CY20_48h <- subcluster51_CY20_48h[1:10, ]
g <- ggplot(
  subcluster51_CY20_48h,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(q-value)")
g <- g + xlab("GO term")
g <- g + labs(title = "subcluster51 CY20 48h GO")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/subcluster51_CY20_48h_GO.png", plot = g)

#subcluster76_CY20_48h
#NumCY20_48h_subcluster
library(ggplot2)
subcluster76_CY20_48h <- EnrichedGO_Summ_CY20_48h_FDR005[EnrichedGO_Summ_CY20_48h_FDR005$CluNum == 76, ]
subcluster76_CY20_48h <- subcluster76_CY20_48h[1:10, ]
g <- ggplot(
  subcluster76_CY20_48h,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(q-value)")
g <- g + xlab("GO term")
g <- g + labs(title = "subcluster76 CY20 48h GO")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/subcluster76_CY20_48h_GO.png", plot = g)

#subcluster139_CY20_48h
#NumCY20_48h_subcluster
library(ggplot2)
subcluster139_CY20_48h <- EnrichedGO_Summ_CY20_48h_FDR005[EnrichedGO_Summ_CY20_48h_FDR005$CluNum == 139, ]
subcluster139_CY20_48h <- subcluster139_CY20_48h[1:10, ]
g <- ggplot(
  subcluster139_CY20_48h,
  aes(
    x = reorder(GO_term, enrichment_score),
    y = enrichment_score
  )
)
g <- g + geom_bar(stat = "identity")
g <- g + coord_flip()
g <- g + ylab("-log2(q-value)")
g <- g + xlab("GO term")
g <- g + labs(title = "subcluster139 CY20 48h GO")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
print(g)
ggsave(file = "~/Nakano_RNAseq/network_analysis/GO_results_Fig/subcluster139_CY20_48h_GO.png", plot = g)


