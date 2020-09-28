####CY15####
n <- 1
CY15_1h_Node <- c()
CY15_1h_annotation <- c()
CY15_1h_expression <- c()
CY15_1h_Nodes <- list()


for(n in 1:10){
  CY15_1h_Node <- unlist(strsplit(as.character(unlist(EnrichedGO_Summ_CY15_1h_FDR005$Node[n])), " "))
  CY15_1h_Node <- CY15_1h_Node[nchar(CY15_1h_Node) == 9]
  temp <- as.character(unlist(allRNASeq$annotation[match(CY15_1h_Node, rownames(allRNASeq))]))
  CY15_1h_annotation <- c(CY15_1h_annotation, temp)
  CY15_1h_expression <- c(CY15_1h_expression, allRNASeq$CY15_1h[match(CY15_1h_Node, rownames(allRNASeq))])
  CY15_1h_Nodes <- c(CY15_1h_Nodes, list(CY15_1h_Node))
  CY15_1h_Node <- c()
  n <- n+1
}
CY15_1h_GOexp <- data.frame(CY15_1h_AGI = unlist(CY15_1h_Nodes),
                            CY15_1h_annotation = CY15_1h_annotation,
                            CY15_1h_expression = CY15_1h_expression)
n <- 1
CY15_3h_Node <- c()
CY15_3h_annotation <- c()
CY15_3h_expression <- c()
CY15_3h_Nodes <- list()

for(n in 1:10){
  CY15_3h_Node <- unlist(strsplit(as.character(unlist(EnrichedGO_Summ_CY15_3h_FDR005$Node[n])), " "))
  CY15_3h_Node <- CY15_3h_Node[nchar(CY15_3h_Node) == 9]
  temp <- as.character(unlist(allRNASeq$annotation[match(CY15_3h_Node, rownames(allRNASeq))]))
  CY15_3h_annotation <- c(CY15_3h_annotation, temp)
  CY15_3h_expression <- c(CY15_3h_expression, allRNASeq$CY15_3h[match(CY15_3h_Node, rownames(allRNASeq))])
  CY15_3h_Nodes <- c(CY15_3h_Nodes, list(CY15_3h_Node))
  CY15_3h_Node <- c()
  n <- n+1
}
CY15_3h_GOexp <- data.frame(CY15_3h_AGI = unlist(CY15_3h_Nodes),
                            CY15_3h_annotation = CY15_3h_annotation,
                            CY15_3h_expression = CY15_3h_expression)
n <- 1
CY15_12h_Node <- c()
CY15_12h_annotation <- c()
CY15_12h_expression <- c()
CY15_12h_Nodes <- list()

for(n in 1:10){
  CY15_12h_Node <- unlist(strsplit(as.character(unlist(EnrichedGO_Summ_CY15_12h_FDR005$Node[n])), " "))
  CY15_12h_Node <- CY15_12h_Node[nchar(CY15_12h_Node) == 9]
  temp <- as.character(unlist(allRNASeq$annotation[match(CY15_12h_Node, rownames(allRNASeq))]))
  CY15_12h_annotation <- c(CY15_12h_annotation, temp)
  CY15_12h_expression <- c(CY15_12h_expression, allRNASeq$CY15_12h[match(CY15_12h_Node, rownames(allRNASeq))])
  CY15_12h_Nodes <- c(CY15_12h_Nodes, list(CY15_12h_Node))
  CY15_12h_Node <- c()
  n <- n+1
}
CY15_12h_GOexp <- data.frame(CY15_12h_AGI = unlist(CY15_12h_Nodes),
                             CY15_12h_annotation = CY15_12h_annotation,
                             CY15_12h_expression = CY15_12h_expression)
n <- 1
CY15_24h_Node <- c()
CY15_24h_annotation <- c()
CY15_24h_expression <- c()
CY15_24h_Nodes <- list()

for(n in 1:10){
  CY15_24h_Node <- unlist(strsplit(as.character(unlist(EnrichedGO_Summ_CY15_24h_FDR005$Node[n])), " "))
  CY15_24h_Node <- CY15_24h_Node[nchar(CY15_24h_Node) == 9]
  temp <- as.character(unlist(allRNASeq$annotation[match(CY15_24h_Node, rownames(allRNASeq))]))
  CY15_24h_annotation <- c(CY15_24h_annotation, temp)
  CY15_24h_expression <- c(CY15_24h_expression, allRNASeq$CY15_24h[match(CY15_24h_Node, rownames(allRNASeq))])
  CY15_24h_Nodes <- c(CY15_24h_Nodes, list(CY15_24h_Node))
  CY15_24h_Node <- c()
  n <- n+1
}
CY15_24h_GOexp <- data.frame(CY15_24h_AGI = unlist(CY15_24h_Nodes),
                             CY15_24h_annotation = CY15_24h_annotation,
                             CY15_24h_expression = CY15_24h_expression)
####CY16####
n <- 1
CY16_1h_Node <- c()
CY16_1h_annotation <- c()
CY16_1h_expression <- c()
CY16_1h_Nodes <- list()

for(n in 1:10){
  CY16_1h_Node <- unlist(strsplit(as.character(unlist(EnrichedGO_Summ_CY16_1h_FDR005$Node[n])), " "))
  CY16_1h_Node <- CY16_1h_Node[nchar(CY16_1h_Node) == 9]
  temp <- as.character(unlist(allRNASeq$annotation[match(CY16_1h_Node, rownames(allRNASeq))]))
  CY16_1h_annotation <- c(CY16_1h_annotation, temp)
  CY16_1h_expression <- c(CY16_1h_expression, allRNASeq$CY16_1h[match(CY16_1h_Node, rownames(allRNASeq))])
  CY16_1h_Nodes <- c(CY16_1h_Nodes, list(CY16_1h_Node))
  CY16_1h_Node <- c()
  n <- n+1
}
CY16_1h_GOexp <- data.frame(CY16_1h_AGI = unlist(CY16_1h_Nodes),
                            CY16_1h_annotation = CY16_1h_annotation,
                            CY16_1h_expression = CY16_1h_expression)
n <- 1
CY16_3h_Node <- c()
CY16_3h_annotation <- c()
CY16_3h_expression <- c()
CY16_3h_Nodes <- list()

for(n in 1:10){
  CY16_3h_Node <- unlist(strsplit(as.character(unlist(EnrichedGO_Summ_CY16_3h_FDR005$Node[n])), " "))
  CY16_3h_Node <- CY16_3h_Node[nchar(CY16_3h_Node) == 9]
  temp <- as.character(unlist(allRNASeq$annotation[match(CY16_3h_Node, rownames(allRNASeq))]))
  CY16_3h_annotation <- c(CY16_3h_annotation, temp)
  CY16_3h_expression <- c(CY16_3h_expression, allRNASeq$CY16_3h[match(CY16_3h_Node, rownames(allRNASeq))])
  CY16_3h_Nodes <- c(CY16_3h_Nodes, list(CY16_3h_Node))
  CY16_3h_Node <- c()
  n <- n+1
}
CY16_3h_GOexp <- data.frame(CY16_3h_AGI = unlist(CY16_3h_Nodes),
                            CY16_3h_annotation = CY16_3h_annotation,
                            CY16_3h_expression = CY16_3h_expression)
n <- 1
CY16_12h_Node <- c()
CY16_12h_annotation <- c()
CY16_12h_expression <- c()
CY16_12h_Nodes <- list()

for(n in 1:10){
  CY16_12h_Node <- unlist(strsplit(as.character(unlist(EnrichedGO_Summ_CY16_12h_FDR005$Node[n])), " "))
  CY16_12h_Node <- CY16_12h_Node[nchar(CY16_12h_Node) == 9]
  temp <- as.character(unlist(allRNASeq$annotation[match(CY16_12h_Node, rownames(allRNASeq))]))
  CY16_12h_annotation <- c(CY16_12h_annotation, temp)
  CY16_12h_expression <- c(CY16_12h_expression, allRNASeq$CY16_12h[match(CY16_12h_Node, rownames(allRNASeq))])
  CY16_12h_Nodes <- c(CY16_12h_Nodes, list(CY16_12h_Node))
  CY16_12h_Node <- c()
  n <- n+1
}
CY16_12h_GOexp <- data.frame(CY16_12h_AGI = unlist(CY16_12h_Nodes),
                             CY16_12h_annotation = CY16_12h_annotation,
                             CY16_12h_expression = CY16_12h_expression)
n <- 1
CY16_24h_Node <- c()
CY16_24h_annotation <- c()
CY16_24h_expression <- c()
CY16_24h_Nodes <- list()

for(n in 1:10){
  CY16_24h_Node <- unlist(strsplit(as.character(unlist(EnrichedGO_Summ_CY16_24h_FDR005$Node[n])), " "))
  CY16_24h_Node <- CY16_24h_Node[nchar(CY16_24h_Node) == 9]
  temp <- as.character(unlist(allRNASeq$annotation[match(CY16_24h_Node, rownames(allRNASeq))]))
  CY16_24h_annotation <- c(CY16_24h_annotation, temp)
  CY16_24h_expression <- c(CY16_24h_expression, allRNASeq$CY16_24h[match(CY16_24h_Node, rownames(allRNASeq))])
  CY16_24h_Nodes <- c(CY16_24h_Nodes, list(CY16_24h_Node))
  CY16_24h_Node <- c()
  n <- n+1
}
CY16_24h_GOexp <- data.frame(CY16_24h_AGI = unlist(CY16_24h_Nodes),
                             CY16_24h_annotation = CY16_24h_annotation,
                             CY16_24h_expression = CY16_24h_expression)
####CY20####
n <- 1
CY20_24h_Node <- c()
CY20_24h_annotation <- c()
CY20_24h_expression <- c()
CY20_24h_Nodes <- list()

for(n in 1:10){
  CY20_24h_Node <- unlist(strsplit(as.character(unlist(EnrichedGO_Summ_CY20_24h_FDR005$Node[n])), " "))
  CY20_24h_Node <- CY20_24h_Node[nchar(CY20_24h_Node) == 9]
  temp <- as.character(unlist(allRNASeq$annotation[match(CY20_24h_Node, rownames(allRNASeq))]))
  CY20_24h_annotation <- c(CY20_24h_annotation, temp)
  CY20_24h_expression <- c(CY20_24h_expression, allRNASeq$CY20_24h[match(CY20_24h_Node, rownames(allRNASeq))])
  CY20_24h_Nodes <- c(CY20_24h_Nodes, list(CY20_24h_Node))
  CY20_24h_Node <- c()
  n <- n+1
}
CY20_24h_GOexp <- data.frame(CY20_24h_AGI = unlist(CY20_24h_Nodes),
                             CY20_24h_annotation = CY20_24h_annotation,
                             CY20_24h_expression = CY20_24h_expression)
n <- 1
CY20_48h_Node <- c()
CY20_48h_annotation <- c()
CY20_48h_expression <- c()
CY20_48h_Nodes <- list()

for(n in 1:10){
  CY20_48h_Node <- unlist(strsplit(as.character(unlist(EnrichedGO_Summ_CY20_48h_FDR005$Node[n])), " "))
  CY20_48h_Node <- CY20_48h_Node[nchar(CY20_48h_Node) == 9]
  temp <- as.character(unlist(allRNASeq$annotation[match(CY20_48h_Node, rownames(allRNASeq))]))
  CY20_48h_annotation <- c(CY20_48h_annotation, temp)
  CY20_48h_expression <- c(CY20_48h_expression, allRNASeq$CY20_48h[match(CY20_48h_Node, rownames(allRNASeq))])
  CY20_48h_Nodes <- c(CY20_48h_Nodes, list(CY20_48h_Node))
  CY20_48h_Node <- c()
  n <- n+1
}
CY20_48h_GOexp <- data.frame(CY20_48h_AGI = unlist(CY20_48h_Nodes),
                             CY20_48h_annotation = CY20_48h_annotation,
                             CY20_48h_expression = CY20_48h_expression)


####並び替え####
CY20_24h_GOexp
CY20_24h_GOexp <- CY20_24h_GOexp[order(CY20_24h_GOexp$CY20_24h_expression, decreasing = T), ]
