"~/Nakano_RNAseq/network_analysis/script/CY151620ExP_CoExpNet/"
#load package----
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
#input base data----
CY.data <- read.table("~/Nakano_RNAseq/network_analysis/base/CY_base/allRNASeq_union.txt", sep = "\t", header = T)
####calculate PCC####
allRNASeq_foldchange <- allRNASeq[, 1:15]
n <- 1
m <- 2
PCC <- c()
PCC_all <- list()
for(n in n:c(nrow(allRNASeq_foldchange)-1)){
  for(m in m:nrow(allRNASeq_foldchange)){
    PCC <- c(PCC, cor(as.numeric(allRNASeq_foldchange[n, ]), as.numeric(allRNASeq_foldchange[m, ]), method = "pearson"))
    m <- m+1
  }
  PCC_all <- c(PCC_all, list(PCC))
  PCC <- c()
  print(n)
  n <- n+1
  m <- n+1
}
####PCC p-value####
allRNASeq_foldchange <- allRNASeq[, 1:15]
n <- 1
m <- 2
PCC_pvalue <- c()
PCC_pvalue_all <- list()
for(n in n:c(nrow(allRNASeq_foldchange)-1)){
  for(m in m:nrow(allRNASeq_foldchange)){
    PCC_pvalue <- c(PCC_pvalue, cor.test(as.numeric(allRNASeq_foldchange[n, ]), as.numeric(allRNASeq_foldchange[m, ]), method = "pearson")$p.value)
    m <- m+1
  }
  PCC_pvalue_all <- c(PCC_pvalue_all, list(PCC_pvalue))
  PCC_pvalue <- c()
  print(n)
  n <- n+1
  m <- n+1
}
####PCC q-value####
PCC_qvalue_all <- p.adjust(unlist(PCC_pvalue_all), method = "BH")
#####Cytoscape_format####
#source genes
m <- 1
source_gene <- c()
source_genes <- list()

for(m in m:c(nrow(allRNASeq)-1)){
  source_gene <- rep(rownames(allRNASeq_foldchange)[m], times = c(nrow(allRNASeq)-m))
  source_genes <- c(source_genes, list(source_gene))
  source_gene <- c()
  print(m)
  m <- m+1
}
#target genes#
n <- 2
target_gene <- c()
target_genes <- list()

for(n in n:nrow(allRNASeq)){
  target_gene <- rownames(allRNASeq)[n:nrow(allRNASeq)]
  target_genes <- c(target_genes, list(target_gene))
  target_gene <- c()
  print(n)
  n <- n+1
}
allRNASeq_foldchange_cytoscape <- data.frame(source_genes = unlist(source_genes), 
                                      interaction_value = 1/abs(unlist(PCC_all)), 
                                      target_genes = unlist(target_genes),
                                      p_value = unlist(PCC_pvalue_all),
                                      q_value = PCC_qvalue_all
)
allRNASeq_cytoscape_th <- allRNASeq_foldchange_cytoscape[allRNASeq_foldchange_cytoscape$q_value < 0.05, ]

####attribute####
arabidopsis_TF <- read.table("~/Nakano_RNAseq/network_analysis/base/TF_arabidopsis.txt", sep = "\t", row.names = 1, header = T)
attrarabidopsis_TF <- data.frame(AGI = rownames(arabidopsis_TF),
                                 TF = rep("TF", times = nrow(arabidopsis_TF))
                                 )
####allCY####
#CY15_FDR005
CY15_FDR005 <- allRNASeq[allRNASeq$CY15_total != 0, ]
#CY16_FDR005
CY16_FDR005 <- allRNASeq[allRNASeq$CY16_total != 0, ]
#CY20_FDR005
CY20_FDR005 <- allRNASeq[allRNASeq$CY20_total != 0, ]
#Venn Diagram#
#allCY#
library(gplots)
CY <- list(CY15_FDR005 = rownames(CY15_FDR005), CY16_FDR005 = rownames(CY16_FDR005), CY20_FDR005 = rownames(CY20_FDR005))
CY_Venn <- venn(CY)
CY_FDR005 <- unlist(attr(CY_Venn,"intersections"))
CY_FDR005 <- data.frame(CY_FDR005)

#CY15_FDR0005
CY15_FDR0005 <- CY15_FDR005 <- allRNASeq[(allRNASeq$CY15_1h_q_value < 0.005 | allRNASeq$CY15_3h_q_value < 0.005 | allRNASeq$CY15_12h_q_value < 0.005 | allRNASeq$CY15_24h_q_value < 0.005 | allRNASeq$CY15_48h_q_value < 0.005), ]
#CY16_FDR0005
CY16_FDR0005 <- CY16_FDR005 <- allRNASeq[(allRNASeq$CY16_1h_q_value < 0.005 | allRNASeq$CY16_3h_q_value < 0.005 | allRNASeq$CY16_12h_q_value < 0.005 | allRNASeq$CY16_24h_q_value < 0.005 | allRNASeq$CY16_48h_q_value < 0.005), ]
#CY20_FDR0005
CY20_FDR0005 <- CY20_FDR005 <- allRNASeq[(allRNASeq$CY20_1h_q_value < 0.005 | allRNASeq$CY20_3h_q_value < 0.005 | allRNASeq$CY20_12h_q_value < 0.005 | allRNASeq$CY20_24h_q_value < 0.005 | allRNASeq$CY20_48h_q_value < 0.005), ]
#Venn Diagram#
#allCY#
library(gplots)
CY <- list(CY15_FDR0005 = rownames(CY15_FDR0005), CY16_FDR0005 = rownames(CY16_FDR0005), CY20_FDR0005 = rownames(CY20_FDR0005))
CY_Venn <- venn(CY)
CY_FDR0005 <- unlist(attr(CY_Venn,"intersections"))
CY_FDR0005 <- data.frame(CY_FDR0005)


####timedata_FDR005_union#####
#CY15_time
CY15_1h_FDR005 <- allRNASeq[allRNASeq$CY15_1h_q_value < 0.05, ]
CY15_3h_FDR005 <- allRNASeq[allRNASeq$CY15_3h_q_value < 0.05, ]
CY15_12h_FDR005 <- allRNASeq[allRNASeq$CY15_12h_q_value < 0.05, ]
CY15_24h_FDR005 <- allRNASeq[allRNASeq$CY15_24h_q_value < 0.05, ]
CY15_48h_FDR005 <- allRNASeq[allRNASeq$CY15_48h_q_value < 0.05, ]
#CY16_time
CY16_1h_FDR005 <- allRNASeq[allRNASeq$CY16_1h_q_value < 0.05, ]
CY16_3h_FDR005 <- allRNASeq[allRNASeq$CY16_3h_q_value < 0.05, ]
CY16_12h_FDR005 <- allRNASeq[allRNASeq$CY16_12h_q_value < 0.05, ]
CY16_24h_FDR005 <- allRNASeq[allRNASeq$CY16_24h_q_value < 0.05, ]
CY16_48h_FDR005 <- allRNASeq[allRNASeq$CY16_48h_q_value < 0.05, ]
#CY20_time
CY20_1h_FDR005 <- allRNASeq[allRNASeq$CY20_1h_q_value < 0.05, ]
CY20_3h_FDR005 <- allRNASeq[allRNASeq$CY20_3h_q_value < 0.05, ]
CY20_12h_FDR005 <- allRNASeq[allRNASeq$CY20_12h_q_value < 0.05, ]
CY20_24h_FDR005 <- allRNASeq[allRNASeq$CY20_24h_q_value < 0.05, ]
CY20_48h_FDR005 <- allRNASeq[allRNASeq$CY20_48h_q_value < 0.05, ]
#Venn Diagram#
#CY15#
library(gplots)
CY15_data <- list(CY15_1h = rownames(CY15_1h_FDR005), CY15_3h = rownames(CY15_3h_FDR005), CY15_12h = rownames(CY15_12h_FDR005), CY15_24h = rownames(CY15_24h_FDR005), CY15_48h = rownames(CY15_48h_FDR005))
CY15_Venn <- venn(CY15_data)
CY15_Venn <- unlist(attr(CY15_Venn,"intersections"))
temp_all <- strsplit(names(CY15_Venn), "h")  #文字列の分割#overlapしている数が複数だとその数字が表示されてしまうので除きたい。例CY15_1h53　この53を除くために以下の作業
CY15_all <- c()
n <- 1
for(n in 1:length(temp_all)){
  temp <- temp_all[[n]]
  temp_nChar <- nchar(temp)     #nchar文字列の長さを調べる    temp_all[[1]]    [1] "CY15_1"  ":CY15_3" "1"   nchar(temp_all[[1]])  [1] 6 7 1
  CY15_all <- c(CY15_all, paste0(paste(temp[temp_nChar >= 6], collapse="h"), "h"))      #>=6     最低の文字列は6（CY15_1）だから
}
names(CY15_Venn) <- CY15_all
#attribute_file#
n <- 1
T_AGI <- c()
T_attribute <- c()
for(n in 1:length(unique(CY15_all))){
  T_AGI <- c(T_AGI, CY15_Venn[names(CY15_Venn) == unique(CY15_all)[n]])
  temp <- switch(unique(CY15_all)[n],
                 "CY15_1h:CY15_3h" = 6,   #[1]
                 "CY15_1h:CY15_12h" = 7,   #[2]
                 "CY15_1h:CY15_24h" = 8,   #[3]
                 "CY15_1h:CY15_48h" = 9,    #[4]
                 "CY15_3h:CY15_12h" = 10,   #[5]
                 "CY15_3h:CY15_24h" = 11,   #[6]
                 "CY15_3h:CY15_48h" = 12,   #[7]
                 "CY15_12h:CY15_24h" = 13,   #[8]
                 "CY15_12h:CY15_48h" = 14,   #[9]
                 "CY15_24h:CY15_48h" = 15,   #[10]
                 "CY15_1h:CY15_3h:CY15_12h" = 16,   #[11]
                 "CY15_1h:CY15_3h:CY15_24h" = 17,   #[12]
                 "CY15_1h:CY15_3h:CY15_48h" = 18,   #[13]
                 "CY15_1h:CY15_12h:CY15_24h" = 19,   #[14]
                 "CY15_3h:CY15_12h:CY15_24h" = 20,   #[15]
                 "CY15_3h:CY15_12h:CY15_48h" = 21,   #[16]
                 "CY15_3h:CY15_24h:CY15_48h" = 22,   #[17]
                 "CY15_12h:CY15_24h:CY15_48h" = 23,   #[18]
                 "CY15_1h:CY15_3h:CY15_12h:CY15_24h" = 24,   #[19]
                 "CY15_1h:CY15_3h:CY15_12h:CY15_48h" = 25,   #[20]
                 "CY15_1h:CY15_12h:CY15_24h:CY15_48h" = 26,   #[21]
                 "CY15_3h:CY15_12h:CY15_24h:CY15_48h" = 27,   #[22]
                 "CY15_1h:CY15_3h:CY15_12h:CY15_24h:CY15_48h" = 28,   #[23]
                 "CY15_1h" = 1,   #[24]
                 "CY15_3h" = 2,   #[25]
                 "CY15_12h" = 3,   #[26]
                 "CY15_24h" = 4,   #[27]
                 "CY15_48h" = 5   #[27]
  )
  T_attribute <- c(T_attribute, rep(temp, times=length(CY15_Venn[names(CY15_Venn) == unique(CY15_all)[n]])))
  print(n)
  n <- n+1
}

CY15_FDR005_time <- data.frame(AGI=T_AGI,
                               CY15_FDR005=T_attribute)
#CY16#
CY16_data <- list(CY16_1h = rownames(CY16_1h_FDR005), CY16_3h = rownames(CY16_3h_FDR005), CY16_12h = rownames(CY16_12h_FDR005), CY16_24h = rownames(CY16_24h_FDR005), CY16_48h = rownames(CY16_48h_FDR005))
CY16_Venn <- venn(CY16_data)
CY16_Venn <- unlist(attr(CY16_Venn,"intersections"))
temp_all <- strsplit(names(CY16_Venn), "h")#文字列の分割#overlapしている数が複数だとその数字が表示されてしまうので除きたい。例CY16_1h53　この53を除くために以下の作業
CY16_all <- c()
n <- 1
for(n in 1:length(temp_all)){
  temp <- temp_all[[n]]
  temp_nChar <- nchar(temp)#nchar文字列の長さを調べる
  CY16_all <- c(CY16_all, paste0(paste(temp[temp_nChar >= 6], collapse="h"), "h"))#>=4の理由はCY16_1hの文字列がヒットするのとoverlapしている数がmax3桁であるから。
}
names(CY16_Venn) <- CY16_all
#attribute_file#
n <- 1
T_AGI <- c()
T_attribute <- c()
for(n in 1:length(unique(CY16_all))){
  T_AGI <- c(T_AGI, CY16_Venn[names(CY16_Venn) == unique(CY16_all)[n]])
  temp <- switch(unique(CY16_all)[n],
                 "CY16_1h:CY16_3h" = 6,          #[1]
                 "CY16_1h:CY16_12h" = 7,          #[2]
                 "CY16_1h:CY16_24h" = 8,          #[3]
                 "CY16_1h:CY16_48h" = 9,          #[4]
                 "CY16_3h:CY16_12h" = 10,          #[5]
                 "CY16_3h:CY16_24h" = 11,          #[6]
                 "CY16_3h:CY16_48h" = 12,          #[7]
                 "CY16_12h:CY16_24h" = 13,          #[8]
                 "CY16_12h:CY16_48h" = 14,          #[9]
                 "CY16_24h:CY16_48h" = 15,          #[10]
                 "CY16_1h:CY16_3h:CY16_12h" = 16,          #[11]
                 "CY16_1h:CY16_3h:CY16_24h" = 17,          #[12]
                 "CY16_1h:CY16_3h:CY16_48h" = 18,          #[13]
                 "CY16_1h:CY16_12h:CY16_24h" = 19,          #[14]
                 "CY16_1h:CY16_24h:CY16_48h" = 20,          #[15]
                 "CY16_3h:CY16_12h:CY16_24h" = 21,          #[16]
                 "CY16_3h:CY16_24h:CY16_48h" = 22,          #[17]
                 "CY16_12h:CY16_24h:CY16_48h" = 23,          #[18]
                 "CY16_1h:CY16_3h:CY16_12h:CY16_24h" = 24,          #[19]
                 "CY16_1h:CY16_12h:CY16_24h:CY16_48h" = 25,          #[20]
                 "CY16_3h:CY16_12h:CY16_24h:CY16_48h" = 26,          #[21]
                 "CY16_1h:CY16_3h:CY16_12h:CY16_24h:CY16_48h" = 27,          #[22]
                 "CY16_1h" = 1,          #[23]
                 "CY16_3h" = 2,          #[24]
                 "CY16_12h" = 3,          #[25]
                 "CY16_24h" = 4,          #[26] 
                 "CY16_48h" = 5          #[27]
)
  T_attribute <- c(T_attribute, rep(temp, times=length(CY16_Venn[names(CY16_Venn) == unique(CY16_all)[n]])))
  print(n)
  n <- n+1
}

CY16_FDR005_time <- data.frame(AGI=T_AGI,
                               CY16_FDR005=T_attribute)
#CY20#
CY20_data <- list(CY20_1h = rownames(CY20_1h_FDR005), CY20_3h = rownames(CY20_3h_FDR005), CY20_12h = rownames(CY20_12h_FDR005), CY20_24h = rownames(CY20_24h_FDR005), CY20_48h = rownames(CY20_48h_FDR005))
CY20_Venn <- venn(CY20_data)
CY20_Venn <- unlist(attr(CY20_Venn,"intersections"))
temp_all <- strsplit(names(CY20_Venn), "h")#文字列の分割#overlapしている数が複数だとその数字が表示されてしまうので除きたい。例CY20_1h53　この53を除くために以下の作業
CY20_all <- c()
n <- 1
for(n in 1:length(temp_all)){
  temp <- temp_all[[n]]
  temp_nChar <- nchar(temp)#nchar文字列の長さを調べる
  CY20_all <- c(CY20_all, paste0(paste(temp[temp_nChar >= 6], collapse="h"), "h"))#>=4の理由はCY20_1hの文字列がヒットするのとoverlapしている数がmax3桁であるから。
}
names(CY20_Venn) <- CY20_all
#attribute_file#
T_AGI <- c()
T_attribute <- c()
for(n in 1:length(unique(CY20_all))){
  T_AGI <- c(T_AGI, CY20_Venn[names(CY20_Venn) == unique(CY20_all)[n]])
  temp <- switch(unique(CY20_all)[n],
                 "CY20_1h:CY20_3h" = 6,          #[1]
                 "CY20_1h:CY20_12h" = 7,          #[2]
                 "CY20_1h:CY20_24h" = 8,          #[3]
                 "CY20_1h:CY20_48h" = 9,          #[4]
                 "CY20_3h:CY20_12h" = 10,          #[5]
                 "CY20_3h:CY20_24h" = 11,          #[6]
                 "CY20_3h:CY20_48h" = 12,          #[7]
                 "CY20_12h:CY20_24h" = 13,          #[8]
                 "CY20_12h:CY20_48h" = 14,          #[9]
                 "CY20_24h:CY20_48h" = 15,          #[10]
                 "CY20_1h:CY20_3h:CY20_12h" = 16,          #[11]
                 "CY20_1h:CY20_3h:CY20_24h" = 17,          #[12]
                 "CY20_1h:CY20_3h:CY20_48h" = 18,          #[13]
                 "CY20_1h:CY20_12h:CY20_24h" = 19,          #[14]
                 "CY20_1h:CY20_24h:CY20_48h" = 20,          #[15]
                 "CY20_3h:CY20_12h:CY20_24h" = 21,          #[16]
                 "CY20_3h:CY20_12h:CY20_48h" = 22,          #[17]
                 "CY20_3h:CY20_24h:CY20_48h" = 23,          #[18]
                 "CY20_12h:CY20_24h:CY20_48h" = 24,          #[19]
                 "CY20_1h:CY20_3h:CY20_12h:CY20_24h" = 25,          #[20]
                 "CY20_1h:CY20_3h:CY20_12h:CY20_48h" = 26,          #[21]
                 "CY20_1h:CY20_3h:CY20_24h:CY20_48h" = 27,          #[22]
                 "CY20_1h:CY20_12h:CY20_24h:CY20_48h" = 28,          #[23]
                 "CY20_3h:CY20_12h:CY20_24h:CY20_48h" = 29,          #[24]
                 "CY20_1h:CY20_3h:CY20_12h:CY20_24h:CY20_48h" = 30,          #[25]
                 "CY20_1h" = 1,          #[26]
                 "CY20_3h" = 2,          #[27]
                 "CY20_12h" = 3,          #[28]
                 "CY20_24h" = 4,          #[29]
                 "CY20_48h" = 5          #[30]
)
  T_attribute <- c(T_attribute, rep(temp, times=length(CY20_Venn[names(CY20_Venn) == unique(CY20_all)[n]])))
}

CY20_FDR005_time <- data.frame(AGI=T_AGI,
                               CY20_FDR005=T_attribute)




####timedata_FDR0005_union####
#CY15_time
CY15_1h_FDR0005 <- allRNASeq[allRNASeq$CY15_1h_q_value < 0.005, ]
CY15_3h_FDR0005 <- allRNASeq[allRNASeq$CY15_3h_q_value < 0.005, ]
CY15_12h_FDR0005 <- allRNASeq[allRNASeq$CY15_12h_q_value < 0.005, ]
CY15_24h_FDR0005 <- allRNASeq[allRNASeq$CY15_24h_q_value < 0.005, ]
CY15_48h_FDR0005 <- allRNASeq[allRNASeq$CY15_48h_q_value < 0.005, ]
#CY16_time
CY16_1h_FDR0005 <- allRNASeq[allRNASeq$CY16_1h_q_value < 0.005, ]
CY16_3h_FDR0005 <- allRNASeq[allRNASeq$CY16_3h_q_value < 0.005, ]
CY16_12h_FDR0005 <- allRNASeq[allRNASeq$CY16_12h_q_value < 0.005, ]
CY16_24h_FDR0005 <- allRNASeq[allRNASeq$CY16_24h_q_value < 0.005, ]
CY16_48h_FDR0005 <- allRNASeq[allRNASeq$CY16_48h_q_value < 0.005, ]
#CY20_time
CY20_1h_FDR0005 <- allRNASeq[allRNASeq$CY20_1h_q_value < 0.005, ]
CY20_3h_FDR0005 <- allRNASeq[allRNASeq$CY20_3h_q_value < 0.005, ]
CY20_12h_FDR0005 <- allRNASeq[allRNASeq$CY20_12h_q_value < 0.005, ]
CY20_24h_FDR0005 <- allRNASeq[allRNASeq$CY20_24h_q_value < 0.005, ]
CY20_48h_FDR0005 <- allRNASeq[allRNASeq$CY20_48h_q_value < 0.005, ]
#Venn Diagram#
#CY15
library(gplots)
CY15_data <- list(CY15_1h = rownames(CY15_1h_FDR0005), CY15_3h = rownames(CY15_3h_FDR0005), CY15_12h = rownames(CY15_12h_FDR0005), CY15_24h = rownames(CY15_24h_FDR0005), CY15_48h = rownames(CY15_48h_FDR0005))
CY15_Venn <- venn(CY15_data)
CY15_Venn <- unlist(attr(CY15_Venn,"intersections"))
temp_all <- strsplit(names(CY15_Venn), "h")  #文字列の分割#overlapしている数が複数だとその数字が表示されてしまうので除きたい。例CY15_1h53　この53を除くために以下の作業
CY15_all <- c()
n <- 1
for(n in 1:length(temp_all)){
  temp <- temp_all[[n]]
  temp_nChar <- nchar(temp)#nchar文字列の長さを調べる
  CY15_all <- c(CY15_all, paste0(paste(temp[temp_nChar >= 6], collapse="h"), "h"))#>=4の理由はCY15_1hの文字列がヒットするのとoverlapしている数がmax3桁であるから。
}
names(CY15_Venn) <- CY15_all
#attribute_file#
n <- 1
T_AGI <- c()
T_attribute <- c()
for(n in 1:length(unique(CY15_all))){
  T_AGI <- c(T_AGI, CY15_Venn[names(CY15_Venn) == unique(CY15_all)[n]])
  temp <- switch(unique(CY15_all)[n],
                 "CY15_1h:CY15_3h" = 6,           #[1]
                 "CY15_1h:CY15_12h" = 7,           #[2]
                 "CY15_1h:CY15_24h" = 8,           #[3]
                 "CY15_1h:CY15_48h" = 9,           #[4]
                 "CY15_3h:CY15_12h" = 10,           #[5]
                 "CY15_3h:CY15_24h" = 11,           #[6]
                 "CY15_3h:CY15_48h" = 12,           #[7]
                 "CY15_12h:CY15_24h" = 13,           #[8]
                 "CY15_12h:CY15_48h" = 14,           #[9]
                 "CY15_24h:CY15_48h" = 15,           #[10]
                 "CY15_1h:CY15_3h:CY15_12h" = 16,           #[11]
                 "CY15_1h:CY15_3h:CY15_24h" = 17,           #[12]
                 "CY15_1h:CY15_3h:CY15_48h" = 18,           #[13]
                 "CY15_1h:CY15_12h:CY15_24h" = 19,           #[14]
                 "CY15_3h:CY15_12h:CY15_24h" = 20,           #[15]
                 "CY15_3h:CY15_12h:CY15_48h" = 21,           #[16]
                 "CY15_3h:CY15_24h:CY15_48h" = 22,           #[17]
                 "CY15_12h:CY15_24h:CY15_48h" = 23,           #[18]
                 "CY15_1h:CY15_3h:CY15_12h:CY15_24h" = 24,           #[19]
                 "CY15_1h:CY15_3h:CY15_12h:CY15_48h" = 25,           #[20]
                 "CY15_1h:CY15_12h:CY15_24h:CY15_48h" = 26,           #[21]
                 "CY15_3h:CY15_12h:CY15_24h:CY15_48h" = 27,           #[22]
                 "CY15_1h:CY15_3h:CY15_12h:CY15_24h:CY15_48h" = 28,           #[23]
                 "CY15_1h" = 1,           #[24]
                 "CY15_3h" = 2,           #[25]
                 "CY15_12h" = 3,           #[26]
                 "CY15_24h" = 4,           #[27]
                 "CY15_48h" = 5           #[28]
  )
  T_attribute <- c(T_attribute, rep(temp, times=length(CY15_Venn[names(CY15_Venn) == unique(CY15_all)[n]])))
  print(n)
  n <- n+1
}
CY15_FDR0005_time <- data.frame(AGI=T_AGI,
                                CY15_FDR0005=T_attribute)
#CY16
CY16_data <- list(CY16_1h = rownames(CY16_1h_FDR0005), CY16_3h = rownames(CY16_3h_FDR0005), CY16_12h = rownames(CY16_12h_FDR0005), CY16_24h = rownames(CY16_24h_FDR0005), CY16_48h = rownames(CY16_48h_FDR0005))
CY16_Venn <- venn(CY16_data)
CY16_Venn <- unlist(attr(CY16_Venn,"intersections"))
temp_all <- strsplit(names(CY16_Venn), "h")#文字列の分割#overlapしている数が複数だとその数字が表示されてしまうので除きたい。例CY16_1h53　この53を除くために以下の作業
CY16_all <- c()
n <- 1
for(n in 1:length(temp_all)){
  temp <- temp_all[[n]]
  temp_nChar <- nchar(temp)#nchar文字列の長さを調べる
  CY16_all <- c(CY16_all, paste0(paste(temp[temp_nChar >= 6], collapse="h"), "h"))#>=4の理由はCY16_1hの文字列がヒットするのとoverlapしている数がmax3桁であるから。
}
names(CY16_Venn) <- CY16_all
#attribute_file#
n <- 1
T_AGI <- c()
T_attribute <- c()
for(n in 1:length(unique(CY16_all))){
  T_AGI <- c(T_AGI, CY16_Venn[names(CY16_Venn) == unique(CY16_all)[n]])
  temp <- switch(unique(CY16_all)[n],
                 "CY16_1h:CY16_3h" = 5,             #[1]
                 "CY16_1h:CY16_12h" = 6,             #[2]
                 "CY16_1h:CY16_24h" = 7,             #[3]
                 "CY16_3h:CY16_12h" = 8,             #[4]
                 "CY16_3h:CY16_24h" = 9,             #[5]
                 "CY16_12h:CY16_24h" = 10,             #[6]
                 "CY16_1h:CY16_3h:CY16_12h" = 11,             #[7]
                 "CY16_1h:CY16_3h:CY16_24h" = 12,             #[8]
                 "CY16_1h:CY16_12h:CY16_24h" = 13,             #[9]
                 "CY16_3h:CY16_12h:CY16_24h" = 14,             #[10]
                 "CY16_1h:CY16_3h:CY16_12h:CY16_24h" = 15,             #[11]
                 "CY16_1h" = 1,          #[12]
                 "CY16_3h" = 2,          #[13]
                 "CY16_12h" = 3,          #[14]
                 "CY16_24h" = 4          #[15]
)
  T_attribute <- c(T_attribute, rep(temp, times=length(CY16_Venn[names(CY16_Venn) == unique(CY16_all)[n]])))
  print(n)
  n <- n+1
}
CY16_FDR0005_time <- data.frame(AGI=T_AGI,
                                CY16_FDR0005=T_attribute)
#CY20
CY20_data <- list(CY20_1h = rownames(CY20_1h_FDR0005), CY20_3h = rownames(CY20_3h_FDR0005), CY20_12h = rownames(CY20_12h_FDR0005), CY20_24h = rownames(CY20_24h_FDR0005), CY20_48h = rownames(CY20_48h_FDR0005))
CY20_Venn <- venn(CY20_data)
CY20_Venn <- unlist(attr(CY20_Venn,"intersections"))
temp_all <- strsplit(names(CY20_Venn), "h")#文字列の分割#overlapしている数が複数だとその数字が表示されてしまうので除きたい。例CY20_1h53　この53を除くために以下の作業
CY20_all <- c()
n <- 1
for(n in 1:length(temp_all)){
  temp <- temp_all[[n]]
  temp_nChar <- nchar(temp)#nchar文字列の長さを調べる
  CY20_all <- c(CY20_all, paste0(paste(temp[temp_nChar >= 6], collapse="h"), "h"))#>=4の理由はCY20_1hの文字列がヒットするのとoverlapしている数がmax3桁であるから。
}
names(CY20_Venn) <- CY20_all
#attribute_file#
T_AGI <- c()
T_attribute <- c()
for(n in 1:length(unique(CY20_all))){
  T_AGI <- c(T_AGI, CY20_Venn[names(CY20_Venn) == unique(CY20_all)[n]])
  temp <- switch(unique(CY20_all)[n],
                 "CY20_24h:CY20_48h" = 3,
                 "CY20_24h" = 1,
                 "CY20_48h" = 2)
  T_attribute <- c(T_attribute, rep(temp, times=length(CY20_Venn[names(CY20_Venn) == unique(CY20_all)[n]])))
}
CY20_FDR0005_time <- data.frame(AGI=T_AGI,
                                CY20_FDR0005=T_attribute)
####up, down regulation####
####FDR005####
#CY15#
#1h
CY15_1h_FDR005_up <- CY15_1h_FDR005[CY15_1h_FDR005$CY15_1h > 0, ]
CY15_1h_FDR005_down <- CY15_1h_FDR005[CY15_1h_FDR005$CY15_1h < 0, ]
CY15_1h_FDR005_up_color <- rep("U", times = nrow(CY15_1h_FDR005_up))
CY15_1h_FDR005_down_color <- rep("D", times = nrow(CY15_1h_FDR005_down))
attrCY15_1h_FDR005 <- data.frame(CY15_1h_FDR005 = c(rownames(CY15_1h_FDR005_up), rownames(CY15_1h_FDR005_down)), 
                                 CY15_1h_FDR005_color = c(CY15_1h_FDR005_up_color, CY15_1h_FDR005_down_color)
                                 )
#3h
CY15_3h_FDR005_up <- CY15_3h_FDR005[CY15_3h_FDR005$CY15_3h > 0, ]
CY15_3h_FDR005_down <- CY15_3h_FDR005[CY15_3h_FDR005$CY15_3h < 0, ]
CY15_3h_FDR005_up_color <- rep("U", times = nrow(CY15_3h_FDR005_up))
CY15_3h_FDR005_down_color <- rep("D", times = nrow(CY15_3h_FDR005_down))
attrCY15_3h_FDR005 <- data.frame(CY15_3h_FDR005 = c(rownames(CY15_3h_FDR005_up), rownames(CY15_3h_FDR005_down)), 
                                 CY15_3h_FDR005_color = c(CY15_3h_FDR005_up_color, CY15_3h_FDR005_down_color)
)
#12h
CY15_12h_FDR005_up <- CY15_12h_FDR005[CY15_12h_FDR005$CY15_12h > 0, ]
CY15_12h_FDR005_down <- CY15_12h_FDR005[CY15_12h_FDR005$CY15_12h < 0, ]
CY15_12h_FDR005_up_color <- rep("U", times = nrow(CY15_12h_FDR005_up))
CY15_12h_FDR005_down_color <- rep("D", times = nrow(CY15_12h_FDR005_down))
attrCY15_12h_FDR005 <- data.frame(CY15_12h_FDR005 = c(rownames(CY15_12h_FDR005_up), rownames(CY15_12h_FDR005_down)), 
                                 CY15_12h_FDR005_color = c(CY15_12h_FDR005_up_color, CY15_12h_FDR005_down_color)
)
#24h
CY15_24h_FDR005_up <- CY15_24h_FDR005[CY15_24h_FDR005$CY15_24h > 0, ]
CY15_24h_FDR005_down <- CY15_24h_FDR005[CY15_24h_FDR005$CY15_24h < 0, ]
CY15_24h_FDR005_up_color <- rep("U", times = nrow(CY15_24h_FDR005_up))
CY15_24h_FDR005_down_color <- rep("D", times = nrow(CY15_24h_FDR005_down))
attrCY15_24h_FDR005 <- data.frame(CY15_24h_FDR005 = c(rownames(CY15_24h_FDR005_up), rownames(CY15_24h_FDR005_down)), 
                                 CY15_24h_FDR005_color = c(CY15_24h_FDR005_up_color, CY15_24h_FDR005_down_color)
)
#48h
CY15_48h_FDR005_up <- CY15_48h_FDR005[CY15_48h_FDR005$CY15_48h > 0, ]
CY15_48h_FDR005_down <- CY15_48h_FDR005[CY15_48h_FDR005$CY15_48h < 0, ]
CY15_48h_FDR005_up_color <- rep("U", times = nrow(CY15_48h_FDR005_up))
CY15_48h_FDR005_down_color <- rep("D", times = nrow(CY15_48h_FDR005_down))
attrCY15_48h_FDR005 <- data.frame(CY15_48h_FDR005 = c(rownames(CY15_48h_FDR005_up), rownames(CY15_48h_FDR005_down)), 
                                 CY15_48h_FDR005_color = c(CY15_48h_FDR005_up_color, CY15_48h_FDR005_down_color)
)

#CY16#
#1h
CY16_1h_FDR005_up <- CY16_1h_FDR005[CY16_1h_FDR005$CY16_1h > 0, ]
CY16_1h_FDR005_down <- CY16_1h_FDR005[CY16_1h_FDR005$CY16_1h < 0, ]
CY16_1h_FDR005_up_color <- rep("U", times = nrow(CY16_1h_FDR005_up))
CY16_1h_FDR005_down_color <- rep("D", times = nrow(CY16_1h_FDR005_down))
attrCY16_1h_FDR005 <- data.frame(CY16_1h_FDR005 = c(rownames(CY16_1h_FDR005_up), rownames(CY16_1h_FDR005_down)), 
                                 CY16_1h_FDR005_color = c(CY16_1h_FDR005_up_color, CY16_1h_FDR005_down_color)
)
#3h
CY16_3h_FDR005_up <- CY16_3h_FDR005[CY16_3h_FDR005$CY16_3h > 0, ]
CY16_3h_FDR005_down <- CY16_3h_FDR005[CY16_3h_FDR005$CY16_3h < 0, ]
CY16_3h_FDR005_up_color <- rep("U", times = nrow(CY16_3h_FDR005_up))
CY16_3h_FDR005_down_color <- rep("D", times = nrow(CY16_3h_FDR005_down))
attrCY16_3h_FDR005 <- data.frame(CY16_3h_FDR005 = c(rownames(CY16_3h_FDR005_up), rownames(CY16_3h_FDR005_down)), 
                                 CY16_3h_FDR005_color = c(CY16_3h_FDR005_up_color, CY16_3h_FDR005_down_color)
)
#12h
CY16_12h_FDR005_up <- CY16_12h_FDR005[CY16_12h_FDR005$CY16_12h > 0, ]
CY16_12h_FDR005_down <- CY16_12h_FDR005[CY16_12h_FDR005$CY16_12h < 0, ]
CY16_12h_FDR005_up_color <- rep("U", times = nrow(CY16_12h_FDR005_up))
CY16_12h_FDR005_down_color <- rep("D", times = nrow(CY16_12h_FDR005_down))
attrCY16_12h_FDR005 <- data.frame(CY16_12h_FDR005 = c(rownames(CY16_12h_FDR005_up), rownames(CY16_12h_FDR005_down)), 
                                  CY16_12h_FDR005_color = c(CY16_12h_FDR005_up_color, CY16_12h_FDR005_down_color)
)
#24h
CY16_24h_FDR005_up <- CY16_24h_FDR005[CY16_24h_FDR005$CY16_24h > 0, ]
CY16_24h_FDR005_down <- CY16_24h_FDR005[CY16_24h_FDR005$CY16_24h < 0, ]
CY16_24h_FDR005_up_color <- rep("U", times = nrow(CY16_24h_FDR005_up))
CY16_24h_FDR005_down_color <- rep("D", times = nrow(CY16_24h_FDR005_down))
attrCY16_24h_FDR005 <- data.frame(CY16_24h_FDR005 = c(rownames(CY16_24h_FDR005_up), rownames(CY16_24h_FDR005_down)), 
                                  CY16_24h_FDR005_color = c(CY16_24h_FDR005_up_color, CY16_24h_FDR005_down_color)
)
#48h
CY16_48h_FDR005_up <- CY16_48h_FDR005[CY16_48h_FDR005$CY16_48h > 0, ]
CY16_48h_FDR005_down <- CY16_48h_FDR005[CY16_48h_FDR005$CY16_48h < 0, ]
CY16_48h_FDR005_up_color <- rep("U", times = nrow(CY16_48h_FDR005_up))
CY16_48h_FDR005_down_color <- rep("D", times = nrow(CY16_48h_FDR005_down))
attrCY16_48h_FDR005 <- data.frame(CY16_48h_FDR005 = c(rownames(CY16_48h_FDR005_up), rownames(CY16_48h_FDR005_down)), 
                                  CY16_48h_FDR005_color = c(CY16_48h_FDR005_up_color, CY16_48h_FDR005_down_color)
)

#CY20#
#1h
CY20_1h_FDR005_up <- CY20_1h_FDR005[CY20_1h_FDR005$CY20_1h > 0, ]
CY20_1h_FDR005_down <- CY20_1h_FDR005[CY20_1h_FDR005$CY20_1h < 0, ]
CY20_1h_FDR005_up_color <- rep("U", times = nrow(CY20_1h_FDR005_up))
CY20_1h_FDR005_down_color <- rep("D", times = nrow(CY20_1h_FDR005_down))
attrCY20_1h_FDR005 <- data.frame(CY20_1h_FDR005 = c(rownames(CY20_1h_FDR005_up), rownames(CY20_1h_FDR005_down)), 
                                 CY20_1h_FDR005_color = c(CY20_1h_FDR005_up_color, CY20_1h_FDR005_down_color)
)
#3h
CY20_3h_FDR005_up <- CY20_3h_FDR005[CY20_3h_FDR005$CY20_3h > 0, ]
CY20_3h_FDR005_down <- CY20_3h_FDR005[CY20_3h_FDR005$CY20_3h < 0, ]
CY20_3h_FDR005_up_color <- rep("U", times = nrow(CY20_3h_FDR005_up))
CY20_3h_FDR005_down_color <- rep("D", times = nrow(CY20_3h_FDR005_down))
attrCY20_3h_FDR005 <- data.frame(CY20_3h_FDR005 = c(rownames(CY20_3h_FDR005_up), rownames(CY20_3h_FDR005_down)), 
                                 CY20_3h_FDR005_color = c(CY20_3h_FDR005_up_color, CY20_3h_FDR005_down_color)
)
#12h
CY20_12h_FDR005_up <- CY20_12h_FDR005[CY20_12h_FDR005$CY20_12h > 0, ]
CY20_12h_FDR005_down <- CY20_12h_FDR005[CY20_12h_FDR005$CY20_12h < 0, ]
CY20_12h_FDR005_up_color <- rep("U", times = nrow(CY20_12h_FDR005_up))
CY20_12h_FDR005_down_color <- rep("D", times = nrow(CY20_12h_FDR005_down))
attrCY20_12h_FDR005 <- data.frame(CY20_12h_FDR005 = c(rownames(CY20_12h_FDR005_up), rownames(CY20_12h_FDR005_down)), 
                                  CY20_12h_FDR005_color = c(CY20_12h_FDR005_up_color, CY20_12h_FDR005_down_color)
)
#24h
CY20_24h_FDR005_up <- CY20_24h_FDR005[CY20_24h_FDR005$CY20_24h > 0, ]
CY20_24h_FDR005_down <- CY20_24h_FDR005[CY20_24h_FDR005$CY20_24h < 0, ]
CY20_24h_FDR005_up_color <- rep("U", times = nrow(CY20_24h_FDR005_up))
CY20_24h_FDR005_down_color <- rep("D", times = nrow(CY20_24h_FDR005_down))
attrCY20_24h_FDR005 <- data.frame(CY20_24h_FDR005 = c(rownames(CY20_24h_FDR005_up), rownames(CY20_24h_FDR005_down)), 
                                  CY20_24h_FDR005_color = c(CY20_24h_FDR005_up_color, CY20_24h_FDR005_down_color)
)
#48h
CY20_48h_FDR005_up <- CY20_48h_FDR005[CY20_48h_FDR005$CY20_48h > 0, ]
CY20_48h_FDR005_down <- CY20_48h_FDR005[CY20_48h_FDR005$CY20_48h < 0, ]
CY20_48h_FDR005_up_color <- rep("U", times = nrow(CY20_48h_FDR005_up))
CY20_48h_FDR005_down_color <- rep("D", times = nrow(CY20_48h_FDR005_down))
attrCY20_48h_FDR005 <- data.frame(CY20_48h_FDR005 = c(rownames(CY20_48h_FDR005_up), rownames(CY20_48h_FDR005_down)), 
                                  CY20_48h_FDR005_color = c(CY20_48h_FDR005_up_color, CY20_48h_FDR005_down_color)
)


####FDR0005####
#CY15#
#1h
CY15_1h_FDR0005_up <- CY15_1h_FDR0005[CY15_1h_FDR0005$CY15_1h > 0, ]
CY15_1h_FDR0005_down <- CY15_1h_FDR0005[CY15_1h_FDR0005$CY15_1h < 0, ]
CY15_1h_FDR0005_up_color <- rep("U", times = nrow(CY15_1h_FDR0005_up))
CY15_1h_FDR0005_down_color <- rep("D", times = nrow(CY15_1h_FDR0005_down))
attrCY15_1h_FDR0005 <- data.frame(CY15_1h_FDR0005 = c(rownames(CY15_1h_FDR0005_up), rownames(CY15_1h_FDR0005_down)), 
                                 CY15_1h_FDR0005_color = c(CY15_1h_FDR0005_up_color, CY15_1h_FDR0005_down_color)
)
#3h
CY15_3h_FDR0005_up <- CY15_3h_FDR0005[CY15_3h_FDR0005$CY15_3h > 0, ]
CY15_3h_FDR0005_down <- CY15_3h_FDR0005[CY15_3h_FDR0005$CY15_3h < 0, ]
CY15_3h_FDR0005_up_color <- rep("U", times = nrow(CY15_3h_FDR0005_up))
CY15_3h_FDR0005_down_color <- rep("D", times = nrow(CY15_3h_FDR0005_down))
attrCY15_3h_FDR0005 <- data.frame(CY15_3h_FDR0005 = c(rownames(CY15_3h_FDR0005_up), rownames(CY15_3h_FDR0005_down)), 
                                 CY15_3h_FDR0005_color = c(CY15_3h_FDR0005_up_color, CY15_3h_FDR0005_down_color)
)
#12h
CY15_12h_FDR0005_up <- CY15_12h_FDR0005[CY15_12h_FDR0005$CY15_12h > 0, ]
CY15_12h_FDR0005_down <- CY15_12h_FDR0005[CY15_12h_FDR0005$CY15_12h < 0, ]
CY15_12h_FDR0005_up_color <- rep("U", times = nrow(CY15_12h_FDR0005_up))
CY15_12h_FDR0005_down_color <- rep("D", times = nrow(CY15_12h_FDR0005_down))
attrCY15_12h_FDR0005 <- data.frame(CY15_12h_FDR0005 = c(rownames(CY15_12h_FDR0005_up), rownames(CY15_12h_FDR0005_down)), 
                                  CY15_12h_FDR0005_color = c(CY15_12h_FDR0005_up_color, CY15_12h_FDR0005_down_color)
)
#24h
CY15_24h_FDR0005_up <- CY15_24h_FDR0005[CY15_24h_FDR0005$CY15_24h > 0, ]
CY15_24h_FDR0005_down <- CY15_24h_FDR0005[CY15_24h_FDR0005$CY15_24h < 0, ]
CY15_24h_FDR0005_up_color <- rep("U", times = nrow(CY15_24h_FDR0005_up))
CY15_24h_FDR0005_down_color <- rep("D", times = nrow(CY15_24h_FDR0005_down))
attrCY15_24h_FDR0005 <- data.frame(CY15_24h_FDR0005 = c(rownames(CY15_24h_FDR0005_up), rownames(CY15_24h_FDR0005_down)), 
                                  CY15_24h_FDR005_color = c(CY15_24h_FDR0005_up_color, CY15_24h_FDR0005_down_color)
)
#48h
CY15_48h_FDR0005_up <- CY15_48h_FDR0005[CY15_48h_FDR0005$CY15_48h > 0, ]
CY15_48h_FDR0005_down <- CY15_48h_FDR0005[CY15_48h_FDR0005$CY15_48h < 0, ]
CY15_48h_FDR0005_up_color <- rep("U", times = nrow(CY15_48h_FDR0005_up))
CY15_48h_FDR0005_down_color <- rep("D", times = nrow(CY15_48h_FDR0005_down))
attrCY15_48h_FDR0005 <- data.frame(CY15_48h_FDR0005 = c(rownames(CY15_48h_FDR0005_up), rownames(CY15_48h_FDR0005_down)), 
                                  CY15_48h_FDR0005_color = c(CY15_48h_FDR0005_up_color, CY15_48h_FDR0005_down_color)
)

#CY16#
#1h
CY16_1h_FDR0005_up <- CY16_1h_FDR0005[CY16_1h_FDR0005$CY16_1h > 0, ]
CY16_1h_FDR0005_down <- CY16_1h_FDR0005[CY16_1h_FDR0005$CY16_1h < 0, ]
CY16_1h_FDR0005_up_color <- rep("U", times = nrow(CY16_1h_FDR0005_up))
CY16_1h_FDR0005_down_color <- rep("D", times = nrow(CY16_1h_FDR0005_down))
attrCY16_1h_FDR0005 <- data.frame(CY16_1h_FDR0005 = c(rownames(CY16_1h_FDR0005_up), rownames(CY16_1h_FDR0005_down)), 
                                 CY16_1h_FDR0005_color = c(CY16_1h_FDR0005_up_color, CY16_1h_FDR0005_down_color)
)
#3h
CY16_3h_FDR0005_up <- CY16_3h_FDR0005[CY16_3h_FDR0005$CY16_3h > 0, ]
CY16_3h_FDR0005_down <- CY16_3h_FDR0005[CY16_3h_FDR0005$CY16_3h < 0, ]
CY16_3h_FDR0005_up_color <- rep("U", times = nrow(CY16_3h_FDR0005_up))
CY16_3h_FDR0005_down_color <- rep("D", times = nrow(CY16_3h_FDR0005_down))
attrCY16_3h_FDR0005 <- data.frame(CY16_3h_FDR0005 = c(rownames(CY16_3h_FDR0005_up), rownames(CY16_3h_FDR0005_down)), 
                                 CY16_3h_FDR0005_color = c(CY16_3h_FDR0005_up_color, CY16_3h_FDR0005_down_color)
)
#12h
CY16_12h_FDR0005_up <- CY16_12h_FDR0005[CY16_12h_FDR0005$CY16_12h > 0, ]
CY16_12h_FDR0005_down <- CY16_12h_FDR0005[CY16_12h_FDR0005$CY16_12h < 0, ]
CY16_12h_FDR0005_up_color <- rep("U", times = nrow(CY16_12h_FDR0005_up))
CY16_12h_FDR0005_down_color <- rep("D", times = nrow(CY16_12h_FDR0005_down))
attrCY16_12h_FDR0005 <- data.frame(CY16_12h_FDR0005 = c(rownames(CY16_12h_FDR0005_up), rownames(CY16_12h_FDR0005_down)), 
                                  CY16_12h_FDR0005_color = c(CY16_12h_FDR0005_up_color, CY16_12h_FDR0005_down_color)
)
#24h
CY16_24h_FDR0005_up <- CY16_24h_FDR0005[CY16_24h_FDR0005$CY16_24h > 0, ]
CY16_24h_FDR0005_down <- CY16_24h_FDR0005[CY16_24h_FDR0005$CY16_24h < 0, ]
CY16_24h_FDR0005_up_color <- rep("U", times = nrow(CY16_24h_FDR0005_up))
CY16_24h_FDR0005_down_color <- rep("D", times = nrow(CY16_24h_FDR0005_down))
attrCY16_24h_FDR0005 <- data.frame(CY16_24h_FDR0005 = c(rownames(CY16_24h_FDR0005_up), rownames(CY16_24h_FDR0005_down)), 
                                  CY16_24h_FDR0005_color = c(CY16_24h_FDR0005_up_color, CY16_24h_FDR0005_down_color)
)
#48h
CY16_48h_FDR0005_up <- CY16_48h_FDR0005[CY16_48h_FDR0005$CY16_48h > 0, ]
CY16_48h_FDR0005_down <- CY16_48h_FDR0005[CY16_48h_FDR0005$CY16_48h < 0, ]
CY16_48h_FDR0005_up_color <- rep("U", times = nrow(CY16_48h_FDR0005_up))
CY16_48h_FDR0005_down_color <- rep("D", times = nrow(CY16_48h_FDR0005_down))
attrCY16_48h_FDR0005 <- data.frame(CY16_48h_FDR0005 = c(rownames(CY16_48h_FDR0005_up), rownames(CY16_48h_FDR0005_down)), 
                                  CY16_48h_FDR0005_color = c(CY16_48h_FDR0005_up_color, CY16_48h_FDR0005_down_color)
)

#CY20#
#1h
CY20_1h_FDR0005_up <- CY20_1h_FDR0005[CY20_1h_FDR0005$CY20_1h > 0, ]
CY20_1h_FDR0005_down <- CY20_1h_FDR0005[CY20_1h_FDR0005$CY20_1h < 0, ]
CY20_1h_FDR0005_up_color <- rep("U", times = nrow(CY20_1h_FDR0005_up))
CY20_1h_FDR0005_down_color <- rep("D", times = nrow(CY20_1h_FDR0005_down))
attrCY20_1h_FDR0005 <- data.frame(CY20_1h_FDR0005 = c(rownames(CY20_1h_FDR0005_up), rownames(CY20_1h_FDR0005_down)), 
                                 CY20_1h_FDR0005_color = c(CY20_1h_FDR0005_up_color, CY20_1h_FDR0005_down_color)
)
#3h
CY20_3h_FDR0005_up <- CY20_3h_FDR0005[CY20_3h_FDR0005$CY20_3h > 0, ]
CY20_3h_FDR0005_down <- CY20_3h_FDR0005[CY20_3h_FDR0005$CY20_3h < 0, ]
CY20_3h_FDR0005_up_color <- rep("U", times = nrow(CY20_3h_FDR0005_up))
CY20_3h_FDR0005_down_color <- rep("D", times = nrow(CY20_3h_FDR0005_down))
attrCY20_3h_FDR0005 <- data.frame(CY20_3h_FDR0005 = c(rownames(CY20_3h_FDR0005_up), rownames(CY20_3h_FDR0005_down)), 
                                 CY20_3h_FDR0005_color = c(CY20_3h_FDR0005_up_color, CY20_3h_FDR0005_down_color)
)
#12h
CY20_12h_FDR0005_up <- CY20_12h_FDR0005[CY20_12h_FDR0005$CY20_12h > 0, ]
CY20_12h_FDR0005_down <- CY20_12h_FDR0005[CY20_12h_FDR0005$CY20_12h < 0, ]
CY20_12h_FDR0005_up_color <- rep("U", times = nrow(CY20_12h_FDR0005_up))
CY20_12h_FDR0005_down_color <- rep("D", times = nrow(CY20_12h_FDR0005_down))
attrCY20_12h_FDR0005 <- data.frame(CY20_12h_FDR0005 = c(rownames(CY20_12h_FDR0005_up), rownames(CY20_12h_FDR0005_down)), 
                                  CY20_12h_FDR0005_color = c(CY20_12h_FDR0005_up_color, CY20_12h_FDR0005_down_color)
)
#24h
CY20_24h_FDR0005_up <- CY20_24h_FDR0005[CY20_24h_FDR0005$CY20_24h > 0, ]
CY20_24h_FDR0005_down <- CY20_24h_FDR0005[CY20_24h_FDR0005$CY20_24h < 0, ]
CY20_24h_FDR0005_up_color <- rep("U", times = nrow(CY20_24h_FDR0005_up))
CY20_24h_FDR0005_down_color <- rep("D", times = nrow(CY20_24h_FDR0005_down))
attrCY20_24h_FDR0005 <- data.frame(CY20_24h_FDR0005 = c(rownames(CY20_24h_FDR0005_up), rownames(CY20_24h_FDR0005_down)), 
                                  CY20_24h_FDR0005_color = c(CY20_24h_FDR0005_up_color, CY20_24h_FDR0005_down_color)
)
#48h
CY20_48h_FDR0005_up <- CY20_48h_FDR0005[CY20_48h_FDR0005$CY20_48h > 0, ]
CY20_48h_FDR0005_down <- CY20_48h_FDR0005[CY20_48h_FDR0005$CY20_48h < 0, ]
CY20_48h_FDR0005_up_color <- rep("U", times = nrow(CY20_48h_FDR0005_up))
CY20_48h_FDR0005_down_color <- rep("D", times = nrow(CY20_48h_FDR0005_down))
attrCY20_48h_FDR0005 <- data.frame(CY20_48h_FDR0005 = c(rownames(CY20_48h_FDR0005_up), rownames(CY20_48h_FDR0005_down)), 
                                  CY20_48h_FDR0005_color = c(CY20_48h_FDR0005_up_color, CY20_48h_FDR0005_down_color)
)


####output####
write.table(allRNASeq_cytoscape_th, file = "~/Nakano_RNAseq/network_analysis/cytoscape/allRNASeq_cytoscape.txt", append=F, quote = F, sep = "\t", row.names = F)
#attribute
#CY all
write.table(CY_FDR005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/union/FDR005/CYFDR005.txt", append=F, quote = F, sep = "\t", row.names = T)
write.table(CY_FDR0005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/union/FDR0005/CYFDR0005.txt", append=F, quote = F, sep = "\t", row.names = T)
#union_FDR005, FDR0005
write.table(CY15_FDR005_time, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/union/FDR005/CY15_FDR005_time.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(CY16_FDR005_time, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/union/FDR005/CY16_FDR005_time.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(CY20_FDR005_time, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/union/FDR005/CY20_FDR005_time.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(CY15_FDR0005_time, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/union/FDR0005/CY15_FDR0005_time.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(CY16_FDR0005_time, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/union/FDR0005/CY16_FDR0005_time.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(CY20_FDR0005_time, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/union/FDR0005/CY20_FDR0005_time.txt", append=F, quote = F, sep = "\t", row.names = F)
#up, down regulation_FDR005
write.table(attrCY15_1h_FDR005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR005/attrCY15_1h_FDR005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY15_3h_FDR005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR005/attrCY15_3h_FDR005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY15_12h_FDR005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR005/attrCY15_12h_FDR005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY15_24h_FDR005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR005/attrCY15_24h_FDR005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY15_48h_FDR005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR005/attrCY15_48h_FDR005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY16_1h_FDR005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR005/attrCY16_1h_FDR005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY16_3h_FDR005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR005/attrCY16_3h_FDR005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY16_12h_FDR005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR005/attrCY16_12h_FDR005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY16_24h_FDR005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR005/attrCY16_24h_FDR005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY16_48h_FDR005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR005/attrCY16_48h_FDR005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY20_1h_FDR005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR005/attrCY20_1h_FDR005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY20_3h_FDR005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR005/attrCY20_3h_FDR005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY20_12h_FDR005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR005/attrCY20_12h_FDR005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY20_24h_FDR005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR005/attrCY20_24h_FDR005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY20_48h_FDR005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR005/attrCY20_48h_FDR005.txt", append=F, quote = F, sep = "\t", row.names = F)
#up, down regulation_FDR0005
write.table(attrCY15_1h_FDR0005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY15_1h_FDR0005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY15_3h_FDR0005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY15_3h_FDR0005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY15_12h_FDR0005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY15_12h_FDR0005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY15_24h_FDR0005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY15_24h_FDR0005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY15_48h_FDR0005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY15_48h_FDR0005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY16_1h_FDR0005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY16_1h_FDR0005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY16_3h_FDR0005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY16_3h_FDR0005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY16_12h_FDR0005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY16_12h_FDR0005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY16_24h_FDR0005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY16_24h_FDR0005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY16_48h_FDR0005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY16_48h_FDR0005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY20_1h_FDR0005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY20_1h_FDR0005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY20_3h_FDR0005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY20_3h_FDR0005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY20_12h_FDR0005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY20_12h_FDR0005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY20_24h_FDR0005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY20_24h_FDR0005.txt", append=F, quote = F, sep = "\t", row.names = F)
write.table(attrCY20_48h_FDR0005, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/up_down/FDR0005/attrCY20_48h_FDR0005.txt", append=F, quote = F, sep = "\t", row.names = F)
#arabidopsis_TF
write.table(attrarabidopsis_TF, "~/Nakano_RNAseq/network_analysis/cytoscape/attribute/attrarabidopsis_TF.txt", append=F, quote = F, sep = "\t", row.names = F)

save.image("~/Nakano_RNAseq/network_analysis/20171030.RData")
