#A_領域のエッジ
A_ATTED_edges <- ATTED_temp[keito3$ATTED_PCC >= (-0.25) & keito3$ATTED_PCC <= (0.25), ]$AGI
A_CY_edges <- CY_temp[keito3$CY_PCC>= 0.75, ]$AGI
#標本のエッジ（A領域のエッジ）A_edges（変える必要ない）
A_edges <- intersect(A_ATTED_edges, A_CY_edges)
#母集団　全エッジ数
Numtotal_edges <- length(keito3$sample)
#全エッジ（母集団）
total_edges <- CY_temp$AGI


total <- c("red", "green", "blue", "black")
#Numsample_edges、目的のエッジ数（ORA59など）
#sample_edges、目的のエッジ(AGI)
i <- 1
pvalue <- c()
for(i in i:length(total)){
  NumTF_edges <- sum(keito3$sample == total[i])
  TF_edges <- CY_temp$AGI[keito3$sample == total[i]]
  TFtotal_edges <- intersect(TF_edges, total_edges)
  ATFtotal_edges <- intersect(A_edges, TFtotal_edges)
  pvalue <- c(pvalue, phyper(c(length(ATFtotal_edges)-1), length(TFtotal_edges), c(Numtotal_edges-length(TFtotal_edges)), length(ATFtotal_edges), lower.tail = F))
  i <- 1
}
