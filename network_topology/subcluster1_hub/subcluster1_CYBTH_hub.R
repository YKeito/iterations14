####pre-data####
subcluster1_node <- read.table("~/Nakano_RNAseq/network_analysis/base/Subcluster_Information/subcluster1nodetable.csv", sep = ",", header = T)
subcluster1_node <- subcluster1_node[order(subcluster1_node$Degree, decreasing=T), ]
subcluster1_hub <- subcluster1_node[1:81, ]$name#hub候補遺伝子のAGI出力

BTHCY151620_FDR005 <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/all_FDR005.txt", sep = "\t", header = T)
#BTHCY151620
BTHCY151620_common <- BTHCY151620_FDR005[BTHCY151620_FDR005$treatment == "CY15:CY16:CY20:BTH", ]$AGI
subcluster1_BTHCY151620_hub <- intersect(BTHCY151620_common, subcluster1_hub)
#CY151620
CY151620_common <- BTHCY151620_FDR005[BTHCY151620_FDR005$treatment == "CY15:CY16:CY20", ]$AGI
subcluster1_CY151620_hub <- intersect(CY151620_common, subcluster1_hub)
#CY1516BTH
CY1516BTH_common <- BTHCY151620_FDR005[BTHCY151620_FDR005$treatment == "CY15:CY16:BTH", ]$AGI
subcluster1_CY1516BTH_hub <- intersect(CY1516BTH_common, subcluster1_hub)
#CY1520BTH
CY1520BTH_common <- BTHCY151620_FDR005[BTHCY151620_FDR005$treatment == "CY15:CY20:BTH", ]$AGI
subcluster1_CY1520BTH_hub <- intersect(CY1520BTH_common, subcluster1_hub)
#CY15BTH
CY15BTH_common <- BTHCY151620_FDR005[BTHCY151620_FDR005$treatment == "CY15:BTH", ]$AGI
subcluster1_CY15BTH_hub <- intersect(CY15BTH_common, subcluster1_hub)
#CY1516
CY1516_common <- BTHCY151620_FDR005[BTHCY151620_FDR005$treatment == "CY15:CY16", ]$AGI
subcluster1_CY1516_hub <- intersect(CY1516_common, subcluster1_hub)
#CY1520
CY1520_common <- BTHCY151620_FDR005[BTHCY151620_FDR005$treatment == "CY15:CY20", ]$AGI
subcluster1_CY1520_hub <- intersect(CY1520_common, subcluster1_hub)
#CY15
CY15_common <- BTHCY151620_FDR005[BTHCY151620_FDR005$treatment == "CY15", ]$AGI
subcluster1_CY15_hub <- intersect(CY15_common, subcluster1_hub)
#CY1620BTH
CY1620BTH_common <- BTHCY151620_FDR005[BTHCY151620_FDR005$treatment == "CY16:CY20:BTH", ]$AGI
subcluster1_CY1620BTH_hub <- intersect(CY1620BTH_common, subcluster1_hub)
#CY16BTH
CY16BTH_common <- BTHCY151620_FDR005[BTHCY151620_FDR005$treatment == "CY16:BTH", ]$AGI
subcluster1_CY16BTH_hub <- intersect(CY16BTH_common, subcluster1_hub)
#CY1620
CY1620_common <- BTHCY151620_FDR005[BTHCY151620_FDR005$treatment == "CY16:CY20", ]$AGI
subcluster1_CY1620_hub <- intersect(CY1620_common, subcluster1_hub)
#CY16
CY16_common <- BTHCY151620_FDR005[BTHCY151620_FDR005$treatment == "CY16", ]$AGI
subcluster1_CY16_hub <- intersect(CY16_common, subcluster1_hub)
#CY20BTH
CY20BTH_common <- BTHCY151620_FDR005[BTHCY151620_FDR005$treatment == "CY20:BTH", ]$AGI
subcluster1_CY20BTH_hub <- intersect(CY20BTH_common, subcluster1_hub)
#CY20
CY20_common <- BTHCY151620_FDR005[BTHCY151620_FDR005$treatment == "CY20", ]$AGI
subcluster1_CY20_hub <- intersect(CY20_common, subcluster1_hub)
#######################################################################
BTHCY_hub <- c(subcluster1_BTHCY151620_hub, subcluster1_CY1516BTH_hub, subcluster1_CY15BTH_hub, subcluster1_CY16BTH_hub, subcluster1_CY20BTH_hub, subcluster1_CY1520BTH_hub, subcluster1_CY1620BTH_hub)
CY_hub <-  c(subcluster1_CY151620_hub, subcluster1_CY15_hub, subcluster1_CY1516_hub, subcluster1_CY1520_hub, subcluster1_CY16_hub, subcluster1_CY1620_hub, subcluster1_CY20_hub)

subcluster1_hub <- data.frame(BTHCY151620 = c(subcluster1_BTHCY151620_hub, rep(" ", times = 3)),
                              CY1516BTH = c(subcluster1_CY1516BTH_hub, rep(" ", times = 15)),
                              CY1520BTH = c(subcluster1_CY1520BTH_hub, rep(" ", times = 18)),
                              CY1620BTH = c(subcluster1_CY1620BTH_hub, rep(" ", times = 21)),
                              CY15BTH = c(subcluster1_CY15BTH_hub, rep(" ", times = 6)),
                              CY16BTH = c(subcluster1_CY16BTH_hub, rep(" ", times = 21)),
                              CY20BTH = c(subcluster1_CY20BTH_hub, rep(" ", times = 21)),
                              CY151620 = subcluster1_CY151620_hub,
                              CY1516 = c(subcluster1_CY1516_hub, rep(" ", times = 15)),
                              CY1520 = c(subcluster1_CY1520_hub, rep(" ", times = 21)),
                              CY1620 = c(subcluster1_CY1620_hub, rep(" ", times = 21)),
                              CY15 = c(subcluster1_CY15_hub, rep(" ", times = 9)),
                              CY16 = c(subcluster1_CY16_hub, rep(" ", times = 21)),
                              CY20 = c(subcluster1_CY20_hub, rep(" ", times = 21))
                              )

write.table(subcluster1_hub, file = "~/Nakano_RNAseq/network_analysis/base/genes_set/subcluster1_hub.txt", append=F, quote = F, sep = "\t", row.names = F)