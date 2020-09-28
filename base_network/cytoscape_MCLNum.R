object <- read.table("~/Nakano_RNAseq/PCC/new/MCL_nodes.csv" ,header=T ,sep="," , row.names=1)
object <- na.omit(object)
allRNASeq_cytoscape_FDR005_MCL <- data.frame(allRNASeq_cytoscape_FDR005,
                                             MCLNum = object[match(allRNASeq_cytoscape_FDR005$source_interaction, object$name), 1])