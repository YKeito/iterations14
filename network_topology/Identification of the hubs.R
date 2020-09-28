####mastercluster####
mastercluster_node <- read.table("~/Nakano_RNAseq/network_analysis/base/Subcluster_Information/mastercluster_node_table.csv", sep = ",", header = T)
mastercluster_node <- mastercluster_node[order(mastercluster_node$Degree, decreasing=T), ]
sum(mastercluster_node$Degree[1:148])/sum(mastercluster_node$Degree)#上位10%でハブ同定
unique(mastercluster_node[1:148, ]$X__mclCluster)#subcluster1のみであった。
mastercluster_hub <- mastercluster_node[1:148, ]$name#hub候補遺伝子のAGI出力
write.table(mastercluster_hub, "~/Nakano_RNAseq/network_analysis/base/Subcluster_Information/mastercluster_hub.txt", append=F, quote = F, sep = "\t", row.names = F)
####subcluster1####
subcluster1_node <- read.table("~/Nakano_RNAseq/network_analysis/base/Subcluster_Information/subcluster1nodetable.csv", sep = ",", header = T)
subcluster1_node <- subcluster1_node[order(subcluster1_node$Degree, decreasing=T), ]
sum(subcluster1_node$Degree[1:81])/sum(subcluster1_node$Degree)#上位10%でハブ同定
unique(subcluster1_node[1:81, ]$X__mclCluster)#subcluster1のみであった。
subcluster1_hub <- subcluster1_node[1:81, ]$name#hub候補遺伝子のAGI出力
write.table(subcluster1_hub, "~/Nakano_RNAseq/network_analysis/base/Subcluster_Information/subcluster1_hub.txt", append=F, quote = F, sep = "\t", row.names = F)

####subcluster2####
subcluster2_node <- read.table("~/Nakano_RNAseq/network_analysis/base/Subcluster_Information/subcluster2nodetable.csv", sep = ",", header = T)
subcluster2_node <- subcluster2_node[order(subcluster2_node$Degree, decreasing=T), ]
sum(subcluster2_node$Degree[1:27])/sum(subcluster2_node$Degree)#上位10%でハブ同定
unique(subcluster2_node[1:27, ]$X__mclCluster)#subcluster2のみであった。
subcluster2_hub <- subcluster2_node[1:27, ]$name#hub候補遺伝子のAGI出力
write.table(subcluster2_hub, "~/Nakano_RNAseq/network_analysis/base/Subcluster_Information/subcluster2_hub.txt", append=F, quote = F, sep = "\t", row.names = F)

####subcluster7####
subcluster7_node <- read.table("~/Nakano_RNAseq/network_analysis/base/Subcluster_Information/subcluster7nodetable.csv", sep = ",", header = T)
subcluster7_node <- subcluster7_node[order(subcluster7_node$Degree, decreasing=T), ]
sum(subcluster7_node$Degree[1])/sum(subcluster7_node$Degree)#上位10%でハブ同定
unique(subcluster7_node[1, ]$X__mclCluster)#subcluster7のみであった。
subcluster7_hub <- subcluster7_node[1, ]$name#hub候補遺伝子のAGI出力
write.table(subcluster7_hub, "~/Nakano_RNAseq/network_analysis/base/Subcluster_Information/subcluster7_hub.txt", append=F, quote = F, sep = "\t", row.names = F)

####subcluster35 hub候補なし####
subcluster35_node <- read.table("~/Nakano_RNAseq/network_analysis/base/Subcluster_Information/subcluster35nodetable.csv", sep = ",", header = T)
subcluster35_node <- subcluster35_node[order(subcluster35_node$Degree, decreasing=T), ]
sum(subcluster35_node$Degree[1])/sum(subcluster35_node$Degree)#上位10%でハブ同定
unique(subcluster35_node[1, ]$X__mclCluster)#subcluster35のみであった。
subcluster35_hub <- subcluster35_node[1, ]$name#hub候補遺伝子のAGI出力
