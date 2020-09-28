#Read hub of all subclusters
filename <- list.files("~/Nakano_RNAseq/network_analysis/base/Subcluster_Information/subcluster_networkanalyzer/", pattern="nodetable.csv", full.names = T)
object_name <- c()
temp <- c(12, 14, 17, 1, 23,  2, 35, 36, 3, 41, 42, 43, 47, 7)
n <- 1
for (n in 1:length(filename)) {
  object_name <- paste0("subcluster", temp[n], "_nodeinfo")
  object_name <- assign(object_name, read.table(filename[n],  header=T, sep=","))
  subcluster_Degree <- object_name[order(object_name$Degree, decreasing=T), ]
  subcluster_BetweennessCentrality <- object_name[order(object_name$BetweennessCentrality, decreasing=T), ]
  itotal <- nrow(subcluster_Degree)
  i <- 1
  for(i in i:itotal){
    temp2 <- sum(subcluster_Degree$Degree[1:i])/sum(subcluster_Degree$Degree)#上位10%でハブ同定
    temp3 <- sum(subcluster_BetweennessCentrality$BetweennessCentrality[1:i])/sum(subcluster_BetweennessCentrality$BetweennessCentrality)#上位10%でハブ同定
    if(temp2 < 0.1){
      hub <- subcluster_Degree[1:i, ] #hub候補遺伝子のAGI出力
    }
    if(temp3 < 0.1){
      bottleneck <- subcluster_BetweennessCentrality[1:i, ] #hub候補遺伝子のAGI出力
    }
    i <- i+1
  }
  File_title <- paste0("~/Nakano_RNAseq/network_analysis/base/Subcluster_Information/hub_list/", paste0("subcluster", temp[n]), "_hub",".txt")
  write.table(hub, File_title, append=F, quote = F, sep = "\t", row.names = F)
  File_title <- paste0("~/Nakano_RNAseq/network_analysis/base/Subcluster_Information/bottleneck_list/", paste0("subcluster", temp[n]), "_bottleneck",".txt")
  write.table(bottleneck, File_title, append=F, quote = F, sep = "\t", row.names = F)
  n <- n + 1
  hub <- c()
  bottleneck <- c()
}
