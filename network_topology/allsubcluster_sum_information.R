filename <- list.files("~/Nakano_RNAseq/network_analysis/base/Subcluster_Information/hublist_5%/", pattern=".txt", full.names = T)
TF_family <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/arabidopsis_TF_family.txt", sep = "\t", header = T, stringsAsFactors = F)
iterations14 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/iteration14_node.txt", sep = "\t", header = T, stringsAsFactors = F)
BTH <- read.table("~/Nakano_RNAseq/network_analysis/base/genes_set/BTH_affected_AGI.txt", sep = "\t", header = T, stringsAsFactors = F)
object_name <- c()
temp_a <- c(1, 2, 3, 4, 5, 6, 7, 9)
n <- 1
for (n in 1:length(filename)){
  object_name <- paste0("subcluster", temp_a[n])
  object_name <- assign(object_name, read.table(filename[n], sep = "\t", header = T, stringsAsFactors = F))
  m <- 1
  treatment <- c()
  AGI <- c()
  AGI_all <- list()
  temp <- list(as.character(unlist(CY15_FDR005)), as.character(unlist(CY16_FDR005)), as.character(unlist(CY20_FDR005)))
  temp2 <- c("CY15", "CY16", "CY20")
  for(m in m:length(temp)){
    AGI <- intersect(unlist(temp[m]), iterations14[iterations14$X__mclCluster == temp_a[n], ]$name)
    AGI_all <- c(AGI_all, list(AGI))
    treatment <- c(treatment, rep(temp2[m], times = length(AGI)))
    AGI <- c()
    m <- m+1
  }
  subcluster_list <- data.frame(treatment = treatment,
                                 AGI = unlist(AGI_all),
                                 TF = rep("NO", times = length(treatment)),
                                 hub = rep("NO", times = length(treatment)), 
                                 BTH = rep("NO", times = length(treatment)),
                                 stringsAsFactors = F
  )
  i <- 1
  for(i in i:length(TF_family$AGI)){
    test <- subcluster_list$AGI == TF_family$AGI[i]
    if(sum(test) >= 1){
      subcluster_list[test, ]$TF <- "YES"
    }
    i <- i+1
  }
  i <- 1
  for(i in i:length(object_name$name)){
    test <- subcluster_list$AGI == object_name$name[i]
    if(sum(test) >= 1){
      subcluster_list[test, ]$hub <- "YES"
    }
    i <- i+1
  }
  i <- 1
  for(i in i:length(BTH$AGI)){
    test <- subcluster_list$AGI == BTH$AGI[i]
    if(sum(test) >= 1){
      subcluster_list[test, ]$BTH <- "YES"
    }
    i <- i+1
  }
  #check ok #length(match(subcluster_hub$name, subcluster_list$AGI))
  File_title <- paste0("~/Nakano_RNAseq/network_analysis/base/Subcluster_Information/Sum_subcluster_information/", paste0("Sum_", "subcluster", temp_a[n]),".txt")
  write.table(subcluster_list, File_title, append=F, quote = F, sep = "\t", row.names = F)
  print(n)
  n <- n + 1
}