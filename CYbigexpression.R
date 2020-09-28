#calculate similarity
a <- 1
i <- 1
m <- 2
itotal <- 4
mtotal <- 5
CY15_all <- c()
CY15 <- c()
CY16_all <- c()
CY16 <- c()
CY20_all <- c()
CY20 <- c()
CY15list <- list()
CY16list <- list()
CY20list <- list()
for(a in a:nrow(allRNASeq)){
  for(i in i:itotal){
    for(m in m:mtotal){
      CY15 <- c(CY15, abs(allRNASeq[a, i] - allRNASeq[a, m]))
      CY16 <- c(CY16, abs(allRNASeq[a, i+5] - allRNASeq[a, m+5]))
      CY20 <- c(CY20, abs(allRNASeq[a, i+10] - allRNASeq[a, m+10]))
      m <- m+1
    }
    CY15_all <- c(CY15_all, CY15)
    CY16_all <- c(CY16_all, CY16)
    CY20_all <- c(CY20_all, CY20)
    CY15 <- c()
    CY16 <- c()
    CY20 <- c()
    i <- i+1
    m <- i+1
  }
  CY15list <- c(CY15list, list(CY15_all))
  CY15_all <- c()
  CY16list <- c(CY16list, list(CY16_all))
  CY16_all <- c()
  CY20list <- c(CY20list, list(CY20_all))
  CY20_all <- c()
  i <- 1
  m <- 2
  print(a)
  a <- a+1
}

NumCluster <- rep("No", times = nrow(allRNASeq)*10)
CY15_bigex <- data.frame(AGI = rep(rownames(allRNASeq), each = 10),
                         timesource = rep(combn(c("1h", "3h", "12h", "24h", "48h"), 2)[1, ], times=nrow(allRNASeq)),
                         bigexpression = unlist(CY15list),
                         timetarget = rep(combn(c("1h", "3h", "12h", "24h", "48h"), 2)[2, ], times=nrow(allRNASeq)),
                         NumMCL = NumCluster,
                         TF = rep("No", times = nrow(allRNASeq)*10),
                         hub = rep("No", times = nrow(allRNASeq)*10),
                         stringsAsFactors = F)
i <- 1
for(i in i:length(iterations14$name)){
  test <- CY15_bigex$AGI == iterations14$name[i]
  if(sum(test) >= 1){
    CY15_bigex[test, ]$NumMCL <- iterations14[i, ]$X__mclCluster
  }
  print(i)
  i <- i+1
}
i <- 1
for(i in i:length(TF_family$AGI)){
  test <- CY15_bigex$AGI == TF_family$AGI[i]
  if(sum(test) >= 1){
    CY15_bigex[test, ]$TF <- "YES"
  }
  print(i)
  i <- i+1
}
filename <- list.files("~/Nakano_RNAseq/network_analysis/base/Subcluster_Information/hublist_5%/", pattern=".txt", full.names = T)
n <- 1
for (n in 1:length(filename)){
  object_name <- paste0("subcluster", temp_a[n])
  object_name <- assign(object_name, read.table(filename[n], sep = "\t", header = T, stringsAsFactors = F))
  i <- 1
  for(i in i:length(object_name$name)){
    test <- CY15_bigex$AGI == object_name$name[i]
    if(sum(test) >= 1){
      CY15_bigex[test, ]$hub <- "YES"
    }
    i <- i+1
  }
  print(n)
}


CY16_bigex <- data.frame(AGI = rep(rownames(allRNASeq), each = 10),
                         timesource = rep(combn(c("1h", "3h", "12h", "24h", "48h"), 2)[1, ], times=nrow(allRNASeq)),
                         bigexpression = unlist(CY16list),
                         timetarget = rep(combn(c("1h", "3h", "12h", "24h", "48h"), 2)[2, ], times=nrow(allRNASeq)),
                         NumMCL = NumCluster,
                         TF = rep("No", times = nrow(allRNASeq)*10),
                         hub = rep("No", times = nrow(allRNASeq)*10),
                         stringsAsFactors = F)
i <- 1
for(i in i:length(iterations14$name)){
  test <- CY16_bigex$AGI == iterations14$name[i]
  if(sum(test) >= 1){
    CY16_bigex[test, ]$NumMCL <- iterations14[i, ]$X__mclCluster
  }
  print(i)
  i <- i+1
}
n <- 1
for (n in 1:length(filename)){
  object_name <- paste0("subcluster", temp_a[n])
  object_name <- assign(object_name, read.table(filename[n], sep = "\t", header = T, stringsAsFactors = F))
  i <- 1
  for(i in i:length(object_name$name)){
    test <- CY16_bigex$AGI == object_name$name[i]
    if(sum(test) >= 1){
      CY16_bigex[test, ]$hub <- "YES"
    }
    i <- i+1
  }
  print(n)
}
i <- 1
for(i in i:length(TF_family$AGI)){
  test <- CY16_bigex$AGI == TF_family$AGI[i]
  if(sum(test) >= 1){
    CY16_bigex[test, ]$TF <- "YES"
  }
  print(i)
  i <- i+1
}

CY20_bigex <- data.frame(AGI = rep(rownames(allRNASeq), each = 10),
                         timesource = rep(combn(c("1h", "3h", "12h", "24h", "48h"), 2)[1, ], times=nrow(allRNASeq)),
                         bigexpression = unlist(CY20list),
                         timetarget = rep(combn(c("1h", "3h", "12h", "24h", "48h"), 2)[2, ], times=nrow(allRNASeq)),
                         NumMCL = NumCluster,
                         TF = rep("No", times = nrow(allRNASeq)*10),
                         hub = rep("No", times = nrow(allRNASeq)*10),
                         stringsAsFactors = F)
i <- 1
for(i in i:length(iterations14$name)){
  test <- CY20_bigex$AGI == iterations14$name[i]
  if(sum(test) >= 1){
    CY20_bigex[test, ]$NumMCL <- iterations14[i, ]$X__mclCluster
  }
  print(i)
  i <- i+1
}
i <- 1
for(i in i:length(TF_family$AGI)){
  test <- CY20_bigex$AGI == TF_family$AGI[i]
  if(sum(test) >= 1){
    CY20_bigex[test, ]$TF <- "YES"
  }
  print(i)
  i <- i+1
}
n <- 1
for (n in 1:length(filename)){
  object_name <- paste0("subcluster", temp_a[n])
  object_name <- assign(object_name, read.table(filename[n], sep = "\t", header = T, stringsAsFactors = F))
  i <- 1
  for(i in i:length(object_name$name)){
    test <- CY20_bigex$AGI == object_name$name[i]
    if(sum(test) >= 1){
      CY20_bigex[test, ]$hub <- "YES"
    }
    i <- i+1
  }
  print(n)
}
