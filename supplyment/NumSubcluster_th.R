#MCL clusering node数の分布
iterations14 <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/iteration14_node.txt", sep = "\t", header = T)
rownames(iterations14) <- iterations14$name

maxiterations14 <- max(iterations14$X__mclCluster)
subcluster_NumNodes <- c()
i <- 1
for(i in i:maxiterations14){
  subcluster_NumNodes <- c(subcluster_NumNodes, sum(iterations14$X__mclCluster == i))
  i <- i+1
}

distribution <- c()
total <- length(unique(subcluster_NumNodes))
n <- 1
for(n in n:total){
  distribution <- c(distribution, sum(subcluster_NumNodes == unique(subcluster_NumNodes)[n]))
  n <- n+1
}
a <- data.frame(subcluster = unique(subcluster_NumNodes),
                distribution = distribution)
hist(a$distribution)

aa <- data.frame(Numsubcluster = c(1:442),
                distribution = rep(a$subcluster, times = a$distribution)
                )
