View(NumMCL)
iteration14 <- na.omit(read.table("~/Nakano_RNAseq/network_analysis/cytoscape/iteration14_node.csv", sep = ",", header = T))

test <- c()
for(i in 1:442){
  test <- c(test, sum(iteration14$X__mclCluster == i))
}
test2 <- c()
for(i in 1:442){
  test2 <- c(test2, sum(NumMCL$X__mclCluster == i))
}
test
test2
sum(test != test2)


test <- c()
for(i in 1:442){
  test <- c(test, list(iteration14$name[iteration14$X__mclCluster == i]))
}
test2 <- c()
for(i in 1:442){
  test2 <- c(test2, list(NumMCL$name[NumMCL$X__mclCluster == i]))
}

test3 <- c()
for(i in 1:442){
  test3 <- c(test3, match(test[i], test2))
}