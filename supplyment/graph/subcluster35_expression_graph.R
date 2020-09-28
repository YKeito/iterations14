library(ggsignif)
library(ggplot2)
library(reshape2)
subcluster35_AGI <- iterations14[iterations14$X__mclCluster == 35, ]$name
data <- allRNASeq[match(subcluster35_AGI, rownames(allRNASeq)), ]
col <- c("#F8766D", "a", "b", "c", "d", "#7CAE00", "e", "f", "g", "h", "#00BFC4")
total <- length(subcluster35_AGI)
i <- 1
for(i in i:total){
  df <- data.frame(time = rep(c("011h", "023h", "0312h", "0424h", "0548h"), time = 3), 
                   expression_change = c(data$CY15_1h[i], data$CY15_3h[i] ,data$CY15_12h[i], data$CY15_24h[i], data$CY15_48h[i],
                                         data$CY16_1h[i], data$CY16_3h[i] ,data$CY16_12h[i], data$CY16_24h[i], data$CY16_48h[i],
                                         data$CY20_1h[i], data$CY20_3h[i] ,data$CY20_12h[i], data$CY20_24h[i], data$CY20_48h[i]),
                   Category = rep(c("CY15", "CY16", "CY20"), each = 5)
  )
  levels(df$time) <- c("1h", "3h", "12h", "24h", "48h")
  
  g <- ggplot(
    df,
    aes (
      x = time,             # 遺伝子別でグルーピング
      y = expression_change,
      group = Category
    )
  )
  g <- g + geom_line(aes(colour=Category))
  g <- g + geom_point(aes(colour=Category))
  title <- paste0(rownames(data)[i], "_", data$genesymbol[i])
  g <- g + ggtitle(title)
  g <- g + xlab("")
  g <- g + ylab("log(treatment/control)")
  xx <- c()
  m <- 1
  for(m in c(19:33)){
    FDR <- data[i, m]
    if(FDR >= 0.05){
      x <- " "
    } else if(FDR < 0.005){
      x <- "**"
    } else{
      x <- "*"
    }
    xx <- c(xx, x)
    m <- m+1
  }
  n <- 1
  for(n in c(1, 6, 11)){
    g <- g + geom_signif(annotations = xx[n:c(n+4)],
                         size = 0,
                         y_position = c(df$expression_change[n] + 0.01, df$expression_change[n+1] + 0.01, df$expression_change[n+2] + 0.01, df$expression_change[n+3] + 0.01, df$expression_change[n+4] + 0.01), 
                         xmin=c(1, 2, 3, 4, 5),
                         xmax = c(1, 2, 3, 4, 5),
                         tip_length = 0,
                         col = col[n]
    )
    n <- n+1
  }
  plot(g)
  Fig_title <- paste0("~/Nakano_RNAseq/network_analysis/results_Fig/subcluster35/", title, ".png")
  ggsave(file = Fig_title, plot = g, dpi = 100)
  i <- i+1
}
