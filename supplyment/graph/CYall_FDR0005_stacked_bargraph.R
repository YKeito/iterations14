#read
common_network <- read.table("~/Nakano_RNAseq/network_analysis/base/Subcluster_Information/subcluster1commonnetworknodetable.csv", sep = ",", header = T)
####FDR0005####
#up
CY15_1h_FDR0005_up <- sum(common_network$CY15_1h_FDR0005_color == "U")
CY15_3h_FDR0005_up <- sum(common_network$CY15_3h_FDR0005_color == "U")
CY15_12h_FDR0005_up <- sum(common_network$CY15_12h_FDR0005_color == "U")
CY15_24h_FDR0005_up <- sum(common_network$CY15_24h_FDR0005_color == "U")
CY15_48h_FDR0005_up <- sum(common_network$CY15_48h_FDR0005_color == "U")
CY16_1h_FDR0005_up <- sum(common_network$CY16_1h_FDR0005_color == "U")
CY16_3h_FDR0005_up <- sum(common_network$CY16_3h_FDR0005_color == "U")
CY16_12h_FDR0005_up <- sum(common_network$CY16_12h_FDR0005_color == "U")
CY16_24h_FDR0005_up <- sum(common_network$CY16_24h_FDR0005_color == "U")
CY16_48h_FDR0005_up <- sum(common_network$CY16_48h_FDR0005_color == "U")
CY20_1h_FDR0005_up <- sum(common_network$CY20_1h_FDR0005_color == "U")
CY20_3h_FDR0005_up <- sum(common_network$CY20_3h_FDR0005_color == "U")
CY20_12h_FDR0005_up <- sum(common_network$CY20_12h_FDR0005_color == "U")
CY20_24h_FDR0005_up <- sum(common_network$CY20_24h_FDR0005_color == "U")
CY20_48h_FDR0005_up <- sum(common_network$CY20_48h_FDR0005_color == "U")
#down
CY15_1h_FDR0005_down <- sum(common_network$CY15_1h_FDR0005_color == "D")
CY15_3h_FDR0005_down <- sum(common_network$CY15_3h_FDR0005_color == "D")
CY15_12h_FDR0005_down <- sum(common_network$CY15_12h_FDR0005_color == "D")
CY15_24h_FDR0005_down <- sum(common_network$CY15_24h_FDR0005_color == "D")
CY15_48h_FDR0005_down <- sum(common_network$CY15_48h_FDR0005_color == "D")
CY16_1h_FDR0005_down <- sum(common_network$CY16_1h_FDR0005_color == "D")
CY16_3h_FDR0005_down <- sum(common_network$CY16_3h_FDR0005_color == "D")
CY16_12h_FDR0005_down <- sum(common_network$CY16_12h_FDR0005_color == "D")
CY16_24h_FDR0005_down <- sum(common_network$CY16_24h_FDR0005_color == "D")
CY16_48h_FDR0005_down <- sum(common_network$CY16_48h_FDR0005_color == "D")
CY20_1h_FDR0005_down <- sum(common_network$CY20_1h_FDR0005_color == "D")
CY20_3h_FDR0005_down <- sum(common_network$CY20_3h_FDR0005_color == "D")
CY20_12h_FDR0005_down <- sum(common_network$CY20_12h_FDR0005_color == "D")
CY20_24h_FDR0005_down <- sum(common_network$CY20_24h_FDR0005_color == "D")
CY20_48h_FDR0005_down <- sum(common_network$CY20_48h_FDR0005_color == "D")

#CYall FDR0005
df <- data.frame(expression_change = rep(c("01up", "02down"), times = 15),
                 Numgenes = c(CY15_1h_FDR0005_up, CY15_1h_FDR0005_down,
                              CY15_3h_FDR0005_up, CY15_3h_FDR0005_down,
                              CY15_12h_FDR0005_up, CY15_12h_FDR0005_down,
                              CY15_24h_FDR0005_up, CY15_24h_FDR0005_down,
                              CY15_48h_FDR0005_up, CY15_48h_FDR0005_down,
                              CY16_1h_FDR0005_up, CY16_1h_FDR0005_down,
                              CY16_3h_FDR0005_up, CY16_3h_FDR0005_down,
                              CY16_12h_FDR0005_up, CY16_12h_FDR0005_down,
                              CY16_24h_FDR0005_up, CY16_24h_FDR0005_down,
                              CY16_48h_FDR0005_up, CY16_48h_FDR0005_down,
                              CY20_1h_FDR0005_up, CY20_1h_FDR0005_down,
                              CY20_3h_FDR0005_up, CY20_3h_FDR0005_down,
                              CY20_12h_FDR0005_up, CY20_12h_FDR0005_down,
                              CY20_24h_FDR0005_up, CY20_24h_FDR0005_down,
                              CY20_48h_FDR0005_up, CY20_48h_FDR0005_down),
                 Category = rep(rep(c("011h", "023h", "0312h", "0424h", "0548h"), each = 2), times = 3),
                 CY = rep(c("CY15", "CY16", "CY20"), each = 10)
)
library(ggplot2)
library(reshape2)
#install.packages("RColorBrewer", dependencies = TRUE)
#library(RColorBrewer)
#display.brewer.all()

levels(df$Category) <- c("1h", "3h", "12h", "24h", "48h")
levels(df$expression_change) <- c("up", "down")
g <- ggplot(
  df,
  aes (
    x = Category,             # 遺伝子別でグルーピング
    y = Numgenes,
    fill = expression_change       # 縦軸を生物種で fill up
  )
)
g <- g + geom_bar(stat = "identity")
g <- g +  theme_bw()
g <- g + scale_colour_manual(values = "Set1")
g <- g + ylab("Number of DEGs")
g <- g + xlab("time-course")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)
CYall <- g + facet_grid(CY ~ .)
plot(CYall)
