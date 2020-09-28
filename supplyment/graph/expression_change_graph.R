#重ねグラフ#
#CY15
df <- data.frame(expression_change = rep(c("up", "down"), times = 5), 
                 Numgenes = c(length(CY15_1h_FDR0005_up_color), length(CY15_1h_FDR0005_down_color), 
                              length(CY15_3h_FDR0005_up_color), length(CY15_3h_FDR0005_down_color), 
                              length(CY15_12h_FDR0005_up_color), length(CY15_12h_FDR0005_down_color), 
                              length(CY15_24h_FDR0005_up_color), length(CY15_24h_FDR0005_down_color), 
                              length(CY15_48h_FDR0005_up_color), length(CY15_48h_FDR0005_down_color)), 
                 Category = rep(c("011h", "023h", "0312h", "0424h", "0548h"), each = 2)
                 )
library(ggplot2)
library(reshape2)
#install.packages("RColorBrewer", dependencies = TRUE)
#library(RColorBrewer)
#display.brewer.all()

levels(df$Category) <- c("1h", "3h", "12h", "24h", "48h")

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
g <- g + geom_text(aes(label = Numgenes), size = 5, hjust = 0.5, vjust = 2.0, position = "stack")
g <- g + scale_colour_manual(values = "Set1")
g <- g + ylab("Number of DEGs")
g <- g + xlab("time-course")
g <- g + labs(title = "CY15 DEGs ( q < 0.005 )")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)

plot(g)
#CY16
df <- data.frame(expression_change = rep(c("up", "down"), times = 5), 
                 Numgenes = c(length(CY16_1h_FDR0005_up_color), length(CY16_1h_FDR0005_down_color), 
                              length(CY16_3h_FDR0005_up_color), length(CY16_3h_FDR0005_down_color), 
                              length(CY16_12h_FDR0005_up_color), length(CY16_12h_FDR0005_down_color), 
                              length(CY16_24h_FDR0005_up_color), length(CY16_24h_FDR0005_down_color), 
                              length(CY16_48h_FDR0005_up_color), length(CY16_48h_FDR0005_down_color)), 
                 Category = rep(c("011h", "023h", "0312h", "0424h", "0548h"), each = 2)
)
library(ggplot2)
library(reshape2)
#install.packages("RColorBrewer", dependencies = TRUE)
#library(RColorBrewer)
#display.brewer.all()

levels(df$Category) <- c("1h", "3h", "12h", "24h", "48h")

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
g <- g + geom_text(aes(label = Numgenes), size = 5, hjust = 0.5, vjust = 1.1, position = "stack")
g <- g + scale_colour_manual(values = "Set1")
g <- g + ylab("Number of DEGs")
g <- g + xlab("time-course")
g <- g + labs(title = "CY16 DEGs ( q < 0.005 )")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)

plot(g)
#CY20
df <- data.frame(expression_change = rep(c("up", "down"), times = 5), 
                 Numgenes = c(length(CY20_1h_FDR0005_up_color), length(CY20_1h_FDR0005_down_color), 
                              length(CY20_3h_FDR0005_up_color), length(CY20_3h_FDR0005_down_color), 
                              length(CY20_12h_FDR0005_up_color), length(CY20_12h_FDR0005_down_color), 
                              length(CY20_24h_FDR0005_up_color), length(CY20_24h_FDR0005_down_color), 
                              length(CY20_48h_FDR0005_up_color), length(CY20_48h_FDR0005_down_color)), 
                 Category = rep(c("011h", "023h", "0312h", "0424h", "0548h"), each = 2)
)
library(ggplot2)
library(reshape2)
#install.packages("RColorBrewer", dependencies = TRUE)
#library(RColorBrewer)
#display.brewer.all()

levels(df$Category) <- c("1h", "3h", "12h", "24h", "48h")

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
g <- g + geom_text(aes(label = Numgenes), size = 5, hjust = 0.5, vjust = 1.1, position = "stack")
g <- g + scale_colour_manual(values = "Set1")
g <- g + ylab("Number of DEGs")
g <- g + xlab("time-course")
g <- g + labs(title = "CY20 DEGs ( q < 0.005 )")
g <- g + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
g <- g + theme_bw(base_size = 20)

plot(g)