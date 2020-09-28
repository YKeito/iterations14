#install.packages("ggsignif")
library("ggsignif")
#View(allRNASeq[rownames(allRNASeq) == "AT1G74710", 19:36])
####AT1G22070 TGA3####
data <- allRNASeq[rownames(allRNASeq) == "AT2G26020" | rownames(allRNASeq) == "AT5G44420" | 
                    rownames(allRNASeq) == "AT5G44430" | rownames(allRNASeq) == "AT2G14610" | 
                    rownames(allRNASeq) == "AT1G64280" | rownames(allRNASeq) == "AT1G22070" | 
                    rownames(allRNASeq) == "AT1G32640", ]
df <- data.frame(time = rep(c("011h", "023h", "0312h", "0424h", "0548h"), time = 3), 
                 expression_change = c(data$CY15_1h[1], data$CY15_3h[1] ,data$CY15_12h[1], data$CY15_24h[1], data$CY15_48h[1],
                                       data$CY16_1h[1], data$CY16_3h[1] ,data$CY16_12h[1], data$CY16_24h[1], data$CY16_48h[1],
                                       data$CY20_1h[1], data$CY20_3h[1] ,data$CY20_12h[1], data$CY20_24h[1], data$CY20_48h[1]),
                 Category = rep(c("CY15", "CY16", "CY20"), each = 5)
)
library(ggplot2)
library(reshape2)
#install.packages("RColorBrewer", dependencies = TRUE)
#library(RColorBrewer)
#display.brewer.all()
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
g <- g + ggtitle("AT1G22070 TGA3")#TGA3
g <- g + xlab("")
g <- g + ylab("log(treatment/control)")
g <- g + geom_signif(annotations = c(" ", "*", "**", "**", " "),
                     size = 0,
                     y_position = c(df$expression_change[1] + 0.01, df$expression_change[2] + 0.01, df$expression_change[3] + 0.01, df$expression_change[4] + 0.01, df$expression_change[5] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#F8766D"
                     )
g <- g + geom_signif(annotations = c(" ", " ", " ", " ", " "),
                     size = 0,
                     y_position = c(df$expression_change[6] + 0.01, df$expression_change[7] + 0.01, df$expression_change[8] + 0.01, df$expression_change[9] + 0.01, df$expression_change[10] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#7CAE00"
)
TGA3 <- g + geom_signif(annotations = c(" ", " ", " ", " ", " "),
                     size = 0,
                     y_position = c(df$expression_change[11] + 0.01, df$expression_change[12] + 0.01, df$expression_change[13] + 0.01, df$expression_change[14] + 0.01, df$expression_change[15] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#00BFC4"
                     )
plot(TGA3)
ggsave(file = "~/Nakano_RNAseq/network_analysis/results_Fig/AT1G22070_TGA3.png", plot = TGA3, dpi = 100)

####AT1G32640 MYC2####
data <- allRNASeq[rownames(allRNASeq) == "AT2G26020" | rownames(allRNASeq) == "AT5G44420" | 
                    rownames(allRNASeq) == "AT5G44430" | rownames(allRNASeq) == "AT2G14610" | 
                    rownames(allRNASeq) == "AT1G64280" | rownames(allRNASeq) == "AT1G22070" | 
                    rownames(allRNASeq) == "AT1G32640", ]
df <- data.frame(time = rep(c("011h", "023h", "0312h", "0424h", "0548h"), time = 3), 
                 expression_change = c(data$CY15_1h[2], data$CY15_3h[2] ,data$CY15_12h[2], data$CY15_24h[2], data$CY15_48h[2],
                                       data$CY16_1h[2], data$CY16_3h[2] ,data$CY16_12h[2], data$CY16_24h[2], data$CY16_48h[2],
                                       data$CY20_1h[2], data$CY20_3h[2] ,data$CY20_12h[2], data$CY20_24h[2], data$CY20_48h[2]),
                 Category = rep(c("CY15", "CY16", "CY20"), each = 5)
)
library(ggplot2)
library(reshape2)
#install.packages("RColorBrewer", dependencies = TRUE)
#library(RColorBrewer)
#display.brewer.all()
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
g <- g + ggtitle("AT1G32640 MYC2")#MYC2
g <- g + xlab("")
g <- g + ylab("log(treatment/control)")
g <- g + geom_signif(annotations = c("**", "**", " ", " ", " "),
                     size = 0,
                     y_position = c(df$expression_change[1] + 0.01, df$expression_change[2] + 0.01, df$expression_change[3] + 0.01, df$expression_change[4] + 0.01, df$expression_change[5] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#F8766D"
)

g <- g + geom_signif(annotations = c(" ", " ", " ", "*", " "),
                     size = 0,
                     y_position = c(df$expression_change[6] + 0.01, df$expression_change[7] + 0.01, df$expression_change[8] + 0.01, df$expression_change[9] + 0.01, df$expression_change[10] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#7CAE00"
)

MYC2 <- g + geom_signif(annotations = c(" ", " ", " ", " ", " "),
                     size = 0,
                     y_position = c(df$expression_change[11] + 0.01, df$expression_change[12] + 0.01, df$expression_change[13] + 0.01, df$expression_change[14] + 0.01, df$expression_change[15] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#00BFC4"
)
plot(MYC2)
ggsave(file = "~/Nakano_RNAseq/network_analysis/results_Fig/AT1G32640_MYC2.png", plot = MYC2, dpi = 100)
####AT1G64280 NPR1####
data <- allRNASeq[rownames(allRNASeq) == "AT2G26020" | rownames(allRNASeq) == "AT5G44420" | 
                    rownames(allRNASeq) == "AT5G44430" | rownames(allRNASeq) == "AT2G14610" | 
                    rownames(allRNASeq) == "AT1G64280" | rownames(allRNASeq) == "AT1G22070" | 
                    rownames(allRNASeq) == "AT1G32640", ]
df <- data.frame(time = rep(c("011h", "023h", "0312h", "0424h", "0548h"), time = 3), 
                 expression_change = c(data$CY15_1h[3], data$CY15_3h[3] ,data$CY15_12h[3], data$CY15_24h[3], data$CY15_48h[3],
                                       data$CY16_1h[3], data$CY16_3h[3] ,data$CY16_12h[3], data$CY16_24h[3], data$CY16_48h[3],
                                       data$CY20_1h[3], data$CY20_3h[3] ,data$CY20_12h[3], data$CY20_24h[3], data$CY20_48h[3]),
                 Category = rep(c("CY15", "CY16", "CY20"), each = 5)
)
library(ggplot2)
library(reshape2)
#install.packages("RColorBrewer", dependencies = TRUE)
#library(RColorBrewer)
#display.brewer.all()
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
g <- g + ggtitle("AT1G64280 NPR1")#NPR1
g <- g + xlab("")
g <- g + ylab("log(treatment/control)")
g <- g + geom_signif(annotations = c(" ", " ", " ", "**", " "),
                     size = 0,
                     y_position = c(df$expression_change[1] + 0.01, df$expression_change[2] + 0.01, df$expression_change[3] + 0.01, df$expression_change[4] + 0.01, df$expression_change[5] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#F8766D"
)
g <- g + geom_signif(annotations = c(" ", " ", " ", " ", " "),
                     size = 0,
                     y_position = c(df$expression_change[6] + 0.01, df$expression_change[7] + 0.01, df$expression_change[8] + 0.01, df$expression_change[9] + 0.01, df$expression_change[10] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#7CAE00"
)
NPR1 <- g + geom_signif(annotations = c(" ", " ", " ", " ", " "),
                     size = 0,
                     y_position = c(df$expression_change[11] + 0.01, df$expression_change[12] + 0.01, df$expression_change[13] + 0.01, df$expression_change[14] + 0.01, df$expression_change[15] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#00BFC4"
)

plot(NPR1)
ggsave(file = "~/Nakano_RNAseq/network_analysis/results_Fig/AT1G64280_NPR1.png", plot = NPR1, dpi = 100)



####AT2G14610 ATPR1_PR 1_PR1__pathogenesis-related gene 1####
data <- allRNASeq[rownames(allRNASeq) == "AT2G26020" | rownames(allRNASeq) == "AT5G44420" | 
                    rownames(allRNASeq) == "AT5G44430" | rownames(allRNASeq) == "AT2G14610" | 
                    rownames(allRNASeq) == "AT1G64280" | rownames(allRNASeq) == "AT1G22070" | 
                    rownames(allRNASeq) == "AT1G32640", ]
df <- data.frame(time = rep(c("011h", "023h", "0312h", "0424h", "0548h"), time = 3), 
                 expression_change = c(data$CY15_1h[4], data$CY15_3h[4] ,data$CY15_12h[4], data$CY15_24h[4], data$CY15_48h[4],
                                       data$CY16_1h[4], data$CY16_3h[4] ,data$CY16_12h[4], data$CY16_24h[4], data$CY16_48h[4],
                                       data$CY20_1h[4], data$CY20_3h[4] ,data$CY20_12h[4], data$CY20_24h[4], data$CY20_48h[4]),
                 Category = rep(c("CY15", "CY16", "CY20"), each = 5)
)
library(ggplot2)
library(reshape2)
#install.packages("RColorBrewer", dependencies = TRUE)
#library(RColorBrewer)
#display.brewer.all()
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
g <- g + ggtitle("AT2G14610 PR1")#ATPR1_PR 1_PR1__pathogenesis-related gene 1
g <- g + xlab("")
g <- g + ylab("log(treatment/control)")
g <- g + geom_signif(annotations = c(" ", " ", "**", "**", " "),
                     size = 0,
                     y_position = c(df$expression_change[1] + 0.01, df$expression_change[2] + 0.01, df$expression_change[3] + 0.01, df$expression_change[4] + 0.01, df$expression_change[5] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#F8766D"
)
g <- g + geom_signif(annotations = c(" ", " ", " ", "**", "*"),
                     size = 0,
                     y_position = c(df$expression_change[6] + 0.01, df$expression_change[7] + 0.01, df$expression_change[8] + 0.01, df$expression_change[9] + 0.01, df$expression_change[10] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#7CAE00"
)
PR1 <- g + geom_signif(annotations = c(" ", " ", " ", "**", "**"),
                     size = 0,
                     y_position = c(df$expression_change[11] + 0.01, df$expression_change[12] + 0.01, df$expression_change[13] + 0.01, df$expression_change[14] + 0.01, df$expression_change[15] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#00BFC4"
)

plot(PR1)
ggsave(file = "~/Nakano_RNAseq/network_analysis/results_Fig/AT2G14610_PR1.png", plot = PR1, dpi = 100)


####AT2G26020 PDF1.2b__plant defensin 1.2b #-0.977224  4.587860####
data <- allRNASeq[rownames(allRNASeq) == "AT2G26020" | rownames(allRNASeq) == "AT5G44420" | 
                    rownames(allRNASeq) == "AT5G44430" | rownames(allRNASeq) == "AT2G14610" | 
                    rownames(allRNASeq) == "AT1G64280" | rownames(allRNASeq) == "AT1G22070" | 
                    rownames(allRNASeq) == "AT1G32640", ]
df <- data.frame(time = rep(c("011h", "023h", "0312h", "0424h", "0548h"), time = 3), 
                 expression_change = c(data$CY15_1h[5], data$CY15_3h[5] ,data$CY15_12h[5], data$CY15_24h[5], data$CY15_48h[5],
                                       data$CY16_1h[5], data$CY16_3h[5] ,data$CY16_12h[5], data$CY16_24h[5], data$CY16_48h[5],
                                       data$CY20_1h[5], data$CY20_3h[5] ,data$CY20_12h[5], data$CY20_24h[5], data$CY20_48h[5]),
                 Category = rep(c("CY15", "CY16", "CY20"), each = 5)
)
library(ggplot2)
library(reshape2)
#install.packages("RColorBrewer", dependencies = TRUE)
#library(RColorBrewer)
#display.brewer.all()
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
g <- g + xlab("")
g <- g + ylab("log(treatment/control)")
g <- g + ggtitle("AT2G26020 PDF1.2b")#PDF1.2b__plant defensin 1.2b
g <- g + geom_signif(annotations = c(" ", "**", "**", "", " "),
                     size = 0,
                     y_position = c(df$expression_change[1] + 0.01, df$expression_change[2] + 0.01, df$expression_change[3] + 0.01, df$expression_change[4] + 0.01, df$expression_change[5] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#F8766D"
)
g <- g + geom_signif(annotations = c(" ", " ", "**", " ", " "),
                     size = 0,
                     y_position = c(df$expression_change[6] + 0.01, df$expression_change[7] + 0.01, df$expression_change[8] + 0.01, df$expression_change[9] + 0.01, df$expression_change[10] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#7CAE00"
)
PDF1.2b <- g + geom_signif(annotations = c(" ", " ", "*", " ", " "),
                     size = 0,
                     y_position = c(df$expression_change[11] + 0.01, df$expression_change[12] + 0.01, df$expression_change[13] + 0.01, df$expression_change[14] + 0.01, df$expression_change[15] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#00BFC4"
)

plot(PDF1.2b)
ggsave(file = "~/Nakano_RNAseq/network_analysis/results_Fig/AT2G26020_PDF1.2b.png", plot = PDF1.2b, dpi = 100)
####AT5G44420 LCR77_PDF1.2_PDF1.2A__plant defensin 1.2 #-1.40608  6.46394####
data <- allRNASeq[rownames(allRNASeq) == "AT2G26020" | rownames(allRNASeq) == "AT5G44420" | 
                    rownames(allRNASeq) == "AT5G44430" | rownames(allRNASeq) == "AT2G14610" | 
                    rownames(allRNASeq) == "AT1G64280" | rownames(allRNASeq) == "AT1G22070" | 
                    rownames(allRNASeq) == "AT1G32640", ]
df <- data.frame(time = rep(c("011h", "023h", "0312h", "0424h", "0548h"), time = 3), 
                 expression_change = c(data$CY15_1h[6], data$CY15_3h[6] ,data$CY15_12h[6], data$CY15_24h[6], data$CY15_48h[6],
                                       data$CY16_1h[6], data$CY16_3h[6] ,data$CY16_12h[6], data$CY16_24h[6], data$CY16_48h[6],
                                       data$CY20_1h[6], data$CY20_3h[6] ,data$CY20_12h[6], data$CY20_24h[6], data$CY20_48h[6]),
                 Category = rep(c("CY15", "CY16", "CY20"), each = 5)
)
library(ggplot2)
library(reshape2)
#install.packages("RColorBrewer", dependencies = TRUE)
#library(RColorBrewer)
#display.brewer.all()
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
g <- g + ggtitle("AT5G44420 PDF1.2a")#LCR77_PDF1.2_PDF1.2A__plant defensin 1.2
g <- g + xlab("")
g <- g + ylab("log(treatment/control)")
g <- g + geom_signif(annotations = c(" ", "**", "**", "**", " "),
                     size = 0,
                     y_position = c(df$expression_change[1] + 0.01, df$expression_change[2] + 0.01, df$expression_change[3] + 0.01, df$expression_change[4] + 0.01, df$expression_change[5] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#F8766D"
)
g <- g + geom_signif(annotations = c(" ", " ", "**", "*", " "),
                     size = 0,
                     y_position = c(df$expression_change[6] + 0.01, df$expression_change[7] + 0.01, df$expression_change[8] + 0.01, df$expression_change[9] + 0.01, df$expression_change[10] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#7CAE00"
)
PDF1.2A <- g + geom_signif(annotations = c(" ", " ", "*", "*", " "),
                     size = 0,
                     y_position = c(df$expression_change[11] + 0.01, df$expression_change[12] + 0.01, df$expression_change[13] + 0.01, df$expression_change[14] + 0.01, df$expression_change[15] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#00BFC4"
)

plot(PDF1.2A)
ggsave(file = "~/Nakano_RNAseq/network_analysis/results_Fig/AT5G44420_PDF1.2a.png", plot = PDF1.2A, dpi = 100)

####AT5G44430 PDF1.2c__plant defensin 1.2C####
data <- allRNASeq[rownames(allRNASeq) == "AT2G26020" | rownames(allRNASeq) == "AT5G44420" | 
                    rownames(allRNASeq) == "AT5G44430" | rownames(allRNASeq) == "AT2G14610" | 
                    rownames(allRNASeq) == "AT1G64280" | rownames(allRNASeq) == "AT1G22070" | 
                    rownames(allRNASeq) == "AT1G32640", ]
df <- data.frame(time = rep(c("011h", "023h", "0312h", "0424h", "0548h"), time = 3), 
                 expression_change = c(data$CY15_1h[7], data$CY15_3h[7] ,data$CY15_12h[7], data$CY15_24h[7], data$CY15_48h[7],
                                       data$CY16_1h[7], data$CY16_3h[7] ,data$CY16_12h[7], data$CY16_24h[7], data$CY16_48h[7],
                                       data$CY20_1h[7], data$CY20_3h[7] ,data$CY20_12h[7], data$CY20_24h[7], data$CY20_48h[7]),
                 Category = rep(c("CY15", "CY16", "CY20"), each = 5)
)
library(ggplot2)
library(reshape2)
#install.packages("RColorBrewer", dependencies = TRUE)
#library(RColorBrewer)
#display.brewer.all()
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
g <- g + ggtitle("AT5G44430 PDF1.2c")#PDF1.2c__plant defensin 1.2C
g <- g + xlab("")
g <- g + ylab("log(treatment/control)")
g <- g + geom_signif(annotations = c(" ", " ", "**", "*", " "),
                     size = 0,
                     y_position = c(df$expression_change[1] + 0.01, df$expression_change[2] + 0.01, df$expression_change[3] + 0.01, df$expression_change[4] + 0.01, df$expression_change[5] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#F8766D"
)
g <- g + geom_signif(annotations = c(" ", " ", "**", " ", "*"),
                     size = 0,
                     y_position = c(df$expression_change[6] + 0.01, df$expression_change[7] + 0.01, df$expression_change[8] + 0.01, df$expression_change[9] + 0.01, df$expression_change[10] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#7CAE00"
)
PDF1.2c <- g + geom_signif(annotations = c(" ", " ", " ", " ", "*"),
                     size = 0,
                     y_position = c(df$expression_change[11] + 0.01, df$expression_change[12] + 0.01, df$expression_change[13] + 0.01, df$expression_change[14] + 0.01, df$expression_change[15] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#00BFC4"
)

plot(PDF1.2c)
ggsave(file = "~/Nakano_RNAseq/network_analysis/results_Fig/AT5G44430_PDF1.2c.png", plot = PDF1.2c, dpi = 100)

####check_cluster####
iterations14[iterations14$name == "AT5G44420", ]$X__mclCluster#pdf1.2b○
iterations14[iterations14$name == "AT2G26020", ]$X__mclCluster#pdf1.2a○
iterations14[iterations14$name == "AT5G44430", ]$X__mclCluster#pdf1.2c○
iterations14[iterations14$name == "AT2G14610", ]$X__mclCluster#PR1○
iterations14[iterations14$name == "AT1G64280", ]$X__mclCluster#NPR1○
iterations14[iterations14$name == "AT1G22070", ]$X__mclCluster#TGA3○
iterations14[iterations14$name == "AT1G32640", ]$X__mclCluster#MYC2○
iterations14[iterations14$name == "AT1G72260", ]$X__mclCluster#Thi2.1×
iterations14[iterations14$name == "AT2G39940", ]$X__mclCluster#COI1×
iterations14[iterations14$name == "AT1G66340", ]$X__mclCluster#ETR1×


####AT3G48100 ARR5####
data <- allRNASeq[rownames(allRNASeq) == "AT3G48100", ]
df <- data.frame(time = rep(c("011h", "023h", "0312h", "0424h", "0548h"), time = 3), 
                 expression_change = c(data$CY15_1h[1], data$CY15_3h[1] ,data$CY15_12h[1], data$CY15_24h[1], data$CY15_48h[1],
                                       data$CY16_1h[1], data$CY16_3h[1] ,data$CY16_12h[1], data$CY16_24h[1], data$CY16_48h[1],
                                       data$CY20_1h[1], data$CY20_3h[1] ,data$CY20_12h[1], data$CY20_24h[1], data$CY20_48h[1]),
                 Category = rep(c("CY15", "CY16", "CY20"), each = 5)
)
library(ggplot2)
library(reshape2)
#install.packages("RColorBrewer", dependencies = TRUE)
#library(RColorBrewer)
#display.brewer.all()
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
g <- g + ggtitle("AT3G48100 ARR5")#ARR5
g <- g + xlab("")
g <- g + ylab("log(treatment/control)")
g <- g + geom_signif(annotations = c("**", " ", " ", "", " "),
                     size = 0,
                     y_position = c(df$expression_change[1] + 0.01, df$expression_change[2] + 0.01, df$expression_change[3] + 0.01, df$expression_change[4] + 0.01, df$expression_change[5] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#F8766D"
)
g <- g + geom_signif(annotations = c("*", " ", " ", " ", " "),
                     size = 0,
                     y_position = c(df$expression_change[6] + 0.01, df$expression_change[7] + 0.01, df$expression_change[8] + 0.01, df$expression_change[9] + 0.01, df$expression_change[10] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#7CAE00"
)
ARR5 <- g + geom_signif(annotations = c(" ", " ", " ", " ", "*"),
                     size = 0,
                     y_position = c(df$expression_change[11] + 0.01, df$expression_change[12] + 0.01, df$expression_change[13] + 0.01, df$expression_change[14] + 0.01, df$expression_change[15] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#00BFC4"
)

plot(ARR5)
ggsave(file = "~/Nakano_RNAseq/network_analysis/results_Fig/AT3G48100_ARR5.png", plot = ARR5, dpi = 100)
####AT3G57040 ARR9####
data <- allRNASeq[rownames(allRNASeq) == "AT3G57040", ]
df <- data.frame(time = rep(c("011h", "023h", "0312h", "0424h", "0548h"), time = 3), 
                 expression_change = c(data$CY15_1h[1], data$CY15_3h[1] ,data$CY15_12h[1], data$CY15_24h[1], data$CY15_48h[1],
                                       data$CY16_1h[1], data$CY16_3h[1] ,data$CY16_12h[1], data$CY16_24h[1], data$CY16_48h[1],
                                       data$CY20_1h[1], data$CY20_3h[1] ,data$CY20_12h[1], data$CY20_24h[1], data$CY20_48h[1]),
                 Category = rep(c("CY15", "CY16", "CY20"), each = 5)
)
library(ggplot2)
library(reshape2)
#install.packages("RColorBrewer", dependencies = TRUE)
#library(RColorBrewer)
#display.brewer.all()
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
g <- g + ggtitle("AT3G57040 ARR9")#ARR9
g <- g + xlab("")
g <- g + ylab("log(treatment/control)")
g <- g + geom_signif(annotations = c(" ", " ", "**", "**", " "),
                     size = 0,
                     y_position = c(df$expression_change[1] + 0.01, df$expression_change[2] + 0.01, df$expression_change[3] + 0.01, df$expression_change[4] + 0.01, df$expression_change[5] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#F8766D"
)
g <- g + geom_signif(annotations = c(" ", " ", " ", " ", " "),
                     size = 0,
                     y_position = c(df$expression_change[6] + 0.01, df$expression_change[7] + 0.01, df$expression_change[8] + 0.01, df$expression_change[9] + 0.01, df$expression_change[10] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#7CAE00"
)
ARR9 <- g + geom_signif(annotations = c(" ", " ", " ", " ", "**"),
                     size = 0,
                     y_position = c(df$expression_change[11] + 0.01, df$expression_change[12] + 0.01, df$expression_change[13] + 0.01, df$expression_change[14] + 0.01, df$expression_change[15] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#00BFC4"
)

plot(ARR9)
ggsave(file = "~/Nakano_RNAseq/network_analysis/results_Fig/AT3G57040_ARR9.png", plot = ARR9, dpi = 100)



####AT4G31920 ARR10####
data <- allRNASeq[rownames(allRNASeq) == "AT4G31920", ]
df <- data.frame(time = rep(c("011h", "023h", "0312h", "0424h", "0548h"), time = 3), 
                 expression_change = c(data$CY15_1h[1], data$CY15_3h[1] ,data$CY15_12h[1], data$CY15_24h[1], data$CY15_48h[1],
                                       data$CY16_1h[1], data$CY16_3h[1] ,data$CY16_12h[1], data$CY16_24h[1], data$CY16_48h[1],
                                       data$CY20_1h[1], data$CY20_3h[1] ,data$CY20_12h[1], data$CY20_24h[1], data$CY20_48h[1]),
                 Category = rep(c("CY15", "CY16", "CY20"), each = 5)
)
library(ggplot2)
library(reshape2)
#install.packages("RColorBrewer", dependencies = TRUE)
#library(RColorBrewer)
#display.brewer.all()
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
g <- g + ggtitle("AT4G31920 ARR10")#ARR10
g <- g + xlab("")
g <- g + ylab("log(treatment/control)")
g <- g + geom_signif(annotations = c("**", "**", " ", " ", " "),
                     size = 0,
                     y_position = c(df$expression_change[1] + 0.01, df$expression_change[2] + 0.01, df$expression_change[3] + 0.01, df$expression_change[4] + 0.01, df$expression_change[5] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#F8766D"
)
g <- g + geom_signif(annotations = c("**", "**", " ", "*", " "),
                     size = 0,
                     y_position = c(df$expression_change[6] + 0.01, df$expression_change[7] + 0.01, df$expression_change[8] + 0.01, df$expression_change[9] + 0.01, df$expression_change[10] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#7CAE00"
)
ARR10 <- g + geom_signif(annotations = c(" ", " ", " ", " ", " "),
                     size = 0,
                     y_position = c(df$expression_change[11] + 0.01, df$expression_change[12] + 0.01, df$expression_change[13] + 0.01, df$expression_change[14] + 0.01, df$expression_change[15] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#00BFC4"
)

plot(ARR10)
ggsave(file = "~/Nakano_RNAseq/network_analysis/results_Fig/AT4G31920_ARR10.png", plot = ARR10, dpi = 100)



####AT2G01830 AHK4####
data <- allRNASeq[rownames(allRNASeq) == "AT2G01830", ]
df <- data.frame(time = rep(c("011h", "023h", "0312h", "0424h", "0548h"), time = 3), 
                 expression_change = c(data$CY15_1h[1], data$CY15_3h[1] ,data$CY15_12h[1], data$CY15_24h[1], data$CY15_48h[1],
                                       data$CY16_1h[1], data$CY16_3h[1] ,data$CY16_12h[1], data$CY16_24h[1], data$CY16_48h[1],
                                       data$CY20_1h[1], data$CY20_3h[1] ,data$CY20_12h[1], data$CY20_24h[1], data$CY20_48h[1]),
                 Category = rep(c("CY15", "CY16", "CY20"), each = 5)
)
library(ggplot2)
library(reshape2)
#install.packages("RColorBrewer", dependencies = TRUE)
#library(RColorBrewer)
#display.brewer.all()
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
g <- g + ggtitle("AT2G01830 AHK4")#AHK4
g <- g + xlab("")
g <- g + ylab("log(treatment/control)")
g <- g + geom_signif(annotations = c(" ", " ", " ", "*", " "),
                     size = 0,
                     y_position = c(df$expression_change[1] + 0.01, df$expression_change[2] + 0.01, df$expression_change[3] + 0.01, df$expression_change[4] + 0.01, df$expression_change[5] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#F8766D"
)
g <- g + geom_signif(annotations = c(" ", " ", " ", " ", " "),
                     size = 0,
                     y_position = c(df$expression_change[6] + 0.01, df$expression_change[7] + 0.01, df$expression_change[8] + 0.01, df$expression_change[9] + 0.01, df$expression_change[10] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#7CAE00"
)
AHK4 <- g + geom_signif(annotations = c(" ", " ", " ", " ", " "),
                     size = 0,
                     y_position = c(df$expression_change[11] + 0.01, df$expression_change[12] + 0.01, df$expression_change[13] + 0.01, df$expression_change[14] + 0.01, df$expression_change[15] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#00BFC4"
)

plot(AHK4)
ggsave(file = "~/Nakano_RNAseq/network_analysis/results_Fig/AT2G01830_AHK4.png", plot = AHK4, dpi = 100)

####AT3G63110 IPT3####
data <- allRNASeq[rownames(allRNASeq) == "AT3G63110", ]
df <- data.frame(time = rep(c("011h", "023h", "0312h", "0424h", "0548h"), time = 3), 
                 expression_change = c(data$CY15_1h[1], data$CY15_3h[1] ,data$CY15_12h[1], data$CY15_24h[1], data$CY15_48h[1],
                                       data$CY16_1h[1], data$CY16_3h[1] ,data$CY16_12h[1], data$CY16_24h[1], data$CY16_48h[1],
                                       data$CY20_1h[1], data$CY20_3h[1] ,data$CY20_12h[1], data$CY20_24h[1], data$CY20_48h[1]),
                 Category = rep(c("CY15", "CY16", "CY20"), each = 5)
)
library(ggplot2)
library(reshape2)
#install.packages("RColorBrewer", dependencies = TRUE)
#library(RColorBrewer)
#display.brewer.all()
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
g <- g + ggtitle("AT3G63110 IPT3")#IPT3
g <- g + xlab("")
g <- g + ylab("log(treatment/control)")
g <- g + geom_signif(annotations = c(" ", "**", "**", "*", " "),
                     size = 0,
                     y_position = c(df$expression_change[1] + 0.01, df$expression_change[2] + 0.01, df$expression_change[3] + 0.01, df$expression_change[4] + 0.01, df$expression_change[5] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#F8766D"
)
g <- g + geom_signif(annotations = c(" ", " ", " ", "*", " "),
                     size = 0,
                     y_position = c(df$expression_change[6] + 0.01, df$expression_change[7] + 0.01, df$expression_change[8] + 0.01, df$expression_change[9] + 0.01, df$expression_change[10] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#7CAE00"
)
IPT3 <- g + geom_signif(annotations = c(" ", " ", " ", " ", " "),
                     size = 0,
                     y_position = c(df$expression_change[11] + 0.01, df$expression_change[12] + 0.01, df$expression_change[13] + 0.01, df$expression_change[14] + 0.01, df$expression_change[15] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#00BFC4"
)

plot(IPT3)
ggsave(file = "~/Nakano_RNAseq/network_analysis/results_Fig/AT3G63110_IPT3.png", plot = IPT3, dpi = 100)

####AT5G19040 IPT5####
data <- allRNASeq[rownames(allRNASeq) == "AT5G19040", ]
df <- data.frame(time = rep(c("011h", "023h", "0312h", "0424h", "0548h"), time = 3), 
                 expression_change = c(data$CY15_1h[1], data$CY15_3h[1] ,data$CY15_12h[1], data$CY15_24h[1], data$CY15_48h[1],
                                       data$CY16_1h[1], data$CY16_3h[1] ,data$CY16_12h[1], data$CY16_24h[1], data$CY16_48h[1],
                                       data$CY20_1h[1], data$CY20_3h[1] ,data$CY20_12h[1], data$CY20_24h[1], data$CY20_48h[1]),
                 Category = rep(c("CY15", "CY16", "CY20"), each = 5)
)
library(ggplot2)
library(reshape2)
#install.packages("RColorBrewer", dependencies = TRUE)
#library(RColorBrewer)
#display.brewer.all()
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
g <- g + ggtitle("AT5G19040 IPT5")#IPT5
g <- g + xlab("")
g <- g + ylab("log(treatment/control)")
g <- g + geom_signif(annotations = c("**", "**", "**", "*", " "),
                     size = 0,
                     y_position = c(df$expression_change[1] + 0.01, df$expression_change[2] + 0.01, df$expression_change[3] + 0.01, df$expression_change[4] + 0.01, df$expression_change[5] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#F8766D"
)
g <- g + geom_signif(annotations = c(" ", "**", " ", "*", " "),
                     size = 0,
                     y_position = c(df$expression_change[6] + 0.01, df$expression_change[7] + 0.01, df$expression_change[8] + 0.01, df$expression_change[9] + 0.01, df$expression_change[10] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#7CAE00"
)
IPT5 <- g + geom_signif(annotations = c(" ", " ", " ", " ", " "),
                     size = 0,
                     y_position = c(df$expression_change[11] + 0.01, df$expression_change[12] + 0.01, df$expression_change[13] + 0.01, df$expression_change[14] + 0.01, df$expression_change[15] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#00BFC4"
)

plot(IPT5)
ggsave(file = "~/Nakano_RNAseq/network_analysis/results_Fig/AT5G19040_IPT5.png", plot = IPT5, dpi = 100)


####AT1G74710 SID2####
data <- allRNASeq[rownames(allRNASeq) == "AT1G74710", ]
df <- data.frame(time = rep(c("011h", "023h", "0312h", "0424h", "0548h"), time = 3), 
                 expression_change = c(data$CY15_1h[1], data$CY15_3h[1] ,data$CY15_12h[1], data$CY15_24h[1], data$CY15_48h[1],
                                       data$CY16_1h[1], data$CY16_3h[1] ,data$CY16_12h[1], data$CY16_24h[1], data$CY16_48h[1],
                                       data$CY20_1h[1], data$CY20_3h[1] ,data$CY20_12h[1], data$CY20_24h[1], data$CY20_48h[1]),
                 Category = rep(c("CY15", "CY16", "CY20"), each = 5)
)
library(ggplot2)
library(reshape2)
#install.packages("RColorBrewer", dependencies = TRUE)
#library(RColorBrewer)
#display.brewer.all()
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
g <- g + ggtitle("AT1G74710 SID2")#SID2
g <- g + xlab("")
g <- g + ylab("log(treatment/control)")
g <- g + geom_signif(annotations = c(" ", "**", "**", "**", " "),
                     size = 0,
                     y_position = c(df$expression_change[1] + 0.01, df$expression_change[2] + 0.01, df$expression_change[3] + 0.01, df$expression_change[4] + 0.01, df$expression_change[5] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#F8766D"
)
g <- g + geom_signif(annotations = c(" ", " ", "*", " ", " "),
                     size = 0,
                     y_position = c(df$expression_change[6] + 0.01, df$expression_change[7] + 0.01, df$expression_change[8] + 0.01, df$expression_change[9] + 0.01, df$expression_change[10] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#7CAE00"
)
SID2 <- g + geom_signif(annotations = c(" ", " ", " ", " ", " "),
                     size = 0,
                     y_position = c(df$expression_change[11] + 0.01, df$expression_change[12] + 0.01, df$expression_change[13] + 0.01, df$expression_change[14] + 0.01, df$expression_change[15] + 0.01), 
                     xmin=c(1, 2, 3, 4, 5),
                     xmax = c(1, 2, 3, 4, 5),
                     tip_length = 0,
                     col = "#00BFC4"
)

plot(SID2)
ggsave(file = "~/Nakano_RNAseq/network_analysis/results_Fig/AT1G74710_SID2.png", plot = SID2, dpi = 100)

#install.packages("gridExtra")
library(gridExtra)
grid.arrange(SID2, NPR1, TGA3, PR1,
             ncol = 2)