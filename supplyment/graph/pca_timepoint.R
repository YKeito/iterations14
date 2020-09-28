test <- t(allRNASeq[, 1:18])
pca_mod <- prcomp(test, scale. = TRUE)
summary(pca_mod)

df_pc <- data.frame(pca_mod$x[, 1:2], 
                    color = rownames(pca_mod$x))


p <- ggplot(df_pc, aes(PC1, PC2, label = color))
p <- p + geom_point()
p <- p + coord_cartesian(xlim = 1.2 * c(min(df_pc$PC1), max(df_pc$PC1)), ylim = 1.2 * c(min(df_pc$PC2), max(df_pc$PC2)))
plot(p)


install.packages("ggrepel") 
library("ggplot2")
library("ggrepel")
p <- ggplot(df_pc, aes(PC1, PC2, label = rownames(df_pc)))
p <- p + geom_point()
p <- p + geom_text_repel()
plot(p)