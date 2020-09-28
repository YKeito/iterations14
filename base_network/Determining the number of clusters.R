determinate_Numsubcluster <- data.frame(Category = rep(c("subcluster"), times = 18), 
                                        iterations = c(1:16, 20, 60), 
                                        Num =  c(rep(1, times = 4), 2, 14, 86, 225, 322, 399, 418, 432, 437, rep(442, times = 5))
                                        )
library(ggplot2)
determinate_Numsubcluster <- 
ggplot(determinate_Numsubcluster, 
       aes(x = iterations,
           y = Num,
           )
       ) + 
  geom_line() +
  ylab("number of subclusters") + 
  ggtitle("Determining the number of subclusters") +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20)) +
  theme_bw(base_size = 20)
ggsave(file = "~/Nakano_RNAseq/network_analysis/cytoscape/image/Determining the number of subclusters.png", plot = determinate_Numsubcluster, dpi = 100)