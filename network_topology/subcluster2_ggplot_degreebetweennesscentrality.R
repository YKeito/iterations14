#subcluster2_node <- read.table("~/Nakano_RNAseq/network_analysis/base/Subcluster_Information/subcluster_networkanalyzer/subcluster2nodetable.csv", sep = ",", header = T)
#subcluster2_node <- subcluster2_node[order(subcluster2_node$Degree, decreasing=T), ]

sample <- rep("grey20", times = length(subcluster2_node$name))
names(sample) <- subcluster2_node$name
sample[match(subcluster2_node$name[subcluster2_node$TF == "TF"], names(sample))] <- "tomato"

g <- data.frame(AGI = subcluster2_node$name,
                degree = subcluster2_node$Degree,
                BetweennessCentrality = subcluster2_node$BetweennessCentrality,
                color = sample,
                stringsAsFactors = F
)

gg <- ggplot(g, aes(x = BetweennessCentrality, y = degree, color = sample))
gg <- gg + geom_point()
gg <- gg +  theme_bw()
gg <- gg +scale_color_identity(guide = "legend")
gg <- gg + theme(legend.position = 'none')
plot(gg)
ggsave(file = "~/Nakano_RNAseq/network_analysis/results_Fig/subcluster2/degree_betweennesscentrality.png", plot = gg)
