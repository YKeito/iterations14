sample <- rep("grey20", times = length(subcluster1_node$name))
names(sample) <- subcluster1_node$name
sample[match(subcluster1_node$name[subcluster1_node$TF == "TF"], names(sample))] <- "tomato"

g <- data.frame(AGI = subcluster1_node$name,
                degree = subcluster1_node$Degree,
                BetweennessCentrality = subcluster1_node$BetweennessCentrality,
                color = sample,
                stringsAsFactors = F
)

gg <- ggplot(g, aes(x = BetweennessCentrality, y = degree, color = sample))
gg <- gg + geom_point()
gg <- gg +  theme_bw()
gg <- gg +scale_color_identity(guide = "legend")
gg <- gg + theme(legend.position = 'none')
plot(gg)
ggsave(file = "~/Nakano_RNAseq/network_analysis/results_Fig/subcluster1/degree_betweennesscentrality.png", plot = gg)
