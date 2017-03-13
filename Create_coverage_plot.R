require(ggplot2)
cov = read.table('coverage.tsv', header = T)
p <- ggplot(cov, aes(Ref.pos, Coverage))
p + geom_point() + ggtitle("Change in Coverage Over Each Position in the Reference Genome") + 
  labs(x = "Position in Reference Genome", y = "Count")
ggsave("coverage.pdf", device = "pdf")