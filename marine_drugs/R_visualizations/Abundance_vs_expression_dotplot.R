library(ggplot2)

data <- read.csv("abundance_vs_expression.txt", header=TRUE, sep="\t", row.names=1)

mean_plot <- ggplot(data, aes(x=Rel_cov, y=Average_expression, group=Sponge, size=4)) + geom_point(aes(shape=Sponge, color=Phylum))

median_plot <- ggplot(data, aes(x=Rel_cov, y=Median_expression, group=Sponge, size=4)) + geom_point(aes(shape=Sponge, color=Phylum))
