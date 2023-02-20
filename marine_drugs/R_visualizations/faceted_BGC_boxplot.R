library(ggplot2)
library(reshape2)

data <- read.csv("All_boxplot_data.txt", header=TRUE, sep="\t", row.names=1)

mean_plot <- ggplot(data, aes(x=Type, y=Average_expression, group=Type)) + geom_boxplot(aes(fill=Type))+ facet_grid(. ~ Sponge_ID) + theme(axis.text.x = element_text(angle = 90))+ coord_cartesian(ylim = c(0, 0.0025))
median_plot <- ggplot(data, aes(x=Type, y=Median_expression, group=Type)) + geom_boxplot(aes(fill=Type))+ facet_grid(. ~ Sponge_ID) + theme(axis.text.x = element_text(angle = 90))
median_plot_zoomed <- ggplot(data, aes(x=Type, y=Median_expression, group=Type)) + geom_boxplot(aes(fill=Type))+ facet_grid(. ~ Sponge_ID) + theme(axis.text.x = element_text(angle = 90))+ coord_cartesian(ylim = c(0, 0.0005))

BGC_phyla_plot <- ggplot(data, aes(x=Phylum, y=Average_expression, group=Phylum)) + geom_boxplot(aes(fill=Phylum)) + facet_grid(. ~ Type) + facet_wrap(. ~ Type,ncol=1) + theme(axis.text.x = element_text(angle = 90))

#Subset out PKS_Other and Saccharide (waste of space - little to show)

BGC_phyla_avg_plot_zoomed <- ggplot(subset(data, Type == "NRPS" | Type == "RiPP"| Type == "Terpene" | Type == "Other"| Type == "PKSI"), aes(x=Phylum, y=Average_expression, group=Phylum)) + geom_boxplot(aes(fill=Phylum)) + facet_grid(. ~ Type) + facet_wrap(. ~ Type,ncol=1) + theme(axis.text.x = element_text(angle = 90))+ coord_cartesian(ylim = c(0, 0.003))
BGC_phyla_med_plot_zoomed <- ggplot(subset(data, Type == "NRPS" | Type == "RiPP"| Type == "Terpene" | Type == "Other"| Type == "PKSI"), aes(x=Phylum, y=Median_expression, group=Phylum)) + geom_boxplot(aes(fill=Phylum)) + facet_grid(. ~ Type) + facet_wrap(. ~ Type,ncol=1) + theme(axis.text.x = element_text(angle = 90))+ coord_cartesian(ylim = c(0, 0.003))
