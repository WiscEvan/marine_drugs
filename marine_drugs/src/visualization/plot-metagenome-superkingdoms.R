library(ggplot2)

# ToColor-blind friendly palettes as discussed here: 
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# scale_fill_manual(values=cbPalette)


setwd("~/marine_drugs/marine_drugs/data/interim/host-assembly/fastas/")

df <- read.table("sponges_superkingdom_counts_stacked.tsv", sep='\t', header = TRUE)
colnames(df)

p <- ggplot(data=df, 
            mapping=aes(x=reorder(sponge, -count), 
                        y=count, 
                        fill=superkingdom)) +
  geom_bar(stat='identity') + 
  theme_bw() +
  coord_flip() +
  labs(
    x="Sponge",
    y='Contig Count', 
    title="Sponge Metagenome Superkingdom counts")
p
ggsave("~/marine_drugs/marine_drugs/reports/figures/sponge-metagenome-superkingdom-counts.pdf", p)  

