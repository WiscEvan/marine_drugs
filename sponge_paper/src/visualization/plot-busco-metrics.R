library(ggplot2)

# ToColor-blind friendly palettes as discussed here: 
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# scale_fill_manual(values=cbPalette)


setwd("~/sponge_paper/sponge_paper/data/interim/host-annotation/busco-marker-assessment/")

df <- read.table("all_analyses_busco_metrics.tsv", sep='\t', header = TRUE)
colnames(df)


p <- ggplot(data=df, 
       mapping=aes(x=reorder(sample, -Percentage.Complete), 
                   y=Percentage.Complete, 
                   fill=busco.mode)) +
  geom_bar(stat='identity', position='dodge') + 
  theme_bw() +
  coord_flip() +
  scale_fill_discrete(name= "BUSCO dataset annotation",
                      labels=c("ab initio from 'eukaryotic' contigs",
                               "ab initio on entire metagenome",
                               "AUGUSTUS proteins for entire metagenome ",
                               "AUGUSTUS proteins from 'eukaryotic' contigs"))+
  labs(
    x="Sample",
    y='Completeness (%)', 
    title="Sponge Marker Completeness Estimation", 
    subtitle = "BUSCO marker assessment")
p
ggsave("~/sponge_paper/sponge_paper/reports/figures/sponge-marker-completeness-estimation.pdf", p)  

