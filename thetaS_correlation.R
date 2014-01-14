## thetaS_correlation.R by Rohan Maddamsetti.

## Do genes with high thetaS tend to have a higher copy number in the
## gut community metagenomics dataset that I found online?

hit.data <- read.csv("/Users/Rohandinho/Desktop/Projects/gene_population_sizes/data/gut_community_hits.csv")
thetaS.data <- read.csv("/Users/Rohandinho/Desktop/Projects/luscombe-reanalysis/data/Martincorena_Maddamsetti_thetaS_estimates.csv")

full.data <- merge(hit.data,thetaS.data)
library(ggplot2)

quartz()
qplot(data=full.data,x=Maddamsetti_thetaS,y=hits)
quartz()
qplot(data=full.data,x=Martincorena_thetaS,y=hits)

no.nas <- na.omit(full.data)
cor(x=no.nas$Martincorena_thetaS,y=no.nas$hits)
