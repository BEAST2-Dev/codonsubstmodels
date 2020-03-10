# plot P(t)

library(tidyverse)
library(ggplot2)
library(reshape2)

WD=paste0("~/WorkSpace/codonsubstmodels/perftest/")
setwd(WD)

pt.all <- read_delim("../p_t_.txt", "\t", comment = "#", col_names = F)
# 3601 is NA
pt.all <- pt.all[colSums(!is.na(pt.all)) > 0]
# equal freq, 0.0167  
tail(pt.all, 10)
nrow(pt.all)

# 4 cols
pt <- pt.all %>% select(1:4,7) 
colnames(pt)[1] <- "distance"
# 		 AAA AAC AAG
# 		  K   N   K
#AAA	K	.   4   1	
colnames(pt)[2:ncol(pt)] <- c("AAA->AAA","AAA->AAC","AAA->AAG","AAA->ACC")
#"AAA->AAC", K -> N, codon changes in more than one codon position
#"AAA->AAG", K -> K, synonymous transition
#"AAA->ACC", K -> T, non-synonymous transversion

pt.melt <- melt(pt, id = c("distance"))

p <- ggplot(pt.melt, aes(distance, value, colour=variable)) + 
  geom_point(shape=1, size=0.5, alpha=.5) +  # geom_point(size=.2, alpha=.5)
  ylim(0,0.5) +
  ggtitle("") + xlab("distance") + ylab("P(d)") + 
  guides(colour=guide_legend(title="i->j")) +
  theme_minimal() 

ggsave("../figures/pt.pdf", p, width = 6, height = 5)  
  
  
