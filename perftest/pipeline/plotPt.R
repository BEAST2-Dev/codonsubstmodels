# plot P(t)

library(tidyverse)
library(ggplot2)
library(reshape2)

WD=paste0("~/WorkSpace/codonsubstmodels/perftest/")
setwd(WD)

pt.all <- read_delim("../pt.txt", "\t", comment = "#", col_names = F)
# 3601 is NA
pt.all <- pt.all[colSums(!is.na(pt.all)) > 0]

step = 0.01

pt <- pt.all %>% select(1,2,3,6) 
# 		 AAA AAC AAG
# 		  K   N   K
#AAA	K	.   4   1	
colnames(pt) <- c("AAA->AAA","AAA->AAC","AAA->AAG","AAA->ACC")
#"AAA->AAC", K -> N, codon changes in more than one codon position
#"AAA->AAG", K -> K, synonymous transition
#"AAA->ACC", K -> T, non-synonymous transversion
pt[["x"]] = 1:nrow(pt.all) * step

pt.melt <- melt(pt, id = c("x"))
colnames(pt.melt) <- c("x","i->j","P(t)")

p <- ggplot(pt.melt, aes(x, value, colour=variable)) + 
  geom_line(size=1, alpha=.5) + ylim(0,0.5) +
  ggtitle("") + xlab("time") + ylab("P(t)") + 
  guides(colour=guide_legend(title="i->j")) +
  theme_minimal() 

ggsave("../figures/pt.pdf", p, width = 6, height = 5)  
  
  
