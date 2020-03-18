# plot P(t)

library(tidyverse)
library(ggplot2)
library(reshape2)

WD=paste0("~/WorkSpace/codonsubstmodels/perftest/")
setwd(WD)

getPt <- function(file) {
	pt.all <- read_delim(file, "\t", comment = "#", col_names = F)
	# 3601 is NA
	pt.all <- pt.all[colSums(!is.na(pt.all)) > 0]
	# equal freq, 0.0167  
	tail(pt.all, 10)
	nrow(pt.all)

	# 4 cols
	pt <- pt.all %>% select(1:4,7,32,(5*60+7)) 
	colnames(pt)[1] <- "distance"
	# 	   AAA AAC AAG
	# 		K   N   K
	#AAA	K	.   4   1	
	colnames(pt)[2:ncol(pt)] <- c("AAA->AAA","AAA->AAC","AAA->AAG","AAA->ACC","AAA->GAA","ACC->ACA")
	#"AAA->AAC", K -> N, 4: non-synonymous transversion
	#"AAA->AAG", K -> K, 1: synonymous transition
	#"AAA->ACC", K -> T, 0: codon changes in more than one codon position
	#"AAA->GAA", K -> E, 3: non-synonymous transition
	#"ACC->ACA", T -> T, 2: synonymous transversion
	return(pt)
}

pt <- getPt("../p_d_0.08_15.txt") 

pt.melt <- melt(pt, id = c("distance"))

p <- ggplot(pt.melt, aes(distance, value, colour=variable)) + 
  geom_point(shape=1, size=0.5, alpha=.9) +  # geom_point(size=.2, alpha=.5)
  ylim(0,0.8) + scale_color_brewer(palette="Set1") +
  ggtitle(paste(nrow(pt), "points")) + xlab("distance") + ylab("P(d)") + 
  guides(colour=guide_legend(title="i->j")) +
  theme_minimal() 

ggsave("../figures/p_d_.pdf", p, width = 6, height = 5)  

pt.omega <- pt[,1:2]
pt.omega$omega = 0.08
for (omega in c(0.01,0.5,1) {
  pt <- getPt(paste0("../p_d_",omega,"_15.txt")) 
  pt.tmp <- pt[,1:2]
  pt.tmp$omega = omega
  pt.omega <- rbind(pt.omega, pt.tmp)
}  
pt.omega$omega <- as.character(pt.omega$omega)

p <- ggplot(pt.omega, aes(distance, `AAA->AAA`, colour=omega)) + 
  geom_point(shape=1, size=0.5, alpha=.9) +  # geom_point(size=.2, alpha=.5)
  scale_color_brewer(palette="Set1") +
  ggtitle("AAA->AAA") + xlab("distance") + ylab("P(d)") + 
  theme_minimal() 

ggsave("../figures/p_d_omega.pdf", p, width = 6, height = 5)  



