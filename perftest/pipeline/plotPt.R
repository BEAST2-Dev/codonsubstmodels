# plot P(t)

library(tidyverse)
library(ggplot2)
library(reshape2)

WD=paste0("~/WorkSpace/codonsubstmodels/perftest/")
setwd(WD)

pt.all <- read_delim("../p_d_0.08_15.txt", "\t", comment = "#", col_names = F)
# 3601 is NA
pt.all <- pt.all[colSums(!is.na(pt.all)) > 0]
# equal freq, 0.0167  
tail(pt.all, 10)
nrow(pt.all)

# 4 cols
pt <- pt.all %>% select(1:4,7,32,(5*60+7)) 
colnames(pt)[1] <- "distance"
# 		 AAA AAC AAG
# 		  K   N   K
#AAA	K	.   4   1	
colnames(pt)[2:ncol(pt)] <- c("AAA->AAA","AAA->AAC","AAA->AAG","AAA->ACC","AAA->GAA","ACC->ACA")
#"AAA->AAC", K -> N, 4: non-synonymous transversion
#"AAA->AAG", K -> K, 1: synonymous transition
#"AAA->ACC", K -> T, 0: codon changes in more than one codon position
#"AAA->GAA", K -> E, 3: non-synonymous transition
#"ACC->ACA", T -> T, 2: synonymous transversion

pt.melt <- melt(pt, id = c("distance"))

p <- ggplot(pt.melt, aes(distance, value, colour=variable)) + 
  geom_point(shape=1, size=0.5, alpha=.9) +  # geom_point(size=.2, alpha=.5)
  ylim(0,0.8) + scale_color_brewer(palette="Set1") +
  ggtitle(paste(nrow(pt), "points")) + xlab("distance") + ylab("P(d)") + 
  guides(colour=guide_legend(title="i->j")) +
  theme_minimal() 

ggsave("../figures/p_d_.pdf", p, width = 6, height = 5)  
  

pt.all <- read_delim("../p_d_0.01_15.txt", "\t", comment = "#", col_names = F)
pt.all <- pt.all[colSums(!is.na(pt.all)) > 0]
pt2 <- pt.all %>% select(1:4,7,32,(5*60+7)) 
colnames(pt2)[1] <- "distance"
colnames(pt2)[2:ncol(pt2)] <- c("AAA->AAA","AAA->AAC","AAA->AAG","AAA->ACC","AAA->GAA","ACC->ACA")

pt.all <- read_delim("../p_d_0.5_15.txt", "\t", comment = "#", col_names = F)
pt.all <- pt.all[colSums(!is.na(pt.all)) > 0]
pt3 <- pt.all %>% select(1:4,7,32,(5*60+7)) 
colnames(pt3)[1] <- "distance"
colnames(pt3)[2:ncol(pt3)] <- c("AAA->AAA","AAA->AAC","AAA->AAG","AAA->ACC","AAA->GAA","ACC->ACA")

pt.all <- read_delim("../p_d_1_15.txt", "\t", comment = "#", col_names = F)
pt.all <- pt.all[colSums(!is.na(pt.all)) > 0]
pt4 <- pt.all %>% select(1:4,7,32,(5*60+7)) 
colnames(pt4)[1] <- "distance"
colnames(pt4)[2:ncol(pt4)] <- c("AAA->AAA","AAA->AAC","AAA->AAG","AAA->ACC","AAA->GAA","ACC->ACA")

pt.omega <- pt[,1:2]
pt.omega$omega = 0.08

pt.omega2 <- pt2[,1:2]
pt.omega2$omega = 0.01

pt.omega3 <- pt3[,1:2]
pt.omega3$omega = 0.5

pt.omega4 <- pt4[,1:2]
pt.omega4$omega = 1

pt.melt <- rbind(pt.omega, pt.omega2, pt.omega3, pt.omega4)
pt.melt$omega <- as.character(pt.melt$omega)

p <- ggplot(pt.melt, aes(distance, `AAA->AAA`, colour=omega)) + 
  geom_point(shape=1, size=0.5, alpha=.9) +  # geom_point(size=.2, alpha=.5)
  scale_color_brewer(palette="Set1") +
  ggtitle("AAA->AAA") + xlab("distance") + ylab("P(d)") + 
  theme_minimal() 

ggsave("../figures/p_d_omega.pdf", p, width = 6, height = 5)  



