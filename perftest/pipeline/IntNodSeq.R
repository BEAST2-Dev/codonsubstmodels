# internal node sequence analysis

library(tidyverse)

nTaxa = 32
WD=paste0("~/WorkSpace/codonsubstmodels/perftest/T",nTaxa)
setwd(WD)

burnin=0.1
# m0.da.ins.txt
ins.log <- file.path(paste0("4t",nTaxa), "m0.da.ins.txt")
# Sample  Node States  
seqs <- read_delim(ins.log, "\t", comment = "#", col_types = cols( States = col_character() ))

# MCMC summary
samples = unique(seqs[["Sample"]])
samples = samples[!is.na(samples)]
cat("Chain length", prettyNum(samples[length(samples)], big.mark=",",scientific=FALSE), 
    ", log every", prettyNum(samples[2], big.mark=",",scientific=FALSE), "samples.\n")

# state is a 2-digit integer from 01 to 60, n.codon = str.len / 2 
str.len = str_length(seqs[["States"]][1])
n.codon = str.len / 2 

node.id = 33
# internal node index starts from nTaxa
sites <- seqs %>% filter(Node==node.id)
# paste0("c", 1:1000), seq(2, 2000, 2)
sites <- sites %>% separate(States, into = paste0("c", 1:n.codon), sep = seq(2, str.len, 2))

stats <- NULL
for (c in 1:n.codon) {
  freq <- table(sites[[paste0("c", c)]])
  freq <- sort(freq, decreasing=T)
  s1 <- tibble(state=names(freq), freq=freq, site=c, order=1:length(freq))
  stats <- bind_rows(stats, s1)
}
stats <- stats %>% mutate(order = as.character(order))

library(ggplot2)

p <- ggplot(stats, aes(site, freq, fill = order)) + 
  geom_bar(stat="identity", alpha = 0.7) + scale_fill_brewer(palette="Set1") +
#  scale_y_log10() +
  ggtitle(paste("Internal Node", node.id)) + ylab("Codon frequency") +
  theme_minimal()
ggsave(paste0("node",node.id,"-hist.pdf"), p, width = 10, height = 5)


#############

#sites %>% select(paste0("c", 1:n.codon)) %>% 

freq <- table(sites[["c6"]])
freq <- sort(freq, decreasing=T)

s1 <- tibble(state=names(freq), freq=freq) %>% mutate(site=1)

# relative br lens to truth
p <- ggplot(s1, aes(site, freq, fill = state)) + 
  geom_bar(stat="identity") +
  #ggtitle(paste(nTaxa, "Taxa", branches, "Branches")) + ylab("ESS / hour") +
  theme_minimal()
