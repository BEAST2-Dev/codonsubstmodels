# internal node sequence analysis

library(tidyverse)

WD="~/WorkSpace/codonsubstmodels/perftest/pipeline"
source(file.path(WD, "TraceUtils.R"))
source(file.path(WD, "Codon.R"))
source(file.path(WD, "EvolerUtils.R"))

n.taxa = 32
WD=paste0("~/WorkSpace/codonsubstmodels/perftest/T",n.taxa)
setwd(WD)

# m0.da.ins.txt
ins.log <- file.path(paste0("4t",n.taxa), "m0.da.ins.txt")

stats.list <- getIntNodeSeqStats(ins.log, burnin=0.1)

### evolver.out
# Note: evolver tree uses the same format as APE tree
nod.states <- getSeqsEvoOut("ancestral.txt", n.taxa=n.taxa, genetic.code="vertebrateMitochondrial")

for (node.id in nodes) {
  
}





node.id = 33

library(ggplot2)

p <- ggplot(stats, aes(site, freq, fill = order)) + 
  geom_bar(stat="identity") + #, alpha = 0.7) + #scale_fill_brewer(palette="Set1") +
  scale_fill_manual(values = c(alpha(c("lightblue"), .2), 
                               alpha(c("red","blue","yellow","purple","orange","brown","green"), .8))) +
#  scale_y_log10() +
  ggtitle(paste("Internal Node", node.id)) + ylab("codon frequency") +
  theme_minimal()
ggsave(paste0("node",node.id,"-freq.pdf"), p, width = 10, height = 5)


### MAP maximum a posteriori
map <- stats %>% filter(order=="1") %>% 
  mutate(state = as.integer(state)) %>% mutate(node = as.character(node.id))   
stopifnot(nrow(map) == n.codon)


# result
map <- map %>% mutate(true.state = as.integer(states.list[[as.character(node.id)]])) %>% 
  mutate(test = true.state == state)

nrow(map[map$test,]) / nrow(map)





#############

#sites %>% select(paste0("c", 1:n.codon)) %>% 

freq <- table(sites[["c6"]])
freq <- sort(freq, decreasing=T)

s1 <- tibble(state=names(freq), freq=freq) %>% mutate(site=1)

# relative br lens to truth
p <- ggplot(s1, aes(site, freq, fill = state)) + 
  geom_bar(stat="identity") +
  #ggtitle(paste(n.taxa, "Taxa", branches, "Branches")) + ylab("ESS / hour") +
  theme_minimal()
