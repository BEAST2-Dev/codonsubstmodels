# internal node sequence analysis

library(tidyverse)

WD="~/WorkSpace/codonsubstmodels/perftest/pipeline"
source(file.path(WD, "TraceUtils.R"))
source(file.path(WD, "Codon.R"))
source(file.path(WD, "EvolerUtils.R"))

n.taxa = 32
# perftest/T32
WD=paste0("~/WorkSpace/codonsubstmodels/perftest/T",n.taxa) 
setwd(WD)

# T32/4t32/m0.da.ins.txt
ins.log <- file.path(paste0("4t",n.taxa), "m0.da.ins.txt")

stats.list <- getIntNodeSeqStats(ins.log, burnin=0.1)
names(stats.list)
#### bug in logging nodeNr not nodeNr+1 ###
#names(stats.list) = as.character(as.integer(names(stats.list)) + 1)
### rm above line after run fixed jar.

### evolver.out
# Note: evolver tree uses the same format as APE tree
nod.states <- getSeqsEvoOut("ancestral.txt", n.taxa=n.taxa, genetic.code="vertebrateMitochondrial")

# node indexes should match
stopifnot(all(names(stats.list) == names(nod.states)))

### make sure the same node index is the same node
# BEAST
edges1 <- stats.list$edges %>% mutate(parent=as.integer(parent), child=as.integer(child)) %>% arrange(child)
# evolver.out
edges2 <- nod.states$edges %>% mutate(parent=as.integer(parent), child=as.integer(child)) %>% arrange(child)
# Note: the tip indexing system must be same
dup1 <- which(duplicated(edges1[["parent"]][1:n.taxa]))
dup2 <- which(duplicated(edges2[["parent"]][1:n.taxa]))
stopifnot(all(dup1 == dup2))

# map edges2 parents to edges1 parents by lineages from tips
edges.map <- edges1 %>% rename(pa1="parent", ch1="child") %>% 
  mutate( pa2=c(edges2[["parent"]][1:n.taxa],rep(NA, nrow(edges2)-n.taxa)),
          ch2=c(edges2[["child"]][1:n.taxa],rep(NA, nrow(edges2)-n.taxa)) )
# pa1   ch1   pa2   ch2
#<int> <int> <int> <int>
#1    35     1    37     1
#2    34     2    38     2
# find row numbers of parents of tips at child1 appearing in parent1,
pa1.idx <- match(edges1[["child"]][(n.taxa+1):nrow(edges1)], edges1[["parent"]][1:n.taxa])
# and use the same row numbers to find edges2 parents of tips, 
ch2 <- edges2[["parent"]][pa1.idx]
stopifnot(length(ch2) == nrow(edges.map)-n.taxa)
# then, fill selecetd internal nodes in child2 according to above mapping
edges.map[["ch2"]][(n.taxa+1):nrow(edges.map)] <- ch2

# fill in parent2 by edges2 mapping
pa2.idx <- match(ch2, edges2[["child"]])
pa2 <- edges2[["parent"]][pa2.idx]
edges.map[["pa2"]][(n.taxa+1):nrow(edges.map)] <- pa2

# repeat to fill in the rest of NA
while (anyNA(edges.map)) {
  na.idx <- which(is.na(edges.map[["pa2"]]))
  
  pa1.idx <- match(edges1[["child"]][na.idx], edges1[["parent"]])
  # use the same row numbers to find edges2 parents of tips  
  ch2 <- edges.map[["pa2"]][pa1.idx]
  stopifnot(length(ch2) == length(na.idx))
  edges.map[["ch2"]][na.idx] <- ch2
  
  # fill in parent2 by edges2 mapping
  pa2.idx <- match(ch2, edges2[["child"]])
  pa2 <- edges2[["parent"]][pa2.idx]
  edges.map[["pa2"]][na.idx] <- pa2
}
#print(edges.map, n=Inf)


### compare to true ancestral states
p.dist = c()
internal.nodes <- suppressWarnings(as.integer(names(stats.list)))
internal.nodes <- internal.nodes[!is.na(internal.nodes)]
for (nod.idx in internal.nodes) {
  # have to map nodes before checking states
  node.true <- edges.map %>% filter(pa1==as.integer(nod.idx)) %>% select(pa2) %>% unlist %>% unique
  stopifnot(length(node.true) == 1)
  
  state.true = as.integer(nod.states[[as.character(node.true)]])
  n.codon = length(state.true)
  
  # MAP maximum a posteriori
  map <- stats.list[[as.character(nod.idx)]] %>% filter(order=="1") %>% 
    mutate(state = as.integer(state)) %>% mutate(node.beast = as.integer(nod.idx)) 
  stopifnot(nrow(map) == n.codon)
  
  map <- map %>% mutate(node.true = as.character(node.true)) %>% 
    mutate(state.true = state.true) %>% 
    mutate(test = state.true == state)
  
  p.dist = c(p.dist, 1 - nrow(map[map[["test"]],]) / nrow(map) )
}
stats <- tibble(node = internal.nodes, p.dist = p.dist)
print(stats, n = Inf)


### compare parsimony ancerstral states in MCMC state0 to true states
#Sample  Node States                                                                                                                                    
#<dbl> <dbl> <chr>                                                                                                                                     
#1      0     0 35..1,34..2,33..3,
#2     NA    33 57574107420152
#3     NA    34 57574107420152
state0 <- read_delim(ins.log, "\t", comment = "#", n_max = n.taxa, 
                     col_types = cols( States = col_character() ))
# last node is root
stopifnot(state0[["Node"]][nrow(state0)]==2*n.taxa-1)
# rm node map in 1st row 
internal.nodes <- state0[["Node"]][state0[["Node"]]>0]
internal.nodes <- state0[["Node"]][state0[["Node"]]>0]
stopifnot(length(internal.nodes)==n.taxa-1)

p.dist = c()
for (nod.idx in internal.nodes) {
  # have to map nodes before checking states
  node.true <- edges.map %>% filter(pa1==as.integer(nod.idx)) %>% select(pa2) %>% unlist %>% unique
  stopifnot(length(node.true) == 1)
  
  state.true = as.integer(nod.states[[as.character(node.true)]])
  n.codon = length(state.true)
  
  # parsimony ancerstral states in MCMC state0 
  str.len = str_length(state0[["States"]][2]) # "States"
  stopifnot(str.len / 2 == n.codon) 
  sites <- state0 %>% filter(Node==nod.idx) %>% 
    separate(States, into = paste0("c", 1:n.codon), sep = seq(2, str.len, 2)) %>% 
    select(paste0("c", 1:n.codon))
  stopifnot(ncol(sites) == n.codon)
  
  state.pars = as.integer(sites[1,])
  
  p.dist = c(p.dist, 1 - length(state.pars[state.pars==state.true]) / length(state.pars) )
}
parsimony <- tibble(node = internal.nodes, p.dist = p.dist)
# compare to MAP
parsimony <- parsimony %>% full_join(stats, by="node") %>% 
  rename(p.dist.0="p.dist.x", p.dist.map="p.dist.y")
print(parsimony, n = Inf)


library(ggplot2)
require(reshape2)
### 95% credible set









######### replaced by using 95% credible set

### plot state0 vs MAP
data.m <- melt(parsimony, id='node')
colnames(data.m)[2] <- "method"
data.m$method <- gsub("p.dist.0", "parsimony", data.m$method)
data.m$method <- gsub("p.dist.map", "MAP", data.m$method)
data.m$method <- factor(data.m$method, levels = unique(data.m$method))

p <- ggplot(data.m, aes(node, value, fill = method)) + 
  geom_bar(position = "dodge", stat="identity") + 
  scale_y_sqrt() +
  ggtitle("Parsimony ASR vs. MAP") + ylab("p distance") +
  theme_minimal()
ggsave(paste0("parsimony-map.pdf"), p, width = 7, height = 5)

### plot MAP maximum a posteriori
nod.idx = 63 # root
codon.freq <- stats.list[[as.character(nod.idx)]] 
#  state  freq  site order
#  <chr> <int> <int> <int>
#1 57       37     1 1    
#2 59       19     1 2    

# top 5 freq
codon.freq <- codon.freq %>% filter(order < 6) %>% mutate(order = as.character(order)) 

p <- ggplot(codon.freq, aes(site, freq, fill = order)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values = c(alpha(c("lightblue"), .2), 
                               alpha(c("red","blue","yellow","purple","orange","brown","green"), .8))) +
  #  scale_y_log10() +
  ggtitle(paste("Internal Node", nod.idx)) + ylab("codon frequency") +
  theme_minimal()
ggsave(paste0("node",nod.idx,"-freq.pdf"), p, width = 10, height = 5)

### one site
site.idx = 1
freq.site <- codon.freq %>% filter(site == site.idx) 
freq.site$state <- factor(freq.site$state, levels = unique(freq.site$state))

p <- ggplot(freq.site, aes(site, freq, fill = state)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = c(alpha(c("lightblue"), .2), 
                               alpha(c("red","blue","yellow","purple","orange","brown","green"), .8))) +
  ggtitle(paste("Internal Node", nod.idx, "Site", site.idx)) + ylab("codon frequency") +
  theme_minimal() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave(paste0("node",nod.idx,"-site",site.idx,"-freq.pdf"), p, width = 5, height = 5)

