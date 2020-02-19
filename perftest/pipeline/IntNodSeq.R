# internal node sequence analysis

library(tidyverse)

WD="~/WorkSpace/codonsubstmodels/perftest/pipeline"
source(file.path(WD, "TraceUtils.R"))
source(file.path(WD, "Codon.R"))
source(file.path(WD, "EvolerUtils.R"))
source(file.path(WD, "TreeUtils.R"))

n.taxa = 32
# perftest/T32
WD=paste0("~/WorkSpace/codonsubstmodels/perftest/T",n.taxa) 
setwd(WD)

# coal or yulelam10
tree.prior = "yulelam10"

### BEAST
# T32/4t32/m0.da.ins.txt
ins.log <- file.path(paste0("t",n.taxa,tree.prior,"DA"),"m0.da.ins.txt")

stats.list <- getIntNodeSeqStats(ins.log, burnin=0.1)
names(stats.list)

### evolver.out
# Note: evolver tree uses the same format as APE tree
nod.states <- getSeqsEvoOut(file.path(tree.prior, "ancestral.txt"), 
                            n.taxa=n.taxa, genetic.code="vertebrateMitochondrial")

stopifnot(length(nod.states[[1]]) == stats.list[["n.codons"]])
# same edges
stopifnot(nrow(nod.states[["edges"]]) == nrow(stats.list[["edges"]]))

### create edges.map to make sure the same node index is the same node 
# stats.list$edges is BEAST, nod.states$edges is evolver.out
edges.map <- mapEdges(stats.list$edges, nod.states$edges)
#print(edges.map, n=Inf)


######  require edges.map before this line ###### 

### 95% credible set
internal.nodes <- suppressWarnings(as.integer(names(stats.list)))
internal.nodes <- internal.nodes[!is.na(internal.nodes)]

### 1. sanity check how true codon falls into x% credible set
cred.check <- tibble(node=internal.nodes) 
for (cred.thre in c(0.05, 0.25, 0.50, 0.75, 0.95)) {
  in.cred = c()
  mean.tot.prob = c()
  for (nod.idx in internal.nodes) {
    # have to map nodes before checking states
    node.true <- edges.map %>% filter(pa1==as.integer(nod.idx)) %>% select(pa2) %>% unlist %>% unique
    stopifnot(length(node.true) == 1)
    
    state.true = as.integer(nod.states[[as.character(node.true)]])
    n.codon = length(state.true)
    
    # 95% credible set at each site (codon)
    freq.table <- stats.list[[as.character(nod.idx)]] %>% group_by(site) 
    out951st <- freq.table %>% filter(cred >= cred.thre) %>% slice(1)
    in95 <-  freq.table %>% filter(cred < cred.thre) %>% bind_rows(out951st) %>% arrange(site)
    stopifnot(length(unique(in95[["site"]])) == n.codon)
    
    cred.set <- in95 %>% group_by(site) %>% mutate(cred95 = paste0(as.integer(state), collapse = ",")) %>% 
      slice(n()) %>% select(site, cred, cred95) %>% arrange(site) %>% as_tibble() 
    stopifnot(nrow(cred.set) == n.codon)
    
    cred.set <- cred.set %>% mutate(state.true = as.character(state.true)) %>%
      mutate(cred95.set = strsplit(cred95, ",")) %>%
      mutate(in.cred = state.true %in% unlist(cred95.set)) 
    
    # convert to percentage
    in.per = nrow(cred.set[cred.set[["in.cred"]],]) / nrow(cred.set) * 100
    in.cred = c(in.cred,  in.per)
    mean.t.p = mean(cred.set[["cred"]])
    mean.tot.prob = c(mean.tot.prob, mean.t.p)
  }
  # all in 95% cred set
  in.cred
  mean.tot.prob
  stopifnot(length(in.cred) == n.taxa - 1 || any(in.cred < cred.thre))
  
  colnm = gsub("0\\.", "cred", as.character(cred.thre))
  colnm2 = paste0(colnm,".prob")
  cred.check <- cred.check %>% add_column(!!colnm := in.cred) %>% add_column(!!colnm2 := mean.tot.prob)
}
print(cred.check, n = Inf)
# crd.prob is the mean of the total probabilities inside the x% credible set
write_delim(cred.check, paste0("t",n.taxa,tree.prior,"-credset.txt"), delim = "\t")


### compare to true ancestral states
p.dist = c()
mean.map.prob = c()
sd.map.prob = c()
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
  mean.map.prob = c(mean.map.prob, mean(map[["prob"]]))
  sd.map.prob = c(sd.map.prob, sd(map[["prob"]]))
}
# p.dist of MAP to true codon, and mean/sd (across sites) of posterior prob that the MAP represents
stats <- tibble(node = internal.nodes, p.dist = p.dist, 
                map.prob.mean = mean.map.prob, map.prob.sd = sd.map.prob)
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

### plot state0 vs MAP
data.m <- parsimony %>% select(node,p.dist.0,p.dist.map) %>% melt(id='node')
colnames(data.m)[2] <- "method"
data.m$method <- gsub("p.dist.0", "parsimony", data.m$method)
data.m$method <- gsub("p.dist.map", "MAP", data.m$method)
data.m$method <- factor(data.m$method, levels = unique(data.m$method))

max.p <- max(data.m$value)
max.p
p <- ggplot(data.m, aes(node, value, fill = method, colour = method)) + 
  geom_bar(position = "dodge", stat="identity") + 
  scale_y_sqrt(breaks=c(0, 0.001, 0.01, 0.05, 0.1)) +
  ggtitle(paste("Parsimony ASR vs. MAP", n.taxa, "Taxa")) + ylab("p distance") +
  theme_minimal()
ggsave(paste0("t",n.taxa,tree.prior,"parsimony-map.pdf"), p, width = 7, height = 5)



### plot MAP maximum a posteriori
nod.idx = 63 # root
codon.freq <- stats.list[[as.character(nod.idx)]] 
#  state  freq  site order
#  <chr> <int> <int> <int>
#1 57       37     1 1    
#2 59       19     1 2    
unique(codon.freq$order)

# merge less freq codons
topn = 5
codon.freq.small <- codon.freq %>% group_by(site) %>% filter(order > topn) %>% 
  summarise(freq=sum(freq), prob=sum(prob), order = "rest") 
# add top 5 freq codons
codon.freq <- codon.freq %>% filter(order < topn+1) %>% select(site,freq,prob,order) %>% 
  mutate(order = as.character(order)) %>% bind_rows(codon.freq.small) 

temp <- codon.freq %>% group_by(site) %>% mutate(cred = cumsum(prob)) 
stopifnot(all(temp[["cred"]] <= 1))

# sort site by MAP prob
codon.prob <- codon.freq %>% mutate(site = as.character(site)) %>% group_by(site) %>% 
  mutate(map.prob = max(prob)) %>% arrange(desc(map.prob), site) 

codon.prob[["site"]] <- factor(codon.prob[["site"]], levels = unique(codon.prob[["site"]]))

# codon prob distribution 
p <- ggplot(codon.prob, aes(site, prob, fill = order)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values = c(alpha(c("lightblue"), .2), 
                               alpha(c("red","blue","yellow","purple","orange","brown","green"), .8))) +
  ggtitle(paste("Internal Node", nod.idx)) + ylab("Probability") +
  #theme_minimal() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave(paste0("t",n.taxa,tree.prior,"-node",nod.idx,"-prob.pdf"), p, width = 7, height = 5)

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
ggsave(paste0("t",n.taxa,tree.prior,"-node",nod.idx,"-site",site.idx,"-freq.pdf"), p, width = 5, height = 5)

