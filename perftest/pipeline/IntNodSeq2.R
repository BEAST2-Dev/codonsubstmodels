

cred.set.check <- list() 
for (cred.thre in c(0.05, 0.25, 0.50, 0.75, 0.95)) {
  #in.cred = c()
  #mean.tot.prob = c()
  n.codon.in.set = list()
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
    
    cred.set <- cred.set %>% mutate(set.n = sapply(cred95.set, length))
    
    n.codon.in.set[[as.character(nod.idx)]] = cred.set[["set.n"]]
    
    # convert to percentage
    #in.per = nrow(cred.set[cred.set[["in.cred"]],]) / nrow(cred.set) * 100
    #in.cred = c(in.cred,  in.per)
    #mean.t.p = mean(cred.set[["cred"]])
    #mean.tot.prob = c(mean.tot.prob, mean.t.p)
  }
  n.codon.in.set = as_tibble(n.codon.in.set)
  stopifnot(ncol(n.codon.in.set) == n.taxa - 1)
  
  thre.txt = gsub("0\\.", "cred", format(cred.thre, nsmall = 2))
  # data frame by threshold
  cred.set.check[[thre.txt]] = n.codon.in.set
}

#nod.idx = 33
nod.idx = 63

cred.set.n = list()
for (cred.thre in c(0.05, 0.25, 0.50, 0.75, 0.95)) {
  thre.txt = gsub("0\\.", "cred", format(cred.thre, nsmall = 2))
  # data frame by threshold
  n.codon.in.set = cred.set.check[[thre.txt]]
  
  cred.set.n[[thre.txt]] = n.codon.in.set[[as.character(nod.idx)]]
}
cred.set.n = as_tibble(cred.set.n) %>% rownames_to_column(var = "site") %>%
  arrange(desc(cred95), desc(cred75), desc(cred50), desc(cred25), desc(cred05))
names(cred.set.n) = gsub("cred", "", names(cred.set.n))
cred.set.n

library(ggplot2)
require(reshape2)

data.m <- melt(cred.set.n, id='site')

p <- ggplot(data.m, aes(reorder(site, -value), value, group = variable)) + 
  geom_line(aes(colour = variable )) +
  labs(color="X% Cred Set") +
  ggtitle(paste("Node", nod.idx)) + 
  ylab("number of codons in X% credible set") + 
  xlab("site index (sorted by codons)") +
  theme_minimal() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave(paste0("node",nod.idx,tree.prior,"-cred-check.pdf"), p, width = 10, height = 5)







