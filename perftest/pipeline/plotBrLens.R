# plot branch lengths

WD="~/WorkSpace/codonsubstmodels/perftest/pipeline"
source(file.path(WD, "TraceUtils.R"))
source(file.path(WD, "TreeUtils.R"))

library(tidyverse)
library(ape)

n.taxa = 32
WD=paste0("~/WorkSpace/codonsubstmodels/perftest/T",n.taxa)
setwd(WD)

# coal or yulelam10
tree.prior = "yulelam10"

### true tree
true.txt <- readLines(paste0("t",n.taxa,tree.prior,".txt"))
true.tre <- read.tree(text=true.txt)
plot(true.tre)
nodelabels()

true.tre$edge.length

#WD=paste0("~/WorkSpace/codonsubstmodels/perftest/T",n.taxa)
#setwd(WD)

burnin=0.1
### standard tree likelihood
# read MCMC trees m0.std.trees
tre.log <- file.path(paste0("t",n.taxa,tree.prior,"STD"),"m0.std.trees")
stats.std <- getBrLensStats(tre.log, burnin)
mean.bl.std <- stats.std$mean
last.bl.std <- stats.std$br.lens[nrow(stats.std$br.lens),]
# by columns
ess.std <- apply(stats.std$br.lens, 2, getESS) 
# total tree length
total.tree.len <- rowSums(stats.std$br.lens)
ess.ttl.std <- getESS(total.tree.len)

### DA tree likelihood
# read MCMC trees m0.da.trees
tre.log <- file.path(paste0("t",n.taxa,tree.prior,"DA"), "m0.da.trees")
stats.da <- getBrLensStats(tre.log, burnin)
mean.bl.da <- stats.da$mean
last.bl.da <- stats.da$br.lens[nrow(stats.da$br.lens),]
ess.da <- apply(stats.da$br.lens, 2, getESS)  
# total tree length
total.tree.len <- rowSums(stats.da$br.lens)
ess.ttl.da <- getESS(total.tree.len)

### plots
branches = 2*n.taxa - 2
stopifnot(length(true.tre$edge.length) == branches)
stopifnot(length(mean.bl.std) == branches)
stopifnot(length(mean.bl.da) == branches)

ess = pmin(ess.std, ess.da)
# cap EES to
#ess[ess > 300] = 300

### create edges.map to make sure the same node index is the same node
# true.tre$edge is true tree, 
true.edges <- as_tibble(true.tre$edge) 
names(true.edges) <- c("parent","child")  
# stats.std$edges is standard likelihood,
# Note: tree is fixed in MCMC
std.edges <- as_tibble(stats.std$edges.list[[1]])
names(std.edges) <- c("parent","child") 
edges.map.std <- mapEdges(true.edges, std.edges)
#print(edges.map.std, n=Inf)
# stats.da$edges is DA
da.edges <- as_tibble(stats.da$edges.list[[1]])
names(da.edges) <- c("parent","child") 
edges.map.da <- mapEdges(true.edges, da.edges)
#print(edges.map.da, n=Inf)

# the indexing at 3 set of trees is same 
all(edges.map.std$ch1 == edges.map.std$ch2)
all(edges.map.da$ch1 == edges.map.da$ch2)

# mean.std the mean of branch lengths sampled from MCMC using BEAST standard tree likelihood,
# mean.da the mean of branch lengths sampled from MCMC using DA
traces <- tibble(branch=1:branches, truth=true.tre$edge.length, mean.std=mean.bl.std, 
                     mean.da=mean.bl.da, ess = ess, ess.std=ess.std, ess.da=ess.da,
                     last.std=unlist(last.bl.std), last.da=unlist(last.bl.da))

# relative mean to truth
#traces.rel = data.frame(branch=1:branches, mean.std=(traces$mean.std-traces$truth), 
#                        mean.da=(traces$mean.da-traces$truth), ess = ess, 
#                        last.std=(traces$last.std-traces$truth), last.da=(traces$last.da-traces$truth))

# summary
#nrow(traces[traces$ess.da < 200,])
#nrow(traces[traces$ess.da < 100,])
#nrow(traces[traces$ess.std < 200,])
#nrow(traces[traces$ess.std < 100,])

# all ESS, branch=0 is total tree length
all.ess <- traces %>% select(branch,ess.std,ess.da) %>% 
  add_row(branch=0,ess.std=ess.ttl.std,ess.da=ess.ttl.da)

library(reshape2)
# branch variable value
data.m <- melt(all.ess, id='branch')
colnames(data.m)[2] <- "simulation"
#comps <- length(unique(data.m$simulation))
data.m[["simulation"]] <- gsub("ess.", "", data.m[["simulation"]])

n.branches = 2*n.taxa - 2

library(ggplot2)
### all ess, branch=0 is total tree length
p <- ggplot(data.m, aes(branch, value)) + 
  geom_point(aes(colour = simulation, shape= simulation), size=.8, alpha=.7) + 
#  geom_hline(yintercept=200, linetype="dashed", color = "red") +
  scale_y_continuous(trans='log10') +
#  annotate("text", x = 0, y = 210, label = "200") +
  ggtitle(paste(n.taxa, "Taxa", n.branches, "Branches")) + 
  ylab("ESS") + xlab("branch index (0 for total tree length)") +
  theme_minimal()
ggsave(paste0("t",n.taxa,tree.prior,"-brlens-ess.pdf"), p, width = 6, height = 5)

### fig to compare 2 likelihoods
# std vs da
min.ess = round(min(traces[["ess"]]) / 10) * 10
stopifnot(min.ess > 200)
p1 <- ggplot(traces, aes(mean.std, mean.da)) +
  geom_point(aes(colour = ess), shape = 1, alpha=.7) +
  scale_colour_gradientn(colors = c("lightblue", "blue")) + # "lightblue", "blue"
  ggtitle(paste(n.taxa, "Taxa", n.branches, "Branches, ESS >", min.ess)) + #"Standard Likelihood vs. Data Augmentation",
  xlab("branch lengths - standard likelihood") + ylab("branch lengths - data augmentation") +
  theme_minimal()
ggsave(paste0("t",n.taxa,tree.prior,"-std-da.pdf"), p1, width = 5, height = 5)


# true br len vs other 2 meothods
traces2 <- traces %>% mutate(std.mean.err = (mean.std-truth)/truth, da.mean.err = (mean.da-truth)/truth) 
print(traces2, n=Inf)
data.m1 <- traces2 %>% select(branch,std.mean.err,da.mean.err) %>% melt(id='branch')
colnames(data.m1)[2] <- "simulation"
data.m1[["simulation"]] <- gsub("\\.mean\\.err", "", data.m1[["simulation"]])
# TODO use 95% hpd
# p1 <- ggplot(data.m1, aes(branch, value, fill = simulation)) + 
#   geom_bar(position = "dodge", stat="identity") + 
#   ggtitle(paste("True Tree vs. Sampled Tree", n.taxa, "Taxa", n.branches, "Branches")) + 
#   ylab("branch lengths") +
#   theme_minimal()
# ggsave(paste0("t",n.taxa,tree.prior,"-truth-mean-err.pdf"), p1, width = 7, height = 5)


# truth vs std,da
data.m2 <- melt(traces[,c("truth","mean.da","mean.std")], id='truth')
colnames(data.m2)[2] <- "simulation"
data.m2[["simulation"]] <- gsub("mean.", "", data.m2[["simulation"]])

p2 <- ggplot(data.m2, aes(truth, value)) + 
  geom_point(aes(colour = simulation, shape= simulation), alpha=.7) + 
  scale_shape_manual(values=3:4) +
  ggtitle(paste("True Tree vs. Sampled Tree at", n.taxa, "Taxa", n.branches, "Branches")) + 
  xlab("true branch lengths") + ylab("mean of sampled branch lengths") +
  theme_minimal()
ggsave(paste0("t",n.taxa,tree.prior,"-truth-mean.pdf"), p2, width = 6, height = 5)

# last
data.m3 <- melt(traces[,c("truth","last.std","last.da")], id='truth')
colnames(data.m3)[2] <- "simulation"
p3 <- ggplot(data.m3, aes(truth, value)) + 
  geom_point(aes(colour = simulation, shape= simulation), alpha=.7) + 
  scale_shape_manual(values=3:4) +
  ggtitle(paste("True Tree vs. Estimated Tree", n.taxa, "Taxa", n.branches, "Branches")) + 
  xlab("true branch lengths") + ylab("last branch lengths") +
  theme_minimal()
ggsave(paste0("t",n.taxa,tree.prior,"-truth-last.pdf"), p3, width = 6, height = 5)


