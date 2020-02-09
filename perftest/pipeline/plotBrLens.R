# plot branch lengths

WD="~/WorkSpace/codonsubstmodels/perftest/pipeline"
source(file.path(WD, "TraceUtils.R"))


library(ape)

n.taxa = 32
WD=paste0("~/WorkSpace/codonsubstmodels/perftest/T",n.taxa)
setwd(WD)

### true tree
true.txt <- readLines(paste0("t",n.taxa,"coal.txt"))
true.tre <- read.tree(text=true.txt)
plot(true.tre)
nodelabels()

true.tre$edge.length

#WD=paste0("~/WorkSpace/codonsubstmodels/perftest/T",n.taxa)
#setwd(WD)

burnin=0.1
### standard tree likelihood
# read MCMC trees m0.std.trees
tre.log <- file.path(paste0("s4t",n.taxa), "m0.std.trees")
stats.std <- getBrLensStats(tre.log, burnin)
mean.bl.std <- stats.std$mean
last.bl.std <- stats.std$br.lens[nrow(stats.std$br.lens),]
# by columns
ess.std <- apply(stats.std$br.lens, 2, getESS)  

### DA tree likelihood
# read MCMC trees m0.da.trees
tre.log <- file.path(paste0("4t",n.taxa), "m0.da.trees")
stats.da <- getBrLensStats(tre.log, burnin)
mean.bl.da <- stats.da$mean
last.bl.da <- stats.da$br.lens[nrow(stats.da$br.lens),]
ess.da <- apply(stats.da$br.lens, 2, getESS)  

### plots
branches = 2*n.taxa - 2
stopifnot(length(true.tre$edge.length) == branches)
stopifnot(length(mean.bl.std) == branches)
stopifnot(length(mean.bl.da) == branches)

ess = pmin(ess.std, ess.da)
# cap EES to
#ess[ess > 300] = 300
traces <- data.frame(branch=1:branches, truth=true.tre$edge.length, mean.std=mean.bl.std, 
                     mean.da=mean.bl.da, ess = ess, ess.std=ess.std, ess.da=ess.da,
                     last.std=unlist(last.bl.std), last.da=unlist(last.bl.da))

# relative mean to truth
traces.rel = data.frame(branch=1:branches, mean.std=(traces$mean.std-traces$truth), 
                        mean.da=(traces$mean.da-traces$truth), ess = ess, 
                        last.std=(traces$last.std-traces$truth), last.da=(traces$last.da-traces$truth))

# summary
nrow(traces[traces$ess.da < 200,])
nrow(traces[traces$ess.da < 100,])
nrow(traces[traces$ess.std < 200,])
nrow(traces[traces$ess.std < 100,])

# ESS

library(reshape2)
# branch variable value
data.m <- melt(traces.rel, id='branch')
colnames(data.m)[2] <- "simulation"
comps <- length(unique(data.m$simulation))

library(ggplot2)
# relative br lens to truth
p <- ggplot(data.m, aes(branch, value)) + 
  geom_point(aes(colour = simulation, shape= simulation), size=.8, alpha=.7) + 
  scale_shape_manual(values=1:comps) + ylab("relative to truth") +
  ggtitle(paste(n.taxa, "Taxa")) +
  theme_minimal()
ggsave(paste0("t",n.taxa,"-relative-brlens.pdf"), p, width = 20, height = 4)

### fig to compare 2 likelihoods
# std vs da
p1 <- ggplot(traces, aes(mean.std, mean.da)) + 
  geom_point(aes(colour = ess), shape = 1, alpha=.7) + 
  scale_colour_gradientn(colors = c("red", "orange", "blue", "blue")) + # "lightblue", "blue"
  ggtitle(paste("Standard Likelihood vs. Data Augmentation", n.taxa, "Taxa")) +
  xlab("standard likelihood branch lengths") + ylab("data augmentation branch lengths") +
  theme_minimal()
ggsave(paste0("t",n.taxa,"-std-da.pdf"), p1, width = 5, height = 5)

# truth vs std,da
data.m2 <- melt(traces[,c("truth","mean.std","mean.da")], id='truth')
colnames(data.m2)[2] <- "simulation"
p2 <- ggplot(data.m2, aes(truth, value)) + 
  geom_point(aes(colour = simulation, shape= simulation), alpha=.7) + 
  scale_shape_manual(values=3:4) +
  ggtitle(paste("True Tree vs. Estimated Tree", n.taxa, "Taxa")) + 
  xlab("true branch lengths") + ylab("mean branch lengths") +
  theme_minimal()
ggsave(paste0("t",n.taxa,"-truth-mean.pdf"), p2, width = 6, height = 5)

# last
data.m3 <- melt(traces[,c("truth","last.std","last.da")], id='truth')
colnames(data.m3)[2] <- "simulation"
p3 <- ggplot(data.m3, aes(truth, value)) + 
  geom_point(aes(colour = simulation, shape= simulation), alpha=.7) + 
  scale_shape_manual(values=3:4) +
  ggtitle(paste("True Tree vs. Estimated Tree", n.taxa, "Taxa")) + 
  xlab("true branch lengths") + ylab("last branch lengths") +
  theme_minimal()
ggsave(paste0("t",n.taxa,"-truth-last.pdf"), p3, width = 6, height = 5)


