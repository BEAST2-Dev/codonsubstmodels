# plot branch lengths


getBrLensStats <- function(tre.log="m0.trees", burnin=0.1, plot.1st.tree=FALSE) {
  require(ape)
  tre <- read.nexus(tre.log)
  # rm burnin
  start = as.integer(burnin * length(tre)) + 1
  cat("Remove burnin ", start-1, " trees from the total of ", length(tre), " trees in ", tre.log, "\n")
  tre <- tre[start:length(tre)]
  
  if (plot.1st.tree) {
    plot(tre[[1]])
    nodelabels()
  }
  # combine the list of phylo$edge.length to a data.frame
  br.lens <- as.data.frame( t(sapply(tre, function(tr) rbind(as.numeric(tr$edge.length)))) )
  # mean
  mean.br.lens <- sapply(br.lens, mean) 
  # standard deviation
  sd.br.lens <- sapply(br.lens, sd) 
  
  list(mean=mean.br.lens, sd=sd.br.lens, br.lens=br.lens)
}

library(ape)

nTaxa = 128
WD=paste0("~/WorkSpace/codonsubstmodels/perftest/","T",nTaxa)
setwd(WD)

# true tree
true.txt <- readLines(paste0("t",nTaxa,"coal.txt"))
true.tre <- read.tree(text=true.txt)
plot(true.tre)
nodelabels()

true.tre$edge.length

burnin=0.1
### standard tree likelihood
# read MCMC trees m0.std.trees
tre.log <- file.path(paste0("s4t",nTaxa), "m0.std.trees")
stats.std <- getBrLensStats(tre.log, burnin)
mean.bl.std <- stats.std$mean
last.bl.std <- stats.std$br.lens[nrow(stats.std$br.lens),]

### DA tree likelihood
# read MCMC trees m0.da.trees
tre.log <- file.path(paste0("4t",nTaxa), "m0.da.trees")
stats.da <- getBrLensStats(tre.log, burnin)
mean.bl.da <- stats.da$mean
last.bl.da <- stats.da$br.lens[nrow(stats.da$br.lens),]

### plots
branches = 2*nTaxa - 2
stopifnot(length(true.tre$edge.length) == branches)
stopifnot(length(mean.bl.std) == branches)
stopifnot(length(mean.bl.da) == branches)

traces <- data.frame(branch=1:branches, truth=true.tre$edge.length, mean.std=mean.bl.std, 
                     mean.da=mean.bl.da, last.std=unlist(last.bl.std), last.da=unlist(last.bl.da))

library(reshape2)
traces.rel = data.frame(branch=1:branches, mean.std=(traces$mean.std-traces$truth), 
                        mean.da=(traces$mean.da-traces$truth), last.std=(traces$last.std-traces$truth), 
                        last.da=(traces$last.da-traces$truth))
# branch variable value
data.m <- melt(traces.rel, id='branch')
colnames(data.m)[2] <- "simulation"
comps <- length(unique(data.m$simulation))

library(ggplot2)
# relative br lens to truth
p <- ggplot(data.m, aes(branch, value)) + 
  geom_point(aes(colour = variable, shape= simulation), size=.8, alpha=.7) + 
  scale_shape_manual(values=1:comps) + ylab("relative to truth") +
  ggtitle(paste(nTaxa, "Taxa")) +
  theme_minimal()
ggsave(paste0("t",nTaxa,"-relative-brlens.pdf"), p, width = 20, height = 4)

# std vs da
p1 <- ggplot(traces, aes(mean.std, mean.da)) + 
  geom_point(shape = 1, alpha=.7) + 
  ggtitle(paste("Standard Likelihood vs. Data Augmentation", nTaxa, "Taxa")) +
  theme_minimal()
ggsave(paste0("t",nTaxa,"-std-da.pdf"), p1, width = 5, height = 5)

# truth vs std,da
data.m2 <- melt(traces[,c("truth","mean.std","mean.da")], id='truth')
colnames(data.m2)[2] <- "simulation"
p2 <- ggplot(data.m2, aes(truth, value)) + 
  geom_point(aes(colour = simulation, shape= simulation), alpha=.7) + 
  scale_shape_manual(values=3:4) +
  ggtitle(paste("Truth vs. Simulations", nTaxa, "Taxa")) + ylab("mean branch lengths") +
  theme_minimal()
ggsave(paste0("t",nTaxa,"-truth-mean.pdf"), p2, width = 6, height = 5)

# last
data.m3 <- melt(traces[,c("truth","last.std","last.da")], id='truth')
colnames(data.m3)[2] <- "simulation"
p3 <- ggplot(data.m3, aes(truth, value)) + 
  geom_point(aes(colour = simulation, shape= simulation), alpha=.7) + 
  scale_shape_manual(values=3:4) +
  ggtitle(paste("Truth vs. Simulations", nTaxa, "Taxa")) + ylab("last branch lengths") +
  theme_minimal()
ggsave(paste0("t",nTaxa,"-truth-last.pdf"), p3, width = 6, height = 5)


