# ESS per hour

WD="~/WorkSpace/codonsubstmodels/perftest/pipeline"
source(file.path(WD, "TraceUtil.R"))

# branch lengths ESS per hour 
getBrLenESSHour <- function(tre.log="m0.trees", screen.log="out.txt", burnin=0.1) {
  require(ape)
  # read MCMC trees m0.std.trees
  stats <- getBrLensStats(tre.log, burnin)
  mean.bl <- stats$mean
  last.bl <- stats$br.lens[nrow(stats$br.lens),]
  # by columns
  ess <- apply(stats$br.lens, 2, getESS)  
  
  # read screen log
  screen.log <- file.path(paste0("s4t",n.taxa), paste0("t",n.taxa,"st.xml_out.txt"))
  screen.info <- scan(screen.log,sep="\n",what="char(0)",skip=110)
  # Total calculation time: 71818.323 seconds
  time <- screen.info[length(screen.info)]
  time <- gsub("^.*time: (.*) seconds.*", "\\1", time)
  time <- as.numeric(time) / 3600 # ESS/hour
  
  ess.per.hour = ess / time
}


n.taxa = 32
WD=paste0("~/WorkSpace/codonsubstmodels/perftest/T",n.taxa)
setwd(WD)

burnin=0.1
### standard tree likelihood
# read MCMC trees m0.std.trees
tre.log <- file.path(paste0("s4t",n.taxa), "m0.std.trees")
# read screen log
screen.log <- file.path(paste0("s4t",n.taxa), paste0("t",n.taxa,"st.xml_out.txt"))
# ESS per hour for all branches
ess.hour.std <- getBrLenESSHour(tre.log, screen.log, burnin)


tre.log <- file.path(paste0("4t",n.taxa), "m0.da.trees")
screen.log <- file.path(paste0("4t",n.taxa), paste0("t",n.taxa,"th4.xml_out.txt"))
ess.hour.da<- getBrLenESSHour(tre.log, screen.log, burnin)

### plot
branches = 2*n.taxa - 2
stopifnot(length(ess.hour.std) == branches)
stopifnot(length(ess.hour.da) == branches)
ess.hour <- data.frame(branch=1:branches, std=ess.hour.std, da=ess.hour.da)

data.m <- melt(ess.hour, id='branch')
colnames(data.m)[2] <- "method"


library(ggplot2)
# relative br lens to truth
p <- ggplot(data.m, aes(method, value)) + 
  geom_boxplot(aes(colour = method)) +
  ggtitle(paste(n.taxa, "Taxa", branches, "Branches")) + ylab("ESS / hour") +
  theme_minimal()
ggsave(paste0("t",n.taxa,"-ess-hour.pdf"), p, width = 5, height = 5)








