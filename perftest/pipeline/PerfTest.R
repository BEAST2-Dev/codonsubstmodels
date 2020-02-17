# ESS per hour

WD="~/WorkSpace/codonsubstmodels/perftest/pipeline"
source(file.path(WD, "TraceUtils.R"))

n.taxa = 32
WD=paste0("~/WorkSpace/codonsubstmodels/perftest/T",n.taxa)
setwd(WD)

# coal or yulelam10
tree.prior = "yulelam10"

burnin=0.1
### standard tree likelihood
# read MCMC trees m0.std.trees
tre.log <- file.path(paste0("t",n.taxa,tree.prior,"STD"),"m0.std.trees")
# read screen log
screen.log <- file.path(paste0("t",n.taxa,tree.prior,"STD"), paste0("t",n.taxa,tree.prior,"STD.xml_out.txt"))
# ESS per hour for all branches
ess.std <- getBrLenESSHour(tre.log, screen.log, burnin)
ess.hour.std <- ess.std[["ess.per.hour.branch.lens"]]
ess.hour.totlen.std <- ess.std[["ess.per.hour.tot.lens"]]

tre.log <- file.path(paste0("t",n.taxa,tree.prior,"DA"), "m0.da.trees")
screen.log <- file.path(paste0("t",n.taxa,tree.prior,"DA"), paste0("t",n.taxa,tree.prior,"DA.xml_out.txt"))
ess.da <- getBrLenESSHour(tre.log, screen.log, burnin)
ess.hour.da <- ess.da[["ess.per.hour.branch.lens"]]
ess.hour.totlen.da <- ess.da[["ess.per.hour.tot.lens"]]


### plot
n.branches = 2*n.taxa - 2
stopifnot(length(ess.hour.std) == n.branches)
stopifnot(length(ess.hour.da) == n.branches)
ess.hour <- data.frame(branch=1:n.branches, std=ess.hour.std, da=ess.hour.da)

require(reshape2)
data.m <- melt(ess.hour, id='branch')
colnames(data.m)[2] <- "method"

library(ggplot2)
# relative br lens to truth
p <- ggplot(data.m, aes(method, value)) + 
  geom_boxplot(aes(colour = method)) +
  ggtitle(paste(n.taxa, "Taxa", n.branches, "Branches")) + ylab("ESS / hour") +
  theme_minimal()
ggsave(paste0("t",n.taxa,tree.prior,"-ess-hour.pdf"), p, width = 5, height = 5)








