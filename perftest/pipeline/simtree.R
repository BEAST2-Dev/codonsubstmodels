# simulate tree

library(ape)

WD="~/WorkSpace/codonsubstmodels/perftest"
setwd(WD)

# c(32, 128, 512)
n.taxa = 32
DIR=paste0("T",n.taxa)
if (!dir.exists(DIR)) {
  dir.create(DIR)
}
# read evolver config
template <- readLines("CodonTemplate.dat")

# working dir
setwd(file.path(WD, DIR))

# simulate coalescent tree
tre = rcoal(n.taxa)
#plot(tre)

### check branch lengths, cannot be too short

# save tree
write.tree(tre, file=paste0("t",n.taxa,"coal.txt"))

# plot tree
pdf(paste0("t",n.taxa,".pdf"))
plot(tre)
nodelabels()
tiplabels(adj = 2)
dev.off()

# replace TAXA and tree in the template
tree <- readLines(paste0("t",n.taxa,"coal.txt"))
evolver <- gsub(pattern = "TAXA", replace = toString(n.taxa), x = template)
evolver <- gsub(pattern = "TREE", replace = tree, x = evolver)
writeLines(evolver, con=paste0("t",n.taxa,".evolver.dat"))

# /Applications/paml4.8/bin/evolver 6 t8.evolver.dat > t8.out.txt
CMD="/Applications/paml4.8/bin/evolver 6"
try(system(paste(CMD, paste0("t",n.taxa,".evolver.dat"), ">", paste0("t",n.taxa,".out.txt")), intern = TRUE))

# load XML

# load codons
