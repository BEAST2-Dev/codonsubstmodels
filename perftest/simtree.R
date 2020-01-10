# simulate tree

library(ape)

WD="~/WorkSpace/codonsubstmodels/perftest"
setwd(WD)

nTaxa = 8
DIR=paste0("T",nTaxa)
if (!dir.exists(DIR)) {
  dir.create(DIR)
}
# read evolver config
template <- readLines("CodonTemplate.dat")

# working dir
setwd(file.path(WD, DIR))

# simulate coalescent tree
tre = rcoal(nTaxa)
#plot(tre)

### check branch lengths, cannot be too short

# save tree
write.tree(tre, file=paste0("t",nTaxa,"coal.txt"))

# plot tree
pdf(paste0("t",nTaxa,".pdf"))
plot(tre)
nodelabels()
tiplabels(adj = 2)
dev.off()

# replace TAXA and tree in the template
tree <- readLines(paste0("t",nTaxa,"coal.txt"))
evolver <- gsub(pattern = "TAXA", replace = toString(nTaxa), x = template)
evolver <- gsub(pattern = "TREE", replace = tree, x = evolver)
writeLines(evolver, con=paste0("t",nTaxa,".evolver.dat"))

# /Applications/paml4.8/bin/evolver 6 t8.evolver.dat > t8.out.txt
CMD="/Applications/paml4.8/bin/evolver 6"
try(system(paste(CMD, paste0("t",nTaxa,".evolver.dat"), ">", paste0("t",nTaxa,".out.txt")), intern = TRUE))

# load XML

# load codons
