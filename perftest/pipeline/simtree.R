# simulate tree

library(ape)

WD="~/WorkSpace/codonsubstmodels/perftest"
setwd(WD)

# read evolver config
template <- readLines("CodonTemplate.dat")

# c(32, 128, 512)
n.taxa = 32
DIR=paste0("T",n.taxa)
if (!dir.exists(DIR)) {
  dir.create(DIR)
}
# working dir
setwd(file.path(WD, DIR))

tree.prior = "yulelamda10"
tree.txt = ""
if (tree.prior == "coal") {
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
  tree.txt <- readLines(paste0("t",n.taxa,tree.prior,".txt"))
  
} else {
  # generate from MASTER YuleTree.xml
  tree.txt <- readLines(paste0("t",n.taxa,tree.prior,".txt"))
  writeLines(tree.txt, con=paste0("t",n.taxa,tree.prior,"-bak.txt"))
  tree.txt
  # add prefix to taxon name
  tree.txt <- gsub("([0-9]+:)", "t\\1", x = tree.txt)
  writeLines(tree.txt, con=paste0("t",n.taxa,tree.prior,".txt"))
}

evolver <- gsub(pattern = "TAXA", replace = toString(n.taxa), x = template)
evolver <- gsub(pattern = "TREE", replace = tree.txt, x = evolver)
writeLines(evolver, con=paste0("t",n.taxa,".evolver.dat"))

# /Applications/paml4.8/bin/evolver 6 t8.evolver.dat > t8.out.txt
CMD="/Applications/paml4.8/bin/evolver 6"
try(system(paste(CMD, paste0("t",n.taxa,".evolver.dat"), ">", paste0("t",n.taxa,".out.txt")), intern = TRUE))

# load XML

# load codons
