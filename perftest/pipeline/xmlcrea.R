# pipeline step 2
# create XML, run simtree before this

WD="~/WorkSpace/codonsubstmodels/perftest/pipeline"
source(file.path(WD, "nexparser.R"))

# c(8, 64, 512)
n.taxa = 32
if (n.taxa == 32) {
  chain.len.da =  "40000000" 
  chain.len.std = "12000000" 
} else if (n.taxa == 64) {
  chain.len.da =  "80000000" 
  chain.len.std = "14000000" 
} else {
  chain.len.da =  "120000000" # 128T 120000000 
  chain.len.std = "15000000" # 128T 15000000 
}
DIR=paste0("T",n.taxa)
# working dir
WD="~/WorkSpace/codonsubstmodels/perftest/"
setwd(file.path(WD, DIR))

# coal or yulelam10
tree.prior = "yulelam10"

###### DA tree likelihood
library(xml2)
# load XML template from ~/WorkSpace/codonsubstmodels/perftest/
if (tree.prior == "coal") {
  template <- read_xml(file.path(WD, "M0DACoalescent.xml"))
} else {
  template <- read_xml(file.path(WD, "M0DAYule.xml"))
}


# load sequences to a 2-column tibble
setwd(file.path(WD, DIR, tree.prior))
cat("Load mc.nex from ", getwd(), "\n")
nex <- readNex("mc.nex", "t\\d+")
# create <data> and save to xml
DATA <- toXMLData(nex, "data.xml", id="alignment")
# back to working dir
setwd(file.path(WD, DIR))

# 1. replace data
node.data <- xml_find_first(template, ".//data")
xml_replace(node.data, DATA)

# 2. replace PI: equal, F1X4, F3X4, F60/F61
PI = "F60/F61"
nodes<-xml_find_all(template, ".//frequencies")
node<-nodes[xml_has_attr(nodes, "pi")]
xml_attr(node, "pi") <- PI

# 3. replace tree
tre.txt <- paste0("t",n.taxa,tree.prior,".txt") # t32yulelam10.txt
TREE <- readLines(tre.txt)
require(ape)
start.tree <- read.tree(text = TREE)
## Note: the tree has to be time tree, be careful when changing branch lengths, instead of node heights. 
# check all branch lengths
stopifnot(all(start.tree$edge.length > 1e-6))
# all branch lengths times 10
start.tree$edge.length <- start.tree$edge.length * 10
# starting tree 
TREE <- write.tree(start.tree)

nodes<-xml_find_all(template, ".//tree")
node<-nodes[xml_has_attr(nodes, "newick")]
xml_attr(node, "newick") <- TREE

# 4. replace MCMC config
nodes<-xml_find_all(template, ".//run")
node<-nodes[xml_has_attr(nodes, "chainLength")]
xml_attr(node, "chainLength") <- chain.len.da

# 5. replace thread
THREAD = 4
nodes<-xml_find_all(template, ".//distribution")
node<-nodes[xml_has_attr(nodes, "threads")]
xml_attr(node, "threads") <- THREAD

# finish XML
write_xml(template, file = paste0("t", n.taxa, tree.prior,"DA.xml"))

###### standard tree likelihood

# load XML template
if (tree.prior == "coal") {
  template <- read_xml(file.path(WD, "M0StandardCoalescent.xml"))
} else {
  template <- read_xml(file.path(WD, "M0StandardYule.xml"))
}


# 1. replace data
node.data <- xml_find_first(template, ".//data")
# use the same from DA
xml_replace(node.data, DATA)

# 2. replace PI: equal, F1X4, F3X4, F60/F61
nodes<-xml_find_all(template, ".//frequencies")
node<-nodes[xml_has_attr(nodes, "pi")]
# use the same from DA
xml_attr(node, "pi") <- PI

# 3. replace tree
nodes<-xml_find_all(template, ".//tree")
node<-nodes[xml_has_attr(nodes, "newick")]
# use the same from DA
xml_attr(node, "newick") <- TREE

# 4. replace MCMC config
nodes<-xml_find_all(template, ".//run")
node<-nodes[xml_has_attr(nodes, "chainLength")]
xml_attr(node, "chainLength") <- chain.len.std

# finish XML
write_xml(template, file = paste0("t", n.taxa, tree.prior,"STD.xml"))

### 
# TREE <- str_replace_all(TREE, ":(\\d+).(\\d+)", ":1.0")
## find all tip branch's indexes 
#tip.edges <- which(start.tree$edge[,2] <= n.taxa)
## and tips height plus 0.1
#start.tree$edge.length[tip.edges] = start.tree$edge.length[tip.edges] + 0.1
###