# create XML, run simtree before this

WD="~/WorkSpace/codonsubstmodels/perftest"
source(file.path(WD, "nexparser.R"))

# c(8, 64, 512)
nTaxa = 128
DIR=paste0("T",nTaxa)
# working dir
setwd(file.path(WD, DIR))

###### DA tree likelihood
library(xml2)
# load XML template
template <- read_xml("../M0DALikelihood.xml")

# load sequences to a 2-column tibble
nex <- readNex("mc.nex", "t\\d+")
# create data
DATA <- toXMLData(nex, "data.xml", id="alignment")
# 1. replace data
node.data <- xml_find_first(template, ".//data")
xml_replace(node.data, DATA)

# 2. replace PI: equal, F1X4, F3X4, F60/F61
PI = "F60/F61"
nodes<-xml_find_all(template, ".//frequencies")
node<-nodes[xml_has_attr(nodes, "pi")]
xml_attr(node, "pi") <- PI

# 3. replace tree
TREE <- readLines(paste0("t",nTaxa,"coal.txt"))
require(ape)
start.tree <- read.tree(text = TREE);
## Note: the tree has to be time tree, be careful when changing branch lengths, instead of node heights.  
# all branch lengths times 10
start.tree$edge.length <- start.tree$edge.length * 10
# starting tree 
TREE <- write.tree(start.tree)

nodes<-xml_find_all(template, ".//tree")
node<-nodes[xml_has_attr(nodes, "newick")]
xml_attr(node, "newick") <- TREE

# 4. replace MCMC config
# TODO

# 5. replace thread
THREAD = 4
nodes<-xml_find_all(template, ".//distribution")
node<-nodes[xml_has_attr(nodes, "threads")]
xml_attr(node, "threads") <- THREAD

# finish XML
write_xml(template, file = paste0("t", nTaxa, "th", THREAD, ".xml"))

###### standard tree likelihood

# load XML template
template <- read_xml("../M0StandardLikelihood.xml")

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
# TODO

# finish XML
write_xml(template, file = paste0("t", nTaxa, "st.xml"))

### 
# TREE <- str_replace_all(TREE, ":(\\d+).(\\d+)", ":1.0")
## find all tip branch's indexes 
#tip.edges <- which(start.tree$edge[,2] <= nTaxa)
## and tips height plus 0.1
#start.tree$edge.length[tip.edges] = start.tree$edge.length[tip.edges] + 0.1
###