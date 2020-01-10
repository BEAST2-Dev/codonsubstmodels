# create XML

WD="~/WorkSpace/codonsubstmodels/perftest"

nTaxa = 8
# working dir
setwd(file.path(WD, DIR))

# load sequences to a 2-column tibble
nex <- readNex("mc.nex", "t\\d+")
# create data
data <- toXMLData(nex, "data.xml", id="alignment")

# load XML template
library(xml2)
template <- read_xml("../M0DALikelihood.xml")


