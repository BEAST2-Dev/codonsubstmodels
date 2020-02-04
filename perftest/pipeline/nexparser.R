# read a paup format (mc.nex) output from evolver

# grep sequences given a pattern of taxa names,
# e.g. "t\\d+" for t1, t2, ...
# restrict to single alignment
# return a 2-column tibble
readNex <- function(file="mc.nex", taxa.pattern="t\\d+") {
  require(tidyverse)
  require(stringr)
  nex <- readLines(file)
  # allow spaces in the pattern
  lins.pattern = paste0("^(\\s*)(", taxa.pattern, ")(\\s*)")
  # grep sequences by a pattern of taxa names.
  lins <- str_subset(nex, lins.pattern)

  # parse to 2 cols
  taxa <- str_extract(lins, taxa.pattern)
  seqs <- str_replace(lins, lins.pattern, "")
  # rm all spaces
  seqs <- str_replace_all(seqs, "\\s", "")

  cat("Load", length(seqs), "sequences having length", str_length(seqs[1]), "from", file, "\n")

  return(tibble(taxa=taxa,seqs=seqs))
}

# convert 2-column data.frame from readNex() to a BEAST XML data section.
# 1st column is taxon names, 2nd is the sequences.
# if file is NULL, then only no XML file output.
toXMLData <- function(nex=tibble(), file="data.xml", id="alignment") {
  require(tidyverse)
  require(xml2)

  doc <- xml_new_root("data")
  xml_set_attrs(doc, c("id" = id, "dataType" = "nucleotide"))

  # add sequences
  for (l in 1:nrow(nex)) {
    xml_add_child(doc, "sequence", taxon=nex[[1]][l], nex[[2]][l])
  }

  # options = 'as_xml' ruins line break
  if (!is.null(file)) write_xml(doc, file = file)

  return(doc)
}


