# evolver

### evolver.out
# Note: evolver tree uses the same format as APE tree
#evo.out <- readLines(paste0("t",n.taxa,".out.txt"))

# return a list of true sequences at internal nodes simulated by evolver.
# list names are the internal node index.
getSeqsEvoOut <- function(out.file="ancestral.txt", n.taxa=NULL, 
                          nod.name="node[0-9]+", genetic.code="universal") {
  evo.anc <- readLines(out.file)
  anc.tru <- evo.anc[grepl(paste0("^",nod.name), evo.anc)] # "^node[0-9]+"
  if (!is.null(n.taxa)) stopifnot(length(anc.tru) == n.taxa-1)
  
  nod.names <- gsub(paste0("^(",nod.name,").*"), "\\1", anc.tru) # "^(node[0-9]+).*"
  seqs <- gsub(paste0("^(",nod.name,")\\s+(.*)"), "\\2", anc.tru) # "^(node[0-9]+)\\s+(.*)"
  # rm spaces
  seqs <- gsub("\\s+", "", seqs)
  stopifnot(str_length(seqs[1]) == n.codon*3)
  
  states <- lapply(seqs, seqToStates, genetic.code=genetic.code)
  # internal nodes index
  nod.idx <- gsub("^node", "", nod.names)
  names(states) <- as.character(nod.idx)
  
  return(states)
}
