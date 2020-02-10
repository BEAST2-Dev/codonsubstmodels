# imported from BEAST
# codons to states, states to codons

require(tidyverse)

# Nucleotides go A, C, G, T - Note: this is not the order used by the Genbank web site
# With the first codon position most significant (i.e. AAA, AAC, AAG, AAT, ACA, etc.).
genetic.code.table = tibble(name=c(
  "universal", "vertebrateMitochondrial", "yeast", "moldProtozoanMitochondrial",
  "mycoplasma", "invertebrateMitochondrial", "ciliate", "echinodermMitochondrial",
  "euplotidNuclear", "bacterial", "alternativeYeast", "ascidianMitochondrial",
  "flatwormMitochondrial", "blepharismaNuclear", "noStops" ),
  code=c(
  # Universal
  "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF",
  # Vertebrate Mitochondrial
  "KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
  # Yeast
  "KNKNTTTTRSRSMIMIQHQHPPPPRRRRTTTTEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
  # Mold Protozoan Mitochondrial
  "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
  # Mycoplasma
  "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
  # Invertebrate Mitochondrial
  "KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
  # Ciliate
  "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF",
  # Echinoderm Mitochondrial
  "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
  # Euplotid Nuclear
  "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF",
  # Bacterial
  "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF",
  # Alternative Yeast
  "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF",
  # Ascidian Mitochondrial
  "KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF",
  # Flatworm Mitochondrial
  "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF",
  # Blepharisma Nuclear
  "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF",
  # No stops
  "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYYQYSSSSWCWCLFLF"
))

# print first n rows of genetic.code.table
printGeneticCodeTable <- function(n = 5) {
  print(genetic.code.table, n = n)
}

# convert one string of triplets into a vector of states as integers
seqToStates <- function(seq.str, genetic.code="universal", 
                        # all 64 triplets
                        triplets = c(
  "AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT",
  "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", # 16
  "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT",
  "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", # 32
  "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT",
  "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", # 48
  "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT",
  "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT" # 64 here
  ) , verbose=FALSE) {
  
  stopifnot(genetic.code %in% genetic.code.table[["name"]]) 
  code = genetic.code.table %>% filter(name == genetic.code) %>% select(code) %>% unlist
  # vertebrateMitochondrial 9 11 49 51
  stop.codon.idx <- str_locate_all(code, "\\*")[[1]][,1]
  stopifnot(genetic.code != "noStops" && length(stop.codon.idx) > 0)
  # rm stop codons
  triplets <- triplets[-stop.codon.idx]
  
  # vaildate seq.str
  if (! (is.character(seq.str) & length(seq.str) == 1) )
    stop("Input seq.str must be one string !")
  stopifnot(str_length(seq.str) %% 3 == 0)
  n.codon <- (str_length(seq.str) / 3)
  # seq.str must a string
  codons <- strsplit(toupper(seq.str), "(?<=.{3})", perl = TRUE)[[1]]
  # vaildate codons
  stopifnot(length(codons) == n.codon)
  
  n.state = length(triplets)
  if (verbose) 
    cat("Convert sequences to", n.codon, "codons with", n.state, "states using", genetic.code, "genetic code.\n")
  
  # map to states, 0-59/60
  states = match(codons, triplets) - 1
  # from 0 to n.state-1
  if ( any(states < 0 | states >= n.state) )
    stop("Incorrect states : ", paste(states[states < 0 | states >= n.state], collapse = ", "))
  return(states)
}


