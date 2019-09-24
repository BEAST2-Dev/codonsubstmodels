# convert frequencies from the fixed order in PAML: 
# TTT, TTC, TTA, TTG, TCT, TCC, ..., GGG    
# to BEAST 2:
# AAA AAC AAG AAT ... TTA TTC TTG TTT

library(tidyverse)
# freqs in PAML order
frequencies=
"0.00983798  0.01745548  0.00222048  0.01443315
0.00844604  0.01498576  0.00190632  0.01239105
0.01064012  0.01887870  0           0
0.00469486  0.00833007  0           0.00688776
0.01592816  0.02826125  0.00359507  0.02336796
0.01367453  0.02426265  0.00308642  0.02006170
0.01722686  0.03056552  0.00388819  0.02527326
0.00760121  0.01348678  0.00171563  0.01115161
0.01574077  0.02792876  0.00355278  0.02309304
0.01351366  0.02397721  0.00305010  0.01982568
0.01702419  0.03020593  0.00384245  0.02497593
0.00751178  0.01332811  0.00169545  0.01102042
0.02525082  0.04480239  0.00569924  0.03704508
0.02167816  0.03846344  0.00489288  0.03180369
0.02730964  0.04845534  0.00616393  0.04006555
0.01205015  0.02138052  0.00271978  0.01767859"

# split string to a char vector 
freqs = gsub("(\\n+)|(\\r+)|(\\s+)", " ", frequencies) %>% strsplit(" ") %>% unlist

# create triplets in PAML order
paml = c("T","C","A","G")
third = rep(paml, 16)
second = rep(c(rep("T", 4), rep("C", 4),rep("A", 4),rep("G", 4)),4)
first = c(rep("T", 16), rep("C", 16),rep("A", 16),rep("G", 16))
tri = paste0(first, second, third)

# create a 2-column (triplet, freqs) mapping between triplets and freqs 
map.paml = tibble(triplet=tri, freqs=freqs)
map.paml 
# print triplets in PAML order
paste(map.paml$triplet, collapse=", ")

# reorder triplets to BEAST order
map.beast = map.paml %>% arrange(triplet)
# print triplets in BEAST order
paste(map.beast$triplet, collapse=", ")
#"AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG","AGT","ATA","ATC","ATG","ATT",
#"CAA","CAC","CAG","CAT","CCA","CCC","CCG","CCT","CGA","CGC","CGG","CGT","CTA","CTC","CTG","CTT",
#"GAA","GAC","GAG","GAT","GCA","GCC","GCG","GCT","GGA","GGC","GGG","GGT","GTA","GTC","GTG","GTT",
#"TAA","TAC","TAG","TAT","TCA","TCC","TCG","TCT","TGA","TGC","TGG","TGT","TTA","TTC","TTG","TTT"

# copy this freqs in BEAST order to XML
paste(map.beast$freqs, collapse=" ")

###### test
filter(map.paml, triplet %in% c("TCT","ACT","GGA")) %>% arrange(triplet) == 
filter(map.beast, triplet %in% c("TCT","ACT","GGA")) %>% arrange(triplet)


