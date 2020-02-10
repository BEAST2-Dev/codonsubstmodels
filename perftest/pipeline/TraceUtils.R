# MCMC post-analysis

# samp.freq is the step size, namely sampling frequency of the samples
getESS <- function(samples, verbose=FALSE) {
  n.samples = length(samples);
  max.lag = min(n.samples, 2000); # Note: lag index starts from 1 NOT 0
  mean.s = mean(samples)
  
  if (verbose) cat("Input", n.samples, "samples, ")
  
  lagged.square.sum=c()
  for (lag in 1:max.lag) {
    lagged.square.sum[lag] = 0
    for (s in 1:(n.samples - lag + 1)) {
      del1 = samples[s] - mean.s;
      # Note: both samples and lag index start from 1 NOT 0
      del2 = samples[s + lag - 1] - mean.s; 
      lagged.square.sum[lag] = lagged.square.sum[lag] + (del1 * del2);
    } # end s loop
    
    lagged.square.sum[lag] = lagged.square.sum[lag] / (n.samples - lag);
    
    if (lag == 1) {
      var.stats = lagged.square.sum[1];
    } else if (lag %% 2 != 0) { # Note: lag index starts from 1 NOT 0
      # fancy stopping criterion :)
      if (lagged.square.sum[lag - 1] + lagged.square.sum[lag] > 0) {
        var.stats = var.stats + 2.0 * (lagged.square.sum[lag - 1] + lagged.square.sum[lag]); 
      } else { # stop for loop
        lagged.square.sum = lagged.square.sum[1:lag]
        break
      }
    }
  }
  
  # auto correlation time
  # if (lagged.square.sum[1] == 0) {
  #   ACT = 0;
  # } else {
  #   ACT = samp.freq * var.stats / lagged.square.sum[1];
  # }
  
  # effective sample size
  if (lagged.square.sum[1] == 0) {
    ESS = 1;
  } else {
    #ESS = (samp.freq * n.samples) / ACT;
    ESS = lagged.square.sum[1] * n.samples / var.stats;
  }
  
  if (verbose) cat("ESS is ", ESS, ".\n")
  
  #list(ESS=ESS, var.stats=var.stats, lagged.square.sum=lagged.square.sum) #ACT=ACT, 
  return(ESS)
}

# Branch lengths stats given a tre.log containing MCMC trees
getBrLensStats <- function(tre.log="m0.trees", burnin=0.1, plot.1st.tree=FALSE) {
  require(ape)
  tre <- read.nexus(tre.log)
  # rm burnin
  start = as.integer(burnin * length(tre)) + 1
  cat("Remove burnin ", start-1, " trees from the total of ", length(tre), " trees in ", tre.log, "\n")
  tre <- tre[start:length(tre)]
  
  if (plot.1st.tree) {
    plot(tre[[1]])
    nodelabels()
  }
  # combine the list of phylo$edge.length to a data.frame
  br.lens <- as.data.frame( t(sapply(tre, function(tr) rbind(as.numeric(tr$edge.length)))) )
  # mean
  mean.br.lens <- sapply(br.lens, mean) 
  # standard deviation
  sd.br.lens <- sapply(br.lens, sd) 
  
  list(mean=mean.br.lens, sd=sd.br.lens, br.lens=br.lens)
}


getIntNodeSeqStats <- function(ins.log="ins.txt", burnin=0.1, 
                               col.names=c("Sample","Node","States")) {
  require(tidyverse)
  # Sample  Node States  
  traces <- read_delim(ins.log, "\t", comment = "#", col_types = cols( States = col_character() ))
  if (! any(colnames(traces) %in% col.names) )
    stop("Incorrect file format in the log file of internal node sequences !\n", 
         "Column names = ", paste(colnames(traces), collapse = ","))
  
  # MCMC summary
  samples = unique(traces[[col.names[1]]]) # "Sample"
  samples = samples[!is.na(samples)]
  nodes = unique(traces[[col.names[2]]]) # "Node"
  nodes = nodes[!is.na(nodes)]
  # state is a 2-digit integer from 00 to 59/60, n.codon = str.len / 2 
  str.len = str_length(traces[[col.names[3]]][1]) # "States"
  n.codon = str.len / 2 
  
  cat("Chain length", prettyNum(samples[length(samples)], big.mark=",",scientific=FALSE), 
      ", log every", prettyNum(samples[2], big.mark=",",scientific=FALSE), "samples, each sample includes",
      length(nodes), "internal nodes [", nodes[1], "-", nodes[length(nodes)], "], ", n.codon, " codons.\n")
  
  stats.list <- list()
  for (node.id in nodes) {
    # internal node index starts from n.taxa, node.id = 33
    node.samples <- traces %>% filter(Node==node.id)
    # rm burnin, +2 to exclude state 0 
    start = as.integer(burnin * nrow(node.samples)) + 2
    cat("Remove burnin ", start-1, " sampled internal node sequences from the total of ", 
        nrow(node.samples), "for node ", node.id, "\n")
    # paste0("c", 1:1000), seq(2, 2000, 2)
    node.samples <- node.samples[start:nrow(node.samples),] %>% 
      separate(States, into = paste0("c", 1:n.codon), sep = seq(2, str.len, 2))
    
    ### posterior distribution
    stats <- NULL
    for (c in 1:n.codon) {
      freq <- table(node.samples[[paste0("c", c)]])
      freq <- sort(freq, decreasing=T)
      s1 <- tibble(state=names(freq), freq=freq, site=c, order=1:length(freq))
      stats <- bind_rows(stats, s1)
    }
    stats <- stats %>% mutate(order = as.character(order))
    
    stats.list[[as.character(node.id)]] <- stats
  }
  
  return(stats.list)
}









