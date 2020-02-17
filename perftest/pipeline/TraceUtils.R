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
  trees <- read.nexus(tre.log)
  # rm burnin
  start = as.integer(burnin * length(trees)) + 1
  cat("Remove burnin ", start-1, " trees from the total of ", length(trees), " trees in ", tre.log, "\n")
  trees <- trees[start:length(trees)]
  
  if (plot.1st.tree) {
    plot(trees[[1]])
    nodelabels()
  }
  # combine the list of phylo$edge.length to a data.frame
  br.lens <- as.data.frame( t(sapply(trees, function(tr) rbind(as.numeric(tr$edge.length)))) )
  # mean
  mean.br.lens <- sapply(br.lens, mean) 
  # standard deviation
  sd.br.lens <- sapply(br.lens, sd) 
  # hpd 95%
  suppressMessages(require(TeachingDemos)) 
  # emp.hpd(trace, conf=0.95)
  hpd95 <- lapply(br.lens, emp.hpd, conf=0.95) 
  
  # all topologies
  edges.list <- lapply(trees, "[[", "edge")
  
  list(mean=mean.br.lens, sd=sd.br.lens, br.lens=br.lens, edges.list=edges.list)
}

# 1st line is node mapping: 50000	0	35..1,34..2,
# 2nd starts sequences: 	33	575741074201...
# internal node index is same to BEAST tree log, where id=nodeNr+1.
getIntNodeSeqStats <- function(ins.log="ins.txt", burnin=0.1, 
                               col.names=c("Sample","Node","States")) {
  require(tidyverse)
  # rm comments # 
  traces <- read_delim(ins.log, "\t", comment = "#", col_types = cols( States = col_character() ))
  if (! any(colnames(traces) %in% col.names) )
    stop("Incorrect file format in the log file of internal node sequences !\n", 
         "Column names = ", paste(colnames(traces), collapse = ","))
  #     Sample  Node States                                                                                                                                    
  #     <dbl> <dbl> <chr>                                                                                                                                     
  #1      0     0 35..1,34..2,33..3,33..4,...
  #2     NA    33 005741074201522219035315...
  #3     NA    34 005741074201522219035315...
  
  # separate nodes map and node states
  nodes.map <- traces %>% filter(Node==0)
  traces <- traces %>% filter(Node!=0)
  # MCMC summary
  samples = unique(nodes.map[[col.names[1]]]) # "Sample"
  in.nodes = unique(traces[[col.names[2]]]) # "Node"
  in.nodes = in.nodes[!is.na(in.nodes)]
  # state is a 2-digit integer from 00 to 59/60, n.codon = str.len / 2 
  str.len = str_length(traces[[col.names[3]]][1]) # "States"
  n.codon = str.len / 2 
  
  cat("Chain length", prettyNum(samples[length(samples)], big.mark=",",scientific=FALSE), 
      ", log every", prettyNum(samples[2], big.mark=",",scientific=FALSE), 
      "samples, each sample includes", length(in.nodes), "internal nodes [", in.nodes[1], "-", 
      in.nodes[length(in.nodes)], "], ", n.codon, " codons.\n")
  
  # rm burnin, +2 to exclude state 0 
  start = as.integer(burnin * length(samples)) + 2
  cat("Remove burnin ", start-1, " sampled internal node sequences from the total of ", 
      length(samples), "\n")
  
  nodes.map <- nodes.map[start:nrow(nodes.map),]
  ### TODO diff topology ?
  # assuming fixed tree
  edges <- nodes.map[["States"]][1] %>% str_split(",") %>% 
    unlist %>% enframe(name = NULL) %>% separate(value, c("parent", "child"))
  edges <- edges %>% drop_na # because of previous logging bug
  # branches == 2 * internal nodes
  stopifnot(nrow(edges) == 2*length(in.nodes))
  ###
  
  # create freq table of states
  freq.tb.list <- list()
  for (node.id in in.nodes) {
    # internal node index starts from n.taxa, node.id = 33
    node.samples <- traces %>% filter(Node==node.id)
    # rm burnin, +2 to exclude state 0 
    node.samples <- node.samples[start:nrow(node.samples),] %>% 
      separate(States, into = paste0("c", 1:n.codon), sep = seq(2, str.len, 2))
    # check
    stopifnot(nrow(node.samples) == nrow(nodes.map))
    
    ### posterior distribution
    freq.tb <- NULL
    for (c in 1:n.codon) {
      freq <- table(node.samples[[paste0("c", c)]])
      freq <- sort(freq, decreasing=T)
      s1 <- tibble(state=names(freq), freq=freq, site=c, order=1:length(freq))
      # calculate prob
      s1 <- s1 %>% mutate(prob = freq / sum(freq)) %>%
        mutate(cred = cumsum(prob)) # used to find 95% credible set
      freq.tb <- bind_rows(freq.tb, s1)
    }
    freq.tb <- freq.tb %>% mutate(order = as.integer(order))
    #  state  freq  site order    prob  cred
    #  <chr> <int> <int> <int>   <dbl> <dbl>
    #1 57      110     1     1 0.692   0.692
    #2 59       44     1     2 0.277   0.969
    #3 27        2     1     3 0.0126  0.981
    freq.tb.list[[as.character(node.id)]] <- freq.tb
  }
  # add edges
  freq.tb.list[["edges"]] <- edges
  freq.tb.list[["n.nodes"]] <- length(in.nodes)
  freq.tb.list[["n.codons"]] <- n.codon
  # names are internal node indexes + "edges"
  return(freq.tb.list) 
}



# branch lengths ESS per hour 
getBrLenESSHour <- function(tre.log="m0.trees", screen.log="out.txt", burnin=0.1) {
  require(ape)
  # read MCMC trees m0.std.trees
  stats <- getBrLensStats(tre.log, burnin)
  mean.bl <- stats$mean
  last.bl <- stats$br.lens[nrow(stats$br.lens),]
  # by columns
  ess <- apply(stats$br.lens, 2, getESS)  
  
  # total tree length
  total.tree.len <- rowSums(stats$br.lens)
  ess.ttl <- getESS(total.tree.len)
  
  # read screen log
  screen.info <- scan(screen.log,sep="\n",what="char(0)",skip=110)
  # Total calculation time: 71818.323 seconds
  time <- screen.info[length(screen.info)]
  time <- gsub("^.*time: (.*) seconds.*", "\\1", time)
  time <- as.numeric(time) / 3600 # ESS/hour
  
  list(ess.per.hour.branch.lens = ess / time, ess.per.hour.tot.lens = ess.ttl / time)
}






