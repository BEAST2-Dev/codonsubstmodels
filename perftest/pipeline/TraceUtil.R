# MCMC post-analysis

# samp.freq is the step size, namely sampling frequency of the samples
getESS <- function(samples) {
  n.samples = length(samples);
  max.lag = min(n.samples, 2000); # Note: lag index starts from 1 NOT 0
  mean.s = mean(samples)
  
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
