# plot result from test.beast.util.BinarySearchSamplingTest.binarySearchOptimalTest()

library(tidyverse)
library(ggplot2)
library(reshape2)

WD=paste0("~/WorkSpace/codonsubstmodels/perftest/")
setwd(WD)

linearBenchmarks <- read_delim("linearBenchmarks.txt", "\t", comment = "#", col_names = F)
colnames(linearBenchmarks)[1] = c("state")
# last col is total
linearBenchmarks = linearBenchmarks[,-ncol(linearBenchmarks)]
linearBenchmarks[["algorithm"]] = "linear"
melt.lin = melt(linearBenchmarks, id = c("state", "algorithm"))

biSeBenchmarks <- read_delim("binarySearchBenchmarks.txt", "\t", comment = "#", col_names = F)
colnames(biSeBenchmarks)[1] = c("state")
# last col is total
biSeBenchmarks = biSeBenchmarks[,-ncol(biSeBenchmarks)]
biSeBenchmarks[["algorithm"]] = "binary-search"
melt.bs = melt(biSeBenchmarks, id = c("state", "algorithm"))

benchmarks <- melt.lin %>% bind_rows(melt.bs)
benchmarks[["state"]] <- factor(benchmarks[["state"]], levels = unique(benchmarks[["state"]]))
benchmarks

p <- ggplot(benchmarks, aes(state, value, colour = algorithm)) + 
  geom_boxplot(width=0.5) + 
  ggtitle("50 tests with 20 million iterations") + ylab("milliseconds") +
  theme_minimal() 
  
ggsave("../figures/linear-binarysearch.pdf", p, width = 6, height = 5)  
  