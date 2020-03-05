# plot result from test.beast.util.BinarySearchSamplingTest.binarySearchOptimalTest()

library(tidyverse)
library(ggplot2)

WD=paste0("~/WorkSpace/codonsubstmodels/perftest/")
setwd(WD)

linearBenchmarks <- read_delim("linearBenchmarks.txt", "\t", comment = "#", col_names = F)
biSeBenchmarks <- read_delim("binarySearchBenchmarks.txt", "\t", comment = "#", col_names = F)
colnames(linearBenchmarks) = c("state", "milliseconds")
colnames(biSeBenchmarks) = c("state", "milliseconds")
linearBenchmarks[[1]] = gsub("\\s:.*", "", linearBenchmarks[[1]])
biSeBenchmarks[[1]] = gsub("\\s:.*", "", biSeBenchmarks[[1]])
linearBenchmarks[["algorithm"]] = "linear"
biSeBenchmarks[["algorithm"]] = "binary-search"

benchmarks <- linearBenchmarks %>% bind_rows(biSeBenchmarks)
benchmarks[["state"]] <- factor(benchmarks[["state"]], levels = unique(benchmarks[["state"]]))
print(benchmarks, n=Inf)

p <- ggplot(benchmarks, aes(state, milliseconds, group = algorithm, colour = algorithm)) + 
  geom_line() + 
  ggtitle("100 million iterations") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = .5)) 
  
ggsave("../figures/linear-binarysearch.pdf", p, width = 12, height = 5)  
  