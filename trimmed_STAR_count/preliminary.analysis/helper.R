if(F){
  "
    Helper functions for the preliminaary analysis.
    Goal:
      1. get df of read count for the top genes for different samples
  "
}

#qqc

library(tidyverse)
library(progress)
library(org.Dm.eg.db)

read.sample <- function(df.path, top.perc = 100){
  # Function for reading in a SINGLE sample and get symbols
  #   1st - Read in original file and only keeping top 10% genes
  #     Did NOT consider the effect of gene length
  df <-
    read_delim(df.path, delim = "\t", 
             col_names = c("gene.id", "count.unstranded",
                           "count.first.read", "count.second.read"),
             col_types = "ciii") %>%
    filter(str_detect(gene.id, "^FBgn")) %>%
    filter(count.unstranded > quantile(count.unstranded, probs = (100-top.perc) / 100))
  #   2nd - Get symbol by gene.id.
  #     It is possible to fail and NA will be the result for those
  id2symbol <- function(id)
    tryCatch(org.Dm.egSYMBOL[[ org.Dm.egFLYBASE2EG[[ id ]] ]],
             error = function(e) NA)
  id2symbol <- Vectorize(id2symbol)
  df %>%
    mutate(symbol = id2symbol(gene.id)) %>%
    relocate(symbol, .after = gene.id)
}

read.samples <- function(named.samples.paths, .progress = T){
  sample.names <- names(named.samples.paths)
  map_dfr(sample.names, function(sample.name){
    if (.progress) message(paste(Sys.time(), 
                                 paste0("Processing ", sample.name),
                                 sep = "\t"))
    sample.path <- named.samples.paths[sample.name]
    read.sample(sample.path) %>%
      mutate(sample = sample.name)
  })
}


