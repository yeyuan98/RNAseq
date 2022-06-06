if(F){
  "
    Analysis functions for the preliminaary analysis.
  "
}

library(magrittr)
library(ggrepel)

# ------ PAIRWISE INTERSECTION ANALYSIS ------
pairwise.intersect.analysis <- function(df, 
                                    what = "symbol",
                                    by = "sample"){
  # Pairwise union analysis
  #   Given df with two columns what and by,
  #   for each `by` pair, how many `what` are shared?
  combn.pairs <- combn(unique(df[[by]]), 2)
  map_dfr(1:dim(combn.pairs)[2], function(idx){
    pair <- combn.pairs[,idx]
    
    df %>%
      filter(df[[by]] %in% pair) %>%  # Get only current pair
      with(intersect(get(what)[get(by) == pair[1]],
                 get(what)[get(by) == pair[2]])) -> shared.genes  # Get union
    
    tibble(combn.1 = pair[1], combn.2 = pair[2],
           count.1 = sum(df[[by]] == pair[1]),
           count.2 = sum(df[[by]] == pair[2]),
           count.shared = length(shared.genes))  # Get result for current pair
  })
}

# ------ TOP GENES WITH LARGEST NUMBER OF READS ------
top.reads.analysis <- function(df, num = 30){
  df %>%
    dplyr::select(symbol, count.unstranded, sample) %>%
    dplyr::rename(count = count.unstranded) %>%
    filter(complete.cases(df)) %>%
    group_by(sample) %>%
    mutate(count.percentage = count / sum(count) * 100) %>%
    slice_max(count.percentage, n = num) %>%
    dplyr::select(-count) %>%
    dplyr::arrange(desc(count.percentage), .by_group = T) %>%
    ungroup()
}

# ------ PAIRWISE FOLDCHANGE PLOTS
foldchange.plot <- function(count.df, s1, s2,
                            fold.change.minmax.values = NA,
                            mean.count.min.value = NA){
  # Plotting fold change and optionally tags the "significant" ones.
  #     s1 / s2 as y
  #     mean(s1, s2) as x
  count.df %<>%
    dplyr::select("symbol", s1, s2)
  count.df <- count.df[complete.cases(count.df),]  # Only genes which are shared
  message(paste0(nrow(count.df), " genes are shared by pair: ", s1, " vs ", s2))
  count.df %<>%
    mutate(fold.change = get(s1) / get(s2),
           mean.count = (get(s1) + get(s2)) / 2) %>%
    arrange(desc(fold.change), desc(mean.count))
  if (!is.na(fold.change.minmax.values)){
    # Remove symbol names for genes whose fold.change is [fc.min, fc.max]
    fc.min <- fold.change.minmax.values[1]
    fc.max <- fold.change.minmax.values[2]
    fc.sig <- count.df$fold.change < fc.min | count.df$fold.change > fc.max
    count.df$symbol[!fc.sig] <- ""
  }
  if (!is.na(mean.count.min.value)){
    # Remove symbol names for genes whose mean.count < min
    count.df$symbol[count.df$mean.count < mean.count.min.value] <- ""
  }
  # Base plot does not include text repel
  count.df %>%
    ggplot(aes(x = mean.count, y = fold.change, label = symbol))+
    scale_x_log10(expand = c(.01,0.1))+
    scale_y_continuous(expand = c(.2,.2))+
    theme_classic(base_size = 24) -> plt
  # If given significant thresholds, add repel
  if (!is.na(fold.change.minmax.values) | !is.na(mean.count.min.value)){
    plt <- plt+
      geom_text_repel(max.overlaps = 100)+
      geom_point(color = ifelse(count.df$symbol == "",
                                "grey50", "red"))
  } else{
    plt <- plt+geom_point()
  }
  plt
}
