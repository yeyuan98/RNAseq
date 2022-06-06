if(F){
  "
    A preliminary analysis script downstream of STAR countPerGene output.
    Goals:
      1. get df of read count for the top genes for different samples
      2. do pairwise comparisons and plot FC ~ averaged.count
  "
}

source("helper.R")
source("analysis.functions.R")

df <- read.samples(c(CT2_RNAi = "../S1.out.tab",
                     CT2_CTR = "../S2.out.tab",
                     CT14_RNAi = "../S3.out.tab",
                     CT14_CTR = "../S4.out.tab"))

# Overall view of the top 10% genes
df %>%
  group_by(sample) %>%
  summarize(num.genes = n(),
            reads.unstranded = sum(count.unstranded),
            .groups = "drop")

# Pairwise shared gene analysis
#   For each sample pair, how many genes are shared?
pairwise.intersect.analysis(df)

# Pairwise fold change analysis
#   For each meaningful sample pair, plot fold changes of the shared genes
#   Getting widened df
df %>% filter(complete.cases(df)) %>% 
  dplyr::select(symbol, count.unstranded, sample) %>% 
  pivot_wider(names_from = sample, values_from = count.unstranded) -> df.fc
#   For each meaningful pair, plot fold change
