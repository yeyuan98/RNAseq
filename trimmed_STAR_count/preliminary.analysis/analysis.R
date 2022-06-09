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

df.fc.2 = as.data.frame(df.fc)
rownames(df.fc.2) = df.fc.2$symbol
df.fc.2 = df.fc.2 %>% dplyr::select(-c("symbol"))
df.fc.2 = as.matrix(df.fc.2)
df.fc.2 = as.matrix(df.fc.2)
df.fc.2[is.na(df.fc.2)] = 0
df.fc.normalized = sweep(df.fc.2,2,colSums(df.fc.2)/1e+6,`/`)
 df.coreclockgene <- df.fc.normalized[c("per","tim","cry","cyc","cwo","slmb","dco","mor"), ]
 df.coreclockgene.dataframe <- as.data.frame(df.coreclockgene)
 df.coreclockgene.dataframe $symbol = rownames(df.coreclockgene.dataframe)
df.fc.normalized.dataframe <- as.data.frame(df.fc.normalized)
df.fc.normalized.dataframe $symbol = rownames(df.fc.normalized.dataframe)

df.filtered = df.fc.normalized.dataframe %>% mutate( mean = (CT2_RNAi+CT2_CTR+CT14_RNAi+CT14_CTR)/4 )
df.filtered = df.filtered %>% filter(mean>50)
df.filtered = df.filtered %>% dplyr::select(-c("mean"))

#foldchange.plot(df.filtered, "CT2_RNAi", "CT2_CTR",
#                fold.change.minmax.values = c(0.5,2),
#                mean.count.min.value = 50) +
#  labs(title="CT2_RNAi/CT2_CTR")+
#  theme(axis.title = element_text(size = 12),
#plot.title = element_text(size = 14))

#foldchange.plot(df.coreclockgene.dataframe, "CT2_RNAi", "CT2_CTR",
#                fold.change.minmax.values = 1,
#                mean.count.min.value = 30) +
#  labs(title="CT2_RNAi/CT2_CTR")+
#  theme(axis.title = element_text(size = 12),
#        plot.title = element_text(size = 14))