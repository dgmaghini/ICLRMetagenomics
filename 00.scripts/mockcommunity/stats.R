library(tidyverse)
library(reshape2)
library(here)

##### READ STATS #####
readstats <- read.table(here("08.mock_stats/readstats.tsv"), header=TRUE)
readstats <- readstats %>% mutate(Method = gsub("nano", "Nano", gsub("short.*", "Short", gsub("inf", "Inf", gsub("_.*", "", file)))))
readstats <- readstats %>% select(Method, num_seqs, sum_len, avg_len, max_len, N50, Q20..., Q30..., GC...)

readstats <- readstats %>% mutate(
  num_seqs = as.numeric(gsub(",", "", num_seqs)) / 1000000, 
  sum_len = as.numeric(gsub(",", "", sum_len)) / 1000000000, 
  avg_len = as.numeric(gsub(",", "", avg_len)) / 1000, 
  max_len = as.numeric(gsub(",", "", max_len)) / 1000, 
  N50 = as.numeric(gsub(",", "", N50)) / 1000
)

names(readstats) <- c("Method", "Number of Reads (M)", "Total Sequencing (Gbp)", "Mean Read Length (Kbp)", 
                      "Maximum Read Length (Kbp)", "Read N50 (kbp)", "Reads ≥ Q20 (%)", "Reads ≥ Q30 (%)", "Read GC Content (%)")


readstats <- t(readstats)

write.table(readstats, here("08.mock_stats/read_table.tsv"), quote=FALSE, col.names=FALSE, sep="\t")


##### ASSEMBLY STATS #####

assemblystats <- read.csv(here("06.mockcommunity/05.subsample_new/1.assemblies/02.quast_new/transposed_report.tsv"), header=TRUE, sep="\t")

# format to have Method, replicate, depth
assemblystats <- assemblystats %>% mutate(Method=gsub("nano", "Nano", gsub("short", "Short", gsub("inf", "Inf", gsub("_.*", "", Assembly))))) %>% 
  mutate(Replicate = gsub(".*_", "", Assembly)) %>%
  mutate(Depth = gsub("_.*", "", gsub(".*subset_", "", Assembly)))

# filter to only have polished nanopore, only highest depth for each method
assemblystats <- assemblystats %>% filter(Method != "unpolished") %>%
  filter(Depth == "10Gb" | Depth == "7Gb")

assemblystats <- assemblystats %>% mutate(N50 = N50/1000) %>% 
  mutate(Total.length.....0.bp. = Total.length.....0.bp./ 1000000)

# group to get summary statistics by method
assembly_groupedstatistics <- assemblystats %>% group_by(Method) %>%
  summarise(
    n = n(), 
    meanassemblylength = mean(Total.length.....0.bp.), 
    assemblylengthse = sd(Total.length.....0.bp.) / sqrt(n), 
    meanN50 = mean(N50), 
    N50se = sd(N50) / sqrt(n)
  )

write.table(assembly_groupedstatistics, here("08.mock_stats/assembly_table.tsv"), quote=FALSE, row.names = FALSE, sep="\t")
##### BINNING STATS #####

binningstats <- read.csv(here("06.mockcommunity/05.subsample_new/2.bins/manual_quast/concatenated_report.tsv"), sep="\t", header=TRUE)

# format to have Method, replicate, depth, organism
binningstats <- binningstats %>% mutate(Method=gsub("nano", "Nano", gsub("short", "Short", gsub("inf", "Inf", gsub("_.*", "", Assembly))))) %>% 
  mutate(Replicate = gsub(".*_", "", gsub("-.*", "", Assembly))) %>%
  mutate(Depth = gsub("_.*", "", gsub(".*subset_", "", Assembly))) %>% 
  mutate(Organism = gsub(".*-", "", Assembly))

# filter to only have polished nanopore, only highest depth for each method
binningstats <- binningstats %>% filter(Method != "unpolished") %>%
  filter(Depth == "10Gb" | Depth == "7Gb")

# add column indicating whether bin is single contig
binningstats <- binningstats %>% mutate(SingleContig = ifelse(X..contigs.....0.bp. == 1, 1, 0))

# group to get summary statistics across organisms for each replicate
binningstats_byrep <- binningstats %>% group_by(Method, Replicate) %>% 
  summarise(
    numbins = n(), 
    singlecontigbins = sum(SingleContig)
  )

# group to get summary statistics across replicates for each method
binningstats_bymethod = binningstats_byrep %>% group_by(Method) %>% 
  summarise(
    n = n(), 
    meannumbins = mean(numbins), 
    numbinsse = sd(numbins) / sqrt(n), 
    meansinglecontigbins = mean(singlecontigbins), 
    singlecontigse = sd(singlecontigbins) / sqrt(n)
  )

write.table(binningstats_bymethod, here("08.mock_stats/binning_table.tsv"), quote=FALSE, row.names = FALSE, sep="\t")
