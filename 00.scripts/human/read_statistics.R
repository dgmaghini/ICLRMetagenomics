library(ggplot2)
library(dplyr)
library(cowplot)
library(here)
library(ggpub)
library(scales)
library(paletteer) 

pal <- paletteer_d("nbapalettes::cavaliers")[c(1,2,3)]
names(pal) <- c("Short", "Infinity", "Nanopore")


#### read statistics ##### 
# read in and format read statistics
readstats <- read.table(here("07.human_stats/readstats.tsv"), header=TRUE)
readstats <- readstats %>% mutate(Method = gsub("_.*", "", file))
readstats <- readstats %>% mutate(Donor=gsub("_.*", "", gsub("-.*", "", gsub("010", "10", gsub("onor", "0",
                                        gsub(".dragen.*", "", gsub(".*D", "D", file)))))))

readstats <- readstats %>% dplyr::select( Method, Donor, num_seqs, sum_len, avg_len, max_len, N50, Q20..., Q30..., GC...)
readstats <- readstats %>% mutate(num_seqs = as.numeric(gsub(",", "", num_seqs)))
readstats <- readstats %>% mutate(sum_len = as.numeric(gsub(",", "", sum_len)))
readstats <- readstats %>% mutate(avg_len = as.numeric(gsub(",", "", avg_len)))
readstats <- readstats %>% mutate(max_len = as.numeric(gsub(",", "", max_len)))
readstats <- readstats %>% mutate(N50 = as.numeric(gsub(",", "", N50)))
readstats <- readstats %>% mutate(Q20... = as.numeric(Q20...))
readstats <- readstats %>% mutate(Q30... = as.numeric(Q30...))
readstats <- readstats %>% mutate(GC... = as.numeric(GC...))

# summarize statistics by method and donor (collate short reads)
grouped_donorstats <- readstats %>%
  group_by(Method, Donor) %>%
  summarise(
    num_seqs_total = sum(num_seqs),
    sum_len_total = sum(sum_len),
    avg_len = sum(sum_len) / sum(num_seqs), 
    max_len = max(max_len), 
    N50 = mean(N50), 
    Q20 = sum((Q20... * num_seqs))/ sum(num_seqs), 
    Q30 = sum((Q30... * num_seqs))/ sum(num_seqs),
    GC = sum((GC... * sum_len))/ sum(sum_len)
  )

# check read max and min 
readmaxes <- grouped_donorstats %>% 
  group_by(Method) %>% 
  summarise(
    readmax = max(sum_len_total), 
    readmin = min(sum_len_total)
  )

readmaxes <- readmaxes %>% mutate(readmax = readmax / 1000000000) %>% 
  mutate(readmin = readmin / 1000000000)

# summarize the statistics by method (collate donors) to average and standard error
grouped_methodstats <- grouped_donorstats %>%
  group_by(Method) %>%
  summarise(
    n = n(),
    num_seqs_avg = mean(num_seqs_total),
    num_seqs_se = sd(num_seqs_total) / sqrt(n), 
    sum_len_avg = mean(sum_len_total),
    sum_len_se = sd(sum_len_total) / sqrt(n), 
    avg_len_avg = mean(avg_len),
    avg_len_se = sd(avg_len) / sqrt(n), 
    max_len_avg = mean(max_len), 
    max_len_se = sd(max_len) / sqrt(n), 
    N50_avg = mean(N50), 
    N50_se = sd(N50) / sqrt(n), 
    Q20_avg = mean(Q20), 
    Q20_se = sd(Q20) / sqrt(n), 
    Q30_avg = mean(Q30), 
    Q30_se = sd(Q30) / sqrt(n), 
    GC_avg = mean(GC), 
    GC_se = sd(GC) / sqrt(n)
  )

grouped_methodstats <- grouped_methodstats %>% mutate(num_seqs_avg = num_seqs_avg / 1000000) %>%
  mutate(num_seqs_se = num_seqs_se / 1000000) %>% 
  mutate(sum_len_avg = sum_len_avg / 1000000000) %>% 
  mutate(sum_len_se = sum_len_se / 1000000000) %>% 
  mutate(avg_len_avg = avg_len_avg / 1000) %>% 
  mutate(avg_len_se= avg_len_se / 1000) %>% 
  mutate(max_len_avg = max_len_avg / 1000) %>% 
  mutate(max_len_se = max_len_se / 1000) %>% 
  mutate(N50_avg = N50_avg / 1000) %>% 
  mutate(N50_se = N50_se / 1000) 
  
  
write.table(grouped_methodstats, here("00.outputs/final/tables/humanmetagenomestats.tsv"), quote=FALSE, row.names = FALSE, sep="\t")



##### assembly statistics ##### 
# Read in QUAST reports
assembly_short <- read.csv(here("04.stat_comparison/shortread_quast_assembly.tsv"), header=TRUE, sep="\t")
assembly_infinity <- read.csv(here("04.stat_comparison/infinity_quast_assembly.tsv"), header=TRUE, sep="\t")
assembly_nanopore <- read.csv(here("04.stat_comparison/nanopore_quast_assembly.tsv"), header=TRUE, sep="\t")

# make donor names same format, add method column
assembly_short <- assembly_short %>% mutate(Method = "Short")

assembly_infinity <- assembly_infinity %>% mutate(Assembly=gsub("D010", "D10", gsub(".asm", "", gsub("Donor", "D0", Assembly)))) %>%
  mutate(Method = "Infinity")

assembly_nanopore <- assembly_nanopore %>% mutate(Assembly=gsub("-NF.*", "", Assembly)) %>% 
  mutate(Method = "Nanopore")

# merge tables
assembly_stats <- rbind(assembly_short, assembly_infinity, assembly_nanopore)
assembly_stats <- assembly_stats %>% mutate(Method = fct_relevel(Method, "Short", "Infinity", "Nanopore")) # reorder rows

assembly_stats <- assembly_stats %>% mutate(N50_kbp = N50 / 1000) %>% 
  mutate(length_Mbp = Total.length.....0.bp. / 1000000)
# summarize statistics by method
grouped_methodassemblystats <- assembly_stats %>%
  group_by(Method) %>% 
  summarise(
    n = n(),
    contigs_avg = mean(X..contigs.....0.bp.),
    contigs_se = sd(X..contigs.....0.bp.) / sqrt(n), 
    length_avg = mean(length_Mbp),
    length_se = sd(length_Mbp) / sqrt(n),
    N50_avg = mean(N50_kbp),
    N50_se = sd(N50_kbp) / sqrt(n)
  )

write.table(grouped_methodassemblystats, here("00.outputs/final/tables/humanmetagenomeassemblystats.tsv"), quote=FALSE, row.names = FALSE, sep="\t")

#####  binning statistics ##### 
short_bins <- read.csv(here("04.stat_comparison/shortread_binning_table_all_full.tsv"), sep="\t", header=TRUE)
short_bins <- mutate(short_bins, Method="Short")
infinity_bins <- read.csv(here("04.stat_comparison/infinity_binning_table_all_full.tsv"), sep="\t", header=TRUE) 
infinity_bins <- mutate(infinity_bins, Method="Infinity")
nanopore_bins <- read.csv(here("04.stat_comparison/nanopore_binning_table_all_full.tsv"), sep="\t", header=FALSE) %>% mutate(Method="Nanopore")
names(nanopore_bins) <- names(infinity_bins)
bins <- rbind(short_bins, infinity_bins, nanopore_bins) %>% filter(Bin != "unbinned")
bins <- bins %>% mutate(Method = fct_relevel(Method, "Short", "Infinity", "Nanopore")) # reorder rows
bins <- bins %>% mutate(Donor = gsub("-.*", "", Sample))
bins <- mutate(bins, Quality = ifelse(Completeness > 90 & Contamination < 5 & tRNA >= 18 & rna.16S > 0 & rna.23S > 0 & rna.5S > 0, "High-quality", ifelse(Completeness >= 50 & Contamination <10, "Medium-quality", "Low-quality")))
bins <- bins %>% mutate(SingleContig = ifelse(X..contigs.....0.bp. == 1, 1, 0))
bins <- bins %>% mutate(MediumQuality = ifelse(Quality == "Medium-quality", 1, 0))
bins <- bins %>% mutate(HighQuality = ifelse(Quality == "High-quality", 1, 0))
bins <- bins %>% mutate(has16S = ifelse(rna.16S >0, 1, 0))
bins <- bins %>% dplyr::select(Method, Donor, MediumQuality, HighQuality, has16S, SingleContig, rna.16S)
names(bins)

# summarize statistics by method and donor (collect info for all bins per donor)
grouped_donorbinstats <- bins %>%
  group_by(Method, Donor) %>%
  summarise(
    num_medquality = sum(MediumQuality),
    num_highquality = sum(HighQuality),
    num_16S = sum(has16S), 
    num_singlecontig = sum(SingleContig), 
    avg_16S = mean(rna.16S)
  )

grouped_methodbinstats <- grouped_donorbinstats %>%
  group_by(Method) %>% 
  summarise(
    n = n(),
    num_medquality_avg = mean(num_medquality),
    num_medquality_se = sd(num_medquality) / sqrt(n), 
    num_highquality_avg = mean(num_highquality),
    num_highquality_se = sd(num_highquality) / sqrt(n), 
    num_16S_avg = mean(num_16S), 
    num_16S_se = sd(num_16S) / sqrt(n), 
    num_singlecontig_avg = mean(num_singlecontig), 
    num_singlecontig_se = sd(num_singlecontig) / sqrt(n), 
    avg16S_avg = mean(avg_16S),
    avg16S_se = sd(avg_16S) / sqrt(n)
  )
  
write.table(grouped_methodbinstats, here("00.outputs/final/tables/humanmetagenomebinstats.tsv"), quote=FALSE, row.names = FALSE, sep="\t")

bins %>% filter(Method == "Infinity") %>% filter(SingleContig == TRUE) %>% select(Total.length, best_species)
