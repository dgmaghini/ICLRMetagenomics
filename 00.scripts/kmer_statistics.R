library(ggplot2)
library(dplyr)
library(tidyverse)
library(cowplot)
library(here)
library(reshape2)
library(paletteer)
library(ggpubr)
library(scales)
library(readr)

##### KMER CLUSTERING FROM SOURMASH #####
pal <- paletteer_d("nbapalettes::cavaliers")[c(1,2,3)]
names(pal) <- c("Short", "Infinity", "Nanopore")
pal <- c(pal, "Mixed"="#919191")

##### K51 
k51 <- read.csv(here("05.kmer_comparison/sourmash/04_sourmash_compare/compare_k51.csv"), header=TRUE)
colnames(k51) <- gsub(".*Donor", "Donor", colnames(k51))
colnames(k51) <- gsub("_con.*", "", colnames(k51))
rownames(k51) <- names(k51)
k51 <- cbind(Sample = rownames(k51), k51)
rownames(k51) <- NULL
data_long <- melt(k51, id.vars=c("Sample"))
names(data_long) <- c("Sample1", "Sample2", "Identity")

# add columns indicating method and donor for each element in comparison
data_long <- data_long %>% mutate(Method1 = ifelse(grepl("nanopore", Sample1), "nanopore", "infinity"))
data_long <- data_long %>% mutate(Method2 = ifelse(grepl("nanopore", Sample2), "nanopore", "infinity"))
data_long <- data_long %>% mutate(Donor1 = gsub("Donor", "", gsub("_.*","",Sample1)))
data_long <- data_long %>% mutate(Donor2 = gsub("Donor", "", gsub("_.*","",Sample2)))

# filter out self comparisons
data_long <- data_long %>% filter(Sample1 != Sample2)

# filter out identical/reversed comparisons
data_long <- data_long %>% mutate(uniqueID=ifelse(as.character(Sample1) > as.character(Sample2), paste(Sample1, Sample2), paste(Sample2, Sample1)))
data_long <- data_long %>% group_by(uniqueID) %>% filter(row_number() == 1)

# add a column indicating whether it's a within method or across method comparison
data_long <- data_long %>% mutate(ComparisonType=ifelse(Method1 == "infinity" && Method2 == "infinity", "Infinity", 
                                                        ifelse(Method1 == "nanopore" && Method2 == "nanopore", "Nanopore", "Mixed")))

data_long <- data_long %>% mutate(ComparisonDonorType=ifelse(Donor1 == Donor2 && ComparisonType == "Mixed", "MixedRelated",
                                                             ifelse(Donor1 != Donor2 && ComparisonType == "Mixed", "MixedUnrelated", 
                                                                    ifelse(Donor1 != Donor2 && ComparisonType == "Nanopore", "NanoporeUnrelated", "InfinityUnrelated"))))
data_long <- data_long %>% mutate(ComparisonDonorType = fct_relevel(ComparisonDonorType, levels=c("MixedRelated", "InfinityUnrelated", "NanoporeUnrelated", "MixedUnrelated"))) %>% 
  arrange(ComparisonDonorType)

data_long$ComparisonDonorType

ggplot(data_long, aes(x = ComparisonDonorType, y=Identity)) + 
  geom_jitter(color="darkgrey", width=0.1) +
  geom_boxplot(outlier.shape=NA, alpha=0.8, aes(fill=ComparisonType)) + 
  theme_bw() +
  scale_fill_manual(values=pal) +
  ggtitle("Similarity (k=51)") +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  scale_x_discrete(limits=c("MixedRelated", "MixedUnrelated", "InfinityUnrelated", "NanoporeUnrelated"),
                   labels = c("Different Methods\nSame Donor", "Different Methods\nDifferent Donor", "Infinity\nDifferent Donor", "Nanopore\nDifferent Donor")) + 
  stat_compare_means(method="wilcox.test", comparisons=list(c("InfinityUnrelated", "NanoporeUnrelated"), c("MixedRelated", "MixedUnrelated"), c("MixedRelated", "InfinityUnrelated"), c("MixedRelated", "NanoporeUnrelated")),
                     label = "p.signif")
 
ggsave(here("00.outputs/kmer_k51.pdf"), dpi=300, w=6, h=5)

##### K21 #
k21 <- read.csv(here("05.kmer_comparison/sourmash/04_sourmash_compare/compare_k21.csv"), header=TRUE)
colnames(k21) <- gsub(".*Donor", "Donor", colnames(k21))
colnames(k21) <- gsub("_con.*", "", colnames(k21))
rownames(k21) <- names(k21)
k21 <- cbind(Sample = rownames(k21), k21)
rownames(k21) <- NULL
data_long <- melt(k21, id.vars=c("Sample"))
names(data_long) <- c("Sample1", "Sample2", "Identity")

# add columns indicating method and donor for each element in comparison
data_long <- data_long %>% mutate(Method1 = ifelse(grepl("nanopore", Sample1), "nanopore", "infinity"))
data_long <- data_long %>% mutate(Method2 = ifelse(grepl("nanopore", Sample2), "nanopore", "infinity"))
data_long <- data_long %>% mutate(Donor1 = gsub("Donor", "", gsub("_.*","",Sample1)))
data_long <- data_long %>% mutate(Donor2 = gsub("Donor", "", gsub("_.*","",Sample2)))

# filter out self comparisons
data_long <- data_long %>% filter(Sample1 != Sample2)

# filter out identical/reversed comparisons
data_long <- data_long %>% mutate(uniqueID=ifelse(as.character(Sample1) > as.character(Sample2), paste(Sample1, Sample2), paste(Sample2, Sample1)))
data_long <- data_long %>% group_by(uniqueID) %>% filter(row_number() == 1)

# add a column indicating whether it's a within method or across method comparison
data_long <- data_long %>% mutate(ComparisonType=ifelse(Method1 == "infinity" && Method2 == "infinity", "Infinity", 
                                                        ifelse(Method1 == "nanopore" && Method2 == "nanopore", "Nanopore", "Mixed")))

data_long <- data_long %>% mutate(ComparisonDonorType=ifelse(Donor1 == Donor2 && ComparisonType == "Mixed", "MixedRelated",
                                                             ifelse(Donor1 != Donor2 && ComparisonType == "Mixed", "MixedUnrelated", 
                                                                    ifelse(Donor1 != Donor2 && ComparisonType == "Nanopore", "NanoporeUnrelated", "InfinityUnrelated"))))
data_long <- data_long %>% mutate(ComparisonDonorType = fct_relevel(ComparisonDonorType, levels=c("MixedRelated", "InfinityUnrelated", "NanoporeUnrelated", "MixedUnrelated"))) %>% 
  arrange(ComparisonDonorType)

data_long$ComparisonDonorType

ggplot(data_long, aes(x = ComparisonDonorType, y=Identity)) + 
  geom_jitter(color="darkgrey", width=0.1) +
  geom_boxplot(outlier.shape=NA, alpha=0.8, aes(fill=ComparisonType)) + 
  theme_bw() +
  scale_fill_manual(values=pal) +
  ggtitle("Similarity (k=21)") +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  scale_x_discrete(limits=c("MixedRelated", "MixedUnrelated", "InfinityUnrelated", "NanoporeUnrelated"),
                   labels = c("Different Methods\nSame Donor", "Different Methods\nDifferent Donor", "Infinity\nDifferent Donor", "Nanopore\nDifferent Donor")) + 
  stat_compare_means(method="wilcox.test", comparisons=list(c("InfinityUnrelated", "NanoporeUnrelated"), c("MixedRelated", "MixedUnrelated"), c("MixedRelated", "InfinityUnrelated"), c("MixedRelated", "NanoporeUnrelated")),
                     label = "p.signif")

ggsave(here("00.outputs/kmer_k21.pdf"), dpi=300, w=6, h=5)



##### KMER JELLYFISH HISTO #####
hist <- read.table(here("05.kmer_comparison/jellyfish/Donor2_infinity/mer_count_hist.tsv"), header=FALSE, sep=" ")
names(hist) <- c("NumberOccurrences", "NumberKmers")

ggplot(hist, aes(x=NumberOccurrences, y=NumberKmers)) +
  geom_vline(xintercept=1000, color="red") +
  geom_col() + 
  ylab("Number of k-mers") + 
  xlab("Number of Occurrences") +
  ggtitle("D02 Infinity k-mer Distribution") +
  scale_y_log10(labels=comma) + 
  scale_x_continuous(breaks=c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)) +
  theme_bw()

ggsave(here("00.outputs/kmer_histogram_example.pdf"), dpi=300, w=8, h=5)
ggsave(here("00.outputs/kmer_histogram_example.jpg"), dpi=300, w=8, h=5)


##### KMER JELLYFISH DESEQ #####
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("genefilter")
# BiocManager::install("DESeq2")
library("DESeq2")
library("genefilter")

## READ IN THE RAW FILE, FORMAT AS LONG, AND WRITE OUT
# counts <- read.table(here("05.kmer_comparison/jellyfish/mer_counts_all.tsv"))
# names(counts) <- c("Sample", "kmer", "count")
# # forwat wide to long
# counts_long <- reshape(counts, idvar="kmer", timevar = "Sample", direction="wide")
# counts_long[is.na(counts_long)] <- 0
# 
# write.table(counts_long, here("05.kmer_comparison/jellyfish/mer_counts_wide.tsv"), quote=FALSE, row.names=FALSE)

# READ IN WIDE FORMAT FILE
counts_long <- read.table(here("05.kmer_comparison/jellyfish/mer_counts_wide.tsv"), header=TRUE, sep=" ")

# turn column names into a metadata table
metadata <- data.frame(colnames(counts_long))
names(metadata) <- c("Sample")
metadata <- metadata %>% filter(Sample != "kmer") %>% filter(Sample != "count.Donor") %>% mutate(Sample=gsub("count.", "", Sample))
metadata <- metadata %>% mutate(Donor = gsub("_.*", "", Sample)) %>% mutate(Method = gsub(".*_", "", Sample))

# clean up column names in the counts table
colnames(counts_long) <- gsub(".*Donor", "Donor", colnames(counts_long))
counts_long <- counts_long %>% select(!Donor)
counts_long <- counts_long %>% remove_rownames %>% column_to_rownames(var="kmer")
counts_long <- counts_long %>% filter(!(row.names(counts_long) %in% c("Kmer")))

dds <- DESeqDataSetFromMatrix(countData = counts_long, colData = metadata, design = ~ Method)
dds <- estimateSizeFactors(dds, type = "poscounts")
idx <- genefilter(counts(dds, normalized = TRUE), pOverA(0.3, 999)) # 20% need to have entries
dds <- dds[idx, ]
dds <- DESeq(dds, parallel = T)
res <- results(dds, name="Method_nanopore_vs_infinity", alpha = 0.05)

resOrdered <- res %>% 
  as_tibble(rownames = "kmer") %>% 
  arrange(pvalue) %>% 
  mutate(group = ifelse(log2FoldChange < 0, "Infinity", "Nanopore"))

resOrdered <- resOrdered %>% mutate(padj_neglog10=ifelse(-log10(padj) > 30, 30, -log10(padj)))
resOrdered <- resOrdered %>% mutate(Enrichment=ifelse(padj < 0.05 & log2FoldChange < -1, "Infinity", ifelse(padj < 0.05 & log2FoldChange > 1, "Nanopore", "Mixed")))

counts_long %>% filter(row.names(counts_long) %in% c("CTATGAGTGTGTGCCTTGTTA"))

# volcano plot?
ggplot(resOrdered, aes(x=log2FoldChange, y=padj_neglog10, color=Enrichment)) +
  geom_hline(yintercept = 1.30103, color="grey", linetype="dashed") +
  geom_vline(xintercept = -1,color="grey", linetype="dotted") + 
  geom_vline(xintercept = 1,color="grey", linetype="dotted") + 
  geom_point(alpha = 0.9) +
  scale_color_manual(values=pal, limits = c("Nanopore", "Infinity")) +
  theme_bw() + 
  xlab("Log2 Fold Change") + 
  ylab("-log10 Adjusted P-value") + 
  theme(text = element_text(size=15))

ggsave(here("00.outputs/kmer_enrichment_volcano.jpg"), dpi = 300, w = 9, h=7)
ggsave(here("00.outputs/kmer_enrichment_volcano.pdf"), dpi = 300, w = 9, h=7)


# next, pull out some individual kmers enriched in nanopore
res_df_filt <- resOrdered %>% 
  filter(padj < 0.05, abs(log2FoldChange) > 15) %>% filter(Enrichment == "Infinity") %>% 
  arrange(log2FoldChange)

outf_infinity <- res_df_filt %>% select(kmer, baseMean, log2FoldChange, padj)

res_df_filt2 <- resOrdered %>% 
  filter(padj < 0.05, abs(log2FoldChange) > 5) %>% filter(Enrichment == "Nanopore") %>% 
  arrange(log2FoldChange)
outf_nanopore <- res_df_filt2 %>% select(kmer, baseMean, log2FoldChange, padj)

write_tsv(outf_infinity, here("05.kmer_comparison/jellyfish/kmers_topinfinity.tsv"), quote="none")
write_tsv(outf_nanopore, here("05.kmer_comparison/jellyfish/kmers_topnanopore.tsv"), quote="none")

names(res_df_filt)
top_kmer <- res_df_filt$kmer[1]
top_kmer2 <- res_df_filt$kmer[2]

counts <- reshape2::melt(as.matrix(counts_long))
names(counts) <- c("kmer", "Sample", "count")
counts <- mutate(counts, Method=gsub(".*_", "", Sample))
counts <- mutate(counts, Method=ifelse(Method == "infinity", "Infinity", "Nanopore"))


a <- ggplot(counts %>% filter(kmer == res_df_filt$kmer[1]), aes(x=Method, y=count)) +
  geom_jitter(width = 0.3, color="darkgrey") +
  geom_boxplot(outlier.shape=NA, alpha=0.8, aes(fill=Method)) + 
  scale_fill_manual(values=pal, limits=c("Infinity", "Nanopore")) + 
  theme_bw() + 
  ylab("k-mer count") +
  ggtitle(res_df_filt$kmer[1]) + 
  theme(legend.position="none", axis.title.x = element_blank(), text = element_text(size=14)) + 
  scale_y_continuous(labels=comma)
b <- ggplot(counts %>% filter(kmer == res_df_filt$kmer[2]), aes(x=Method, y=count)) +
  geom_jitter(width = 0.3, color="darkgrey") +
  geom_boxplot(outlier.shape=NA, alpha=0.8, aes(fill=Method)) + 
  scale_fill_manual(values=pal, limits=c("Infinity", "Nanopore")) + 
  theme_bw() + 
  ylab("k-mer count") +
  ggtitle(res_df_filt$kmer[2]) + 
  theme(legend.position="none", axis.title.x = element_blank(), text = element_text(size=14)) + 
  scale_y_continuous(labels=comma)
c <- ggplot(counts %>% filter(kmer == res_df_filt$kmer[3]), aes(x=Method, y=count)) +
  geom_jitter(width = 0.3, color="darkgrey") +
  geom_boxplot(outlier.shape=NA, alpha=0.8, aes(fill=Method)) + 
  scale_fill_manual(values=pal, limits=c("Infinity", "Nanopore")) + 
  theme_bw() + 
  ylab("k-mer count") +
  ggtitle(res_df_filt$kmer[3]) + 
  theme(legend.position="none", axis.title.x = element_blank(), text = element_text(size=14)) + 
  scale_y_continuous(labels=comma)
d <- ggplot(counts %>% filter(kmer == res_df_filt$kmer[4]), aes(x=Method, y=count)) +
  geom_jitter(width = 0.3, color="darkgrey") +
  geom_boxplot(outlier.shape=NA, alpha=0.8, aes(fill=Method)) + 
  scale_fill_manual(values=pal, limits=c("Infinity", "Nanopore")) + 
  theme_bw() + 
  ylab("k-mer count") +
  ggtitle(res_df_filt$kmer[4]) + 
  theme(legend.position="none", axis.title.x = element_blank(), text = element_text(size=14)) + 
  scale_y_continuous(labels=comma)

plot_grid(a,b,c,d, nrow=2, ncol=2, scale=0.9)
ggsave(here("00.outputs/kmer_nanopore_enriched.jpg"), w = 9, h=6, dpi=300)


# next, pull out some individual kmers enriched in infinity
res_df_filt <- resOrdered %>% 
  filter(padj < 0.05, abs(log2FoldChange) > 20) %>% filter(Enrichment == "Infinity") %>% 
  arrange(log2FoldChange)


a <- ggplot(counts %>% filter(kmer == res_df_filt$kmer[1]), aes(x=Method, y=count)) +
  geom_jitter(width = 0.3, color="darkgrey") +
  geom_boxplot(outlier.shape=NA, alpha=0.8, aes(fill=Method)) + 
  scale_fill_manual(values=pal, limits=c("Infinity", "Nanopore")) + 
  theme_bw() + 
  ylab("k-mer count") +
  ggtitle(res_df_filt$kmer[1]) + 
  theme(legend.position="none", axis.title.x = element_blank(), text = element_text(size=14))
b <- ggplot(counts %>% filter(kmer == res_df_filt$kmer[2]), aes(x=Method, y=count)) +
  geom_jitter(width = 0.3, color="darkgrey") +
  geom_boxplot(outlier.shape=NA, alpha=0.8, aes(fill=Method)) + 
  scale_fill_manual(values=pal, limits=c("Infinity", "Nanopore")) + 
  theme_bw() + 
  ylab("k-mer count") +
  ggtitle(res_df_filt$kmer[2]) + 
  theme(legend.position="none", axis.title.x = element_blank(), text = element_text(size=14))
c <- ggplot(counts %>% filter(kmer == res_df_filt$kmer[3]), aes(x=Method, y=count)) +
  geom_jitter(width = 0.3, color="darkgrey") +
  geom_boxplot(outlier.shape=NA, alpha=0.8, aes(fill=Method)) + 
  scale_fill_manual(values=pal, limits=c("Infinity", "Nanopore")) + 
  theme_bw() + 
  ylab("k-mer count") +
  ggtitle(res_df_filt$kmer[3]) + 
  theme(legend.position="none", axis.title.x = element_blank(), text = element_text(size=14))
d <- ggplot(counts %>% filter(kmer == res_df_filt$kmer[4]), aes(x=Method, y=count)) +
  geom_jitter(width = 0.3, color="darkgrey") +
  geom_boxplot(outlier.shape=NA, alpha=0.8, aes(fill=Method)) + 
  scale_fill_manual(values=pal, limits=c("Infinity", "Nanopore")) + 
  theme_bw() + 
  ylab("k-mer count") +
  ggtitle(res_df_filt$kmer[4]) + 
  theme(legend.position="none", axis.title.x = element_blank(), text = element_text(size=14))

plot_grid(a,b,c,d, nrow=2, ncol=2, scale=0.9)
ggsave(here("00.outputs/kmer_infinity_enriched.jpg"), w = 9, h=6, dpi=300)

