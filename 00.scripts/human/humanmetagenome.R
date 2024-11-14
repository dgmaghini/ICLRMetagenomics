library(ggplot2)
library(dplyr)
library(cowplot)
library(here)
library(ggpubr)
library(scales)
library(paletteer)
library(ggExtra)
library(ggbeeswarm)

pal <- c("#8C1515", "#E98300", "#007C92") 
names(pal) <- c("Short", "ICLR", "ONT")

# Read in QUAST reports
assembly_short <- read.csv(here("04.stat_comparison/shortread_quast_assembly.tsv"), header=TRUE, sep="\t")
assembly_infinity <- read.csv(here("04.stat_comparison/infinity_quast_assembly.tsv"), header=TRUE, sep="\t")
assembly_nanopore <- read.csv(here("04.stat_comparison/nanopore_quast_assembly.tsv"), header=TRUE, sep="\t")

# make donor names same format, add method column
assembly_short <- assembly_short %>% mutate(Method = "Short")
assembly_infinity <- assembly_infinity %>% mutate(Assembly=gsub("D010", "D10", gsub(".asm", "", gsub("Donor", "D0", Assembly)))) %>%
  mutate(Method = "ICLR")
assembly_nanopore <- assembly_nanopore %>% mutate(Assembly=gsub("-NF.*", "", Assembly)) %>% 
  mutate(Method = "ONT")

# merge tables
assembly_stats <- rbind(assembly_short, assembly_infinity, assembly_nanopore)
assembly_stats <- assembly_stats %>% mutate(Method = fct_relevel(Method, "Short", "ICLR", "ONT")) # reorder rows
assembly_stats <- assembly_stats %>% mutate(TotalLengthMb=Total.length/1000000)
assembly_stats <- assembly_stats %>% mutate(N50_Mb=N50/1000000)
assembly_stats <- assembly_stats %>% mutate(LargestContigMb=Largest.contig/1000000)

mytheme <- theme(panel.grid = element_blank(), text = element_text(size = 10), axis.title.y = element_text(size=10), 
                 axis.title.x = element_blank(), 
                 legend.position = "none")

plot_dots <- function(dataframe, ycol, ytitle, yscale, lims, breakpoints) {
  ggplot(dataframe, aes(x=Method, y=get(ycol), color=Method)) + 
    geom_line(aes(group=Assembly), color="lightgrey", alpha=0.5) + 
    geom_point(size=2.5) + 
    theme_bw() + 
    ylab(ytitle) + 
    scale_y_continuous(limits = lims, breaks=breakpoints) + 
    scale_color_manual(values=pal) + 
    mytheme +  
    stat_compare_means(test="wilcox.test", 
                       comparisons=list(c("Short", "ICLR"), c("ICLR", "ONT"), c("Short", "ONT")), 
                       paired=TRUE, group.by="Assembly", label="p.signif", tip.length=0, size=3, label.y=yscale)
  
}

assemblyLength <- plot_dots(assembly_stats, "TotalLengthMb", "Assembly Length (Mb)", c(640, 675, 740), c(0, 800), c(0, 200, 400, 600, 800))
assemblyN50 <- plot_dots(assembly_stats, "N50_Mb", "Assembly N50 (Mb)", c(0.185, 0.202, 0.215), c(0, 0.23), c(0, 0.05, 0.10, 0.15, 0.20))
assemblyLargestContig <- plot_dots(assembly_stats, "LargestContigMb", "Largest Contig (Mb)", c(4.8, 5.2, 5.6), c(0,6), c(0, 2, 4, 6))

### Binning plots
# read in short read binning table
short_bins <- read.csv(here("04.stat_comparison/shortread_binning_table_all_full.tsv"), sep="\t", header=TRUE)
short_bins <- mutate(short_bins, Method="Short")

# read in infinity binning table
infinity_bins <- read.csv(here("04.stat_comparison/infinity_binning_table_all_full.tsv"), sep="\t", header=TRUE) 
infinity_bins <- mutate(infinity_bins, Method="ICLR")

# read in long read binning table
nanopore_bins <- read.csv(here("04.stat_comparison/nanopore_binning_table_all_full.tsv"), sep="\t", header=FALSE) %>% mutate(Method="ONT")
names(nanopore_bins) <- names(infinity_bins)

bins <- rbind(short_bins, infinity_bins, nanopore_bins) %>% filter(Bin != "unbinned")
bins <- bins %>% mutate(Method = fct_relevel(Method, "Short", "ICLR", "ONT")) # reorder rows

bins <- bins %>% mutate(ContaminationTransform = ifelse(Contamination == 0, 0.001, Contamination))
bins <- bins %>% mutate(N50_Mb = N50/1000000)

plot_box_beeswarm_log <- function(dataframe, ycol, ytitle, yscale) {
  ggplot(dataframe, aes(x=Method, y=get(ycol), fill=Method)) + 
    geom_beeswarm(color="darkgrey", size=0.8, cex=0.1) + 
    geom_boxplot(outlier.shape=NA, alpha = 0.7) + 
    theme_bw() + 
    scale_fill_manual(values=pal) +
    scale_y_log10(breaks=c(0.01, 1, 100), labels=c(0.01, 1, 100)) + 
    stat_compare_means(test = "wilcox.test", comparisons = list( c("Short", "ICLR"), c("ICLR", "ONT"), c("Short", "ONT")), 
                       label="p.signif", tip.length=0, label.y = yscale, size=3) + 
    ylab(ytitle) +
    mytheme
}

plot_box_beeswarm <- function(dataframe, ycol, ytitle, yscale, breakpoints) {
  ggplot(dataframe, aes(x=Method, y=get(ycol), fill=Method)) + 
    geom_beeswarm(color="darkgrey", size=0.8, cex=0.25) + 
    geom_boxplot(outlier.shape=NA, alpha = 0.7) + 
    theme_bw() + 
    scale_fill_manual(values=pal) +
    theme(legend.position = "none") + 
    stat_compare_means(test = "wilcox.test", comparisons = list( c("Short", "ICLR"), c("ICLR", "ONT"), c("Short", "ONT")), 
                       label="p.signif", tip.length=0, label.y = yscale, size=3) + 
    ylab(ytitle) +
    scale_y_continuous(expand = expansion(mult=0.1), breaks=breakpoints) + 
    mytheme
}

binContamination <- plot_box_beeswarm_log(bins, "ContaminationTransform", "Bin Contamination (%)", c(2, 2.1, 2.5))
binCompleteness <- plot_box_beeswarm(bins, "Completeness", "Bin Completeness (%)", c(100, 102, 111), c(0, 25, 50, 75, 100))
binN50 <- plot_box_beeswarm(bins, "N50_Mb", "Bin N50 (Mb)", c(4, 4.2, 4.5), c(0, 1, 2, 3, 4, 5))

# gene length plotting
nanopore_gene_length <- read.csv(here("03.bin_comparison/01.prokka/gene_lengths_nanopore.tsv"), header=FALSE, sep="\t")
infinity_gene_length <- read.csv(here("03.bin_comparison/01.prokka/gene_lengths_infinity.tsv"), header=FALSE, sep="\t")
short_gene_length <- read.csv(here("03.bin_comparison/01.prokka/gene_lengths_short.tsv"), header=FALSE, sep="\t")

names(nanopore_gene_length) <- c("Sample", "Length")
names(infinity_gene_length) <- c("Sample", "Length")
names(short_gene_length) <- c("Sample", "Length")

gene_lengths <- rbind(nanopore_gene_length, infinity_gene_length, short_gene_length)
gene_lengths <- gene_lengths %>% 
  mutate(Method = gsub("nanopore", "ONT", gsub("short", "Short", gsub("infinity", "ICLR", gsub("_.*", "", Sample))))) %>% 
  mutate(Donor = gsub("-N.*", "", gsub(".*D", "D", Sample)))

gene_lengths <- gene_lengths %>% mutate(Method = factor(Method, levels=c("Short", "ICLR", "ONT")))

binGeneLength <- plot_box_beeswarm(gene_lengths, "Length", "Mean Gene Length (bp)", c(1175, 1190, 1250), c(600, 800, 1000, 1200))

#### PLOT FIGURE
top_row <- plot_grid(assemblyLength, assemblyN50, assemblyLargestContig, nrow=1, labels=c("a", "b", "c"))
top_row

bottom_row <- plot_grid(binCompleteness, binContamination, binN50, binGeneLength, nrow=1, labels=c("d", "e", "f", "g"), rel_widths = c(1,1.02,0.95,1.049))
bottom_row

plot_grid(top_row, bottom_row, nrow=2, rel_heights = c(1, 1.2))
ggsave(here("00.outputs/final/human.pdf"), w=8, h=5, dpi=300)
ggsave(here("00.outputs/final/human.jpg"), w=8, h=5, dpi=300)



