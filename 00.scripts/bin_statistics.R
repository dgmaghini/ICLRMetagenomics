library(ggplot2)
library(tidyverse)
library(cowplot)
library(here)
library(ggpubr)
library(scales)
library(paletteer) 
library(ggExtra)


##### ALL BY ALL BIN COMPARISON ####
pal <- paletteer_d("nbapalettes::cavaliers")[c(1,2,3)]
names(pal) <- c("Short", "Infinity", "Nanopore")
# read in short read binning table
# read in infinity binning table
# read in long read binning table
short_bins <- read.csv(here("04.stat_comparison/shortread_binning_table_all_full.tsv"), sep="\t", header=TRUE)
short_bins <- mutate(short_bins, Method="Short")

infinity_bins <- read.csv(here("04.stat_comparison/infinity_binning_table_all_full.tsv"), sep="\t", header=TRUE) 
infinity_bins <- mutate(infinity_bins, Method="Infinity")

nanopore_bins <- read.csv(here("04.stat_comparison/nanopore_binning_table_all_full.tsv"), sep="\t", header=FALSE) %>% mutate(Method="Nanopore")
names(nanopore_bins) <- names(infinity_bins)

bins <- rbind(short_bins, infinity_bins, nanopore_bins) %>% filter(Bin != "unbinned")
bins <- bins %>% mutate(Method = fct_relevel(Method, "Short", "Infinity", "Nanopore")) # reorder rows

bins <- bins %>% mutate(ContaminationTransform = ifelse(Contamination == 0, 0.001, Contamination))
bins <- bins %>% mutate(N50_kb = N50/1000)
main <- ggplot(bins, aes(x=Completeness, y=ContaminationTransform, color=Method)) + 
  geom_point(alpha=0.5) + 
  theme_classic() + 
  scale_color_manual(values=pal) + 
  xlim(c(0,100)) +
  scale_y_log10() +
  xlab("Bin Completeness") + 
  ylab("Bin Contamination") + 
  theme(legend.position = c(0.1, 0.9)) + 
  theme(plot.margin = unit(c(0,0,0,0), "cm")) + 
  theme(panel.border = element_rect(color = "black", fill=NA, size=1), panel.grid.major = element_line())
#ggMarginal(p, type="density", margins = "y")

contam <- ggplot(bins, aes(y=ContaminationTransform, fill=Method)) + 
  geom_density(alpha=0.3, color=NA) + 
  theme_classic() +
  scale_y_log10() + 
  scale_fill_manual(values=pal) +
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), axis.title.y = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

contam

complete <- ggplot(bins, aes(y=Completeness, fill=Method)) + 
  geom_density(alpha = 0.3, color=NA) + 
  theme_classic() +
  scale_fill_manual(values=pal) +
  ylim(c(0, 100)) + 
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y=element_blank(), axis.ticks.y = element_blank()) + 
  coord_flip() + 
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

complete

plot_grid(complete, NULL, main, contam, ncol=2, nrow=2, align="hv", rel_heights = c(0.1, 1), rel_widths = c(1, 0.1))
#ggsave(here("00.outputs/completeness_by_contamination.pdf"), w=6, h=6, dpi=300)


# redo

cont <- ggplot(bins, aes(x=Method, y=ContaminationTransform, fill=Method)) + 
  geom_jitter(width=0.1, color="darkgrey") + 
  geom_boxplot(outlier.shape=NA, alpha = 0.7) + 
  theme_bw() + 
  scale_fill_manual(values=pal) +
  theme(legend.position = "none") + 
  scale_y_log10(labels = comma) + 
  stat_compare_means(test = "wilcox.test", comparisons = list( c("Short", "Infinity"), c("Infinity", "Nanopore"), c("Short", "Nanopore")), 
                     label="p.signif", tip.length=0, label.y=c(2, 2.05, 2.3)) + 
  ylab("Bin Contamination (%)") +
  theme(axis.title.x = element_blank(), text = element_text(size=15))

comp <- ggplot(bins, aes(x=Method, y=Completeness, fill=Method)) + 
  geom_jitter(width=0.1,  color="darkgrey") + 
  geom_boxplot(outlier.shape=NA, alpha = 0.7) + 
  theme_bw() + 
  scale_fill_manual(values=pal) +
  theme(legend.position = "none") + 
  stat_compare_means(test = "wilcox.test", comparisons = list( c("Short", "Infinity"), c("Infinity", "Nanopore"), c("Short", "Nanopore")), 
                     label="p.signif", tip.length=0, label.y = c(100, 101, 108)) + 
  ylab("Bin Completeness (%)") +
  theme(axis.title.x = element_blank(), text = element_text(size=15)) + 
  scale_y_continuous(limits=c(0, 115), breaks=c(25, 50, 75, 100))

n50 <- ggplot(bins, aes(x=Method, y=N50_kb, fill=Method)) + 
  geom_jitter(width=0.1, color="darkgrey") + 
  geom_boxplot(outlier.shape=NA, alpha = 0.7) + 
  theme_bw() + 
  scale_fill_manual(values=pal) + 
  theme(legend.position="none") + 
  ylab("Bin N50 (kb)") + 
  stat_compare_means(test = "wilcox.test", comparisons = list( c("Short", "Infinity"), c("Infinity", "Nanopore"), c("Short", "Nanopore")), 
                     label="p.signif", tip.length=0, label.y=c(4000, 4200, 4500)) + 
  theme(axis.title.x = element_blank(), text = element_text(size=15)) 
  
n50 
plot_grid(comp, cont, n50, nrow=1, ncol=3, align = "hv", scale = 0.9)

# plots
#    number of bins by sample by method
#    number of bins by quality by sample by method

binsnew <- bins %>% mutate(Quality = ifelse(bin.quality.call %in% c("4) high quality Bowers"), "High", 
                                                   ifelse(bin.quality.call %in% c("0) really bad", "1) low quality"), "Low", "Medium"))) 
                           
                           
counts <- binsnew %>% group_by(Sample, Method, Quality) %>% summarize(count = n())

counts <- counts %>% mutate(Quality = fct_relevel(Quality, "Low", "Medium", "High")) # reorder rows

## FIXME would be great to add stats to this
countplot <- ggplot(counts, aes(x=Quality, y=count, fill=Method)) + 
  geom_point(color = "darkgrey", position=position_dodge(width=.75)) + 
  geom_boxplot(outlier.shape=NA, alpha=0.7) +
  theme_bw() + 
  scale_fill_manual(values = pal) + 
  scale_color_manual(values = pal) + 
  #theme(legend.position="none") + 
  ylab("Number of Bins") + 
  theme(axis.title.x = element_blank(), text = element_text(size =15)) + 
  theme(legend.direction = "horizontal", legend.title = element_blank())
countplot
#ggsave(here("00.outputs/bincounts.pdf"), w=5, h=5, dpi=300)

legend <- get_legend(countplot + theme(legend.box.margin = margin(0, 0, 0, 12)))
a <- plot_grid(countplot + theme(legend.position="none"), comp, cont, n50, nrow=1, ncol=4, scale = 0.9, rel_widths=c(1.5, 1, 1, 1))
plot_grid(a, legend, nrow=2, ncol = 1, rel_heights = c(5, .4))
#ggsave(here("00.outputs/bin_statistics.pdf"), w=15, h=6, dpi=300)


#### MATCHING BIN COMPARISON #####
matching_bins <- read.csv(here("03.bin_comparison/matching_bins_paths.tsv"), header=TRUE, sep="\t")

infinity_bins <- read.csv(here("04.stat_comparison/infinity_binning_table_all_full.tsv"), sep="\t", header=TRUE) 
headernames <- names(infinity_bins)
colnames(infinity_bins) <- paste(colnames(infinity_bins), "infinity", sep="_")
# add Query identifier to infinity bins
infinity_bins <- mutate(infinity_bins, Query=paste("infinity_", Sample_infinity, "_", Bin_infinity, ".fa", sep=""))


nanopore_bins <- read.csv(here("04.stat_comparison/nanopore_binning_table_all_full.tsv"), sep="\t", header=FALSE) %>% mutate(Method="Nanopore")
names(nanopore_bins) <- headernames
colnames(nanopore_bins) <- paste(colnames(nanopore_bins), "nanopore", sep="_")

# add Reference identifier to nanopore bins
nanopore_bins <- mutate(nanopore_bins, Reference=paste("nanopore_", Sample_nanopore, "_", Bin_nanopore, ".fa", sep=""))


# merge all the bin statistics into a single table
matching_bins <- matching_bins %>% select(Query, Reference)
matching_bins <- merge(matching_bins, infinity_bins, by="Query")
matching_bins <- merge(matching_bins, nanopore_bins, by="Reference")


# calculate fold change in completeness, contamination, N50, size
#matching_bins <- mutate(matching_bins, changeN50=log2(N50_infinity/N50_nanopore))
matching_bins <- mutate(matching_bins, changeN50=log2(N50_infinity) - log2(N50_nanopore))
matching_bins <- mutate(matching_bins, changeSize=log2(Size.Mb_infinity) - log2(Size.Mb_nanopore))
matching_bins <- mutate(matching_bins, changeCompleteness=log2(Completeness_infinity) - log2(Completeness_nanopore))
matching_bins <- mutate(matching_bins, changeContamination=log2(ifelse(Contamination_infinity == 0, 0.1, Contamination_infinity)) - log2(ifelse(Contamination_nanopore == 0, 0.1, Contamination_nanopore)))


matching_bins %>% select(Size.Mb_infinity, Size.Mb_nanopore, changeSize)

p <- ggplot(matching_bins, aes(x=changeSize, y=changeN50)) + 
  geom_hline(yintercept=0, color="lightgrey") +
  geom_vline(xintercept=0, color="lightgrey") + 
  geom_point() + 
  theme_bw() + 
  xlab("Genome Size Change\nlog2(Infinity/Nanopore)") + 
  ylab("Genome N50 Change\nlog2(Infinity/Nanopore)") + 
  scale_x_continuous(labels=comma, limits=c(-4, 4)) + 
  ylim(c(-6.5, 6.5)) + 
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(), text=element_text(size=15)) 

p1 <- ggMarginal(p, type="density", fill="grey", color="white", size=10)
p1
#ggsave(here("00.outputs/paired_bins_size_scatter.pdf"), p1, w=6, h=6, dpi=300)

ggplot(matching_bins, aes(x=changeContamination, y=changeCompleteness)) + 
  geom_hline(yintercept=0, color="lightgrey") +
  geom_vline(xintercept=0, color="lightgrey") + 
  geom_point() + 
  theme_bw() + 
  xlab("Contamination Log2 Fold Change (Infinity/Nanopore)") + 
  ylab("Completeness N50 Log2 Fold Change (Infinity/Nanopore)") + 
  scale_x_continuous(labels=comma, limits=c(-10, 10)) + 
  scale_y_continuous(limits=c(-4, 4)) +
  theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank())

matching_bins %>% arrange(changeN50) %>% select(Sample_infinity, Bin_infinity, Sample_nanopore, Bin_nanopore, Size.Mb_nanopore, Size.Mb_infinity, lca_species_nanopore, changeSize, changeN50) %>% head()
