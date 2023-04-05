library(ggplot2)
library(tidyverse)
library(cowplot)
library(here)
library(ggpubr)
library(scales)
library(paletteer)
pal <- paletteer_d("nbapalettes::cavaliers")[c(1,2,3)]
names(pal) <- c("Short", "Infinity", "Nanopore")

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

# get median and standard deviation 
s <- assembly_stats %>% filter(Method == "Short")
i <- assembly_stats %>% filter(Method == "Infinity") 
n <- assembly_stats %>% filter(Method == "Nanopore")
median(n$N50)
sd(n$N50)

# assembly N50 plot with connecting lines
n50 <- ggplot(assembly_stats, aes(x=Method, y=N50, color=Method)) +
  geom_line(aes(group=Assembly), color="lightgrey", alpha=0.5) + 
  geom_point(size=2.5) + 
  theme_bw() + 
  ylab("Assembly N50 (bp)") +   
  scale_y_continuous(labels=comma) + 
  theme(text=element_text(size=15), axis.title.x = element_blank(), legend.position="none") +
  scale_color_manual(values=pal) + 
  stat_compare_means(test="wilcox.test", comparisons=list(c("Short", "Infinity"), c("Infinity", "Nanopore")), paired=TRUE, group.by="Assembly", label="p.signif", tip.length=0)
n50
# assembly contigs > 1000bp plot with connecting lines
num_contigs <- ggplot(assembly_stats, aes(x=Method, y=X..contigs.....1000.bp., color=Method)) +
  geom_line(aes(group=Assembly), color="lightgrey", alpha=0.5) + 
  geom_point(size=2.5) + 
  theme_bw() + 
  ylab("Contigs > 1,000 bp") + 
  scale_y_continuous(labels=comma) + 
  theme(text=element_text(size=15), axis.title.x = element_blank(), legend.position="none") + 
  scale_color_manual(values=pal) + 
  stat_compare_means(test="wilcox.test", comparisons=list(c("Short", "Infinity"), c("Infinity", "Nanopore")), paired=TRUE, group.by="Assembly", label="p.signif", tip.length=0)

# assembly N90 plot with connecting lines
n90 <- ggplot(assembly_stats, aes(x=Method, y=N90, color=Method)) +
  geom_line(aes(group=Assembly), color="lightgrey", alpha=0.5) + 
  geom_point(size=2.5) + 
  theme_bw() + 
  ylab("N90 (bp)") + 
  scale_y_continuous(labels=comma) + 
  scale_color_manual(values=pal) + 
  theme(text=element_text(size=15), axis.title.x = element_blank(), legend.position="none") + 
  stat_compare_means(test="wilcox.test", comparisons=list(c("Short", "Infinity"), c("Infinity", "Nanopore")), paired=TRUE, group.by="Assembly", label="p.signif", tip.length=0)

# assembly length plot with connecting lines
assembly_len <- ggplot(assembly_stats, aes(x=Method, y=Total.length, color=Method)) +
  geom_line(aes(group=Assembly), color="lightgrey", alpha=0.5) + 
  geom_point(size=2.5) + 
  theme_bw() + 
  ylab("Total Assembly Length (bp)") + 
  scale_y_continuous(labels=comma) + 
  scale_color_manual(values=pal) + 
  theme(text=element_text(size=15), axis.title.x = element_blank(), legend.position="none") + 
  stat_compare_means(test="wilcox.test", comparisons=list(c("Short", "Infinity"), c("Infinity", "Nanopore")), paired=TRUE, group.by="Assembly", label="p.signif", tip.length=0)

# assembly largest contig plot with connecting lines
largest_contig <- ggplot(assembly_stats, aes(x=Method, y=Largest.contig, color=Method)) +
  geom_line(aes(group=Assembly), color="lightgrey", alpha=0.5) + 
  geom_point(size=2.5) + 
  theme_bw() + 
  ylab("Largest Contig (bp)") + 
  scale_y_continuous(labels=comma) +   
  scale_color_manual(values=pal) + 
  theme(text=element_text(size=15), axis.title.x = element_blank(), legend.position="none") + 
  stat_compare_means(test="wilcox.test", comparisons=list(c("Short", "Infinity"), c("Infinity", "Nanopore")), paired=TRUE, group.by="Assembly", label="p.signif", tip.length=0)

plot_grid(n50, assembly_len, num_contigs, largest_contig, nrow=2, ncol=2, rel_widths = c(1, 1, 1, 1), align="hv", scale=0.9)
ggsave(here("00.outputs/assembly_comparison.pdf"), dpi=300, w=10, h=8)
