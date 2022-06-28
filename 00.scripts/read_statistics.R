library(ggplot2)
library(tidyverse)
library(cowplot)
library(here)
library(ggpub)
library(scales)
library(paletteer) 

pal <- paletteer_d("nbapalettes::cavaliers")[c(1,2,3)]
names(pal) <- c("Short", "Infinity", "Nanopore")

# to do 
# number of reads per sample
# number of bases per sample
# read length distribution
nanopore_report <- read.csv(here("04.stat_comparison/nanopore_reads/report.tsv"), sep="\t", header=FALSE)
names(nanopore_report) <- c("Sample", "Reads", "Bases")
nanopore_report <- mutate(nanopore_report, Method = "Nanopore")

short_report <- read.csv(here("04.stat_comparison/shortread_reads/report.tsv"), sep="\t", header=TRUE)
names(short_report) <- c("Sample", "Reads", "Bases")
short_report <- mutate(short_report, Method = "Short")

infinity_report <- read.csv(here("04.stat_comparison/infinity_reads/report.tsv"), sep="\t", header=FALSE)
names(infinity_report) <- c("Sample", "Reads", "Bases")
infinity_report <- mutate(infinity_report, Method="Infinity")

report <- rbind(nanopore_report, short_report, infinity_report)
report <- report %>% mutate(Method = fct_relevel(Method, "Short", "Infinity", "Nanopore")) # reorder rows


ggplot(report, aes(x=Reads, y=Bases, color=Method)) + 
  geom_point() + 
  theme_bw() + 
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma, limits=c(0,35000000000)) 


bases <- ggplot(report, aes(x=Method, y=Bases, fill=Method)) + 
  geom_jitter(width=0.05, color="darkgrey") + 
  geom_boxplot(alpha = 0.7, outlier.shape=NA) +
  theme_bw() + 
  scale_y_continuous(labels=comma, limits=c(0, NA)) + 
  scale_fill_manual(values=pal) + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        text = element_text(size=15))

reads <- ggplot(report, aes(x=Method, y=Reads, fill=Method)) + 
  geom_jitter(width=0.05, color="darkgrey")+
  geom_boxplot(alpha=0.7, outlier.shape=NA) + 
  scale_fill_manual(values=pal) +
  scale_y_log10(labels=comma, breaks=c(1000000, 10000000,100000000), limits=c(999999, NA)) + 
  theme_bw() + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        text = element_text(size=15)) 
reads

plot_grid(reads, bases, nrow=1, ncol=2, align="hv")
ggsave(here("00.outputs/readcounts.pdf"), dpi=300, w=10, h=4)


nanopore_dist <- read.csv(here("04.stat_comparison/nanopore_reads/read_dist.tsv"), sep="\t", header=TRUE)
names(nanopore_dist) <- c("Sample", "Length", "Count")
nanopore_dist <- mutate(nanopore_dist, Method="Nanopore")

infinity_dist <- read.csv(here("04.stat_comparison/infinity_reads/read_dist.tsv"), sep="\t", header=FALSE)
names(infinity_dist) <- c("Sample", "Length", "Count")
infinity_dist <- mutate(infinity_dist, Method="Infinity")

dist <- rbind(nanopore_dist, infinity_dist) 

ggplot(dist %>% filter(Sample == "D08-NF-R1"), aes(x=Length, y=Count, fill=Method)) + 
  geom_col(alpha=0.5, position="dodge") + 
  theme_bw()  + 
  scale_x_continuous(limits = c(NA, 20000), labels = comma) + 
  scale_y_continuous(labels = comma) +
  ylab("Number of Reads") + 
  xlab("Read Length (bp)") + 
  scale_fill_manual(values=pal[2:3]) + 
  theme(text = element_text(size=15), legend.position=c(0.9, 0.9)) + 
  ggtitle("Donor 08 Read Lengths")

ggsave(here("00.outputs/readlengths.pdf"), dpi=300, h=6, w=11)
