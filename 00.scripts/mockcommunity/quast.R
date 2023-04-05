library(here)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(paletteer)
library(cowplot)
library(raster)
library(xts)

##### ORIGINAL QUAST (Not subsampled) #####
pal <- paletteer_d("nbapalettes::cavaliers")[c(1,2,3)]
pal <- c(pal, "#C9A964FF")
names(pal) <- c("Short", "Infinity", "Nanopore\n(Polished)", "Nanopore\n(Unpolished)")

report <- read.csv(here("06.mockcommunity/03.genomequast/transposed_report.tsv"), sep="\t", header=TRUE)

report <- report %>% mutate(Method=ifelse(gsub("_.*", "", Assembly) == "03.alignreference", "Infinity", ifelse(
  gsub("_.*", "", Assembly) == "unpolished", "Nanopore\n(Unpolished)", "Nanopore\n(Polished)")))

report <- report %>% mutate(Genome = gsub(".*_", "", Assembly))

report <- report %>% mutate(Method = fct_relevel(Method, "Infinity", "Nanopore\n(Unpolished)", "Nanopore\n(Polished)")) 

a <- ggplot(report, aes(x=Method, y=X..indels.per.100.kbp)) + 
  geom_line(aes(group=Genome), color="lightgrey", alpha=0.5) +
  geom_point(aes(color=Method), size=3) + 
  scale_color_manual(values=pal[c(2,3,4)])+
  ylab("Indels per 100 kbp") + 
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.title.x = element_blank(), 
        legend.position = "none") + 
  stat_compare_means(test="wilcox.test", comparisons=list(c("Infinity", "Nanopore\n(Unpolished)"), c("Nanopore\n(Unpolished)", "Nanopore\n(Polished)"), c("Infinity", "Nanopore\n(Polished)")), 
                     paired=TRUE, group.by="Genome", label="p.signif", tip.length=0, label.y = c(165, 160, 175)) +
  theme(text=element_text(size=14))
a
ggsave(here("00.outputs/mockcommunity/quast_indels.pdf"), w=5, h=5, dpi=300)

b <- ggplot(report, aes(x=Method, y=X..mismatches.per.100.kbp)) + 
  geom_line(aes(group=Genome), color="lightgrey", alpha=0.5) +
  geom_point(aes(color=Method), size=3) + 
  scale_color_manual(values=pal[c(2,3,4)])+
  ylab("Mismatches per 100 kbp") + 
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.title.x = element_blank(), 
        legend.position = "none") + 
  stat_compare_means(test="wilcox.test", comparisons=list(c("Infinity", "Nanopore\n(Unpolished)"), c("Nanopore\n(Unpolished)", "Nanopore\n(Polished)"), c("Infinity", "Nanopore\n(Polished)")), 
                     paired=TRUE, group.by="Genome", label="p.signif", tip.length=0, label.y = c(315, 305, 335)) +
  theme(text=element_text(size=14))
b
ggsave(here("00.outputs/mockcommunity/quast_mismatches.pdf"), w=5, h=5, dpi=300)


plot_grid(a,b, nrow=1)
ggsave(here("00.outputs/mockcommunity/quast_accuracy.pdf"), w=7, h=3.5, dpi=300)

ggplot(report, aes(x=Method, y=Duplication.ratio)) + 
  geom_hline(yintercept=1, color="darkgrey") +
  geom_line(aes(group=Genome), color="lightgrey", alpha=0.5) +
  geom_point(aes(color=Method), size=3) +
  #geom_text(inherit.aes = FALSE, data = report %>% filter(Method == "Nanopore\n(Polished)") %>% filter(Duplication.ratio > 1.2), aes(x=Method, y=Duplication.ratio, label=Genome, vjust=-1)) +
  scale_color_manual(values=pal[c(2,3,4)])+
  ylim(0.998, 1.022) +
  ylab("Duplication Ratio") + 
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.title.x = element_blank(), 
        legend.position = "none") + 
  stat_compare_means(test="wilcox.test", comparisons=list(c("Infinity", "Nanopore\n(Unpolished)"), c("Nanopore\n(Unpolished)", "Nanopore\n(Polished)"), c("Infinity", "Nanopore\n(Polished)")), 
                     paired=TRUE, group.by="Genome", label="p.signif", tip.length=0, label.y = c(1.015, 1.018, 1.021)) +
  theme(text=element_text(size=14))

ggsave(here("00.outputs/mockcommunity/quast_duplicationratio_updated.pdf"), w=5, h=5, dpi=300)


ggplot(report, aes(x=Method, y=Genome.fraction....)) + 
  geom_line(aes(group=Genome), color="lightgrey", alpha=0.5) +
  geom_point(aes(color=Method), size=3) + 
  geom_text(inherit.aes = FALSE, data = report %>% filter(Method == "Nanopore\n(Polished)") %>% filter(Genome.fraction.... < 95), aes(x=Method, y=Genome.fraction...., label=Genome), vjust=-1) +
  scale_color_manual(values=pal[c(2,3,4)])+
  ylab("Genome Fraction") + 
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.title.x = element_blank(), 
        legend.position = "none") + 
  ylim(85,105) +
  stat_compare_means(test="wilcox.test", comparisons=list(c("Infinity", "Nanopore\n(Unpolished)"), c("Nanopore\n(Unpolished)", "Nanopore\n(Polished)"), c("Infinity", "Nanopore\n(Polished)")), 
                     paired=TRUE, group.by="Genome", label="p.signif", tip.length=0, label.y = c(102, 103, 104)) +
  theme(text=element_text(size=14))

ggsave(here("00.outputs/mockcommunity/quast_genomefraction.pdf"), w=5, h=5, dpi=300)


ggplot(report, aes(x=Method, y=X..contigs.....0.bp.)) + 
  geom_line(aes(group=Genome), color="lightgrey", alpha=0.5) +
  geom_point(aes(color=Method), size=3) + 
  geom_text(inherit.aes = FALSE, data = report %>% filter(Method == "Nanopore\n(Polished)") %>% filter(X..contigs.....0.bp. > 50), aes(x=Method, y=X..contigs.....0.bp., label=Genome), vjust=-1) +
  scale_color_manual(values=pal[c(2,3,4)])+
  ylab("Number of Contigs") + 
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.title.x = element_blank(), 
        legend.position = "none") + 
  stat_compare_means(test="wilcox.test", comparisons=list(c("Infinity", "Nanopore\n(Unpolished)"), c("Nanopore\n(Unpolished)", "Nanopore\n(Polished)"), c("Infinity", "Nanopore\n(Polished)")), 
                     paired=TRUE, group.by="Genome", label="p.signif", tip.length=0, label.y = c(241, 253, 265)) +
  theme(text=element_text(size=14))

ggsave(here("00.outputs/mockcommunity/quast_contignumber.pdf"), w=5, h=5, dpi=300)

ggplot(report, aes(x=Method, y=X..contigs.....0.bp.)) + 
  geom_line(aes(group=Genome), color="lightgrey", alpha=0.5) +
  geom_point(aes(color=Method), size=3) + 
  geom_text(inherit.aes = FALSE, data = report %>% filter(Method == "Nanopore\n(Polished)") %>% filter(Genome == "Escherichia" | Genome == "Salmonella" | Genome == "Listeria"), aes(x=Method, y=X..contigs.....0.bp., label=Genome), hjust=-0.5) +
  scale_color_manual(values=pal[c(2,3,4)])+
  ylab("Number of Contigs") + 
  theme_bw() + 
  ylim(0,23) +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.title.x = element_blank(), 
        legend.position = "none") + 
  theme(text=element_text(size=14))

ggsave(here("00.outputs/mockcommunity/quast_contignumberzoom.pdf"), w=7, h=5, dpi=300)


ggplot(report, aes(x=Method, y=GC.... - Reference.GC....)) + 
  geom_hline(yintercept = 0, color="darkgrey") +
  geom_line(aes(group=Genome), color="lightgrey", alpha=0.5) +
  geom_point(aes(color=Method), size=3) + 
  scale_color_manual(values=pal[c(2,3,4)])+
  ylab("Increase in GC Relative to Reference") + 
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.title.x = element_blank(), 
        legend.position = "none") + 
  stat_compare_means(test="wilcox.test", comparisons=list(c("Infinity", "Nanopore\n(Unpolished)"), c("Nanopore\n(Unpolished)", "Nanopore\n(Polished)"), c("Infinity", "Nanopore\n(Polished)")), 
                     paired=TRUE, group.by="Genome", label="p.signif", tip.length=0, label.y = c(0.25, 0.27, 0.29)) +
  theme(text=element_text(size=14))

ggsave(here("00.outputs/mockcommunity/quast_gcchange.pdf"), w=5, h=5, dpi=300)


ggplot(report, aes(x=Reference.GC...., y=GC.... - Reference.GC...., color=Method)) + 
  geom_hline(yintercept = 0, color="darkgrey") +
  geom_point(aes(color=Method), size=3) + 
  scale_color_manual(values=pal[c(2,å3,4)])+
  xlab("Reference %GC") + 
  ylab("Change in %GC Relative to Reference") + 
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank()) +
  theme(text=element_text(size=14))

ggsave(here("00.outputs/mockcommunity/quast_gcchangebyref.pdf"), w=7, h=5, dpi=300)


##### QUAST ON SUBSAMPLED DATA, AGGREGATED, DEPRECATED #####
pal <- paletteer_d("nbapalettes::cavaliers")[c(1,2,3)]
pal <- c(pal, "#C9A964FF", "#064396FF")
names(pal) <- c("Short", "Infinity1", "Nanopore\n(Polished)", "Nanopore\n(Unpolished)", "Infinity2")
report <- read.csv(here("06.mockcommunity/04.subsample/01.quast/allgenomes/transposed_report.tsv"), sep="\t", header=TRUE)

report <- report %>% mutate(Method = gsub("_subset_.*assembly", "", Assembly))
report <- report %>% mutate(Method = gsub("infinity", "Infinity", gsub("_polished", "\n(Polished)", gsub("_unpolished", "\n(Unpolished)", gsub("nano", "Nano", Method)))))
                            

report <- report %>% mutate(Depth = gsub(".*subset_", "", gsub("_assembly.*", "", Assembly)))
report <- report %>% mutate(DepthNumeric = as.numeric(gsub("Mb", "", gsub("Gb", "", Depth))))
report <- report %>% mutate(DepthNumeric = ifelse(DepthNumeric == 500, 0.5, DepthNumeric))
report <- report %>% mutate(Method = fct_relevel(Method, "Infinity1", "Infinity2", "Nanopore\n(Unpolished)", "Nanopore\n(Polished)")) 

a <- ggplot(report, aes(x=DepthNumeric, y=N50)) + 
  geom_line(aes(color=Method)) + 
  geom_point(aes(color=Method)) + 
  scale_color_manual(values=pal[c(2,3,4,5)])+
  xlab("Depth (Gb)") + 
  ylab("Contig N50") +
  theme_bw() + 
  theme(legend.position = "bottom")

b <- ggplot(report, aes(x=DepthNumeric, y=Total.length.....10000.bp.)) + 
  geom_line(aes(color=Method)) + 
  geom_point(aes(color=Method)) + 
  scale_color_manual(values=pal[c(2,3,4,5)])+
  xlab("Depth (Gb)") + 
  ylab("Assembly Length (>10 kbp)") +
  theme_bw() 

c <- ggplot(report, aes(x=DepthNumeric, y=X..mismatches.per.100.kbp)) + 
  geom_line(aes(color=Method)) + 
  geom_point(aes(color=Method)) + 
  scale_color_manual(values=pal[c(2,3,4,5)])+
  xlab("Depth (Gb)") + 
  ylab("Mismatches per 100 Kbp") +
  theme_bw()

d <- ggplot(report, aes(x=DepthNumeric, y=X..indels.per.100.kbp)) + 
  geom_line(aes(color=Method)) + 
  geom_point(aes(color=Method)) + 
  scale_color_manual(values=pal[c(2,3,4,5)])+
  xlab("Depth (Gb)") + 
  ylab("Indels per 100 Kbp") +
  theme_bw()

z <- get_legend(a)

e <- plot_grid(a + theme(legend.position = "none"), b + theme(legend.position = "none"), 
          c + theme(legend.position = "none"), d + theme(legend.position = "none"), align = "hv")

plot_grid(e, z, nrow=2, ncol=1, rel_heights = c(1, 0.1))
ggsave(here("00.outputs/mockcommunity/quast_depthreport.pdf"), w=5, h=5, dpi=300)

a <- ggplot(report, aes(x=DepthNumeric, y=Duplication.ratio)) + 
  geom_line(aes(color=Method)) + 
  geom_point(aes(color=Method)) + 
  scale_color_manual(values=pal[c(2,3,4,5)])+
  ylim(1,1.007) + 
  xlab("Depth (Gb)") + 
  ylab("Duplication Ratio") +
  theme_bw() + 
  theme(legend.position = "bottom")

b <- ggplot(report, aes(x=DepthNumeric, y=Genome.fraction....)) + 
  geom_line(aes(color=Method)) + 
  geom_point(aes(color=Method)) + 
  scale_color_manual(values=pal[c(2,3,4,5)])+
  ylim(0,100) + 
  xlab("Depth (Gb)") + 
  ylab("Genome Fraction") +
  theme_bw()

z <- get_legend(a)

temp <- plot_grid(a + theme(legend.position = "none"), 
          b + theme(legend.position = "none"),  
          rel_widths = c(1,1))
plot_grid(temp, z, nrow=2, ncol=1, rel_heights = c(1, 0.1))

ggsave(here("00.outputs/mockcommunity/quast_depthduplicationreport.pdf"), w=4.5, h=2.8, dpi=300)


#### QUAST ON SUBSAMPLED DATA, PER GENOME, DEPRECATED ####

pal <- paletteer_d("nbapalettes::cavaliers")[c(1,2,3)]
pal <- c(pal, "#C9A964FF", "#064396FF")
names(pal) <- c("Short", "Infinity1", "Nanopore\n(Polished)", "Nanopore\n(Unpolished)", "Infinity2")

report <- read.csv(here("06.mockcommunity/04.subsample/01.quast/transposed_report.tsv"), sep="\t", header=TRUE)

report <- report %>% mutate(Method=gsub("genomes_", "", gsub("_subset.*", "", Assembly)))
report <- report %>% mutate(Method=gsub("infinity", "Infinity", gsub("^polished", "Nanopore\n(Polished)", 
                                                                     gsub("^unpolished", "Nanopore\n(Unpolished)", Method))))

report <- report %>% mutate(Genome = gsub(".*_", "", Assembly))
report <- report %>% mutate(Depth = gsub("_.*", "", gsub(".*subset_", "", Assembly)))
report <- report %>% mutate(DepthNumeric = as.numeric(gsub("Mb", "", gsub("Gb", "", Depth))))
report <- report %>% mutate(DepthNumeric = ifelse(DepthNumeric == 500, 0.5, DepthNumeric))

report <- report %>% mutate(Method = fct_relevel(Method, "Infinity1", "Infinity2", "Nanopore\n(Unpolished)", "Nanopore\n(Polished)")) 

ggplot(report, aes(x=DepthNumeric, y=N50)) + 
  geom_line(aes(color=Method)) + 
  geom_point(aes(color=Method)) + 
  scale_color_manual(values=pal[c(2,3,4,5)])+
  xlab("Depth (Gb)") + 
  ylab("Contig N50") +
  theme_bw() + 
  facet_wrap(~Genome)

ggsave(here("00.outputs/mockcommunity/quastdepth_N50.pdf"), w=6, h=5, dpi=300)

ggplot(report, aes(x=DepthNumeric, y=L50)) + 
  geom_line(aes(color=Method)) + 
  geom_point(aes(color=Method)) + 
  scale_color_manual(values=pal[c(2,3,4,5)])+
  xlab("Depth (Gb)") + 
  ylab("L50") +
  theme_bw() + 
  facet_wrap(~Genome, scales="free_y")

ggsave(here("00.outputs/mockcommunity/quastdepth_L50.pdf"), w=6.5, h=5, dpi=300)

ggplot(report, aes(x=DepthNumeric, y=L50)) + 
  geom_line(aes(color=Method)) + 
  geom_point(aes(color=Method)) + 
  scale_color_manual(values=pal[c(2,3,4,5)])+
  xlab("Depth (Gb)") + 
  ylab("L50") +
  theme_bw() + 
  facet_wrap(~Genome, scales="free_y")

ggplot(report, aes(x=DepthNumeric, y=Total.length.....10000.bp.)) + 
  geom_line(aes(color=Method)) + 
  geom_point(aes(color=Method)) + 
  scale_color_manual(values=pal[c(2,3,4,5)])+
  xlab("Depth (Gb)") + 
  ylab("Total.length.....10000.bp.") +
  theme_bw() + 
  facet_wrap(~Genome, scales="free_y")

ggplot(report, aes(x=DepthNumeric, y=X..mismatches.per.100.kbp)) + 
  geom_line(aes(color=Method)) + 
  geom_point(aes(color=Method)) + 
  scale_color_manual(values=pal[c(2,3,4,5)])+
  xlab("Depth (Gb)") + 
  ylab("X..mismatches.per.100.kbp") +
  theme_bw() + 
  facet_wrap(~Genome, scales="free_y")

ggplot(report, aes(x=as.factor(DepthNumeric), y=X..indels.per.100.kbp, color=Method, fill=Method)) + 
  geom_point(position=position_dodge(width = 0.75)) + 
  geom_boxplot(aes(fill=Method), alpha=0.8, outlier.shape=NA) + 
  scale_fill_manual(values=pal[c(2,3,4,5)])+ 
  scale_color_manual(values=pal[c(2,3,4,5)]) + 
  xlab("Depth (Gb)") + 
  ylab("X..indels.per.100.kbp") + 
  theme_bw()

ggplot(report, aes(x=as.factor(DepthNumeric), y=N50, color=Method, fill=Method)) + 
  geom_point(position=position_dodge(width = 0.75)) + 
  geom_boxplot(aes(fill=Method), alpha=0.8, outlier.shape=NA) + 
  scale_fill_manual(values=pal[c(2,3,4,5)])+ 
  scale_color_manual(values=pal[c(2,3,4,5)]) + 
  xlab("Depth (Gb)") + 
  ylab("N50 (bp)") + 
  theme_bw()

ggsave(here("00.outputs/mockcommunity/quastdepth_N50_boxplot.pdf"), w=7, h=4, dpi=300)



##### QUAST ON SUBSAMPLED DATA, PER GENOME WITH BINNING ASSIGNMENTS, DEPRECATED #####


pal <- paletteer_d("nbapalettes::cavaliers")[c(1,2,3)]
pal <- c(pal, "#C9A964FF", "#064396FF")
names(pal) <- c("Short", "Infinity1", "Nanopore\n(Polished)", "Nanopore\n(Unpolished)", "Infinity2")

report <- read.csv(here("06.mockcommunity/04.subsample/01.quast/transposed_report.tsv"), sep="\t", header=TRUE)
bin_assignments <- read.csv(here("06.mockcommunity/04.subsample/00.bins/bin_assignments.tsv"), sep="\t", header=TRUE)
names(bin_assignments) <- c("Assembly", "Species")
bin_assignments <- bin_assignments %>% mutate(Assembly=gsub(".fa", "", Assembly))
report <- merge(report, bin_assignments, by="Assembly")
report <- report %>% mutate(Method = gsub(".*_", "", gsub("-.*", "", Assembly)))
report <- report %>% mutate(Method=gsub("infinity", "Infinity", gsub("^polished", "Nanopore\n(Polished)", 
                                                                     gsub("^unpolished", "Nanopore\n(Unpolished)", Method))))

report <- report %>% mutate(Depth = gsub("_.*", "", gsub(".*subset_", "", Assembly)))
report <- report %>% mutate(DepthNumeric = as.numeric(gsub("Mb", "", gsub("Gb", "", Depth))))
report <- report %>% mutate(DepthNumeric = ifelse(DepthNumeric == 500, 0.5, DepthNumeric))
report <- report %>% mutate(LengthRelativeReference = Total.length.....0.bp. / Reference.length)

report <- report %>% mutate(Method = fct_relevel(Method, "Infinity1", "Infinity2", "Nanopore\n(Unpolished)", "Nanopore\n(Polished)")) 

ggplot(report, aes(x=DepthNumeric, y=N50)) + 
  geom_line(aes(color=Method)) + 
  geom_point(aes(color=Method)) + 
  scale_color_manual(values=pal[c(2,3,4,5)])+
  xlab("Depth (Gb)") + 
  ylab("Contig N50") +
  theme_bw() + 
  facet_wrap(~Species)

ggsave(here("00.outputs/mockcommunity/quastdepthbinning_N50.pdf"), w=6, h=5, dpi=300)


ggplot(report, aes(x=DepthNumeric, y=LengthRelativeReference)) + 
  geom_line(aes(color=Method)) + 
  geom_point(aes(color=Method)) + 
  scale_color_manual(values=pal[c(2,3,4,5)])+
  xlab("Depth (Gb)") + 
  ylab("Length Relative to Reference") +
  theme_bw() + 
  facet_wrap(~Species)

ggsave(here("00.outputs/mockcommunity/quastdepthbinning_length.pdf"), w=6, h=5, dpi=300)


ggplot(report, aes(x=DepthNumeric, y=Largest.contig)) + 
  geom_line(aes(color=Method)) + 
  geom_point(aes(color=Method)) + 
  scale_color_manual(values=pal[c(2,3,4,5)])+
  xlab("Depth (Gb)") + 
  ylab("Largest Contig") +
  theme_bw() + 
  facet_wrap(~Species)

ggsave(here("00.outputs/mockcommunity/quastdepthbinning_largest_contig.pdf"), w=6, h=5, dpi=300)


ggplot(report, aes(x=DepthNumeric, y=Genome.fraction....)) + 
  geom_line(aes(color=Method)) + 
  geom_point(aes(color=Method)) + 
  scale_color_manual(values=pal[c(2,3,4,5)])+
  xlab("Depth (Gb)") + 
  ylab("Genome Fraction") +
  theme_bw() + 
  facet_wrap(~Species)
ggsave(here("00.outputs/mockcommunity/quastdepthbinning_genomefraction.pdf"), w=6, h=5, dpi=300)

ggplot(report, aes(x=DepthNumeric, y=X..indels.per.100.kbp)) + 
  geom_line(aes(color=Method)) + 
  geom_point(aes(color=Method)) + 
  scale_color_manual(values=pal[c(2,3,4,5)])+
  xlab("Depth (Gb)") + 
  ylab("Indels per 100 Kbp") +
  theme_bw() + 
  facet_wrap(~Species)
ggsave(here("00.outputs/mockcommunity/quastdepthbinning_indels.pdf"), w=6, h=5, dpi=300)

ggplot(report %>% filter(Species != "Saccharomyces"), aes(x=DepthNumeric, y=X..indels.per.100.kbp)) + 
  geom_line(aes(color=Method)) + 
  geom_point(aes(color=Method)) + 
  scale_color_manual(values=pal[c(2,3,4,5)])+
  ylim(0,60) +
  xlab("Depth (Gb)") + 
  ylab("Indels per 100 Kbp") +
  theme_bw() + 
  facet_wrap(~Species)
ggsave(here("00.outputs/mockcommunity/quastdepthbinning_indelszoom.pdf"), w=6, h=5, dpi=300)

ggplot(report, aes(x=DepthNumeric, y=Duplication.ratio)) + 
  geom_line(aes(color=Method)) + 
  geom_point(aes(color=Method)) + 
  scale_color_manual(values=pal[c(2,3,4,5)])+
  xlab("Depth (Gb)") + 
  ylab("Duplication Ratio") +
  theme_bw() + 
  facet_wrap(~Species)

ggsave(here("00.outputs/mockcommunity/quastdepthbinning_duplicationratio.pdf"), w=6, h=5, dpi=300)

ggplot(report, aes(x=as.factor(DepthNumeric), y=N50, fill=Method, color=Method)) + 
  geom_point(aes(color=Method), position=position_dodge(width=0.75)) +
  geom_boxplot(outlier.shape=NA, alpha=0.8) + 
  scale_fill_manual(values=pal[c(2,3,4,5)])+
  scale_color_manual(values=pal[c(2,3,4,5)])+
  xlab("Depth (Gb)") + 
  ylab("N50 (bp)") +
  theme_bw()

ggsave(here("00.outputs/mockcommunity/quastdepthbinning_N50boxplot.pdf"), w=6, h=4, dpi=300)


#### QUAST ON 10x SUBSAMPLED DATA, ASSEMBLY #####

mytheme <- theme(panel.grid = element_blank(), text = element_text(size = 12))
# pal <- paletteer_d("nbapalettes::cavaliers")[c(1,2,3)]
# pal <- c(pal, "#C9A964FF", "#064396FF")
pal <- c("#8C1515", "#FFAF47", "#007C92", "#004552", "#E98300")
names(pal) <- c("Short", "Infinity1", "Nanopore\n(Polished)", "Nanopore\n(Unpolished)", "Infinity2")
report <- read.csv(here("06.mockcommunity/05.subsample_new/1.assemblies/02.quast_new/transposed_report.tsv"), sep="\t", header=TRUE)
report <- report %>% mutate(Method=gsub("short.*", "Short", gsub("nanopore.*", "Nanopore\n(Polished)", gsub("unpolished_nanopore.*", "Nanopore\n(Unpolished)", gsub("infinity1.*", "Infinity1", gsub("infinity2.*", "Infinity2", Assembly))))))
report <- report %>% mutate(Depth=as.numeric(gsub(".*Mb", "0.5", gsub("Gb", "", gsub("_.*", "", gsub(".*subset_", "", Assembly))))))
report <- report %>% mutate(Replicate=gsub(".*_", "", Assembly))
report <- report %>% mutate(Method = fct_relevel(Method, "Short", "Infinity1", "Infinity2", "Nanopore\n(Unpolished)", "Nanopore\n(Polished)")) 
report <- report %>% mutate(Total.length.....10000.bp.=Total.length.....10000.bp./1000000)
report <- report %>% mutate(N50=N50/1000000)

statsN50table <- report %>%
  group_by(Method, Depth) %>%
  summarise(
    MeanN50 = mean(N50),
    N50sd = sd(N50),
    n = n(),
    N50se = N50sd / sqrt(n), 
    MeanAssemblyLength = mean(Total.length.....10000.bp.), 
    AssemblyLengthsd = sd(Total.length.....10000.bp.),
    AssemblyLengthse = AssemblyLengthsd / sqrt(n), 
    MeanContigs = mean(X..contigs), 
    Contigssd = sd(X..contigs), 
    Contigsse = Contigssd / sqrt(n), 
    MeanGenomeFraction = mean(Genome.fraction....), 
    GenomeFractionsd = sd(Genome.fraction....), 
    GenomeFractionse = GenomeFractionsd / sqrt(n), 
    MeanDuplicationRatio = mean(Duplication.ratio), 
    DuplicationRatiosd = sd(Duplication.ratio), 
    DuplicationRatiose = DuplicationRatiosd / sqrt(n), 
    MeanMismatches = mean(X..mismatches.per.100.kbp), 
    Mismatchessd = sd(X..mismatches.per.100.kbp), 
    Mismatchesse = Mismatchessd / sqrt(n), 
    MeanIndels = mean(X..indels.per.100.kbp), 
    Indelssd = sd(X..indels.per.100.kbp), 
    Indelsse = Indelssd / sqrt(n)
  )

# plot_ribbon <- function(dataframe, ycolumn, title) {
#   ggplot(dataframe, aes(x=Depth, y=get(ycolumn))) + 
#     geom_smooth(stat='summary', alpha=0.1, aes(fill=Method, color=Method), fun.data = median_hilow, fun.args = list(conf.int = 0.5)) + 
#     scale_fill_manual(values=pal) + 
#     scale_color_manual(values=pal) + 
#     xlab("Depth (Gb)") + 
#     ylab(title) + 
#     theme_bw() + 
#     mytheme
# }
# 
# plot_points <- function(dataframe, ycolumn, title) {
#   ggplot(dataframe, aes(x=Depth, y=get(ycolumn))) + 
#     geom_point(aes(color=Method, fill=Method)) + 
#     scale_fill_manual(values=pal) + 
#     scale_color_manual(values=pal) + 
#     xlab("Depth (Gb)") + 
#     ylab(title) + 
#     theme_bw() + 
#     mytheme
# }

plot_curve <- function(dataframe, meancolumn, secolumn, ytitle) {
  ggplot(dataframe, aes(x=Depth, y=get(meancolumn), fill=Method, color=Method)) + 
    geom_errorbar(aes(ymin=get(meancolumn)-get(secolumn), ymax=get(meancolumn)+get(secolumn)), width=.2,
                  position=position_dodge(0.05)) + 
    geom_line() + 
    geom_point(size=3) +
    scale_fill_manual(values=pal) + 
    scale_color_manual(values=pal) +
    xlab("Depth (Gb)") +
    ylab(ytitle) +
    theme_bw() + 
    mytheme
}

#smallpal <- paletteer_d("nbapalettes::cavaliers")[c(1,2,3)]
smallpal <- c("#8C1515", "#007C92", "#E98300")
names(smallpal) <- c("Short", "Nanopore", "Infinity")

# filter stats table to only one infinity method, one nanopore
#statsN50table <- statsN50table %>% filter(Method != "Nanopore\n(Unpolished" & Method != "Infinity1") %>% mutate(Method=gsub("Infinity2", "Infinity", Method)) %>% mutate(Method=gsub("Nano.*", "Nanopore", Method))
a <- plot_curve(statsN50table, meancolumn = "MeanN50", secolumn = "N50se", "Contig N50 (Mb)")
a
ggsave(here("00.outputs/mockcommunity/quast_subsampled_N50.pdf"), w=4, h=3, dpi=300)

d <- plot_curve(statsN50table, meancolumn = "MeanAssemblyLength", secolumn = "AssemblyLengthse", "Assembly Length (Mb)")
d
# b <- plot_curve(statsN50table, meancolumn = "MeanContigs", secolumn = "Contigsse", "Number of Contigs")
# b
# ggsave(here("00.outputs/mockcommunity/quast_subsampled_contigs.pdf"), w=4, h=3, dpi=300)
# 
# e <- plot_curve(statsN50table, meancolumn = "MeanMismatches", secolumn = "Mismatchesse", "Mismatches per 100 kb")
# e
# ggsave(here("00.outputs/mockcommunity/quast_subsampled_mismatches.pdf"), w=4, h=3, dpi=300)

c <- plot_curve(statsN50table, meancolumn = "MeanIndels", secolumn = "Indelsse", "Indels per 100 kb")
c
ggsave(here("00.outputs/mockcommunity/quast_subsampled_indels.pdf"), w=4, h=3, dpi=300)

# testleg <- get_legend(c)
# plot_grid(d + theme(legend.position = "none"), 
#           b + theme(legend.position = "none"),
#           a + theme(legend.position = "none"),
#           c + theme(legend.position = "none"),
#           nrow=1, ncol=4,
#           rel_widths = c(1,1,1,1))
# 
# ggsave(here("00.outputs/mockcommunity/quast_subsampled_summary.pdf"), w=9.5, h=3, dpi=300)
# a + theme(legend.position = "bottom") 
# ggsave(here("00.outputs/mockcommunity/temp_for_legend.pdf"), w=9.5, h=3, dpi=300)

mockfig_top <- plot_grid(d + theme(legend.position = "none"),
                         a + theme(legend.position = "none"), 
                         c + theme(legend.position = "none"), 
                         nrow=1, labels = c("a", "b", "c"))
mockfig_top


#report <- report %>% filter(Method != "Nanopore\n(Unpolished" & Method != "Infinity1") %>% mutate(Method=gsub("Infinity2", "Infinity", Method)) %>% mutate(Method=gsub("Nano.*", "Nanopore", Method))

# plot_ribbon(report, ycolumn = "N50", "Contig N50 (bp)")
# plot_ribbon(report, ycolumn = "Total.length.....10000.bp.", "Assembly Length (bp)")
# plot_ribbon(report, ycolumn = "X..contigs", "Number of Contigs")
# plot_ribbon(report, ycolumn = "Duplication.ratio", "Duplication Ratio")
# plot_ribbon(report, ycolumn = "X..mismatches.per.100.kbp", "Mismatches per 100 kbp")
# plot_ribbon(report, ycolumn = "X..indels.per.100.kbp", "Indels per 100 kbp")
# 
# plot_points(report, ycolumn = "N50", "Contig N50 (bp)")
# a <- plot_points(report, ycolumn = "X..contigs.....0.bp.", "Number of Contigs >  0 bp")
# b <- plot_points(report, ycolumn = "X..contigs.....1000.bp.", "Number of Contigs >  1000 bp")
# c <- plot_points(report, ycolumn = "X..contigs.....5000.bp.", "Number of Contigs >  5000 bp")
# plot_grid(a + theme(legend.position = "none"), b + theme(legend.position = "none"), 
#           c + theme(legend.position = "none"), nrow = 1, ncol = 3)
# 
# 
# plot_points(report, ycolumn = "X..contigs", "Number of Contigs")
# plot_points(report, ycolumn = "X..mismatches.per.100.kbp", "Mismatches per 100 kbp")

#### QUAST ON 10x SUBSAMPLED DATA, BINNING #####
mytheme <- theme(panel.grid = element_blank(), text = element_text(size = 12))
mytheme_facet <- theme(panel.grid = element_blank(), text = element_text(size = 12), legend.position = "bottom", 
                       strip.background = element_rect(color = "white", fill="white"), strip.text = element_text(face = "italic", size=10))

# pal <- paletteer_d("nbapalettes::cavaliers")[c(1,2,3)]
# pal <- c(pal, "#C9A964FF", "#064396FF")
# names(pal) <- c("Short", "Infinity1", "Nanopore\n(Polished)", "Nanopore\n(Unpolished)", "Infinity2")
report <- read.csv(here("06.mockcommunity/05.subsample_new/2.bins/manual_quast/concatenated_report.tsv"), sep="\t", header=TRUE)
report <- report %>% mutate(Method=gsub("short.*", "Short", gsub("nanopore.*", "Nanopore\n(Polished)", gsub("unpolished_nanopore.*", "Nanopore\n(Unpolished)", gsub("infinity1.*", "Infinity1", gsub("infinity2.*", "Infinity2", Assembly))))))
report <- report %>% mutate(Depth=as.numeric(gsub(".*Mb", "0.5", gsub("Gb", "", gsub("_.*", "", gsub(".*subset_", "", Assembly))))))
report <- report %>% mutate(Replicate=gsub(".*_", "", gsub("-.*", "", Assembly)))
report <- report %>% mutate(Method = fct_relevel(Method, "Short", "Infinity1", "Infinity2", "Nanopore\n(Unpolished)", "Nanopore\n(Polished)")) 
report <- report %>% mutate(Species = gsub("_", " ", gsub("_complete", "", gsub(".*-", "", Assembly))))
report <- report %>% mutate(Grouping = paste(Method, Depth, Species, sep="-"))
report <- report %>% mutate(N50 = N50/1000000)
# facet and plot line with error bars
species_abbreviations <- data.frame(Species=c("Bacillus subtilis", "Enterococcus faecalis", "Escherichia coli", 
                                        "Listeria monocytogenes", "Pseudomonas aeruginosa", "Saccharomyces cerevisiae", 
                                        "Salmonella enterica", "Staphylococcus aureus"), 
                            Abbreviation=c("B. subtilis", "E. faecalis", "E. coli", "L. monocytogenes", 
                                           "P. aeruginosa", "S. cerevisiae", "S. enterica", "S. aureus"),
                            ExpectedLength=c(4.045, 2.845, 4.875, 2.992, 6.792, 12.1, 4.760, 2.730))

report <- merge(report, species_abbreviations, all.x = TRUE, by=c("Species"))
statstable <- report %>%
  group_by(Method, Depth, Abbreviation) %>%
  summarise(
    MeanN50 = mean(N50),
    N50sd = sd(N50),
    n = n(),
    N50se = N50sd / sqrt(n), 
    MeanAssemblyLength = mean(Total.length.....10000.bp.), 
    AssemblyLengthsd = sd(Total.length.....10000.bp.),
    AssemblyLengthse = AssemblyLengthsd / sqrt(n), 
    MeanContigs = mean(X..contigs), 
    Contigssd = sd(X..contigs), 
    Contigsse = Contigssd / sqrt(n), 
    MeanGenomeFraction = mean(Genome.fraction....), 
    GenomeFractionsd = sd(Genome.fraction....), 
    GenomeFractionse = GenomeFractionsd / sqrt(n), 
    MeanDuplicationRatio = mean(Duplication.ratio), 
    DuplicationRatiosd = sd(Duplication.ratio), 
    DuplicationRatiose = DuplicationRatiosd / sqrt(n), 
    MeanMismatches = mean(X..mismatches.per.100.kbp), 
    Mismatchessd = sd(X..mismatches.per.100.kbp), 
    Mismatchesse = Mismatchessd / sqrt(n), 
    MeanIndels = mean(X..indels.per.100.kbp), 
    Indelssd = sd(X..indels.per.100.kbp), 
    Indelsse = Indelssd / sqrt(n), 
    ExpectedLength = mean(ExpectedLength)
  )



plot_curve_facet <- function(dataframe, meancolumn, secolumn, ytitle) {
  ggplot(dataframe, aes(x=Depth, y=get(meancolumn), fill=Method, color=Method)) +
    geom_hline(aes(yintercept=ExpectedLength), color="grey", linetype="dashed") + 
    geom_errorbar(aes(ymin=get(meancolumn)-get(secolumn), ymax=get(meancolumn)+get(secolumn)), width=.2,
                  position=position_dodge(0.05)) + 
    geom_line() + 
    geom_point(size=3) +
    scale_fill_manual(values=pal) + 
    scale_color_manual(values=pal) +
    xlab("Depth (Gb)") +
    ylab(ytitle) +
    theme_bw() + 
    facet_wrap(~Abbreviation, nrow = 2) +
    scale_y_continuous(limits = c(0,7)) +
    scale_x_continuous(limits = c(0, 11), breaks=c(0,2.5, 5, 7.5, 10)) +
    mytheme_facet
}

plot_curve_facet8 <- function(dataframe, meancolumn, secolumn, ytitle) {
  ggplot(dataframe, aes(x=Depth, y=get(meancolumn), fill=Method, color=Method)) + 
    geom_errorbar(aes(ymin=get(meancolumn)-get(secolumn), ymax=get(meancolumn)+get(secolumn)), width=.2,
                  position=position_dodge(0.05)) + 
    geom_line() + 
    geom_point(size=3) +
    scale_fill_manual(values=pal) + 
    scale_color_manual(values=pal) +
    xlab("Depth (Gb)") +
    ylab(ytitle) +
    theme_bw() + 
    facet_wrap(~Abbreviation, nrow = 1) +
    scale_x_continuous(limits = c(0, 11), breaks=c(0,2.5, 5, 7.5, 10)) +
    mytheme_facet
}

# plot_curve <- function(dataframe, meancolumn, secolumn, ytitle) {
#   ggplot(dataframe, aes(x=Depth, y=get(meancolumn), fill=Method, color=Method)) + 
#     geom_line() + 
#     geom_point(size=3) +
#     scale_fill_manual(values=pal) + 
#     scale_color_manual(values=pal) +
#     geom_errorbar(aes(ymin=get(meancolumn)-get(secolumn), ymax=get(meancolumn)+get(secolumn)), width=.2,
#                   position=position_dodge(0.05)) + 
#     xlab("Depth (Gb)") +
#     ylab(ytitle) +
#     theme_bw() + 
#     scale_x_continuous(limits = c(0, 11), breaks=c(0,2.5, 5, 7.5, 10)) +
#     mytheme + 
#     theme(legend.position = "none")
# }

mockfig_bottom <- plot_curve_facet(statstable, "MeanN50", "N50se", "Contig N50 (Mb)")
mockfig_bottom

# mockfig_bottom8 <- plot_curve_facet8(statstable, "MeanN50", "N50se", "Contig N50 (Mb)")
# mockfig_bottom8

plot_grid(mockfig_top, mockfig_bottom, labels = c("", "d"), nrow=2, ncol=1, rel_heights = c(1, 1.6))
ggsave(here("00.outputs/final/mockcommunity_justN50.pdf"), w=8, h=7, dpi=300)
ggsave(here("00.outputs/final/mockcommunity_justN50.jpg"), w=8, h=7, dpi=300)


# plot_curve_facet(statstable, "MeanGenomeFraction", "GenomeFractionse", "Genome Fraction")
# ggsave(here("00.outputs/mockcommunity/binning/quast_genome_fraction.pdf"), w=9, h=5, dpi=300)
# 
# plot_curve_facet(statstable, "MeanAssemblyLength", "AssemblyLengthse", "Assembly Length (bp)")
# ggsave(here("00.outputs/mockcommunity/binning/quast_assembly_length.pdf"), w=9, h=5, dpi=300)
# 
# plot_curve_facet(statstable, "MeanIndels", "Indelsse", "Indels per 100 kb")
# ggsave(here("00.outputs/mockcommunity/binning/quast_indels.pdf"), w=9, h=5, dpi=300)
# 
# plot_curve_facet(statstable, "MeanMismatches", "Mismatchesse", "Mismatches per 100 kb")
# ggsave(here("00.outputs/mockcommunity/binning/quast_mistmatches.pdf"), w=9, h=5, dpi=300)
# 
# plot_curve_facet(statstable, "MeanContigs", "Contigsse", "Number of Contigs")
# ggsave(here("00.outputs/mockcommunity/binning/quast_contig_number.pdf"), w=9, h=5, dpi=300)
# 
# 
# plot_curve(statstable %>% filter(Species == "Saccharomyces cerevisiae"), "MeanN50", "N50se", "Contig N50 (bp)")
# ggsave(here("00.outputs/mockcommunity/binning/quast_contig_N50_yeast.pdf"), w=4, h=3.2, dpi=300)
# 
# # plotting number of bins per species per method
# ggplot(statstable, aes(x=Depth, y=n, fill=Method, color=Method)) + 
#     geom_line() + 
#     geom_point(size=3) +
#     scale_fill_manual(values=pal) + 
#     scale_color_manual(values=pal) +
#     xlab("Depth (Gb)") +
#     ylab("Number of Genomes") +
#     theme_bw() + 
#     facet_wrap(~Species, nrow = 2) +
#     scale_x_continuous(limits = c(0, 11), breaks=c(0,2.5, 5, 7.5, 10)) +
#     mytheme_facet
# 
# ggsave(here("00.outputs/mockcommunity/binning/quast_genome_number.pdf"), w=9, h=5, dpi=300)
# 
