library(here)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(paletteer)
library(cowplot)
library(raster)
library(xts)
library(forcats)
library(scales)

#### QUAST ON 10x SUBSAMPLED DATA, ASSEMBLY #####

mytheme <- theme(panel.grid = element_blank(), text = element_text(size = 12))
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

smallpal <- c("#8C1515", "#007C92", "#E98300")
names(smallpal) <- c("Short", "Nanopore", "Infinity")


a <- plot_curve(statsN50table %>% filter(Method != "Nanopore\n(Unpolished)", Method != "Infinity1"), meancolumn = "MeanN50", secolumn = "N50se", "Contig N50 (Mb)")
a
d <- plot_curve(statsN50table%>% filter(Method != "Nanopore\n(Unpolished)", Method != "Infinity1"), meancolumn = "MeanAssemblyLength", secolumn = "AssemblyLengthse", "Assembly Length (Mb)")
d
c <- plot_curve(statsN50table%>% filter(Method != "Infinity1"), meancolumn = "MeanIndels", secolumn = "Indelsse", "Indels per 100 kb")
c

mockfig_top <- plot_grid(d + theme(legend.position = "none"),
                         a + theme(legend.position = "none"), 
                         c + theme(legend.position = "none"), 
                         nrow=1, labels = c("a", "b", "c"))
mockfig_top

#### QUAST ON 10x SUBSAMPLED DATA, BINNING #####
mytheme <- theme(panel.grid = element_blank(), text = element_text(size = 12))
mytheme_facet <- theme(panel.grid = element_blank(), text = element_text(size = 12), legend.position = "bottom", legend.title=element_blank(),
                       strip.background = element_rect(color = "white", fill="white"), strip.text = element_text(face = "italic", size=10))

pal <- c("#8C1515", "#007C92", "#004552", "#E98300")
names(pal) <- c("Short", "ONT\n(Polished)", "ONT\n(Unpolished)", "ICLR")

report <- read.csv(here("06.mockcommunity/05.subsample_new/2.bins/manual_quast/concatenated_report.tsv"), sep="\t", header=TRUE)
report <- report %>% mutate(Method=gsub("short.*", "Short", gsub("nanopore.*", "ONT\n(Polished)", gsub("unpolished_nanopore.*", "ONT\n(Unpolished)", gsub("infinity1.*", "ICLR1", gsub("infinity2.*", "ICLR", Assembly))))))
report <- report %>% mutate(Depth=as.numeric(gsub(".*Mb", "0.5", gsub("Gb", "", gsub("_.*", "", gsub(".*subset_", "", Assembly))))))
report <- report %>% mutate(Replicate=gsub(".*_", "", gsub("-.*", "", Assembly)))
report <- report %>% mutate(Method = fct_relevel(Method, "Short", "ICLR1", "ICLR", "ONT\n(Unpolished)", "ONT\n(Polished)")) 
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

statstable$Abbreviation <- factor(statstable$Abbreviation, levels=c("B. subtilis", "E. coli", "E. faecalis", "L. monocytogenes", "P. aeruginosa", "S. aureus", "S. enterica", "S. cerevisiae"))                           # Replicate data
statstable$Method <- factor(statstable$Method, levels=c("Short", "ONT\n(Polished)", "ONT\n(Unpolished)", "ICLR", "ICLR1"))

plot_curve_facet <- function(dataframe, meancolumn, secolumn, ytitle) {
  ggplot(dataframe, aes(x=Depth, y=get(meancolumn), fill=Method, color=Method)) +
    geom_hline(aes(yintercept=ExpectedLength, linetype="Expected\nGenome Size"), color="grey") + 
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
    #scale_y_log10(labels = comma) + 
    scale_x_continuous(limits = c(0, 11), breaks=c(0,2.5, 5, 7.5, 10)) +
    mytheme_facet + 
    scale_linetype_manual(name = "", values="dashed") 
}

mockfig_bottom <- plot_curve_facet(statstable%>% filter(Method != "ONT\n(Unpolished)", Method != "ICLR1"), "MeanN50", "N50se", "Contig N50 (Mb)")
mockfig_bottom

plot_grid(mockfig_top, mockfig_bottom, labels = c("", "d"), nrow=2, ncol=1, rel_heights = c(1, 1.6))



ggsave(here("00.outputs/final/mockcommunity_justN50.pdf"), w=8, h=7, dpi=300)
ggsave(here("00.outputs/final/mockcommunity_justN50.jpg"), w=8, h=7, dpi=300)



