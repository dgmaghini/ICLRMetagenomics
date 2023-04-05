library(here)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(paletteer)
library(cowplot)

# Read in manual binning via alignment table
manual_df <- read.table(here("06.mockcommunity/05.subsample_new/2.bins/binning_manual/contig_lengths_all.tsv"), header=FALSE, sep="\t")
names(manual_df) <- c("Contig", "ContigAssignment", "ContigLength", "Sample")
manual_df <- manual_df %>% dplyr::select(Contig, ContigAssignment, Sample)
# Read in binning from DAS_tool table
das_df <- read.table(here("06.mockcommunity/05.subsample_new/2.bins/binning_das_tool/contig_assignments_all.tsv"), header=FALSE, sep="\t")
names(das_df) <- c("DASToolBin", "Contig", "Sample")

# Read in the contig lengths table (generated from assembly .fa.fai files)
lengths_df <- read_table(here("06.mockcommunity/05.subsample_new/1.assemblies/contig_lengths.tsv"), col_names=FALSE)
names(lengths_df) <- c("Contig", "ContigLength", "Sample")

# merge dataframes
df <- merge(das_df, manual_df, by=c("Contig", "Sample"), all.x = TRUE)
df <- df %>% mutate(ContigAssignment = ifelse(is.na(ContigAssignment), "None", ContigAssignment))
df <- merge(df, lengths_df, by=c("Contig", "Sample"), all.x = TRUE)
names(df)
df <- df %>% filter(!grepl("unbinned", DASToolBin))

# Calculate total length per bin
binbreakdown <- df %>%
  group_by(Sample, DASToolBin) %>%
  summarise(
    TotalBinLength = sum(ContigLength)
  )

# Calculate total length per organism per bin
binbreakdown2 <- df %>%
  group_by(Sample, DASToolBin, ContigAssignment) %>%
  summarise(
    TotalTaxonLength = sum(ContigLength)
  )

# merge length dataframes
binbreakdown <- merge(binbreakdown, binbreakdown2, all.y=TRUE)
binbreakdown <- binbreakdown %>% mutate(TaxonPortion=TotalTaxonLength/TotalBinLength)

# filter bin table to only have top taxon per bin
bin_toptaxon <- binbreakdown %>% group_by(Sample, DASToolBin) %>%
  filter(TaxonPortion == max(TaxonPortion))

# shows the number of bins representing each organism for each sample
bin_check <- bin_toptaxon %>% group_by(Sample, ContigAssignment) %>%
  summarise(count=n())

bin_toptaxon <- bin_toptaxon %>% mutate(MisassignedPortion = (1-TaxonPortion))

pal <- paletteer_d("nbapalettes::cavaliers")[c(1,2,3)]
pal <- c(pal, "#C9A964FF", "#064396FF")
names(pal) <- c("Short", "Infinity1", "Nanopore\n(Polished)", "Nanopore\n(Unpolished)", "Infinity2")
bin_toptaxon <- bin_toptaxon %>% mutate(Method=gsub("short.*", "Short", gsub("nanopore.*", "Nanopore\n(Polished)", gsub("unpolished_nanopore.*", "Nanopore\n(Unpolished)", gsub("infinity1.*", "Infinity1", gsub("infinity2.*", "Infinity2", Sample))))))
bin_toptaxon <- bin_toptaxon %>% mutate(Depth=as.numeric(gsub(".*Mb", "0.5", gsub("Gb", "", gsub("_.*", "", gsub(".*subset_", "", Sample))))))
bin_toptaxon <- bin_toptaxon %>% mutate(Replicate=gsub(".*_", "", gsub("-.*", "", Sample)))
bin_toptaxon <- bin_toptaxon %>% mutate(Method = fct_relevel(Method, "Short", "Infinity1", "Infinity2", "Nanopore\n(Unpolished)", "Nanopore\n(Polished)")) 
bin_toptaxon <- bin_toptaxon %>% mutate(Species = gsub("_", " ", gsub("_complete", "", gsub(".*-", "", ContigAssignment))))

statstable <- bin_toptaxon %>%
  group_by(Method, Depth, Species) %>%
  summarise(
    MeanMisassigned = mean(MisassignedPortion),
    Misassignedsd = sd(MisassignedPortion),
    n = n(),
    Misassignedse = Misassignedsd / sqrt(n),
  )

plot_curve_facet <- function(dataframe, meancolumn, secolumn, ytitle) {
  ggplot(dataframe, aes(x=Depth, y=get(meancolumn), fill=Method, color=Method)) +
    geom_line() +
    geom_point(size=3) +
    scale_fill_manual(values=pal) +
    scale_color_manual(values=pal) +
    geom_errorbar(aes(ymin=get(meancolumn)-get(secolumn), ymax=get(meancolumn)+get(secolumn)), width=.2,
                  position=position_dodge(0.05)) +
    xlab("Depth (Gb)") +
    ylab(ytitle) +
    theme_bw() +
    facet_wrap(~Species, nrow = 2) +
    scale_x_continuous(limits = c(0, 11), breaks=c(0,2.5, 5, 7.5, 10)) +
    mytheme_facet
}

plot_curve_facet(statstable, "MeanMisassigned", "Misassignedse", "non-Majority Species Portion per Bin")
ggsave(here("00.outputs/mockcommunity/binning/das_tool_misassignment_fraction.pdf"), w=9, h=5, dpi=300)
