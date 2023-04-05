library(here)
library(tidyverse)
library(reshape2)

nanopore_length <- read.table(here("08.mock_stats/reads/nanopore_length.tsv"), header=FALSE, sep=" ")
infinity1_length <- read.table(here("08.mock_stats/reads/infinity1_length.tsv"), header=FALSE, sep=" ")
infinity2_length <- read.table(here("08.mock_stats/reads/infinity2_length.tsv"), header=FALSE, sep=" ")

names(nanopore_length) <- c("Length", "Count")
names(infinity1_length) <- c("Length", "Count")
names(infinity2_length) <- c("Length", "Count")

nanopore_length <- nanopore_length %>% mutate(Method = "Nanopore")
infinity1_length <- infinity1_length %>% mutate(Method = "Infinity")
infinity2_length <- infinity2_length %>% mutate(Method = "Infinity2")

lengthdf <- rbind(nanopore_length, infinity1_length, infinity2_length)

lengthdf <- lengthdf %>% mutate(bin = floor(Length/500) * 500)

lengthdf <- lengthdf %>% mutate(totalBases = Count * Length)



pal <- paletteer_d("nbapalettes::cavaliers")[c(1,2,3)]
pal <- c(pal, "#C9A964FF", "#064396FF")
names(pal) <- c("Short", "Infinity", "Nanopore", "NanoporeUnpolished", "Infinity2")


lengthdf <- lengthdf %>% group_by(Method, bin) %>% 
  summarise(
    ReadCount = sum(Count), 
    BaseCount = sum(totalBases)
  )

readcount <- ggplot(lengthdf %>% filter(Method != "Infinity2"), aes(x=bin, y=ReadCount, fill=Method, color=Method)) + 
  geom_col(position=position_dodge2(width=80, preserve="single")) + 
  theme_bw() + 
  scale_x_continuous(limits = c(0,160000)) + 
  scale_fill_manual(values=pal, limits = force) + 
  scale_color_manual(values=pal, limits = force) + 
  xlab("Read Length (bp)") + 
  ylab("Number of Reads") + 
  scale_y_log10() + 
  theme(panel.grid = element_blank(), 
        text = element_text(size=9), legend.title=element_blank())
readcount

basecount <- ggplot(lengthdf %>% filter(Method != "Infinity2"), aes(x=bin, y=BaseCount, fill=Method, color=Method)) + 
  geom_col(position=position_dodge2(width=80, preserve="single")) + 
  theme_bw() + 
  scale_x_continuous(limits = c(0,160000)) + 
  scale_fill_manual(values=pal, limits = force) + 
  scale_color_manual(values=pal, limits = force) + 
  xlab("Read Length (bp)") + 
  ylab("Number of Bases") + 
  scale_y_log10() + 
  theme(panel.grid = element_blank(), 
        text = element_text(size=9), legend.title=element_blank())
basecount

leg <- get_legend(readcount)
main <- plot_grid(readcount + theme(legend.position = "NA"), 
          basecount + theme(legend.position = "NA"),
          leg,
          nrow = 1, ncol=3, labels = c("a", "b", ""), rel_widths = c(1,1,0.3))
main


ggsave(here("00.outputs/final/suppfigures/mockcommunity_readlengths.pdf"), dpi=300, w=6, h=2.7)
ggsave(here("00.outputs/final/suppfigures/mockcommunity_readlengths.jpg"), dpi=300, w=6, h=2.7)

