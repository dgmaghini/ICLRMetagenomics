library(here)
library(tidyverse)
library(reshape2)

nanopore_length <- read.table(here("07.human_stats/reads/nanopore_D1_length.tsv"), header=FALSE, sep=" ")
infinity_length <- read.table(here("07.human_stats/reads/infinity_D1_length.tsv"), header=FALSE, sep=" ")

names(nanopore_length) <- c("Length", "Count")
names(infinity_length) <- c("Length", "Count")

nanopore_length <- nanopore_length %>% mutate(Method = "Nanopore")
infinity_length <- infinity_length %>% mutate(Method = "Infinity")


lengthdf <- rbind(nanopore_length, infinity_length)

lengthdf <- lengthdf %>% mutate(bin = floor(Length/500) * 500)

lengthdf <- lengthdf %>% mutate(totalBases = Count * Length)



pal <- paletteer_d("nbapalettes::cavaliers")[c(1,2,3)]
names(pal) <- c("Short", "Infinity", "Nanopore")


lengthdf <- lengthdf %>% group_by(Method, bin) %>% 
  summarise(
    ReadCount = sum(Count), 
    BaseCount = sum(totalBases)
  )

readcount <- ggplot(lengthdf, aes(x=bin, y=ReadCount, fill=Method, color=Method)) + 
  geom_col(position=position_dodge2(width=80, preserve="single")) + 
  theme_bw() + 
  scale_x_continuous(limits = c(0,122000)) + 
  scale_fill_manual(values=pal, limits = force) + 
  scale_color_manual(values=pal, limits = force) + 
  xlab("Read Length (bp)") + 
  ylab("Number of Reads") + 
  scale_y_log10() + 
  theme(panel.grid = element_blank(), 
        text = element_text(size=9), legend.title=element_blank())
readcount

basecount <- ggplot(lengthdf, aes(x=bin, y=BaseCount, fill=Method, color=Method)) + 
  geom_col(position=position_dodge2(width=80, preserve="single")) + 
  theme_bw() + 
  scale_x_continuous(limits = c(0,122000)) + 
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


ggsave(here("00.outputs/final/suppfigures/human_readlengths.pdf"), dpi=300, w=6, h=2.7)
ggsave(here("00.outputs/final/suppfigures/human_readlengths_notlog.jpg"), dpi=300, w=6, h=2.7)

