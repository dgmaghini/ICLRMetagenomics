library(here)
library(dplyr)
library(reshape2)
library(ggplot2)
library(cowplot)

nanopore_length <- read.table(here("08.mock_stats/reads/nanopore_length.tsv"), header=FALSE, sep=" ")
infinity2_length <- read.table(here("08.mock_stats/reads/infinity2_length.tsv"), header=FALSE, sep=" ")

names(nanopore_length) <- c("Length", "Count")
names(infinity2_length) <- c("Length", "Count")

nanopore_length <- nanopore_length %>% mutate(Method = "ONT")
infinity2_length <- infinity2_length %>% mutate(Method = "ICLR")

lengthdf <- rbind(nanopore_length, infinity2_length)
lengthdf <- lengthdf %>% mutate(bin = floor(Length/500) * 500)
lengthdf <- lengthdf %>% mutate(totalBases = Count * Length)

pal <- c("#8C1515", "#E98300", "#007C92") 
names(pal) <- c("Short", "ICLR", "ONT")

lengthdf <- lengthdf %>% group_by(Method, bin) %>% 
  summarise(
    ReadCount = sum(Count), 
    BaseCount = sum(totalBases)
  )

basecount2 <- ggplot(lengthdf, aes(x=bin/1000, y=BaseCount, fill=Method, color=Method)) + 
  geom_area(alpha = 0.5, position = "identity")+
  #geom_line()+
  theme_bw() + 
  scale_x_continuous(limits = c(0,100)) + 
  scale_fill_manual(values=pal, limits = force) + 
  scale_color_manual(values=pal, limits = force) + 
  xlab("Read Length (kbp)") + 
  ylab("Number of Bases") + 
  scale_y_log10() + 
  theme(panel.grid = element_blank(), 
        text = element_text(size=9), legend.title=element_blank())
basecount2

readcount2 <- ggplot(lengthdf, aes(x=bin/1000, y=ReadCount, fill=Method, color=Method)) + 
  geom_area(alpha = 0.5, position = "identity")+
  theme_bw() + 
  scale_x_continuous(limits = c(0,100)) + 
  scale_fill_manual(values=pal, limits = force) + 
  scale_color_manual(values=pal, limits = force) + 
  xlab("Read Length (kbp)") + 
  ylab("Number of Bases") + 
  scale_y_log10() + 
  theme(panel.grid = element_blank(), 
        text = element_text(size=9), legend.title=element_blank())
readcount2

leg <- get_legend(readcount2)
main <- plot_grid(readcount2 + theme(legend.position = "NA", panel.background = element_rect(color = "black", linewidth = 0.8), 
                                     axis.text = element_text(color = "black", size = 8), axis.ticks = element_line(color = "black"), 
                                     axis.title = element_text(color = "black")), 
                  basecount2 + theme(legend.position = "NA", panel.background = element_rect(color = "black", linewidth = 0.8), 
                                     axis.text = element_text(color = "black", size = 8), axis.ticks = element_line(color = "black"), 
                                     axis.title = element_text(color = "black")),
                  leg,
                  nrow = 1, ncol=3, labels = c("a", "b", ""), rel_widths = c(1,1,0.3))
main

ggsave(here("00.outputs/final/suppfigures/mockcommunity_readlengths_density.pdf"), dpi=300, w=5.5, h=2.4)
ggsave(here("00.outputs/final/suppfigures/mockcommunity_readlengths_density.jpg"), dpi=300, w=5.5, h=2.4)


