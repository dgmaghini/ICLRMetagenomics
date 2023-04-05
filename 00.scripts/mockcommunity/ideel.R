library(ggplot2)
library(here)
library(tidyverse)
library(paletteer)

pal <- paletteer_d("nbapalettes::cavaliers")[c(1,2,3)]
pal <- c(pal, "#C9A964FF")
names(pal) <- c("Short", "Infinity", "Nanopore (Polished)", "Nanopore (Unpolished)")

infinity <- read.table(here("06.mockcommunity/01.ideel/lengths/infinity.data"), sep="\t")
nano_polished <- read.table(here("06.mockcommunity/01.ideel/lengths/nanopore.data"), sep="\t")
nano_unpolished <- read.table(here("06.mockcommunity/01.ideel/lengths/nanopore_unpolished.data"), sep="\t")
names(infinity) <- c("qlen", "slen")
names(nano_polished) <- c("qlen", "slen")
names(nano_unpolished) <- c("qlen", "slen")

infinity <- mutate(infinity, Method="Infinity")
nano_polished <- mutate(nano_polished, Method="Nanopore (Polished)")
nano_unpolished <- mutate(nano_unpolished, Method="Nanopore (Unpolished)")

all <- rbind(infinity, nano_polished, nano_unpolished)

pseudogenes <- sum(all$qlen / all$slen < 0.9)
print(paste0('Encountered genes < 0.9 reference length: ', pseudogenes))

all <- mutate(all, QueryLengthPortion=qlen/slen)

p <- ggplot(all, aes(x=qlen/slen)) + 
  geom_histogram(aes(fill=Method), color='grey25', bins=20, position="dodge") +
  xlab('query length / hit length') +
  ylab('frequency') +
  # scale_y_log10() +
  scale_x_continuous(limits=c(0, 1.3)) +
  theme_minimal() + 
  scale_fill_manual(values=pal[c(2,3,4)]) + 
  theme(text=element_text(size=12))
p
ggsave(here("00.outputs/mockcommunity/ideel.pdf"), w=8, h=3, dpi=300)
ggsave(here("00.outputs/mockcommunity/ideel.jpg"), w=8, h=3, dpi=300)


ggplot(all, aes(x=qlen/slen)) + 
  geom_density(aes(fill=Method, color=Method), alpha=0.3) +
  xlab('query length / hit length') +
  #scale_y_log10() +
  scale_x_continuous(limits=c(0, 1.3)) +
  theme_minimal()


ggplot(all, aes(x=Method, y=qlen/slen)) + 
  geom_violin(aes(fill=Method, color=Method), alpha=0.3) +
  xlab('query length / hit length') +
  # scale_y_log10() +
  scale_y_continuous(limits=c(0.975, 1.025)) +
  theme_minimal()

ggplot(all, aes(x=Method, y=qlen/slen)) + 
  geom_boxplot(aes(fill=Method, color=Method), alpha=0.3, outlier.shape=NA) +
  xlab('query length / hit length') +
  # scale_y_log10() +
  scale_y_continuous(limits=c(0.975, 1.025)) +
  theme_minimal()



inf <- all %>% filter(Method == "Infinity")
median(inf$QueryLengthPortion)

nanou <- all %>% filter(Method == "Nanopore (Unpolished)")
median(nanou$QueryLengthPortion)

nanop <- all %>% filter(Method == "Nanopore (Polished)") 
median(nanop$QueryLengthPortion)
