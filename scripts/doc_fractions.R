#!/usr/bin/Rscript

library(ggplot2)
library(plyr)

doc_filename <- commandArgs(trailingOnly=T)[1]
if (is.na(doc_filename) || !file.exists(doc_filename)) {
  print('Usage: Rscript <script.R> <doc file>')
  quit('no', 1)
}

doc <- read.table(doc_filename, header=T, sep='\t')
doc <- ddply(doc, .(sample), transform, cumulative=1-cumsum(fraction))

p <- qplot(count, cumulative, data=doc, color=sample, geom='line') +
  geom_vline(xintercept = 20, color='red', alpha=.5) +
  geom_hline(yintercept = .80, color='red', alpha=.5) +
  scale_x_log10(breaks=c(seq(1, 9), seq(10, 30, 5), seq(40, 60, 10), 100, 200)) +
  scale_y_continuous(breaks=seq(0, 1, .1)) +
  labs(title='Cumulative target region coverage', x='Depth', y='Fraction of target bases >= depth') +
  guides(col = guide_legend(nrow = 16))
ggsave(plot=p, filename=paste0(doc_filename, '.png'), width=15, height=10, dpi=75)
