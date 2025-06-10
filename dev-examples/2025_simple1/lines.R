# Plot repo sizes

library(data.table)
library(ggplot2)

dat <- fread('lines.csv')

dat[, clines := cumsum(lines), by = ver]
dat[, rank := 1:.N, by = ver]

ggplot(dat, aes(rank, clines, colour = ver)) +
  geom_point() +
  geom_line(alpha = 0.3) +
  theme_bw() +
  ylim(0, NA) +
  theme(legend.position = 'top') +
  labs(colour = 'ABM version', x = 'File size rank', y = 'Cumulative lines of code')
ggsave('lines.png', height = 8/2, width = 16/2)
