# Plot repo sizes

library(data.table)
library(ggplot2)

dat <- fread('lines.csv')

dat[, clines := cumsum(lines), by = ver]
dat[, rank := 1:.N, by = ver]

ggplot(dat, aes(

