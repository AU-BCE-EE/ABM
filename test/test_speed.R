rm(list = ls())

ff <- list.files('../R', full.names = T)
for (i in ff) source(i)

library(profvis)

profvis({
  
  replicate(abm(10000), 10)
  
  
})
