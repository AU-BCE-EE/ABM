# Runs all Rmd files in dev-tests

library(rmarkdown)

# Find Rmd files
ff <- list.files(pattern = '\\.Rmd$', recursive = TRUE, full.names = TRUE)

# Set sink to a log file
sink(file('log.txt', 'wt'), type = 'message')

# Render all of them
for (i in ff) {
  render(i)
}

sink()
