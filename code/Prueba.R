install.packages('knitr',lib='~/software')

myPaths <- .libPaths()

myPaths <- c(myPaths, '~/software')

.libPaths(myPaths)  # add new path