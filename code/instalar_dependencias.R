dir.create(paste0(Sys.getenv("R_LIBS_USER"),"/discalculia"), recursive = TRUE)  # create personal library
.libPaths(paste0(Sys.getenv("R_LIBS_USER"),"/discalculia"))  # add to the path

path_discalculia <- paste0(Sys.getenv("R_LIBS_USER"),"/discalculia")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = path_discalculia)


if (!require("igraph", quietly = TRUE))
  install.packages("igraph", lib = path_discalculia)

if (!require("xtable", quietly = TRUE))
  install.packages("xtable", lib = path_discalculia)


if (!require("STRINGdb", quietly = TRUE))
  BiocManager::install("STRINGdb", lib = path_discalculia)

if (!require("linkcomm", quietly = TRUE))
  BiocManager::install("linkcomm", lib = path_discalculia)

if (!require("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db", lib = path_discalculia)


if (!require("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler", lib = path_discalculia)


