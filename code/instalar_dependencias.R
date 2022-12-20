#.libPaths(c(.libPaths(), temp <- tempdir()))
path_discal <- "~/deps_R"

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = path_discal)

install.packages("igraph")
install.packages("xtable")

BiocManager::install("STRINGdb", lib = path_discal)
BiocManager::install("linkcomm", lib = path_discal)
BiocManager::install("clusterProfiler", lib = path_discal)
BiocManager::install("org.Hs.eg.db", lib = path_discal)

library(igraph, lib.loc = path_discal)
library(STRINGdb, lib.loc = path_discal)
library(linkcomm, lib.loc = path_discal)
library(clusterProfiler, lib.loc = path_discal)
library(org.Hs.eg.db, lib.loc = path_discal)
library(xtable, lib.loc = path_discal)
