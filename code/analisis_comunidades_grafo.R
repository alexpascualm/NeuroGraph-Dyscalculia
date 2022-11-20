library(igraph)

################# Grafos con iGraph.

## Limpiar terminal de R
cat("\014")

## Borrar variables de R
rm(list=ls())

nodes <- read.table("./string_node_degrees.tsv", sep = '\t',header = TRUE, as.is=T)


links <- read.table("./string_interactions.tsv", sep = '\t',header = TRUE, as.is=T)

## Crea un nuevo grafo a partir de los nodos y de las aristas proporcionadas. La red será dirigida.
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)


lo <- layout_with_kk(net) # create a layout
lo <- norm_coords(lo, ymin=-1, ymax=1, xmin=-1, xmax=1)

## Representa graficamente
plot(net, edge.arrow.size=0.2, edge.curved=0.1, vertex.size=15, vertex.color="white", vertex.frame.color="black", vertex.label.color="black",vertex.label.cex=0.4, layout=lo)


# Detección de comunidades por clustering
community <- cluster_edge_betweenness(net)


dendPlot(community)

plot(community, net)



length(community)
membership(community)

##Detección de comunidades mediante optimización voraz.

cfg <- cluster_fast_greedy(as.undirected(net))

dendPlot(cfg)

plot(cfg, as.undirected(net))

length(cfg)
membership(cfg)


#Detección de comunidades mediante propagación de etiquetas.
clp <- cluster_label_prop(net)
dendPlot(clp)
plot(clp, net)

length(clp)
membership(clp)


