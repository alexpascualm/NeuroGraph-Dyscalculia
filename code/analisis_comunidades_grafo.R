library(igraph)
library(STRINGdb)
library(linkcomm)
library(clusterProfiler)
library(org.Hs.eg.db)
library(xtable)

## Limpiar terminal de R
cat("\014")

## Borrar variables de R
rm(list=ls())

nodes <- read.table("./string_node_degrees.tsv", sep = '\t',header = TRUE, as.is=T)


links <- read.table("./string_interactions.tsv", sep = '\t',header = TRUE, as.is=T)


##### Grafos con iGraph ####

## Crea un nuevo grafo a partir de los nodos y de las aristas proporcionadas (La red será dirigida)

net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)


lo <- layout_with_kk(net) # Creamos la disposición (layout) de los nodos
lo <- norm_coords(lo, ymin=-1, ymax=1, xmin=-1, xmax=1)

## Representa graficamente
plot(net, edge.arrow.size=0.2, edge.curved=0.1, vertex.size=15, vertex.color="white", vertex.frame.color="black", vertex.label.color="black",vertex.label.cex=0.4, layout=lo)


## Detección de comunidades por clustering ##
community <- cluster_edge_betweenness(net)

dendPlot(community) 

plot(community, net, edge.arrow.size=0.05, edge.curved=0.1, vertex.size=13, vertex.color="white", vertex.frame.color="black", vertex.label.color="black",vertex.label.cex=0.4, layout=lo)

length(community)
membership(community)


## Detección de comunidades mediante optimización voraz ##

cfg <- cluster_fast_greedy(as.undirected(net))

dendPlot(cfg)

plot(cfg, as.undirected(net), vertex.size=13, vertex.color="white", vertex.frame.color="black", vertex.label.color="black",vertex.label.cex=0.4, layout=lo)

length(cfg)
membership(cfg)


#### NetworkPropagation de la red ####
string_db <- STRINGdb$new( version="11", species=9606, score_threshold=400, input_directory="" )
string.network <- string_db$get_graph()

hits <- nodes$X9606.ENSP00000442954 # Nos quedamos con los String_ID de nuestro conjunto de genes
hits.network <- string_db$get_subnetwork(hits) # Creamos una network usando StringDb

first.neigh <- (neighbors(graph = string.network, v = V(hits.network)$name, mode = "all"))$name # Encontramos una serie de nodos vecinos
hits.network <- string_db$get_subnetwork(unique(c(V(hits.network)$name, first.neigh))) # Unimos la red de nodos originales con la red de vecinos

DFNetwork <- igraph::as_data_frame(hits.network)


## Analisis por linkcomm ##
DC_lc <- getLinkCommunities(DFNetwork,hcmethod = "single") #Comunidades por LinkComm

plot(DC_lc,
     type = "graph",
     vsize = 20,
     vshape = "circle",
     vlabel = FALSE,
     layout = layout.fruchterman.reingold)



DC_lc$clustsizes



Muestra <- function(Clustnumber){
  
  Nodes <- getNodesIn(DC_lc,clusterids = c(Clustnumber)) # Obtenemos los genes de la comunidad seleccionada
  Nodes<-gsub("9606.","",Nodes) # Procesamos su formato string
  genes = bitr(Nodes, fromType="ENSEMBLPROT", toType="ENTREZID", OrgDb="org.Hs.eg.db") # Los pasamos a tipo ENTREZID
  genes$ENTREZID=as.numeric(genes$ENTREZID) # Pasamos estos ENTREZID a numérico
  return(genes)
  
}

Comm1 <- Muestra(42)




# Enriquecimiento funcional

ego <- enrichGO(gene          = Comm1$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

df_enrich = as.data.frame(ego@result)

df_enrich$qvalue = NULL
rownames(df_enrich) = NULL

print(xtable(df_enrich, type = "latex"), file = "Tabla_Encriquecimiento_Funcional.tex")
