.libPaths(paste0(Sys.getenv("R_LIBS_USER"),"/discalculia"))  # add to the path

suppressMessages(library(igraph))
suppressMessages(library(STRINGdb))
suppressMessages(library(linkcomm))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(xtable))

## Limpiar terminal de R ##
cat("\014")

## Borrar variables de R ##
rm(list=ls())

nodes <- read.table("./string_node_degrees.tsv", sep = '\t',header = TRUE, as.is=T)


links <- read.table("./string_interactions.tsv", sep = '\t',header = TRUE, as.is=T)


##### Grafos con iGraph ####

## Crea un nuevo grafo a partir de los nodos y de las aristas proporcionadas (La red será dirigida) ##

net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)


lo <- layout_with_kk(net) # Creamos la disposición (layout) de los nodos
lo <- norm_coords(lo, ymin=-1, ymax=1, xmin=-1, xmax=1)

## Representa graficamente ##

pdf(file="../results/Red_de_genes_Original.pdf")
plot(net, edge.arrow.size=0.2, edge.curved=0.1, vertex.size=15, vertex.color="white", vertex.frame.color="black", vertex.label.color="black",vertex.label.cex=0.4, layout=lo)
dev.off()

## Detección de comunidades por proximidad ##

community <- cluster_edge_betweenness(net)

dendPlot(community) 

pdf(file="../results/comunidades_por_proximidad.pdf")
plot(community, net, edge.arrow.size=0.05, edge.curved=0.1, vertex.size=13, vertex.color="white", vertex.frame.color="black", vertex.label.color="black",vertex.label.cex=0.4, layout=lo)
dev.off()

length(community)
membership(community)


## Detección de comunidades mediante optimización voraz ##

cfg <- cluster_fast_greedy(as.undirected(net))

dendPlot(cfg)

pdf(file="../results/comunidades_por_voraz.pdf")
plot(cfg, as.undirected(net), vertex.size=13, vertex.color="white", vertex.frame.color="black", vertex.label.color="black",vertex.label.cex=0.4, layout=lo)
dev.off()

length(cfg)
membership(cfg)


#### NetworkPropagation de la red ####

string_db <- STRINGdb$new( version="11", species=9606, score_threshold=400, input_directory="" )
string.network <- string_db$get_graph() # Obtenemos una red de nodos de string del organismo Homo Sapiens

hits <- nodes$identifier # Nos quedamos con los identificadores de String de nuestro conjunto de genes
hits.network <- string_db$get_subnetwork(hits)  # Creamos una network de nuestros nodos usando StringDb

first.neigh <- (neighbors(graph = string.network, v = V(hits.network)$name, mode = "all"))$name   # Encontramos una serie de nodos vecinos
hits.network <- string_db$get_subnetwork(unique(c(V(hits.network)$name, first.neigh))) # Unimos la red de nodos originales con la red de vecinos

DFNetwork <- igraph::as_data_frame(hits.network)





## Analisis por linkcomm ##

pdf(file="../results/dendograma_por_linkcomm.pdf")
DC_lc <- getLinkCommunities(DFNetwork,hcmethod = "single") # Comunidades por LinkComm
dev.off()

pdf(file="../results/comunidades_por_linkcomm.pdf")
plot(DC_lc,
     type = "graph",
     vsize = 20,
     vshape = "circle",
     vlabel = FALSE,
     layout = layout.fruchterman.reingold)
dev.off()


max(DC_lc$numclusters)

pdf(file="../results/comunidades_por_linkcomm_barplot.pdf")
barplot(DC_lc$clustsize, xlab="Id comunity", ylab="Size of cluster") # Tamanos de los diferentes clusters
dev.off()

#### Analisis de nuestra comunidad de genes originales ####

## Funcion para traducir a ENTREZID los genes de una comunidad dada ##

ID_to_EntrezID <- function(Clustnumber){
  Nodes <- getNodesIn(DC_lc,clusterids = c(Clustnumber))           # Obtenemos los genes de la comunidad seleccionada
  Nodes<-gsub("9606.","",Nodes)                                                        # Procesamos su formato string
  genes = bitr(Nodes, fromType="ENSEMBLPROT", toType="ENTREZID", OrgDb="org.Hs.eg.db") # Los pasamos a tipo ENTREZID
  genes$ENTREZID=as.numeric(genes$ENTREZID)                                            # Pasamos estos ENTREZID a numérico
  return(genes)
}

### Analisis de nuestra comunidad de genes originales ###

originales <- gsub("9606.","",nodes$identifier)
ororiginales <- bitr(originales, fromType="ENSEMBLPROT", toType="ENTREZID", OrgDb="org.Hs.eg.db") # Los pasamos a tipo ENTREZID
originales <- as.numeric(originales$ENTREZID)
ego <- enrichGO(gene          = originales,
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

head(ego)

# A continuación se procede a guardar toda la tabla en un formato legible por latex

df_enrich = as.data.frame(ego@result)

df_enrich$qvalue = NULL
rownames(df_enrich) = NULL

df_enrich2 <- df_enrich

df_enrich$BgRatio = NULL
df_enrich$pvalue = NULL
df_enrich$p.adjust = NULL
df_enrich$geneID = NULL
df_enrich$Count = NULL

df_enrich2$Description = NULL
df_enrich2$GeneRatio = NULL

dir.create(file.path("../results/", "/comunidad_Original"))
print(xtable(df_enrich, type = "latex"), file = "../results/comunidad_Original/Tabla_Encriquecimiento_Funcional_parte1.tex") # Se guardara en el file especificado
print(xtable(df_enrich2, type = "latex"), file = "../results/comunidad_Original/Tabla_Encriquecimiento_Funcional_parte2.tex") # Se guardara en el file especificado


## Enriquecimiento funcional comunidad 64 ##

Comm <- ID_to_EntrezID(1) # Id de la comunidad del cluster que se quiera estudiar

ego <- enrichGO(gene          = Comm$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

head(ego)

# Guardamos las tablas en un formato legible por latex #

df_enrich = as.data.frame(ego@result)

df_enrich$qvalue = NULL
rownames(df_enrich) = NULL

df_enrich2 <- df_enrich
df_enrich3 <- df_enrich

df_enrich$BgRatio = NULL
df_enrich$pvalue = NULL
df_enrich$p.adjust = NULL
df_enrich$geneID = NULL
df_enrich$Count = NULL

df_enrich2$Description = NULL
df_enrich2$GeneRatio = NULL
df_enrich2$geneID = NULL

df_enrich3$Description = NULL
df_enrich3$GeneRatio = NULL
df_enrich3$BgRatio = NULL
df_enrich3$pvalue = NULL
df_enrich3$p.adjust = NULL
df_enrich3$Count = NULL

dir.create(file.path("../results/", "/comunidad_64"))
print(xtable(df_enrich, type = "latex"), file = "../results/comunidad_64/Tabla_Encriquecimiento_Funcional_parte1.tex") # Se guardara en el archivo especificado
print(xtable(df_enrich2, type = "latex"), file = "../results/comunidad_64/Tabla_Encriquecimiento_Funcional_parte2.tex") # Se guardara en el archivo especificado
print(xtable(df_enrich3, type = "latex"), file = "../results/comunidad_64/Tabla_Encriquecimiento_Funcional_parte3.tex") # Se guardara en el archivo especificado



## Enriquecimiento funcional comunidad 70 ##

Comm <- ID_to_EntrezID(70) # Id de la comunidad del cluster que se quiera estudiar

ego <- enrichGO(gene          = Comm$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

head(ego)

# Guardamos las tablas en un formato legible por latex #

df_enrich = as.data.frame(ego@result)

df_enrich$qvalue = NULL
rownames(df_enrich) = NULL

df_enrich2 <- df_enrich
df_enrich3 <- df_enrich

df_enrich$BgRatio = NULL
df_enrich$pvalue = NULL
df_enrich$p.adjust = NULL
df_enrich$geneID = NULL
df_enrich$Count = NULL

df_enrich2$Description = NULL
df_enrich2$GeneRatio = NULL
df_enrich2$geneID = NULL

df_enrich3$Description = NULL
df_enrich3$GeneRatio = NULL
df_enrich3$BgRatio = NULL
df_enrich3$pvalue = NULL
df_enrich3$p.adjust = NULL
df_enrich3$Count = NULL

dir.create(file.path("../results/", "/comunidad_70"))
print(xtable(df_enrich, type = "latex"), file = "../results/comunidad_70/Tabla_Encriquecimiento_Funcional_parte1.tex") # Se guardara en el archivo especificado
print(xtable(df_enrich2, type = "latex"), file = "../results/comunidad_70/Tabla_Encriquecimiento_Funcional_parte2.tex") # Se guardara en el archivo especificado
print(xtable(df_enrich3, type = "latex"), file = "../results/comunidad_70/Tabla_Encriquecimiento_Funcional_parte3.tex") # Se guardara en el archivo especificado



## Enriquecimiento funcional comunidad 33 ##

Comm <- ID_to_EntrezID(33) # Id de la comunidad del cluster que se quiera estudiar

ego <- enrichGO(gene          = Comm$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

head(ego)

# Guardamos las tablas en un formato legible por latex #

df_enrich = as.data.frame(ego@result)

df_enrich$qvalue = NULL
rownames(df_enrich) = NULL

df_enrich2 <- df_enrich
df_enrich3 <- df_enrich

df_enrich$BgRatio = NULL
df_enrich$pvalue = NULL
df_enrich$p.adjust = NULL
df_enrich$geneID = NULL
df_enrich$Count = NULL

df_enrich2$Description = NULL
df_enrich2$GeneRatio = NULL
df_enrich2$geneID = NULL

df_enrich3$Description = NULL
df_enrich3$GeneRatio = NULL
df_enrich3$BgRatio = NULL
df_enrich3$pvalue = NULL
df_enrich3$p.adjust = NULL
df_enrich3$Count = NULL

dir.create(file.path("../results/", "/comunidad_33"))
print(xtable(df_enrich, type = "latex"), file = "../results/comunidad_33/Tabla_Encriquecimiento_Funcional_parte1.tex") # Se guardara en el archivo especificado
print(xtable(df_enrich2, type = "latex"), file = "../results/comunidad_33/Tabla_Encriquecimiento_Funcional_parte2.tex") # Se guardara en el archivo especificado
print(xtable(df_enrich3, type = "latex"), file = "../results/comunidad_33/Tabla_Encriquecimiento_Funcional_parte3.tex") # Se guardara en el archivo especificado
