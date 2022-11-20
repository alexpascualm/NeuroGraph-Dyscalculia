library(igraph)

################# Grafos con iGraph.

## Limpiar terminal de R
cat("\014")

## Borrar variables de R
rm(list=ls())

nodes <- read.table("./string_node_degrees.tsv", sep = '\t',header = TRUE, as.is=T)


links <- read.table("./string_interactions.tsv", sep = '\t',header = TRUE, as.is=T)

## Crea un nuevo grafo a partir de los nodos y de las aristas proporcionadas. La red será dirigida.
g <- graph_from_data_frame(d=links, vertices=nodes, directed=T)


## Functiones ##
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

sequential.attacks.targeted <- function(grafo){
  q = seq(from=0,to=1,by=0.01)
  g = grafo
  S = max(components(grafo)$csize)/vcount(grafo)
  contador = S
  removalset = NULL
  v = vcount(grafo)
  s = max(components(grafo)$csize)
  for(i in q){
    if(max(components(g)$csize)/vcount(grafo) >0.05){
      removalset <- names(sort(degree(g),decreasing = T)[1:(i*vcount(g))])
      g <- delete.vertices(graph = g, v = removalset)
      S = c(S, max(components(g)$csize)/vcount(grafo))
      v <- c(v, vcount(g))
      s = c(s, max(components(g)$csize))
      contador = max(components(g)$csize)/vcount(grafo)
    }
  }
  S.vs.q <- tbl_df(data.frame(cbind(q[1:length(S)],S,s,v)))
  names(S.vs.q) <- c("q", "S", "s", "v")
  return(S.vs.q)
}


sequential.attacks.random <- function(grafo){
  q = seq(from=0,to=1,by=0.01)
  g = grafo
  S = max(components(grafo)$csize)/vcount(grafo)
  contador = S
  removalset = NULL
  v = vcount(grafo)
  s = max(components(grafo)$csize)
  for(i in q){
    if(contador > 0.05 & length(removalset) < vcount(g)/2){
      removalset <- sample(x = V(g)$name, size = 10, replace = F)
      g <- delete.vertices(graph = g, v = removalset)
      S = c(S, max(components(g)$csize)/vcount(grafo))
      v <- c(v, vcount(g))
      s = c(s, max(components(g)$csize))
      contador = max(components(g)$csize)/vcount(grafo)
    }
  }
  S.vs.q <- tbl_df(data.frame(cbind(q[1:length(S)],S,s,v)))
  names(S.vs.q) <- c("q", "S", "s", "v")
  return(S.vs.q)
}

ataques_dirigidos <- sequential.attacks.targeted(g)
ataques_dirigidos

#qplot(data=ataques_dirigidos, x=ataques_dirigidos$q, y=ataques_dirigidos$S, xlab = "Iteracion de eliminacion",
#      ylab = "Nodos que siguen conectados a la red frente a los totales", xlim = c(0, 1))


ataques_aleatorios <- sequential.attacks.random(g)
ataques_aleatorios

qplot(data=ataques_aleatorios, x=ataques_aleatorios$q, y=ataques_aleatorios$S, xlab = "Iteracion de eliminacion",
      ylab = "Nodos que siguen conectados a la red frente a los totales", xlim = c(0, 1), ylim= c(0, 1))
