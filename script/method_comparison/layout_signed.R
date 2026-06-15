#' Signed-Weighted Graph Layout
#'
#' @description It elaborates the coordinates for the representation of
#' the vertices of the graph considering only the links with a positive sign.

# igraph: used to create and analyze network graphs
# https://cran.r-project.org/web/packages/igraph/index.html
library(igraph)

LAYOUT_SIGNED <- function(g){
  g.sub <- subgraph.edges(graph = g,
                                  eids = which(E(g)$weight>0),
                                  delete.vertices = FALSE)
  
  layout <- layout.fruchterman.reingold(g.sub)
  
  return(layout)
}
