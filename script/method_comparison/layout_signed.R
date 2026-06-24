#' Signed-Weighted Graph Layout
#'
#' @description It elaborates the coordinates for the representation of
#' the vertices of the graph considering only the links with a positive sign.

# igraph: used to create and analyze network graphs
# https://cran.r-project.org/web/packages/igraph/index.html
library(igraph)

LAYOUT_SIGNED <- function(g){
  
  # Extract a subgraph containing only edges with a positive weight
  g.sub <- subgraph.edges(graph = g,
                          eids = which(E(g)$weight > 0),
                          delete.vertices = FALSE)
  
  # Compute vertex coordinates using the Fruchterman-Reingold force-directed algorithm
  layout <- layout.fruchterman.reingold(g.sub)
  
  # Return the matrix of (x, y) coordinates for each vertex
  return(layout)
}