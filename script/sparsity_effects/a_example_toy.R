# a_example_toy.R
# 
# Purpose: 
#   This script generates a minimal network illustration to visually 
#   explain the concept of a single edge between two labeled nodes in a larger graph.
#
# Inputs:
#   - No external files required. All inputs are defined inline:
#     * adj <- matrix(0, nrow = 40, ncol = 40)
#     * adj[5, 25] <- adj[25, 5] <- 1
# Outputs:
#   - a single .png image showing a graph with the connection of two nodes
#
# extrafont: used to use fonts other than the standard PostScript fonts
# https://cran.r-project.org/web/packages/extrafont/index.html
library(extrafont)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# igraph: used to create and analyze network graphs
# https://cran.r-project.org/web/packages/igraph/index.html
library(igraph)

# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)

# import and load fonts for postscript devices
font_import()
loadfonts(device = "postscript")

# generate a 40x40 adjacency matrix filled with zeros (no edges)
adj <- matrix(0, nrow = 40, ncol = 40)

# add a single undirected edge between node 5 and node 25
adj[5, 25] <- adj[25, 5] <- 1

# create an undirected graph from the adjacency matrix
g <- graph_from_adjacency_matrix(adj, mode = "undirected")

# style the edge: dark gray, width 2, dashed line
E(g)$color <- "darkgray"
E(g)$width <- 2
E(g)$lty   <- 2

# set seed for reproducibility of random node sizes
set.seed(4)

# assign a random size to each of the 40 nodes
V(g)$size <- runif(n = 40, min = 3, max = 12.5)

# open a PNG device (1200x1200 px, 300 dpi)
png(filename = here("Plots", "example_netw.png"), width = 1200, height = 1200, res = 300)

# remove margins
par(mar = c(0,0,0,0))

# plot the graph in a circular layout, without node labels
plot(g, layout = layout_in_circle(g), vertex.label = NA)

# label node 5 as "I" and node 25 as "J" using Times New Roman
text( 0.925,  0.7, "I", family = "Times New Roman", cex = 2)
text(-0.925, -0.7, "J", family = "Times New Roman", cex = 2)

# close the PNG device and save the file
dev.off()