# Converts base R plots and grobs into ggplot2 objects (as.grob)
# https://cran.r-project.org/web/packages/ggplotify/index.html
library(ggplotify)

# ggplot2 extensions for publication-ready graphics (scatter, stat_cor)
# https://cran.r-project.org/web/packages/ggpubr/index.html
library(ggpubr)

# Low-level graphics: viewport and grob management
# https://cran.r-project.org/web/packages/grid/index.html
library(grid)

# Arranges multiple plots in grids (grid.arrange, tableGrob)
# https://cran.r-project.org/web/packages/gridExtra/index.html
library(gridExtra)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# Sparse microbial network estimation (SPIEC-EASI: MB and GLASSO)
# https://github.com/zdk123/SpiecEasi  ← Bioconductor/GitHub
library(SpiecEasi)

# Collection of packages for data wrangling and visualization
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)

# -------- create output directories --------

dir.create(here("script", "tmp_files"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("outputs"),             recursive = TRUE, showWarnings = FALSE)


# -------- load custom functions --------

# Centered Log-Ratio transformation
source(here("script", "method_comparison", "CLR.R"))

# Upper triangle extraction from a matrix
source(here("script", "method_comparison", "TRIU.R"))

# Network layout separating positive and negative edges
source(here("script", "method_comparison", "layout_signed.R"))


# -------- read and filter data --------

# Load OTU abundance matrix (samples x OTUs)
otu  <- readRDS(here("data", "otu_HMP2.rds" ))
# Load sample metadata
meta <- readRDS(here("data", "meta_HMP2.rds"))
# Load OTU taxonomic annotations
taxa <- readRDS(here("data", "taxonomy.rds" ))

# Select samples from subject 69-001 in healthy status
otu.69001.H <- otu[meta$SubjectID == "69-001" & meta$CL4_2 == "Healthy", ]

# Filter rare OTUs:
# prevalence >= 33% of samples
otu.filt  <- otu.69001.H[, colSums(otu.69001.H > 0) / nrow(otu.69001.H) >= 0.33]
# median of non-zero values >= 5 reads
otu.filt  <- otu.filt[, apply(otu.filt, 2, function(x) median(x[x > 0]) >= 5)]
taxa.filt <- taxa[colnames(otu.filt), ]


# -------- Correlation Methods --------


# -------- SPIEC-EASI GLASSO --------

# Sparse network estimation via Graphical LASSO on compositional data
# ncores = 1: parallelization disabled to avoid workers not finding SpiecEasi
res.gl <- spiec.easi(data = otu.filt, 
                     method = 'glasso',
                     lambda.max = 0.75, 
                     lambda.min.ratio = 0.5,
                     pulsar.params = list(ncores = 1, thresh = 0.05))

# Adjacency matrix: binarized partial correlations (+1/-1)
# getOptCov() is valid for GLASSO, which estimates a covariance matrix
adj.gl <- cov2cor(as.matrix(getOptCov(res.gl)))
adj.gl <- adj.gl * as.matrix(getRefit(res.gl))
colnames(res.gl$est$data) -> colnames(adj.gl) -> rownames(adj.gl)
adj.gl[adj.gl >  0] <-  1
adj.gl[adj.gl <  0] <- -1


# -------- SPIEC-EASI MB --------

# Sparse network estimation via Meinshausen-Bühlmann neighborhood selection
res.mb <- spiec.easi(data = otu.filt, 
                     method = 'mb',
                     lambda.max = 0.75, 
                     lambda.min.ratio = 0.5,
                     pulsar.params = list(ncores = 1, thresh = 0.05))

# MB estimates directional regression coefficients (getOptBeta), not a covariance matrix.
# Symmetrize by averaging i->j and j->i coefficients to obtain an undirected adjacency matrix
beta.mb          <- as.matrix(getOptBeta(res.mb))
adj.mb           <- (beta.mb + t(beta.mb)) / 2
colnames(adj.mb) <- rownames(adj.mb) <- colnames(otu.filt)
adj.mb[adj.mb >  0] <-  1
adj.mb[adj.mb <  0] <- -1


# -------- SparCC --------

# Log-ratio-based correlations, robust to compositionality
res.cc <- sparcc(otu.filt)$Cor
colnames(res.cc) <- rownames(res.cc) <- colnames(otu.filt)


# -------- Rho (propr) --------

# Symmetric proportionality measure between OTU pairs
res.rho <- propr::propr(counts = otu.filt, metric = "rho")@matrix


# -------- Pearson + CLR --------

# Pearson correlation after CLR transformation + Bonferroni p-value correction
res.clr <- cor(CLR(otu.filt), method = "pearson")
diag(res.clr) <- 0

# https://cran.r-project.org/web/packages/psych/index.html
# psych::corr.p used for significance testing
# Compute p-values for the correlation matrix with Bonferroni correction
p.adjust <- psych::corr.p(r = res.clr,
                          n = nrow(otu.filt),
                          adjust = "bonferroni", 
                          ci = FALSE)$p
# Mirror upper triangle into lower triangle to make the matrix symmetric
p.adjust[lower.tri(p.adjust)] <- t(p.adjust)[lower.tri(p.adjust)]
diag(p.adjust) <- 1
# Keep only significant correlations
adj.clr <- res.clr * (p.adjust <= 0.05)


# -------- prepare data for plots --------

# CLR correlations filtered by GLASSO mask (selected edges only)
res.gl.clr <- TRIU(res.clr * abs(adj.gl)) %>% .[. != 0]

# Dataframe for histogram: full CLR distribution vs. GLASSO subset
df.gl <- tibble("method" = rep("clr", length(TRIU(res.clr))),
                "value"  = TRIU(res.clr)) %>%
   rbind(tibble("method" = rep("gl", length(res.gl.clr)),
                "value"  = res.gl.clr))


# -------- Plots --------


# -------- panel A: Scatter Pearson+CLR vs SparCC --------

p1 <- ggscatter(
  data.frame("PearsonCLR" = TRIU(res.clr), "SparCC" = TRIU(res.cc)),
  x   = "PearsonCLR", 
  y   = "SparCC",
  add = "reg.line", 
  conf.int   = TRUE,
  add.params = list(color = "red", fill = "lightgray")) +
  
  stat_cor(aes(label = after_stat(r.label)),
                   label.x =  0.45, 
                   label.y = -0.25, 
                   size    =  6.00) +
  
  theme_bw() +
  
  xlab("Pearson+CLR") +
  
  theme(plot.title = element_text(hjust = 0.5))


# -------- panel B: Scatter Pearson+CLR vs Rho --------

p2 <- ggscatter(
  data.frame("PearsonCLR" = TRIU(res.clr), "Rho" = TRIU(res.rho)),
  x   = "PearsonCLR", 
  y   = "Rho",
  add = "reg.line", 
  conf.int   = TRUE,
  add.params = list(color = "red", fill = "lightgray")) +
  
  stat_cor(aes(label = after_stat(r.label)),
                   label.x =  0.45,
                   label.y = -0.25, 
                   size    =  6.00) +
  
  theme_bw() + xlab("Pearson+CLR") +
  
  theme(plot.title = element_text(hjust = 0.5))


# -------- panel C: Histogram of CLR vs GLASSO correlation distributions --------

# Dashed line at y = 200 as visual reference threshold
p4.all <- ggplot(df.gl, aes(x = value, fill = method, color = method)) +
  
  geom_histogram(position = "identity", alpha = 0.5, breaks = seq(-1, 1, by = 0.1)) +
  
  geom_hline(yintercept  = 200, linetype = "twodash", color = "black", linewidth = 1) +
  
  xlim(c(-1, 1)) + 
  
  theme_bw() +
  
  scale_color_manual(values = c("steelblue1", "forestgreen")) +
  
  scale_fill_manual(values  = c("steelblue1", "forestgreen")) +
  
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right", legend.direction = "vertical")

# zoomed inset of the histogram tail (y <= 200)
p4.zoom <- ggplot(df.gl, aes(x = value, fill = method, color = method)) +
  
  geom_histogram(position   = "identity", alpha = 0.5, breaks = seq(-1, 1, by = .1)) +
  
  coord_cartesian(ylim      = c(0, 200)) +
  
  scale_color_manual(values = c("steelblue1", "forestgreen")) +
  
  scale_fill_manual(values  = c("steelblue1", "forestgreen")) +
  
  theme_void() + theme(legend.position = "none") +
  
  geom_rect(aes(xmin  = -Inf, 
                xmax  =  Inf, 
                ymin  = -Inf, 
                ymax  = +Inf),
                col   = "black", 
                alpha = 0, 
            linewidth = 1)

# Extract the actual max count from the histogram to position the inset dynamically
# This avoids hardcoding y coordinates that may not match the data range
y_max <- max(ggplot_build(p4.all)$data[[1]]$count)

# Composition: main histogram + inset in the upper left
# ymin/ymax are set proportionally to the actual Y range for robust positioning
p4 <- p4.all +
  
  annotation_custom(ggplotGrob(p4.zoom),
                    xmin = -1.1,
                    xmax = -0.3,
                    ymin = y_max * 0.70,
                    ymax = y_max * 0.95)


# -------- build graphs --------

# igraph: library for network analysis and visualization
# https://cran.r-project.org/web/packages/igraph/index.html
g.mb  <- igraph::graph_from_adjacency_matrix(adj.mb,  mode = "undirected", weighted = TRUE)
g.gl  <- igraph::graph_from_adjacency_matrix(adj.gl,  mode = "undirected", weighted = TRUE)
g.clr <- igraph::graph_from_adjacency_matrix(adj.clr, mode = "undirected", weighted = TRUE)


# -------- network visualization --------

# qualpalr: generates qualitatively distinct color palettes
# https://cran.r-project.org/web/packages/qualpalr/index.html
# set seed for reproducibility
set.seed(42)
colpal <- rownames(qualpalr::qualpal(n = length(unique(taxa.filt[, "family"])))$RGB)
names(colpal) <- unique(taxa.filt[, "family"])

# Node size proportional to mean CLR abundance
vertex.size <- colMeans(CLR(otu.filt) - min(CLR(otu.filt)))


# -------- panel D: CLR network --------

# blue edges = positive correlations, red edges = negative correlations
set.seed(42)
p.graph.CLR <- as.grob(~plot(g.clr, vertex.label = NA,
                             vertex.color = colpal[taxa.filt[, "family"]],
                             vertex.size  = vertex.size,
                             edge.color   = ifelse(igraph::E(g.clr)$weight > 0,
                                                   rgb(0, 0, 1), 
                                                   rgb(1, 0, 0)),
                             edge.width   = 0.5,
                             layout       = LAYOUT_SIGNED(g.clr)))


# -------- panel E: GLASSO network --------

# Same layout as CLR network for direct visual comparison
set.seed(42)
p.graph.glasso <- as.grob(~plot(g.gl, vertex.label = NA,
                                vertex.color = colpal[taxa.filt[, "family"]],
                                vertex.size  = vertex.size,
                                edge.color   = ifelse(igraph::E(g.gl)$weight > 0,
                                                      rgb(0, 0, 1), 
                                                      rgb(1, 0, 0)),
                                edge.width   = 0.5,
                                layout       = LAYOUT_SIGNED(g.clr)))


# -------- panel F: Venn diagram of shared edges between GLASSO and CLR --------
# ggVennDiagram: Venn diagrams with ggplot2
# https://cran.r-project.org/web/packages/ggVennDiagram/index.html
p.ven <- ggVennDiagram::ggVennDiagram(
  x = list(
    "GLASSO" = paste(igraph::as_edgelist(g.gl) [, 1], "-",
                     igraph::as_edgelist(g.gl) [, 2], sep = ""),
    "CLR"    = paste(igraph::as_edgelist(g.clr)[, 1], "-",
                     igraph::as_edgelist(g.clr)[, 2], sep = "")
  ), label_alpha = 0) +
  
  ggplot2::scale_fill_gradient("Shared \n Links", low = "white", high = "red") +
  
  ggplot2::scale_color_manual(values = c("gray20", "gray20", "gray20")) +
  
  theme(legend.position = "right")


# -------- save final outputs --------

# -------- Main figure: 6 panels (A-F), 3 columns x 2 rows --------

png(here("outputs", "Methods_comparison.png"), width = 3600, height = 2400, res = 300)
ggarrange(
  plotlist = list(p1, p2, p4,
                  p.graph.CLR, 
                  p.graph.glasso,
                  p.ven + 
                    
                  theme(legend.position = "bottom")),
  labels = c("A", "B", "C", "D", "E", "F"),
  ncol   = 3, 
  nrow   = 2
)
dev.off()


# -------- save network plots to disk --------

# These PNG files are required by cowplot::draw_image() in the next section
# MB network
p.mb <- png(filename = here("script", "tmp_files", "graph_mb.png"), width = 1200, height = 1200, res = 200)
set.seed(42)
plot(g.mb, vertex.label = NA,
     vertex.color = colpal[taxa.filt[, "family"]],
     vertex.size  = vertex.size,
     edge.color   = ifelse(igraph::E(g.mb)$weight > 0, 
                           rgb(0, 0, 1), 
                           rgb(1, 0, 0)),
     edge.width   = 0.5,
     layout       = LAYOUT_SIGNED(g.mb))
dev.off()

# GLASSO network
p.gl <- png(filename = here("script", "tmp_files", "graph_gl.png"), width = 1200, height = 1200, res = 200)
set.seed(42)
plot(g.gl, vertex.label = NA,
     vertex.color = colpal[taxa.filt[, "family"]],
     vertex.size  = vertex.size,
     edge.color   = ifelse(igraph::E(g.gl)$weight > 0, 
                           rgb(0, 0, 1), 
                           rgb(1, 0, 0)),
     edge.width   = 0.5,
     layout       = LAYOUT_SIGNED(g.gl))
dev.off()

# CLR network
p.clr <- png(filename = here("script", "tmp_files", "graph_clr.png"), width = 1200, height = 1200, res = 200)
set.seed(42)
plot(g.clr, vertex.label = NA,
     vertex.color = colpal[taxa.filt[, "family"]],
     vertex.size  = vertex.size,
     edge.color   = ifelse(igraph::E(g.clr)$weight > 0, 
                           rgb(0, 0, 1), 
                           rgb(1, 0, 0)),
     edge.width   = 0.5,
     layout       = LAYOUT_SIGNED(g.clr))
dev.off()

# Venn diagram
png(filename = here("script", "tmp_files", "ggvenn.png"), width = 1200, height = 1200, res = 200)
print(p.ven)
dev.off()


# -------- taxonomy color legend (family level) --------

p.clr <- png(filename = here("script", "tmp_files", "colorLegend.png"), width = 1200, height = 1200, res = 200)
par(mar = c(2, 0, 2, 0))
plot.new()
legend("center", legend = names(colpal), fill = colpal, title = "Family", ncol = 2)
dev.off()


# -------- Summary table: edge count and % positive/negative edges per method --------

# Summarize edge count and proportion of positive/negative edges for MB, GLASSO, 
# and CLR networks in a a data frame 

netw.info <- data.frame(
  "Edges"    = sapply(list(g.mb, g.gl, g.clr), igraph::ecount),
  "Positive" = paste(100 * sapply(list(g.mb, g.gl, g.clr),
                                  function(x) round(sum(igraph::E(x)$weight > 0) /
                                                      igraph::ecount(x), 2)), "%", sep = ""),
  "Negative" = paste(100 * sapply(list(g.mb, g.gl, g.clr),
                                  function(x) round(sum(igraph::E(x)$weight < 0) /
                                                      igraph::ecount(x), 2)), "%", sep = "")
)
rownames(netw.info) <- c("MB", "GLASSO", "CLR")

# generate a .png that allows to visualize the results
p.inf <- png(filename = here("script", "tmp_files", "table.png"), width = 1200, height = 1200, res = 300)
print(
  ggplot() +
  
  theme(axis.line        = element_blank(), axis.text.x      = element_blank(),
        axis.text.y      = element_blank(), axis.ticks       = element_blank(),
        axis.title.x     = element_blank(), axis.title.y     = element_blank(),
        legend.position  = "none",          panel.background = element_blank(),
        panel.border     = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.background  = element_blank()) +
  
  annotation_custom(gridExtra::tableGrob(netw.info),
                    xmin = -Inf, 
                    xmax =  Inf, 
                    ymin = -Inf, 
                    ymax =  Inf)
  )
dev.off()


# -------- Alternative composition using pre-saved images from disk --------
# cowplot: advanced composition of ggplot2 graphics
# https://cran.r-project.org/web/packages/cowplot/index.html

# Load pre-saved network plots and supplementary figures as ggdraw objects
p.mb  <- cowplot::ggdraw() + cowplot::draw_image(here("script", "tmp_files", "graph_mb.png"))
p.gl  <- cowplot::ggdraw() + cowplot::draw_image(here("script", "tmp_files", "graph_gl.png"))
p.clr <- cowplot::ggdraw() + cowplot::draw_image(here("script", "tmp_files", "graph_clr.png"))
p.leg <- cowplot::ggdraw() + cowplot::draw_image(here("script", "tmp_files", "colorLegend.png"))
p.ven <- cowplot::ggdraw() + cowplot::draw_image(here("script", "tmp_files", "ggvenn.png"))
p.inf <- cowplot::ggdraw() + cowplot::draw_image(here("script", "tmp_files", "table.png"))

# Arrange all panels into a 3-column grid and export as a single PNG
png(filename = here("outputs", "Graph_Comparison.png"), width = 1600, height = 1000, res = 300)
grid.arrange(p.mb, p.gl, p.clr, p.ven, p.leg, p.inf, ncol = 3)
dev.off()