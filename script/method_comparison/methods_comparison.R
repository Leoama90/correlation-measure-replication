# -------- LIBRARIES --------

# Collection of packages for data wrangling and visualization
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)

# Sparse microbial network estimation (SPIEC-EASI: MB and GLASSO)
# https://github.com/zdk123/SpiecEasi  ← Bioconductor/GitHub
library(SpiecEasi)

# Proportionality analysis for compositional data (rho, phi)
# https://cran.r-project.org/web/packages/propr/index.html
library(propr)

# ggplot2 extensions for publication-ready graphics (scatter, stat_cor)
# https://cran.r-project.org/web/packages/ggpubr/index.html
library(ggpubr)

# Low-level graphics: viewport and grob management
# https://cran.r-project.org/web/packages/grid/index.html
library(grid)

# Converts base R plots and grobs into ggplot2 objects (as.grob)
# https://cran.r-project.org/web/packages/ggplotify/index.html
library(ggplotify)

# Arranges multiple plots in grids (grid.arrange, tableGrob)
# https://cran.r-project.org/web/packages/gridExtra/index.html
library(gridExtra)


# -------- CUSTOM FUNCTIONS --------

source("CLR.R")           # Centered Log-Ratio transformation
source("TRIU.R")          # Upper triangle extraction from a matrix
source("LAYOUT_SIGNED.R") # Network layout separating positive and negative edges


# -------- READ AND FILTER DATA --------

otu  <- readRDS("../../data/otu_HMP2.rds")
meta <- readRDS("../../data/meta_HMP2.rds")
taxa <- readRDS("../../data/taxonomy.rds")

# Select samples from subject 69-001 in healthy status
otu.69001.H <- otu[meta$SubjectID == "69-001" & meta$CL4_2 == "Healthy", ]

# Filter rare OTUs:
#   - prevalence >= 33% of samples
#   - median of non-zero values >= 5 reads
otu.filt  <- otu.69001.H[, colSums(otu.69001.H > 0) / nrow(otu.69001.H) >= .33]
otu.filt  <- otu.filt[, apply(otu.filt, 2, function(x) median(x[x > 0]) >= 5)]
taxa.filt <- taxa[colnames(otu.filt), ]


# -------- CORRELATION METHODS --------

# -- SPIEC-EASI GLASSO --
# Sparse network estimation via Graphical LASSO on compositional data
res.gl <- spiec.easi(data = otu.filt, method = 'glasso',
                     lambda.max = .75, lambda.min.ratio = .5,
                     pulsar.params = list(ncores = 6, thresh = 0.05))

# Adjacency matrix: binarized partial correlations (+1/-1)
adj.gl <- cov2cor(as.matrix(getOptCov(res.gl)))
adj.gl <- adj.gl * as.matrix(getRefit(res.gl))
colnames(res.gl$est$data) -> colnames(adj.gl) -> rownames(adj.gl)
adj.gl[adj.gl >  0] <-  1
adj.gl[adj.gl <  0] <- -1

# -- SparCC --
# Log-ratio-based correlations, robust to compositionality
res.cc <- sparcc(otu.filt)$Cor
colnames(res.cc) <- rownames(res.cc) <- colnames(otu.filt)

# -- Rho (propr) --
# Symmetric proportionality measure between OTU pairs
res.rho <- propr::propr(counts = otu.filt, metric = "rho")@matrix

# -- Pearson + CLR --
# Pearson correlation after CLR transformation + Bonferroni p-value correction
# psych::corr.p used for significance testing
# https://cran.r-project.org/web/packages/psych/index.html
res.clr <- cor(CLR(otu.filt), method = "pearson")
diag(res.clr) <- 0

p.adjust <- psych::corr.p(r = res.clr, n = nrow(otu.filt),
                          adjust = "bonferroni", ci = FALSE)$p
p.adjust[lower.tri(p.adjust)] <- t(p.adjust)[lower.tri(p.adjust)]
diag(p.adjust) <- 1
adj.clr <- res.clr * (p.adjust <= .05)  # Keep only significant correlations


# -------- PREPARE DATA FOR PLOTS --------

# CLR correlations filtered by GLASSO mask (selected edges only)
res.gl.clr <- TRIU(res.clr * abs(adj.gl)) %>% .[. != 0]

# Dataframe for histogram: full CLR distribution vs. GLASSO subset
df.gl <- tibble("method" = rep("clr", length(TRIU(res.clr))),
                "value"  = TRIU(res.clr)) %>%
  rbind(tibble("method" = rep("gl", length(res.gl.clr)),
               "value"  = res.gl.clr))


# -------- PLOTS --------

# -- Panel A: Scatter Pearson+CLR vs SparCC --
p1 <- ggpubr::ggscatter(
  data.frame("PearsonCLR" = TRIU(res.clr), "SparCC" = TRIU(res.cc)),
  x = "PearsonCLR", y = "SparCC",
  add = "reg.line", conf.int = TRUE,
  add.params = list(color = "red", fill = "lightgray")) +
  ggpubr::stat_cor(aes(label = after_stat(r.label)),
                   label.x = .45, label.y = -.25, size = 6) +
  theme_bw() + xlab("Pearson+CLR") +
  theme(plot.title = element_text(hjust = 0.5))

# -- Panel B: Scatter Pearson+CLR vs Rho --
p2 <- ggpubr::ggscatter(
  data.frame("PearsonCLR" = TRIU(res.clr), "Rho" = TRIU(res.rho)),
  x = "PearsonCLR", y = "Rho",
  add = "reg.line", conf.int = TRUE,
  add.params = list(color = "red", fill = "lightgray")) +
  ggpubr::stat_cor(aes(label = after_stat(r.label)),
                   label.x = .45, label.y = -.25, size = 6) +
  theme_bw() + xlab("Pearson+CLR") +
  theme(plot.title = element_text(hjust = 0.5))

# -- Panel C: Histogram of CLR vs GLASSO correlation distributions --
# Dashed line at y=200 as visual reference threshold
p4.all <- ggplot(df.gl, aes(x = value, fill = method, color = method)) +
  geom_histogram(position = "identity", alpha = .5, breaks = seq(-1, 1, by = .1)) +
  geom_hline(yintercept = 200, linetype = "twodash", color = "black", linewidth = 1) +
  xlim(c(-1, 1)) + theme_bw() +
  scale_color_manual(values = c("steelblue1", "forestgreen")) +
  scale_fill_manual(values  = c("steelblue1", "forestgreen")) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right", legend.direction = "vertical")

# Zoomed inset of the histogram tail (y <= 200)
p4.zoom <- ggplot(df.gl, aes(x = value, fill = method, color = method)) +
  geom_histogram(position = "identity", alpha = .5, breaks = seq(-1, 1, by = .1)) +
  coord_cartesian(ylim = c(0, 200)) +
  scale_color_manual(values = c("steelblue1", "forestgreen")) +
  scale_fill_manual(values  = c("steelblue1", "forestgreen")) +
  theme_void() + theme(legend.position = "none") +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = +Inf),
            col = "black", alpha = 0, linewidth = 1)

# Composition: main histogram + inset in the upper left
p4 <- p4.all +
  annotation_custom(ggplotGrob(p4.zoom), xmin = -1.1, xmax = -.3, ymin = 1400, ymax = 1900)


# -------- BUILD GRAPHS --------

# igraph: library for network analysis and visualization
# https://cran.r-project.org/web/packages/igraph/index.html
g.gl  <- igraph::graph_from_adjacency_matrix(adj.gl,  mode = "undirected", weighted = TRUE)
g.clr <- igraph::graph_from_adjacency_matrix(adj.clr, mode = "undirected", weighted = TRUE)

# -- Panel F: Venn diagram of shared edges between GLASSO and CLR --
# ggVennDiagram: Venn diagrams with ggplot2
# https://cran.r-project.org/web/packages/ggVennDiagram/index.html
p.venn <- ggVennDiagram::ggVennDiagram(
  x = list(
    "GLASSO" = paste(igraph::as_edgelist(g.gl)[, 1],  "-",
                     igraph::as_edgelist(g.gl)[, 2],  sep = ""),
    "CLR"    = paste(igraph::as_edgelist(g.clr)[, 1], "-",
                     igraph::as_edgelist(g.clr)[, 2], sep = "")
  ), label_alpha = 0) +
  ggplot2::scale_fill_gradient("Shared \n Links", low = "white", high = "red") +
  ggplot2::scale_color_manual(values = c("gray20", "gray20", "gray20")) +
  theme(legend.position = "right")


# -------- NETWORK VISUALIZATION --------

# qualpalr: generates qualitatively distinct color palettes
# https://cran.r-project.org/web/packages/qualpalr/index.html
set.seed(42)
colpal <- rownames(qualpalr::qualpal(n = length(unique(taxa.filt[, "family"])))$RGB)
names(colpal) <- unique(taxa.filt[, "family"])

# Node size proportional to mean CLR abundance
vertex.size <- colMeans(CLR(otu.filt) - min(CLR(otu.filt)))

# -- Panel D: CLR network --
# Blue edges = positive correlations, red edges = negative correlations
set.seed(42)
p.graph.CLR <- as.grob(~plot(g.clr, vertex.label = NA,
                             vertex.color = colpal[taxa.filt[, "family"]],
                             vertex.size  = vertex.size,
                             edge.color   = ifelse(igraph::E(g.clr)$weight > 0,
                                                   rgb(0, 0, 1), rgb(1, 0, 0)),
                             edge.width   = .5,
                             layout       = LAYOUT_SIGNED(g.clr)))

# -- Panel E: GLASSO network --
# Same layout as CLR network for direct visual comparison
set.seed(42)
p.graph.glasso <- as.grob(~plot(g.gl, vertex.label = NA,
                                vertex.color = colpal[taxa.filt[, "family"]],
                                vertex.size  = vertex.size,
                                edge.color   = ifelse(igraph::E(g.gl)$weight > 0,
                                                      rgb(0, 0, 1), rgb(1, 0, 0)),
                                edge.width   = .5,
                                layout       = LAYOUT_SIGNED(g.clr)))


# -------- SAVE FINAL OUTPUT --------

# -- Main figure: 6 panels (A-F), 3 columns x 2 rows --
png("../Plots/Methods_comparison.png", width = 3600, height = 2400, res = 300)
ggarrange(
  plotlist = list(p1, p2, p4,
                  p.graph.CLR, p.graph.glasso,
                  p.venn + theme(legend.position = "bottom")),
  labels = c("A", "B", "C", "D", "E", "F"),
  ncol = 3, nrow = 2
)
dev.off()

# -- Taxonomy color legend (family level) --
png(filename = "scripts/tmp_files/colorLegend.png", width = 1200, height = 1200, res = 200)
par(mar = c(2, 0, 2, 0))
plot.new()
legend("center", legend = names(colpal), fill = colpal, title = "Family", ncol = 2)
dev.off()

# -- Summary table: edge count and % positive/negative edges per method --
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

png(filename = "scripts/tmp_files/table.png", width = 1200, height = 1200, res = 300)
ggplot() +
  theme(axis.line = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        legend.position = "none", panel.background = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), plot.background = element_blank()) +
  annotation_custom(gridExtra::tableGrob(netw.info),
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
dev.off()

# -- Alternative composition using pre-saved images from disk --
# cowplot: advanced composition of ggplot2 graphics
# https://cran.r-project.org/web/packages/cowplot/index.html
p.mb  <- cowplot::ggdraw() + cowplot::draw_image("scripts/tmp_files/graph_mb.png")
p.gl  <- cowplot::ggdraw() + cowplot::draw_image("scripts/tmp_files/graph_gl.png")
p.clr <- cowplot::ggdraw() + cowplot::draw_image("scripts/tmp_files/graph_clr.png")
p.leg <- cowplot::ggdraw() + cowplot::draw_image("scripts/tmp_files/colorLegend.png")
p.ven <- cowplot::ggdraw() + cowplot::draw_image("scripts/tmp_files/ggvenn.png")
p.inf <- cowplot::ggdraw() + cowplot::draw_image("scripts/tmp_files/table.png")

png(filename = "outputs/Graph_Comparison.png", width = 1600, height = 1000, res = 300)
gridExtra::grid.arrange(p.mb, p.gl, p.clr, p.ven, p.leg, p.inf, ncol = 3)
dev.off()