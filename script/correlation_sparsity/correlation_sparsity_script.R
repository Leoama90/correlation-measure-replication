# correlation_sparsity.R
#
# Purpose:
#   Investigate how the structure of the true correlation matrix affects
#   CLR compositional bias, comparing two scenarios: block-structured
#   correlations and random Gaussian correlations.
#
# Input:
#   - No external data files.
#   - The script generates all required inputs internally by:
#     * constructing block correlation matrices over a grid of parameters
#     * generating dense random Gaussian correlation matrices
#     * simulating compositional data with ToyModel
#
#   In the block experiment, the input parameter grid is:
#   - Dimension = 500
#   - taxa_connected = 200 down to 0
#   - correlation_value = 0.5 to 0.9 by 0.1
#
#   In the Gaussian experiment, the script generates 100 random dense
#   correlation matrices from symmetric Gaussian random matrices.
#
# Outputs:
#   - results_block.rds
#   - results_gauss.rds
#   - correlation_density_biases.png
#
#   The final figure combines:
#   - Panel A: CLR MAE distribution across Gaussian correlation matrices
#   - Panel B: CLR MAE vs edge density for block-structured correlations
# doSNOW: parallel backend for foreach, supports progress bars via snow clusters
# https://cran.r-project.org/web/packages/doSNOW/index.html
library(doSNOW)

# foreach: provides the %dopar% operator for parallel for-loops
# https://cran.r-project.org/web/packages/foreach/index.html
library(foreach)

# here: project-oriented file paths
# [https://cran.r-project.org/package=here](https://cran.r-project.org/package=here)
library(here)

# mvtnorm: generates multivariate normal and t distributions
# https://cran.r-project.org/web/packages/mvtnorm/index.html
library(mvtnorm)

# Matrix: tools for working with dense and sparse matrices, including nearPD()
# https://cran.r-project.org/web/packages/Matrix/index.html
library(Matrix)

# parallel: base R package for creating and managing parallel clusters
# https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf
library(parallel)

# tidyverse: collection of packages for data manipulation and visualization (ggplot2, dplyr, ...)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)

# ToyModel: simulates microbiome-like compositional data with known correlation structure
# https://github.com/MaStatLab/ToyModel
library(ToyModel)


# -------- different densities --------

# set up the parallel backend with doSNOW
num_cores <- 4
cl <- makeCluster(num_cores)
registerDoSNOW(cl)

# define the parameters
Dimension <- 500
Connected_Taxa_Sequence <- seq(200, 0, - 1)
Correlation_Sequence    <- seq(0.5, 0.9, 0.1)

# total possible edges
total_possible_edges    <- Dimension * (Dimension - 1) / 2

# set up progress bar
total_iterations <- length(Correlation_Sequence) * length(Connected_Taxa_Sequence)
pb <- txtProgressBar(max = total_iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# perform parallel computation with progress tracking
results_block <-  foreach(corr_val = Correlation_Sequence, .combine = "rbind") %:%
  foreach(n_taxa_connected = Connected_Taxa_Sequence, .combine = "rbind", .packages = c("mvtnorm", "Matrix", "ToyModel"), .options.snow = opts) %dopar% {
    
    # start with a diagonal matrix of ones
    corM <- diag(1, Dimension, Dimension)
    
    # get the current upper triangle values excluding the diagonal
    corM[1:n_taxa_connected, 1:n_taxa_connected] <- corr_val
    
    # number of edges currently connected
    num_edges <- sum(corM[lower.tri(corM)] != 0)
    edge_density <- num_edges / total_possible_edges
    
    # ensure the matrix is positive semidefinite
    corM_PD <- nearPD(corM, corr = TRUE, keepDiag = TRUE)$mat
    corM_PD <- as.matrix(corM_PD)
    
    toy <- toy_model(n      = 10^4, 
                     cor    = corM_PD, 
                     M      = 1,
                     qdist  = qnorm, 
                     param  = c(mean = 0, sd = 1),
                     method = "pearson",
                     force.positive = TRUE)
    
    data.frame("d"                 = Dimension, 
               "taxa_connected"    = n_taxa_connected,
               "edge_density"      = edge_density,
               "correlation_value" = corr_val,
               "ERR_CLR"           = mean(abs(toy$cor_NorTA - toy$cor_CLR)))
  }

# close progress bar and stop the cluster
close(pb)
stopCluster(cl)

saveRDS(results_block, here("script", "correlation_sparsity", "results_block.rds"))

# -------- gaussian distributed correlation values --------

# set up the parallel backend with doSNOW
num_cores <- 4
cl <- makeCluster(num_cores)
registerDoSNOW(cl)

# set up progress bar
total_iterations <- 100
pb       <- txtProgressBar(max = total_iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts     <- list(progress = progress)

# -------- gaussian random values --------

n <- Dimension
results_gauss <- foreach(i = 1:100, 
                         .combine  = "rbind", 
                         .packages = c("mvtnorm", "Matrix", "ToyModel"), 
                         .options.snow = opts) %dopar% {
  
  # step 1: generate a random matrix with gaussian entries
  random_matrix <- matrix(rnorm(n * n), n, n)
  
  # step 2: make the matrix symmetric
  symmetric_matrix <- (random_matrix + t(random_matrix)) / 2
  
  # step 3: standardize the diagonal to make it a correlation matrix
  diag(symmetric_matrix) <- 1
  
  # step 4: ensure the matrix is positive semidefinite
  eigen_values <- eigen(symmetric_matrix)
  # set any negative eigenvalues to a small positive value (e.g., 1e-10)
  eigen_values$values[eigen_values$values < 0] <- 1e-10
  # reconstruct the matrix
  correlation_matrix <- eigen_values$vectors %*% diag(eigen_values$values) %*% t(eigen_values$vectors)
  
  # verify the matrix is now a valid correlation matrix
  correlation_matrix <- nearPD(correlation_matrix, corr = TRUE, keepDiag = TRUE)$mat
  correlation_matrix <- as.matrix(correlation_matrix)
  
  toy <- toy_model(n      = 10^4, 
                   cor    = correlation_matrix, 
                   M      = 1,
                   qdist  = qnorm, 
                   param  = c(mean = 0, sd = 1),
                   method = "pearson",
                   force.positive = TRUE)
  
  data.frame("d"       = n,
             "i"       = i,
             "ERR_CLR" = mean(abs(toy$cor_NorTA - toy$cor_CLR)))
  
}

# save results
saveRDS(results_gauss, here("script", "correlation_sparsity", "results_gauss.rds"))


# create the plot with improved aesthetics
p_gauss <- results_gauss %>%
  ggplot(aes(x = factor(1), 
             y = ERR_CLR)) +
  
  geom_violin(fill   = "#69b3a2", 
              color  = "#1b4f4a", 
              alpha  = 0.5) +
  
  geom_boxplot(width = 0.1, fill = "white", 
               color = "#2a2a2a", 
               alpha = 0.9, outlier.shape = NA) +  
  geom_jitter(width  = 0.15, color = "#fc9272", alpha = 0.6, size = 2) +
  
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.x  = element_blank(),
        plot.title   = element_text(hjust = 0.5),
        legend.position = "none") +
  
  labs(y = "MAE") +
  
  geom_hline(yintercept = median(results_gauss$ERR_CLR), 
             linetype   = "dotdash", 
             color      = "darkgoldenrod1", 
             linewidth  = 0.75)
#annotate("text", x = .2, y = Inf, label = "A", hjust = -0.1, vjust = 1.1, size = 6, color = "black")
print(p_gauss)

p_block <- results_block %>%
  filter(edge_density <=.10) %>%
  filter(correlation_value %in% c(0.5, 0.7, 0.9)) %>%
  mutate(edge_density = 100*edge_density) %>%
  ggplot(aes(x     = edge_density, 
             y     = ERR_CLR, 
             color = as.factor(correlation_value))) +
  
  geom_point() +
  
  ylab("MAE") +
  
  xlab("Correlation Density") +
  
  theme_bw() +
  
  theme(legend.position = "right") +
  
  labs(color = "Correlation Value") +
  
  xlab("Edge Density") +
  #annotate("text", x = 1, y = Inf, label = "B", hjust = 1.1, vjust = 1.1, size = 6, color = "black") +
  
  scale_x_continuous(labels = scales::label_percent(scale = 1)) +
  
  geom_hline(yintercept = median(results_gauss$ERR_CLR), 
             linetype   = "dotdash", 
             color      = "darkgoldenrod1", 
             linewidth  = 1.5) +
  
  annotate("text", 
           x     =   10, 
           y     =   median(results_gauss$ERR_CLR), 
           label =   sprintf("median mae from gaussian-weighted\ncorrelation matrices ≈ %.0e", median(results_gauss$ERR_CLR)), 
           hjust =   1.1, 
           vjust = - 0.1, 
           size  =   4, 
           color =   "darkgoldenrod4") +
  
  labs(y     = "MAE", 
       x     = "Edge Density (%)", 
       color = "Correlation Value")
print(p_block)

# arrange together multiple ggplots (p_gauss and p_block in this case)
pall <- ggpubr::ggarrange(
  
  p_gauss, p_block, widths = c(0.3, 0.7), labels = c("A", "B")
  
)

# absolute path to Plots/ relative to project root
plots_dir <- here("Plots")

# create the folder if it doesn't exist yet
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# generate the plot
png(filename = here("Plots", "correlation_density_biases.png"), 
    width    = 2800, 
    height   = 1500, 
    res      = 300)  
print(pall)

# UI reminder where to search the generated plot
cat("the plot has been saved in the 'Plots' folder with the name 'correlation_density_biases.png' ")


dev.off()