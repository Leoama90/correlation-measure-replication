# biases_parameters_example.R
#
# Purpose:
#   Visually demonstrate how compositional bias affects correlation
#   estimation under different normalization methods.
#
#   All synthetic datasets are generated with no true correlation
#   (cor = identity matrix), so any correlation pattern visible in the
#   plots is an artifact introduced by the normalization procedure.
#
#   The figure compares three settings:
#   - 5 dimensions with high evenness (~100% Pielou)
#   - 5 dimensions with medium evenness (~50% Pielou)
#   - 30 dimensions with medium evenness (~50% Pielou)
#
#   Each dataset is shown under three representations:
#   Normal, L1, and CLR.
#
# Outputs:
#   - High-resolution PNG figure saved to outputs/.
#
# here: easy locate files using top-level file directories
# https://cran.r-project.org/web/packages/here/vignettes/here.html
library(here)

# ToyModel: generates compositional toy datasets
# https://github.com/nome-repo/ToyModel
library(ToyModel)


# vegan: community ecology analysis
# https://cran.r-project.org/package=vegan
library(vegan)

# defines a function that computes the Pielou mean across all samples
mean_Pielou <- function(x){
  # divides Shannon diversity by log of number of species, then averages across samples
  P <- mean(apply(x, 1, vegan::diversity) / log(ncol(x)))
  # returns the result rounded to 2 decimal places
  return(round(P, 2))
}

# creates 5x5, 30x30 and 100x100 identity matrices to use as correlation structures
cor_D5   <- diag(5)
cor_D30  <- diag(30)
cor_D100 <- diag(100)

# sets the seed for reproducibility
set.seed(10)

# generates a toy dataset with 5 dimensions and high evenness (Pielou ~ 100%)
toy_D5_P100 <- toy_model(n      = 10^4, 
                         cor    = cor_D5, 
                         M      = 1, 
                         qdist  = qnorm,
                         param  = c("mean" = 0, "sd" = 1), force.positive = T, 
                         method = "pearson")

# generates a toy dataset with 5 dimensions and medium evenness (Pielou ~ 50%)
toy_D5_P50  <- toy_model(n      = 10^4, 
                         cor    = cor_D5, 
                         M      = 15.5, 
                         qdist  = qnorm,
                         param  = c("mean" = 0, "sd" = 1), force.positive = T, 
                         method = "pearson")

# generates a toy dataset with 30 dimensions and medium evenness (Pielou ~ 50%)
toy_D30_P50 <- toy_model(n      = 10^4, 
                         cor    = cor_D30, 
                         M      = 64, 
                         qdist  = qnorm,
                         param  = c("mean" = 0, "sd" = 1), force.positive = T, 
                         method = "pearson")

# creates the outputs folder at the project root if it does not already exist
dir.create(here("outputs"), showWarnings = FALSE)

# opens a PNG graphics device and sets the output file path, size and resolution
png(filename = here("Plots", "biases_example.png"), 
    width    = 8000, 
    height   = 4400, 
    res      = 600)

# arranges the plot area into a 4x3 grid, with the first row used as a header
layout(matrix(c(1:12), 
              nrow    = 4, 
              byrow   = TRUE), 
              heights = c(0.1, 0.3, 0.3, 0.3))
# sets plot margins
par(mar = c(0, 4, 1, 4))

# adds a blank plot with the column header "Generated"
plot.new()
text(x    = 0.5, 
     y    = 0.5, 
     "Generated", 
     cex  = 2, 
     font = 2)

# adds a blank plot with the column header "L1"
plot.new()
text(x    = 0.5, 
     y    = 0.5, 
     "L1", 
     cex  = 2, 
     font = 2)

# adds a blank plot with the column header "CLR"
plot.new()
text(x    = 0.5, 
     y    = 0.5, 
     "CLR", 
     cex  = 2, 
     font = 2)

# plots row A: 5 dimensions, high evenness, under Normal, L1 and CLR normalization
plot(toy_D5_P100, "Normal", vertex.label = NA)
text(x   = -1, 
     y   = -1, 
     "A", 
     cex =  2)
plot(toy_D5_P100, "L1",     vertex.label = NA)
plot(toy_D5_P100, "CLR",    vertex.label = NA)

# plots row B: 5 dimensions, medium evenness, under Normal, L1 and CLR normalization
plot(toy_D5_P50, "Normal",  vertex.label = NA)
text(x   = -1, 
     y   = -1, 
     "B", 
     cex =  2)
plot(toy_D5_P50, "L1",      vertex.label = NA)
plot(toy_D5_P50, "CLR",     vertex.label = NA)

# plots row C: 30 dimensions, medium evenness, under Normal, L1 and CLR normalization
plot(toy_D30_P50, "Normal", vertex.label = NA)
text(x   = -1, 
     y   = -1,
     "C", 
     cex = 2)
plot(toy_D30_P50, "L1",     vertex.label = NA)
plot(toy_D30_P50, "CLR",    vertex.label = NA)

# allows drawing outside the plot area
par(xpd  = NA)

# draws horizontal lines to visually separate the rows
abline(h = 6.5 , col = "darkgreen", lty = 3)
abline(h = 3.75, col = "darkgreen", lty = 3)
abline(h = 1.25, col = "darkgreen", lty = 3)

# UI reminder where to search the generated plot
cat("the plot has been saved in the 'Plots' folder with the name 'biases_example.png' ")

# closes the graphics device and saves the PNG file to disk
dev.off()

# opens the output PNG in the default image viewer once the script completes
browseURL(list.files(
          path       = here(),
          pattern    = "biases_example.png",
          full.names = TRUE,
          recursive  = TRUE
          ))