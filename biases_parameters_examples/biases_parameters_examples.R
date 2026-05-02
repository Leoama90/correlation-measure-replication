# load required libraries
library(ToyModel)
library(vegan)
library(here)

# defines a function that computes the Pielou mean across all samples
mean_Pielou <- function(x){
  # divides Shannon diversity by log of number of species, then averages across samples
  P <- mean( apply(x,1,vegan::diversity) / log(ncol(x)))
  # returns the result rounded to 2 decimal places
  return(round(P,2))
}

# creates 5x5, 30x30 and 100x100 identity matrices to use as correlation structures
cor_D5 <- diag(5)
cor_D30 <- diag(30)
cor_D100 <- diag(100)

# sets the seed
set.seed(10)

# generates a toy dataset with 5 dimensions and high evenness (Pielou ~ 100%)
toy_D5_P100 <- ToyModel::toy_model(n=10^4, cor=cor_D5, M=1, qdist=qnorm,
                                   param=c("mean"=0, "sd"=1), force.positive=T, 
                                   method="pearson")

# generates a toy dataset with 5 dimensions and medium evenness (Pielou ~ 50%)
toy_D5_P50 <- ToyModel::toy_model(n=10^4, cor=cor_D5, M=15.5, qdist=qnorm,
                                  param=c("mean"=0, "sd"=1), force.positive=T, 
                                  method="pearson")

# generates a toy dataset with 30 dimensions and medium evenness (Pielou ~ 50%)
toy_D30_P50 <- ToyModel::toy_model(n=10^4, cor=cor_D30, M=64, qdist=qnorm,
                                   param=c("mean"=0, "sd"=1), force.positive=T, 
                                   method="pearson")

# creates the outputs folder at the project root if it does not already exist
dir.create(here("outputs"), showWarnings = FALSE)

# opens a PNG graphics device and sets the output file path, size and resolution
png(filename=here("outputs", "biases_example.png"), width=8000, height=4400, res=600)

# arranges the plot area into a 4x3 grid, with the first row used as a header
layout(matrix(c(1:12), nrow=4, byrow=T), heights=c(.1,.3,.3,.3))
# sets plot margins
par(mar=c(0,4,1,4))

# adds a blank plot with the column header "Generated"
plot.new()
text(x=.5, y=.5, "Generated", cex=2, font=2)

# adds a blank plot with the column header "L1"
plot.new()
text(x=.5, y=.5, "L1", cex=2, font=2)

# adds a blank plot with the column header "CLR"
plot.new()
text(x=.5, y=.5, "CLR", cex=2, font=2)

# plots row A: 5 dimensions, high evenness, under Normal, L1 and CLR normalization
plot(toy_D5_P100, "Normal", vertex.label=NA)
text(x=-1,y=-1,"A", cex=2)
plot(toy_D5_P100, "L1", vertex.label=NA)
plot(toy_D5_P100, "CLR", vertex.label=NA)

# plots row B: 5 dimensions, medium evenness, under Normal, L1 and CLR normalization
plot(toy_D5_P50, "Normal", vertex.label=NA)
text(x=-1,y=-1,"B", cex=2)
plot(toy_D5_P50, "L1", vertex.label=NA)
plot(toy_D5_P50, "CLR", vertex.label=NA)

# plots row C: 30 dimensions, medium evenness, under Normal, L1 and CLR normalization
plot(toy_D30_P50, "Normal", vertex.label=NA)
text(x=-1,y=-1,"C", cex=2)
plot(toy_D30_P50, "L1", vertex.label=NA)
plot(toy_D30_P50, "CLR", vertex.label=NA)

# allows drawing outside the plot area
par(xpd = NA)
# draws horizontal lines to visually separate the rows
abline