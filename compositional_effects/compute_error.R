library(tidyverse)
library(foreach)
library(doParallel)
library(vegan)
library(VGAM)
library(MASS)
library(here)

# set Parameters, D=Dimensionality, M=Magnification Factor
D <- seq(5,200,by=50)
#M <- c(1:10,seq(12,98,by=2),seq(100,500,by=25),seq(550,1500,by=50),seq(2000,10000,by=500))
M <- c(1,2)

# elaborates Transformations Effects
#------------------------------------------------------------------------------#
# create a cluster, selecting all cores minus one, to allow parallel execution
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl) 

# set seed for reproducibility
set.seed(42)

# creates a data frame where to store the result of the simulation loop
df <-
  # for each (d, m), simulate data and compute error/pielou metrics in parallel
  foreach(d=D, .combine="rbind") %:%
  foreach(m=M, .combine="rbind", .packages=c("ToyModel")) %dopar% {
    
    # simulate ToyModel data with NorTA structure and given parameters
    toy <- toy_model(n=10^4, cor=diag(d), M=m,
                     qdist=qnorm, 
                     param=c(mean=0, sd=1),
                     method="pearson",
                     force.positive=TRUE)
    # compute correlation errors and mean Pielou evenness
    data.frame("d"=d, 
               "m"=m, 
               "ERR_L1" = mean(abs(toy$cor_NorTA - toy$cor_L1)),
               "ERR_CLR"= mean(abs(toy$cor_NorTA - toy$cor_CLR)),
               "pielou" = mean(apply(toy$NorTA,1,vegan::diversity) / log(d)))
  }

# stop the cluster to reduce resource use
stopCluster(cl)

# create the output folder if it does not exist
dir.create(here("compositional_effects"), showWarnings = FALSE)
# save the results to the crated (or already existing) output folder
saveRDS(df, here("compositional_effects", "Compositional_Effects_01.rds"))
