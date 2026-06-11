# doParallel: parallel backend for foreach
# [https://cran.r-project.org/package=doParallel](https://cran.r-project.org/package=doParallel)
library(doParallel)

# foreach: parallel foreach loops
# [https://cran.r-project.org/package=foreach](https://cran.r-project.org/package=foreach)
library(foreach)

# here: project-oriented file paths
# [https://cran.r-project.org/package=here](https://cran.r-project.org/package=here)
library(here)

# MASS: statistical modeling and analysis
# [https://cran.r-project.org/package=MASS](https://cran.r-project.org/package=MASS)
library(MASS)

# tidyverse: data science toolkit (dplyr, ggplot2, tidyr, etc.)
# [https://tidyverse.org](https://tidyverse.org)
library(tidyverse)

# VGAM: vector generalized linear models
# [https://cran.r-project.org/package=VGAM](https://cran.r-project.org/package=VGAM)
library(VGAM)


# -------- set parameters --------

dimensions    <- seq(5, 200, by = 50)
magnification <- c(1, 2)


# -------- elaborates transformations effects --------

# create a cluster, selecting all cores minus one, to allow parallel execution
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# set seed for reproducibility
set.seed(42)

# creates a data frame where to store the result of the simulation loop
simulation_results <-
  # for each (d, m), simulate data and compute error/pielou metrics in parallel
  foreach(d = dimensions,    .combine = "rbind") %:%
  foreach(m = magnification, .combine = "rbind", .packages = c("ToyModel")) %dopar% {
    # simulate ToyModel data with NorTA structure and given parameters
    toy <- toy_model(n      = 10^4,
                     cor    = diag(d),
                     M      = m,
                     qdist  = qnorm,
                     param  = c(mean = 0, sd = 1),
                     method = "pearson",
                     force.positive = TRUE)
    # compute correlation errors and mean Pielou evenness
    data.frame("d" = d,
               "m" = m,
               "ERR_L1"  = mean(abs(toy$cor_NorTA - toy$cor_L1 )),
               "ERR_CLR" = mean(abs(toy$cor_NorTA - toy$cor_CLR)),
               "Pielou"  = mean(apply(toy$NorTA, 1, vegan::diversity) / log(d)))
  }

# stop the cluster to reduce resource use
stopCluster(cl)


# -------- save results --------

# create output folder
dir.create(here("compositional_effects"), showWarnings = FALSE)

# save the results in output folder
saveRDS(simulation_results, here("compositional_effects", "compositional_effects_01.rds"))
