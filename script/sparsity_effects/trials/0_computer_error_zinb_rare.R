# doSNOW: parallel backend for foreach, supports progress bars via snow clusters
# https://cran.r-project.org/web/packages/doSNOW/index.html
library(doSNOW)

# foreach: parallel foreach loops
# https://cran.r-project.org/package=foreach
library(foreach)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# progress: displays text progress bars for long-running loops
# https://cran.r-project.org/web/packages/progress/index.html
library(progress)

# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)

# ToyModel: lightweight package with necessary function to generate metagenomics data
# https://github.com/Fuschi/ToyModel
library(ToyModel)

# VGAM: vector generalized linear models
# https://cran.r-project.org/package=VGAM
library(VGAM)

# SpiecEasi: package with spiec.easi and sparCC methods
# https://github.com/zdk123/SpiecEasi
library(SpiecEasi)


# -------- read data --------

# Load OTU abundance matrix (samples x OTUs)
otu <- readRDS(list.files(
  path       = here(),
  pattern    = "otu_HMP2.rds",
  full.names = TRUE,
  recursive  = TRUE
))
# Load sample metadata
meta <- readRDS(list.files(
  path       = here(),
  pattern    = "meta_HMP2.rds",
  full.names = TRUE,
  recursive  = TRUE
))


# -------- filter data --------

# Subset samples belonging to subject 69-001 in the healthy state
otu_69001_H <- otu[meta$SubjectID == "69-001" & meta$CL4_2 == "Healthy", ]

# Remove OTUs present in fewer than 33% of samples (low prevalence)
otu_filt <- otu_69001_H[, colSums(otu_69001_H > 0)/nrow(otu_69001_H)   >= 0.33]
# Further remove OTUs whose median non-zero count is below 5 (too rare)
otu_filt <- otu_filt[, apply(otu_filt, 2, function(x) median(x[x > 0]) >= 5)]
# Free memory by removing objects no longer needed
rm(otu, otu_69001_H, meta)


# Fit ZINB params from real data for each filtered OTUs and
# save the first and ninen-th decil of the distribution
HMP2_quantile_params <- apply(otu_filt, 2, function(x){
    fitdistr(as.numeric(x), "zinegbin")$par
}) %>% apply(1, function(x) quantile(x, probs = c(0.1, 0.9)))

# Create the progression bar
nIteration <- 100
pb <- progress_bar$new(
  format = "[:bar] :elapsed | eta: :eta",
  total  = nIteration * 7600,
  width  = 60)
progress <- function(n){pb$tick()}

set.seed(42)
result <- data.frame()
for(iter in 1:nIteration){
  
  
# -------- START OUTER LOOP --------
  
# generation of the underline dataset
params_random_HMP2 <- data.frame(
  "munb"  = runif(n = 200,
                min = HMP2_quantile_params["10%", "munb"],
                max = HMP2_quantile_params["90%", "munb"]),
  "size"  = runif(n = 200,
                min = HMP2_quantile_params["10%", "size"],
                max = HMP2_quantile_params["90%", "size"]),
  "pstr0" = runif(n = 200,
                min = HMP2_quantile_params["10%", "pstr0"],
                max = HMP2_quantile_params["90%", "pstr0"]))

# generate a dummy dataset
random_HMP2 <- toy_model(n     = 10^4, 
                         cor   = diag(200), 
                         M     = 1,
                         qdist = VGAM::qzinegbin,
                         param = params_random_HMP2)
random_cor0_HMP2 <- random_HMP2$cor_normal

# create parameters of the supervised couple of variables
params_set <- expand_grid(
  "pstr0_1" = seq( 0  , 0.95, by = 0.05),
  "pstr0_2" = seq( 0  , 0.95, by = 0.05),
  "cor"     = seq(-0.9, 0.9 , by = 0.1 )
) %>%
  mutate("munb_1" = runif(n = n(),
                        min = HMP2_quantile_params["10%", "munb"],
                        max = HMP2_quantile_params["90%", "munb"]),
         "munb_2" = runif(n = n(),
                        min = HMP2_quantile_params["10%", "munb"],
                        max = HMP2_quantile_params["90%", "munb"]),
         "size_1" = runif(n = n(),
                        min = HMP2_quantile_params["10%", "size"],
                        max = HMP2_quantile_params["90%", "size"]),
         "size_2" = runif(n = n(),
                        min = HMP2_quantile_params["10%", "size"],
                        max = HMP2_quantile_params["90%", "size"])) %>%
  as.data.frame()

# create the cluster for parallel execution
cl <- makeCluster(6)
registerDoSNOW(cl)
  
  
  # -------- START INNER LOOP --------
  
  df <- foreach(i = 1:nrow(params_set), 
                .combine  = "rbind",
                .packages = c("ToyModel"),
                .options.snow = list(progress = progress)) %dopar% {
                 
                 # Generate another dummy dataset with the specified ZINB params and correlation  
                 couple <- toy_model(n       = 10^4, 
                                     cor     = params_set[i, "cor"], 
                                     M       = 1,
                                     qdist   = VGAM::qzinegbin,
                                     param   = data.frame(
                                     "munb"  = c(params_set[i, "munb_1"], params_set[i, "munb_2"]),
                                     "size"  = c(params_set[i, "size_1"], params_set[i, "size_2"]),
                                     "pstr0" = c(params_set[i, "pstr0_1"], params_set[i, "pstr0_2"])
                                     ))
                  
                  random_HMP2_NorTA_i        <- random_HMP2$NorTA
                  random_HMP2_NorTA_i[, 25]  <- couple$NorTA[, 1]
                  random_HMP2_NorTA_i[, 125] <- couple$NorTA[, 2]
                  
                  cor_PCLR <- random_HMP2_NorTA_i %>% ToyModel::clr() %>% cor
                  
                  data.frame("iteration" = iter,
                             "cor_input" = params_set[i, "cor"],
                             "cor_normal" = couple$cor_normal[1, 2],
                             "munb_1" = params_set[i, "munb_1"],
                             "munb_2" = params_set[i, "munb_2"],
                             "size_1" = params_set[i, "size_1"],
                             "size_2" = params_set[i, "size_2"],
                             "pstr0_1" = params_set[i, "pstr0_1"],
                             "pstr0_2" = params_set[i, "pstr0_2"],
                             "cor_NorTA_PCLR" = cor_PCLR[25, 125])
                }
  
  # -------- END INNER LOOP --------
  
# close parallel cluster to release resources
stopCluster(cl)
# Append this iteration's results to the already existing data frame
result <- bind_rows(result, df)


# -------- END OUTER LOOP --------

}

# save results in the appropriate folder
saveRDS(result, here("script", "sparsity_effects", "trials", "Sparsity_Effects_zinbin.rds"))

# This script investigates how sparsity (excess zeros) in metagenomic data distorts 
# correlation estimates between OTUs. 
# Using real HMP2 data to calibrate realistic ZINB parameters, 
# it repeatedly generates a 200-OTU synthetic background 
# community and embeds a focal pair of OTUs with known correlation and 
# varying zero-inflation levels. 
# For each scenario, it computes the CLR-based correlation of the focal pair within the full community, allowing a 
# direct comparison between the true input correlation and the one recovered after compositional transformation.