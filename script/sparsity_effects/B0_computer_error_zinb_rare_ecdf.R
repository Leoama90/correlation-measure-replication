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

# SpiecEasi: package with spiec.easi and sparCC methods
# https://github.com/zdk123/SpiecEasi
library(SpiecEasi)

# tidyverse: useful to manage data (dplyr) and make nice plots (ggplot2)
# https://cran.r-project.org/web/packages/tidyverse/index.html
library(tidyverse)

# ToyModel: lightweight package with necessary function to generate metagenomics data
# https://github.com/Fuschi/ToyModel
library(ToyModel)

# VGAM: vector generalized linear models
# https://cran.r-project.org/package=VGAM
library(VGAM)


# -------- read data --------

# load OTU abundance matrix (samples x OTUs)
otu  <- readRDS(here("data", "otu_HMP2.rds" ))
# load sample metadata
meta <- readRDS(here("data", "meta_HMP2.rds"))

# select samples belonging to subject 69-001 in healthy status
otu_69001_H <- otu[meta$SubjectID == "69-001" & meta$CL4_2 == "Healthy", ]

# remove rarest OTUs: keep only OTUs present in at least 33% of samples (prevalence filter)
otu_filt <- otu_69001_H[, colSums(otu_69001_H > 0) / nrow(otu_69001_H) >= .33]
# further filter: keep only OTUs whose median non-zero abundance is >= 5 (abundance filter)
otu_filt <- otu_filt[, apply(otu_filt, 2, function(x) median(x[x > 0]) >= 5)]
# rarefy: rescale all samples to the same total read count (minimum library size)
otu_filt <- round(otu_filt / rowSums(otu_filt) * min(rowSums(otu_filt)))
# remove objects no longer needed to free memory
rm(otu, otu_69001_H, meta)

# fit ZINB (zero-inflated negative binomial) parameters from real data for each filtered OTU;
# returns a matrix with one row per OTU and columns: munb, size, pstr0
HMP2_params <- apply(otu_filt, 2, function(x) {
  SpiecEasi::fitdistr(as.numeric(x), "zinegbin")$par
}) %>%
  t() %>%
  as.data.frame()

# compute the 10th and 90th percentile of each ZINB parameter across OTUs,
# used to trim extreme parameter values before simulation
HMP2_quantile_params <- HMP2_params %>%
  apply(2, function(x) quantile(x, probs = c(.1, .9)))

# retain only OTUs whose fitted parameters fall within the [10%, 90%] range
# for all three parameters (munb, size, pstr0), removing outlier OTUs
HMP2_params_filt <- HMP2_params %>%
  filter(
    munb  >= HMP2_quantile_params["10%", "munb" ],
    munb  <= HMP2_quantile_params["90%", "munb" ],
    size  >= HMP2_quantile_params["10%", "size" ],
    size  <= HMP2_quantile_params["90%", "size" ],
    pstr0 >= HMP2_quantile_params["10%", "pstr0"],
    pstr0 <= HMP2_quantile_params["90%", "pstr0"]
  )

# initialize a progress bar tracking total iterations across inner and outer loops
nIteration <- 100
pb <- progress_bar$new(
  format = "[:bar] :elapsed | eta: :eta",
  total  = nIteration * 380,  
  width  = 60
)
# callback passed to doSNOW so each parallel worker advances the progress bar by one tick
progress <- function(n) { pb$tick() }

set.seed(42)
result <- data.frame()

for (iter in 1:nIteration) {
  
  # -------- START OUTER LOOP --------
  
  # sample 200 random parameter combinations from the filtered HMP2 ZINB distributions
  # using uniform random quantiles, then simulate a 10,000-sample background dataset
  # with no correlation structure (identity correlation matrix)
  params_random_HMP2 <- data.frame(
    "munb"  = quantile(HMP2_params_filt$munb,  probs = runif(runif(200))),
    "size"  = quantile(HMP2_params_filt$size,  probs = runif(runif(200))),
    "pstr0" = quantile(HMP2_params_filt$pstr0, probs = runif(runif(200)))
  )
  
  random_HMP2 <- ToyModel::toy_model(
    n     = 10^4,   
    cor   = diag(200),  
    M     = 1,
    qdist = VGAM::qzinegbin,
    param = params_random_HMP2
  )
  # store the underlying normal correlation matrix for later reference
  random_cor0_HMP2 <- random_HMP2$cor_normal
  
  
  # -------- create parameters of the supervised couple of variables --------
  
  # build a grid of all combinations of input correlation and zero-inflation (pstr0),
  # then randomly draw munb and size parameters from the filtered HMP2 empirical distribution
  params_set <- expand_grid(
    "cor"     = seq(-.9, .9, by = .1),   
    "pstr0_1" = seq(0, .95, by = .05)    
  ) %>%
    mutate(
      pstr0_2 = pstr0_1, .after = pstr0_1,        
      munb_1  = quantile(HMP2_params_filt$munb, probs = runif(n())),
      munb_2  = quantile(HMP2_params_filt$munb, probs = runif(n())),
      size_1  = quantile(HMP2_params_filt$size, probs = runif(n())),
      size_2  = quantile(HMP2_params_filt$size, probs = runif(n()))
    ) %>%
    as.data.frame()
  
  # create a SNOW cluster with 6 workers for parallel execution of the inner loop
  cl <- makeCluster(6)
  registerDoSNOW(cl)
  
  
  # -------- START INNER LOOP --------
  
  
  # iterate over all parameter combinations in parallel;
  # for each combination, simulate a supervised OTU pair, embed it into the
  # background dataset, and compute the CLR-transformed Pearson correlation
  df <- foreach(
    i             = 1:nrow(params_set),
    # stack results into a single data frame
    .combine      = "rbind",      
    # packages required on each worker
    .packages     = c("ToyModel"),  
    .options.snow = list(progress = progress)
  ) %dopar% {
    
    # simulate a pair of OTUs with the specified correlation and ZINB parameters
    couple <- ToyModel::toy_model(
      n     = 10^4,
      cor   = params_set[i, "cor"],
      M     = 1,
      qdist = VGAM::qzinegbin,
      param = data.frame(
        "munb"  = c(params_set[i,  "munb_1"], params_set[i, "munb_2" ]),
        "size"  = c(params_set[i,  "size_1"], params_set[i, "size_2" ]),
        "pstr0" = c(params_set[i, "pstr0_1"], params_set[i, "pstr0_2"])
      )
    )
    
    # embed the simulated OTU pair into the background dataset at positions 25 and 125,
    # replacing the original OTUs with the supervised couple
    random_HMP2_NorTA_i         <- random_HMP2$NorTA
    random_HMP2_NorTA_i[, 25]   <- couple$NorTA[, 1]
    random_HMP2_NorTA_i[, 125]  <- couple$NorTA[, 2]
    
    # apply CLR (centred log-ratio) transformation and compute the full correlation matrix;
    # the entry [25, 125] gives the PCLR correlation between the supervised OTU pair
    cor_PCLR <- random_HMP2_NorTA_i %>%
      ToyModel::clr() %>%
      cor()
    
    # return one row per parameter combination with inputs and the estimated correlation
    data.frame(
      "iteration"      = iter,
      "cor_input"      = params_set[i, "cor"    ],  
      "cor_normal"     = couple$cor_normal[1, 2 ],  
      "munb_1"         = params_set[i, "munb_1" ],
      "munb_2"         = params_set[i, "munb_2" ],
      "size_1"         = params_set[i, "size_1" ],
      "size_2"         = params_set[i, "size_2" ],
      "pstr0_1"        = params_set[i, "pstr0_1"],
      "pstr0_2"        = params_set[i, "pstr0_2"],
      "cor_NorTA_PCLR" = cor_PCLR[25, 125]          
    )
  }
  
  
  # -------- END INNER LOOP --------
  
  
  # shut down the parallel cluster to release resources
  stopCluster(cl)
  # accumulate results from this outer iteration
  result <- bind_rows(result, df)
  
  
  # -------- END OUTER LOOP --------
  
}

# save the full simulation results to disk for downstream analysis
saveRDS(result, here("script", "sparsity_effects", "Sparsity_Effects_zinbin_rare_ecdf.rds"))