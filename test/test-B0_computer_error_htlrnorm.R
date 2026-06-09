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

# -------- test the function that import datas --------
otu  <- readRDS(here("data", "otu_HMP2.rds" ))
meta <- readRDS(here("data", "meta_HMP2.rds"))

test_that("loaded object is not correct (otu should be a matrix, meta should be a list)!", {
  expect_true(is.matrix(otu))
  expect_true(is.list(meta ))
})

# -------- test the filtered data from the script --------

# select samples belonging to subject 69-001 in healthy status
otu.69001.H <- otu[meta$SubjectID == "69-001" & meta$CL4_2 == "Healthy", ]
# remove rarest OTUs: keep only those present in at least 33% of samples
otu.filt    <- otu.69001.H[, colSums(otu.69001.H > 0)/nrow(otu.69001.H) >= 0.33]

# test that the filter reduced the column number of data
test_that("Otu.filt has more columns than otu.69001.H!", {
  expect_gte(ncol(otu.69001.H), ncol(otu.filt))
})

# keep only OTUs whose median non-zero abundance is at least 5
otu.filt    <- otu.filt[, apply(otu.filt, 2, function(x) median(x[x > 0]) >= 5)]

# test that all matrix entries have median >= 5
test_that("Some OTU has median non-zero abundance < 5!", {
  medians   <- apply(otu.filt, 2, function(x) median(x[x > 0]))
  expect_true(all(medians >= 5))
})

