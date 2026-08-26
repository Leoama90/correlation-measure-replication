# demo_clr_pearson_old.R
#
# Purpose:
#   Demonstrates the usage of the function "clr_on_data()" (from the script 
#   "clr_pearson"), loading a real metagenomic
#   OTU count table and running the full filtering + pseudocount + CLR +
#   Pearson correlation pipeline on it.
#
# Inputs:
#   - clr_pearson.R (sourced below)
#   - otu_HMP2.rds: a metagenomic OTU count table (samples x OTUs)
#
# Outputs (saved in the "outputs/" folder at the project root):
#   - outputs/samp_filt.rds: the filtered OTU table (tibble)
#   - outputs/y_clr.rds: the CLR-transformed OTU table (matrix)
#   - outputs/cor_matrix.rds: the OTU x OTU Pearson correlation matrix (matrix)
#
# Used scripts:
#   - clr_pearson.R


# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# tidyverse: data science toolkit (dplyr, ggplot2, tidyr, etc.)
# https://tidyverse.org
library(tidyverse)


# -------- load the clr_pearson.R script --------

# load clr_on_data() and its dependencies (datasum, filt_data, pseudocount)
source(
  list.files(
    path = here(),
    pattern = "^clr_pearson\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- read raw data and run the pipeline --------

# Load the OTU abundance matrix (samples x OTUs) from its .rds file
otu <- readRDS(
  list.files(
    path = here(),
    pattern = "^otu_HMP2\\.rds$",
    full.names = TRUE,
    recursive = TRUE
  )
)
# run the full pipeline: filtering, pseudocount, CLR, Pearson correlation
result <- clr_on_data(otu)


# -------- save each element of the result list as a separate .rds file --------

# create the outputs folder once, before the loop, if it does not exist yet
dir.create(here("outputs"), showWarnings = FALSE)
# each file is named after the element itself (e.g. samp_filt.rds, y_clr.rds, ...)
for (name in names(result)) {
  saveRDS(result[[name]], here("outputs", paste0(name, ".rds")))
}


# -------- report to the user which files were saved and where --------

# the list of saved files is built dynamically from names(result), so this
# message stays correct even if clr_on_data() changes in the future
cat("\n")
cat("The pipeline result was saved as", length(result), "separate .rds files ")
cat("in the 'outputs/' folder at the project root:\n")
for (name in names(result)) {
  cat("  - outputs/", paste0(name, ".rds"), "\n", sep = "")
}


# -------- remind the user what each file contains --------
cat("\n")
cat("Reminder of what each file contains:\n")
cat("  - samp_filt: the filtered OTU table (tibble)\n")
cat("  - y_clr: the CLR-transformed OTU table (matrix)\n")
cat("  - cor_matrix: the OTU x OTU Pearson correlation matrix (matrix)\n")


# -------- reload the saved .rds files into the environment --------

# even though result already holds these objects in memory, reloading
# them from disk confirms the files were saved correctly and gives the
# user ready-to-view objects under clear, dedicated names. Since we just
# wrote these files ourselves, we read them back directly by path,
# instead of searching the whole project for them.
samp_filt_view <- readRDS(here("outputs", "samp_filt.rds"))
y_clr_view <- readRDS(here("outputs", "y_clr.rds"))
cor_matrix_view <- readRDS(here("outputs", "cor_matrix.rds"))


# -------- let the user know the objects are ready to inspect --------

cat("\n")
cat("The saved files were reloaded into the environment as:\n")
cat("  - samp_filt_view\n")
cat("  - y_clr_view\n")
cat("  - cor_matrix_view\n")
cat("Use View(samp_filt_view), View(y_clr_view), or View(cor_matrix_view) ")
cat("to inspect them.\n")
