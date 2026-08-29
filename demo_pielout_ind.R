# demo_pielou_ind.R
#
# Purpose:
#   The aim of this script is to show how the script pielou_ind.R works,
#   printing on screen information regarding the function's parameters
#   and output, checking:
#       - what the Pielou index measures and how to interpret it
#       - how to call pielou_ind() to get a single dataset-wide value
#       - how to call pielou_ind() with per_sample = TRUE to get one
#         value per sample instead
#       - how the function behaves on real (sparse) simulated data,
#         including samples that contain zeroes
#
# Inputs:
#   - pielou_ind.R (sourced below)
#   - data_sim_ph_driven.R (sourced below, used only to generate a
#     small demo dataset)
#
# Outputs:
#   - a series of information and computed diversity values printed on
#     screen

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)


# -------- load the function to demonstrate --------

source(
  list.files(
    path = here(),
    pattern = "^pielou_ind\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- load a function that can generate simulated metagenomic data --------

source(
  list.files(
    path = here(),
    pattern = "^data_sim_ph_driven\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- data generation for the demo --------

# data generation with data_sim_ph_driven.R script
# since this is a demo script, the number of columns and rows was kept relatively low
new_data <- data_sim_ph_driven(n = 10, N = 25, ph = 6.0, n_groups = 3, seed = 42)

# extract the data to analyze from the final output of the previous function
# (for more info check the function data_sim_ph_driven.R or the README.md)
piel_demo_data <- new_data$sim_counts


# -------- descriptive part --------

cat("\n#------------------------------------------------------------------------------#\n")
cat("\nPielou index is a parameter that describes the evenness of a biological community\n")
cat("A value close to 1 means all taxa are equally abundant within a sample;\n")
cat("a value close to 0 means abundance is concentrated in very few taxa.\n")
cat("#------------------------------------------------------------------------------#\n")


# -------- parameters explanation --------

cat("\npielou_ind() takes two arguments:\n
      - x: a matrix/data.frame of non-negative values, samples on rows,\n
        taxa on columns (zeroes are allowed, see below);\n
      - per_sample: if FALSE (default), returns a single value averaged\n
        across all samples; if TRUE, returns one value per sample instead.\n")


# -------- dataset-wide Pielou index --------

cat("\n#------------------------------------------------------------------------------#\n")
cat("\nComputing the Pielou index averaged across all samples ",
    "(per_sample = FALSE, the default):\n\n")

mean_pielou <- pielou_ind(piel_demo_data)
cat("Mean Pielou index for this dataset:", round(mean_pielou, 3), "\n")


# -------- per-sample Pielou index --------

cat("\n#------------------------------------------------------------------------------#\n")
cat("\nComputing the Pielou index for each sample individually ",
    "(per_sample = TRUE):\n\n")

per_sample_pielou <- pielou_ind(piel_demo_data, per_sample = TRUE)
cat("First 5 per-sample Pielou values:\n")
print(round(head(per_sample_pielou, 5), 3))
cat("\nRange across all", length(per_sample_pielou), "samples: [",
    round(min(per_sample_pielou), 3), ",", round(max(per_sample_pielou), 3), "]\n")


# -------- zero handling --------

cat("\n#------------------------------------------------------------------------------#\n")
cat("\nThis dataset is sparse (it comes from data_sim_ph_driven.R, which\n")
cat("introduces zero-inflation): about",
    round(100 * mean(piel_demo_data == 0), 1),
    "% of its entries are zero.\n")
cat("pielou_ind() handles this directly: zero entries contribute 0 to\n")
cat("the entropy sum (the standard convention 0*log(0) = 0), so no\n")
cat("pseudocount or zero-replacement step is needed before calling it.\n")
cat("#------------------------------------------------------------------------------#\n")