# demo_data_sim_ph_driven.R
#
# Purpose:
#   The aim of this script is to show how the script data_sim_ph_driven.R
#   works, printing on screen information regarding the script's
#   parameters and output, then generating a demo dataset and showing
#   its actual results (instead of asking the user to inspect it
#   manually), to help the user understand how it works.
#
# Inputs:
#   - data_sim_ph_driven.R (sourced below)
#
# Outputs:
#   - a series of information printed on screen, and a demo dataset
#     (demo_data) left in the environment for further inspection

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)


# -------- load the function to demonstrate --------

source(
  list.files(
    path = here(),
    pattern = "^data_sim_ph_driven\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- initial explanation --------

cat("\n#------------------------------------------------------------------#\n")
cat("\nThis script is a demonstration of how data_sim_ph_driven.R works\n")
cat("\nThe final output of the function is a list of 7 elements:\n
    - mat: the correlation matrix used to create the dataset;\n
    - groups: a vector with one entry per taxon, giving which of the
      n_groups latent correlation groups that taxon belongs to;\n
    - ph_optima: the optimal pH of each taxon;\n
    - phi_per_taxon: determines how likely each taxon is to be absent
      (zero) depending on its match with the environmental pH;\n
    - sim_data: intermediate continuous data generated according to the
      correlation structure;\n
    - sim_counts: the final simulated count data;\n
    - total_zero_rate: the overall proportion of zero values in the
      simulated count matrix.\n")
cat("The data we are looking for is sim_counts\n")
cat("\n#------------------------------------------------------------------#\n")


# -------- explanation of the function's parameters --------

cat("\nHere is a brief explanation of the parameters needed to run the data simulation.\n
      The first parameters we need are:\n
            - n = number of simulated taxa (this will give the number of columns);\n
            - N = number of samples to simulate (this will give the number of rows);\n
            - ph = the pH value shared by all the taxa;\n
            - n_groups = number of groups which contain correlated taxa.\n")
cat("\nThere are other parameters with default values, which are:\n
            - ph_min: minimum pH of the dataset (default 5.5);\n
            - ph_max: maximum pH of the dataset (default 7.5);\n
            - sigma_min: minimum pH tolerance (default 0.3);\n
            - sigma_max: maximum pH tolerance (default 0.8);\n
            - phi_max: upper bound on zero-inflation probability
              (0 to 1): even a taxon maximally mismatched to the environmental pH
              is capped at this probability, so it is never certain to be absent
              from every sample (default 0.9);\n
            - mu = mean parameter of the negative binomial component
              of the ZINB distribution, shared by all simulated taxa
              (default 20);\n
            - size = dispersion parameter of the negative binomial component
              (default 30);\n
            - seed = seed for the rng (default NULL).\n")


# -------- generate the demo dataset --------

cat("\n#------------------------------------------------------------------#\n")
cat("\nFor a simple data generation we will use:
            - n = 20,
            - N = 50,
            - ph = 6.5,
            - n_groups = 5,
            - seed = 42.\n")

# generate the demonstrative data using the values mentioned above
demo_data <- data_sim_ph_driven(n = 20, N = 50, ph = 6.5, n_groups = 5, seed = 42)

cat(
  "\nDataset generated.
Dimensions of sim_counts:\n",
  "- number of rows:", nrow(demo_data$sim_counts), "\n",
  "- number of columns:", ncol(demo_data$sim_counts), "\n",
  "\n"
)

cat("Overall zero rate: ", round(demo_data$total_zero_rate, 3), "\n")


# -------- show which taxa are best/worst matched to this pH --------

cat("\n#------------------------------------------------------------------#\n")
cat("\nAt pH = 6.5, the taxon best matched to the environment (lowest\n")
cat("zero-inflation probability) and the worst matched (highest) are:\n\n")

best_taxon <- which.min(demo_data$phi_per_taxon)
worst_taxon <- which.max(demo_data$phi_per_taxon)

cat(
  "Best match:  ", colnames(demo_data$sim_counts)[best_taxon],
  "- pH optimum:", round(demo_data$ph_optima[best_taxon], 2),
  "- phi:", round(demo_data$phi_per_taxon[best_taxon], 3), "\n"
)
cat(
  "Worst match: ", colnames(demo_data$sim_counts)[worst_taxon],
  "- pH optimum:", round(demo_data$ph_optima[worst_taxon], 2),
  "- phi:", round(demo_data$phi_per_taxon[worst_taxon], 3), "\n"
)


# -------- preview the final simulated count matrix --------

cat("\n#------------------------------------------------------------------#\n")
cat("\nFirst 5 samples (rows) and 5 taxa (columns) of sim_counts:\n\n")
print(demo_data$sim_counts[1:5, 1:5])

cat("\nThe full dataset is available in demo_data$sim_counts.\n")
cat("Use View(demo_data$sim_counts) to inspect it interactively.\n")
