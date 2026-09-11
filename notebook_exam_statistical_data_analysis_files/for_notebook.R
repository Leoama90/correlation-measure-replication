# for_notebook.R
#
# Purpose:
#   here are stored all the data generations and/or the usage of other scripts.
#   Few scripts in the 02_new_scripts/ folder have interactive user parts (where the user is
#   asked to give a number between 0 and 1). Running those scripts in the notebook
#   might lead to notebook freezing.
#   The outputs of the scripts will be used in the notebook code chunks.
#
# Inputs:
#   a series of hardcoded instruction put here only to show how scripts work
#
# Outputs:
#   various outputs of various scripts
#
# Note:
# to make the reading of this script (slightly) more intuitive, the scripts
# that need to be used will be sourced only in their region (e.g.: the filt_data.R
# script will be sourced under its commented title).
#
# Note 2:
# to maintain a logical flow every script is sourced in the same order as the notebook.
# Althought, since some need data to transform, the data_sim_ph_driven.R is sourced first,
# to allow the generation of some dummy data to process.
#
# Note 3:
# Rstudio can turn code lines into comments by pressing ctrl + shift + c.
# Every code chunk was run by commenting the others.

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# brings into scope the data_sim_ph_driven.R (the only one dummy data generator)
source(
  list.files(
    path = here(),
    pattern = "^data_sim_ph_driven\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- filt_data.R chunk --------

# # brings into scope the filt_data.R script
# source(
#   list.files(
#     path = here(),
#     pattern = "^filt_data\\.R$",
#     full.names = TRUE,
#     recursive = TRUE
#   )
# )
#
# # generate data for the filt_data.R explanation code chunk for the notebook
# data_for_notebook <- data_sim_ph_driven(n_taxa = 9,
#                                         N_sample = 8,
#                                         n_groups = 3,
#                                         ph = 5.6,
#                                         sigma_min = 0.01,
#                                         sigma_max = 0.2,
#                                         seed = 42)
#
# # take the simulated data from the list and put it into another variable
# data_for_notebook_count <- data_for_notebook$sim_counts
#
# # show the generated data
# print(data_for_notebook_count)
#
# # separation line
# cat("\n#---------------------------------------------------------------------------#")
#
# # show and impervious message to the user to remember that this is the filt_data chunk
# cat("\nTHIS IS FOR THE FILT_DATA CHUNK! READ MEEEEEEE!\n")
#
# # apply the filt_data() function
# filtered_data_for_notebook_count <- filt_data(data_for_notebook_count)
#
# # show the newly filtered data
# print(filtered_data_for_notebook_count)
#
# # show and impervious message to the user to remember that this is the end of the filt_data chunk
# cat("\nTHIS. IS. THE. END. OF THE FILT_DATA CHUNK!!\n")
#
# # write Data in .rds file format
# saveRDS(as.matrix(data_for_notebook_count), here("notebook_exam_statistical_data_analysis_files", "data_for_notebook_count.rds"))
# saveRDS(as.matrix(filtered_data_for_notebook_count$samp_filt), here("notebook_exam_statistical_data_analysis_files", "filtered_data_for_notebook_count.rds"))
# # separation line
# cat("\n#---------------------------------------------------------------------------#")


# -------- generate_matrix_factor.R code chunk --------

# # brings into scope the filt_data.R script
# source(
#   list.files(
#     path = here(),
#     pattern = "^generate_matrix_factors\\.R$",
#     full.names = TRUE,
#     recursive = TRUE
#   )
# )
#
# # generate a dummy matrix with n = 10, n_groups = 2, seed = 42
# corr_mat_notebook <- generate_matrix_factors(10, 3, seed = 42)
#
# # print the actual matrix
# print(corr_mat_notebook$mat)

# --------- clr_pearson.R code chunk ---------

# brings into scope the clr_pearson.R (the only one dummy data generator)
source(
  list.files(
    path = here(),
    pattern = "^clr_pearson\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)

# brings the dummy pre-generated data into scope
data_for_notebook_count <- readRDS(
  list.files(
    path = here(),
    pattern = "^data_for_notebook_count\\.rds$",
    full.names = TRUE,
    recursive = TRUE
  )
)

# apply the clr_on_data() function
processed_data <- clr_on_data(data_for_notebook_count)

# print the filtered samples
print(processed_data$samp_filt)

# print the log of every value
print(processed_data$y_clr)

# print correlation matrix
print(processed_data$cor_matrix)
