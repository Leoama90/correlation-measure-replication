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

# -------- library & sourced data_sim_ph_driven code area --------

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
# # print title for correlation matrix
# cat("\ndata_for_notebook\n")
# 
# # show the generated data
# print(data_for_notebook_count)
# 
# # print title for correlation matrix
# cat("\nCorrelation matrix of data_for_notebook\n")
# 
# # print correlation matrix of the newly generated data
# print(data_for_notebook$mat)
#
# # separation line
# cat("\n#---------------------------------------------------------------------------#")
#
# # show and impervious message to the user to remember that this is the filt_data chunk
# cat("\nTHIS IS FOR THE FILT_DATA CHUNK! READ MEEEEEEE!\n")
#
# # apply the filt_data() function using 0.33 as prevalence threshold
# filtered_data_for_notebook_count <- filt_data(data_for_notebook_count, prevalence_threshold = 0.33)
# #
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






# -------- pseudocount.R code chunk --------

# brings into scope the filt_data.R script
source(
  list.files(
    path = here(),
    pattern = "^pseudocount\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)



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




# -------- clr_pearson.R code chunk -------- ---------

# # brings into scope the clr_pearson.R (the only one dummy data generator)
# source(
#   list.files(
#     path = here(),
#     pattern = "^clr_pearson\\.R$",
#     full.names = TRUE,
#     recursive = TRUE
#   )
# )
# 
# # brings the pre-generated dummy data into scope
# data_for_notebook_count <- readRDS(
#   list.files(
#     path = here(),
#     pattern = "^data_for_notebook_count\\.rds$",
#     full.names = TRUE,
#     recursive = TRUE
#   )
# )
# 
# # print the pre-generated data
# print(data_for_notebook_count)
# 
# # apply the clr_on_data() function
# processed_data <- clr_on_data(data_for_notebook_count)
# 
# # separation line
# cat("\n#-----------------------------------------------#\n")
# 
# # print the filtered samples
# print(processed_data$samp_filt)
# 
# # separation line
# cat("\n#-----------------------------------------------#\n")
# 
# # points out the transformed matrix
# cat("\nThis is the centered-log transformed matrix:\n")
# 
# # print the log of every value
# print(processed_data$y_clr, digits = 1)
# 
# # point down the rowsums
# cat("\n...................................................\n")
# 
# print(rowSums(processed_data$y_clr), digits = 2)
# 
# # point the row sums
# cat("\n^^^^^^^^^^^^^^^^^^^^^^rowSums^^^^^^^^^^^^^^^^^^^^^^\n")
# 
# # separation line
# cat("\n#-----------------------------------------------#\n")
# 
# # cat to indicate correlation matrix
# cat("\nThis is the correlation matrix:\n")
# 
# # print correlation matrix
# print(processed_data$cor_matrix, digits = 2)
# 
# # runs a datasum check to find eigenvalues
# datasum(processed_data$cor_matrix)
# 
# # --- save clr_on_data outputs in .rds ---
# 
# # save the filtered data as .rds file
# saveRDS(as.matrix(processed_data$samp_filt), here("notebook_exam_statistical_data_analysis_files", "clr_on_data_output", "processed_data_samp_filt.rds"))
# # save the transformed data as .rds file
# saveRDS(as.matrix(processed_data$y_clr), here("notebook_exam_statistical_data_analysis_files", "clr_on_data_output", "processed_data_y_clr.rds"))
# # save the correlation matrix as .rds file
# saveRDS(as.matrix(processed_data$cor_matrix), here("notebook_exam_statistical_data_analysis_files", "clr_on_data_output", "processed_data_cor_matrix.rds"))

# brings the clr_pearson script into scope
source(
  list.files(
    path = here(),
    pattern = "^clr_pearson\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)

# generate data for the example
data_for_clr_pearson <- data_sim_ph_driven(
  n_taxa = 15,
  N_sample = 25,
  n_groups = 7,
  ph = 6.0,
  sigma_min = 0.01,
  sigma_max = 0.2,
  seed = 42
)

# isolate the count from the newly generated data
data_for_clr_pearson_counts <- data_for_clr_pearson$sim_counts

# print the counts
print(data_for_clr_pearson_counts)

# apply clr_on_data() function to the newly generated data
transf_data <- clr_on_data(data_for_clr_pearson_counts,
                           prevalence_threshold = 0.3,
                           threshold_pct = 0.5
)

# print the filtered data
print(round(transf_data$samp_filt, 2))

# print the log-transformed data
print(round(transf_data$y_clr, 2))

# print the sums of the rows
row_sums <- rowSums(transf_data$y_clr)
# set a tolerance for very small values
tol <- 1e-14
# impose value = 0 for values < tolerance
row_sums[abs(row_sums) < tol] <- 0
# row sums
print(row_sums)
# -------- NorTa_simulation.R code chunk --------

# # brings into scope the NorTa_simulation.R script
# source(
#   list.files(
#     path = here(),
#     pattern = "^NorTa_simulation\\.R$",
#     full.names = TRUE,
#     recursive = TRUE
#   )
# )
# 
# # explain the user how the datas are generated with which parameters
# cat("\nThe following data are generated with the norta_simulation()
# function defined in the NorTa_simulation.R script, using the
# following parameters:
#       - n = 9;\n
#       - n_groups = 3;\n
#       - N = 11;\n
#       - seed = 42;\n")
# 
# # additional space 
# cat("\n")
# 
# # generate data with n = 9, n_groups = 3, N = 11, seed = 42
# norta_data_01 <- norta_simulation(n = 9, n_groups = 3, N = 11, seed = 42) 
# 
# # title of the data
# cat("\nData generated with n_groups = 3 and n = 9\n")
# # print generated data sim_counts
# print(norta_data_01$sim_counts)
# 
# # title of the first correlation matrix
# cat("\nCorrelation matrix with n_groups = 3 and n = 9\n")
# # print generated data sim_counts
# print(norta_data_01$R_true)
# 
# # explain the user how the datas are generated with which parameters
# cat("\nThe following data are generated with the norta_simulation()
# function defined in the NorTa_simulation.R script, using the
# following parameters:
#       - n = 10;\n
#       - n_groups = 10;\n
#       - N = 10;\n
#       - seed = 42;\n
# With a number of groups (n_groups) equal to the number of taxa 
# (number of taxa is the parameter n), it will resolve in having
# zero taxa correlated to each other, leading to a unitary correlation matrix.\n")
# 
# # additional space 
# cat("\n")
# 
# # generate data with n = 10, n_groups = 10, N = 10, seed = 42 (gives an identity matrix)
# norta_data_02 <- norta_simulation(n = 10, n_groups = 10, N = 10, seed = 42) 
# 
# # title of the data
# cat("\nData generated with n_groups = n\n")
# # print generated data sim_counts
# print(norta_data_02$sim_counts)
# 
# # title of the second correlation matrix
# cat("\nCorrelation matrix with n_groups = n\n")
# # print generated data sim_counts
# print(norta_data_02$R_true)


# # write Data in .rds file format
# saveRDS(as.matrix(norta_data$sim_counts), here("notebook_exam_statistical_data_analysis_files", "norta_data_counts.rds"))

# -------- data_sim_ph_driven.R code chunk --------




