# pielou_confrontation.R
#
# Purpose:
#   the aim of this script is to use the data_sim_ph_driven.R script
#   to generate different artificial metagnomic datasets (each one with different
#   ph values) to test how the Pielou index of the simulated biological community
#   varies.
#
# Inputs:
#   - data_sim_ph_driven.R script
#   - pielou_index.R script
# Outputs:
#   - a histogram plot showing how the Pielou index varies 
#

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# tidyverse: data science toolkit (dplyr, ggplot2, tidyr, etc.)
# [https://tidyverse.org](https://tidyverse.org)
library(tidyverse)


# -------- load the functions necessary to make the script happen --------

source(
  list.files(
    path = here(),
    pattern = "^data_sim_ph_driven\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)

source(
  list.files(
    path = here(),
    pattern = "^pielou_ind\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- generate datasets across the pH gradient with a loop --------

ph_values <- seq(5.6, 7.4, by = 0.1)

# empty list to collect one row per pH value
pielou_list <- list()

for (i in seq_along(ph_values)) {
  ph <- ph_values[i]

  # generate the simulated dataset for this pH
  sim <- data_sim_ph_driven(n = 40, N = 50, ph = ph, n_groups = 10, seed = 42)

  # compute its Pielou index and store it alongside the pH value
  pielou_list[[i]] <- data.frame(
    ph = ph,
    pielou = pielou_ind(sim$sim_counts)
  )
}

# combine all rows into a single dataframe

pielou_data <- bind_rows(pielou_list)

# plot Pielou index as a function of pH, as a bar chart

p <-ggplot(pielou_data, aes(x = ph, y = pielou)) +
    geom_col(fill = "steelblue") +
    ggtitle("Pielou Index variation with pH") +
    xlab("pH") +
    ylab("Pielou Index") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 15)
    )
# print on screen the histogram, call needed because ggplots diagrams are not 
# print automatically when sourced
print(p)

# prints on screen the calculated Pielou values along with their pH
cat("\n#----------------------------------------------------------------------#\n")
cat("\npH and Pielou index calculated are:\n")
print(pielou_data)
cat("\n#----------------------------------------------------------------------#\n")


# -------- generate datasets across the pH (with ph_min = 4.0 and ph_max = 9.0) gradient with a loop --------

ph_values <- seq(5.6, 7.4, by = 0.1)

# empty list to collect one row per pH value
pielou_list_01 <- list()

for (i in seq_along(ph_values)) {
  ph <- ph_values[i]
  
  # generate the simulated dataset for this pH
  sim <- data_sim_ph_driven(n = 40, N = 50, ph = ph, n_groups = 10, mu = 5, size = 40, seed = 42)
  
  # compute its Pielou index and store it alongside the pH value
  pielou_list_01[[i]] <- data.frame(
    ph = ph,
    pielou_01 = pielou_ind(sim$sim_counts)
  )
}

# combine all rows into a single dataframe

pielou_data_01 <- bind_rows(pielou_list)

# plot Pielou index as a function of pH, as a bar chart

p_01 <-ggplot(pielou_data, aes(x = ph, y = pielou)) +
  geom_col(fill = "steelblue") +
  ggtitle("Pielou Index variation with pH (with ph range from 4.0 to 9.0)") +
  xlab("pH") +
  ylab("Pielou Index") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15)
  )
# print on screen the histogram, call needed because ggplots diagrams are not 
# print automatically when sourced
print(p_01)

# prints on screen the calculated Pielou values along with their pH
cat("\n#----------------------------------------------------------------------#\n")
cat("\npH and Pielou index calculated are:\n")
print(pielou_data_01)
cat("\n#----------------------------------------------------------------------#\n")

# UI reminder where to search the generated plot
# cat("\nthe plots have been saved in the 'Plots' folder with the names:\n
#      - ph_pielou_relation_00.png\n
#      - ph_pielou_relation_01.png;")

