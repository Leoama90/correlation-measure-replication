# ph_pielou_confrontation.R
#
# Purpose:
#   Uses data_sim_ph_driven.R to generate five families of simulated
#   metagenomic datasets, each covering the same pH gradient
#   (3.6 to 10.4, step 0.4), but with a different pH tolerance range
#   (sigma_min/sigma_max) for the simulated taxa: 0.1-0.2, 0.3-0.4,
#   0.5-0.6, 0.7-0.8, and 0.9-1.0. 
#   All other parameters (n, N, n_groups, phi_max, ph_min, ph_max, seed)
#   are kept fixed across the five scenarios, isolating pH tolerance
#   as the only variable of interest.
#   For each scenario, computes the Pielou index at every pH value, to
#   study how community evenness responds to an environmental gradient
#   as a function of how ecologically specialized (narrow tolerance) or
#   generalist (wide tolerance) the simulated taxa are.
#
# Inputs:
#   - data_sim_ph_driven.R script
#   - pielou_ind.R script
#
# Outputs:
#   - five bar charts printed to screen (one per tolerance scenario),
#     each showing Pielou index vs pH
#   - Plots/ph_pielou_confrontation.png: the five charts combined into
#     a single 3x2 figure, with a caption panel in the last cell

# gridExtra: arrange multiple plots (ggplot or grob) in a grid
# https://cran.r-project.org/web/packages/gridExtra/index.html
library(gridExtra)
# grid: low-level graphics, needed here for textGrob()
# https://cran.r-project.org/web/packages/grid/index.html
library(grid)

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


# -------- define the ph range for this script --------

ph_values <- seq(3.6, 10.4, by = 0.4)


# -------- 00 generate datasets across the pH gradient with a loop --------


# empty list to collect one row per pH value
pielou_list <- list()

for (i in seq_along(ph_values)) {
  ph_sim <- ph_values[i]

  # generate the simulated dataset for this pH
  sim <- data_sim_ph_driven(
    n = 100, N = 120,
    ph = ph_sim, ph_min = 3.5, ph_max = 10.5,
    n_groups = 10, sigma_min = 0.1, sigma_max = 0.2,
    phi_max = 0.9, seed = 42
  )

  # compute its Pielou index and store it alongside the pH value
  pielou_list[[i]] <- data.frame(
    ph = ph_sim,
    pielou = round(pielou_ind(sim$sim_counts), 2)
  )
}

# combine all rows into a single dataframe
pielou_data <- bind_rows(pielou_list)

# plot Pielou index as a function of pH, as a bar chart
p <- ggplot(pielou_data, aes(x = ph, y = pielou)) +
  geom_col(fill = "steelblue") +
  ggtitle("Pielou Index variation with pH (tolerance 0.1 - 0.2)") +
  xlab("pH") +
  ylab("Pielou Index") +
  coord_cartesian(xlim = c(3.0, 11.0), ylim = c(0, 1)) +
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
# cat("\nColumn number of pielou_data_01 is:", ncol(pielou_data), "\n")
# cat("\nRow number of pielou_data_01 is:", nrow(pielou_data), "\n")
cat("\n#----------------------------------------------------------------------#\n")


# -------- 01 generate a dataset which differs from the previous one --------

# empty list to collect one row per pH value
pielou_list_01 <- list()

for (i in seq_along(ph_values)) {
  ph_sim <- ph_values[i]

  # generate the simulated dataset for this pH, with different parameters from the previous simulation
  sim_01 <- data_sim_ph_driven(
    n = 100, N = 120,
    ph = ph_sim, ph_min = 3.5, ph_max = 10.5,
    sigma_min = 0.3, sigma_max = 0.4, phi_max = 0.9,
    n_groups = 10, seed = 42
  )

  # compute its Pielou index and store it alongside the pH value
  pielou_list_01[[i]] <- data.frame(
    ph = ph_sim,
    pielou_01 = round(pielou_ind(sim_01$sim_counts), 2)
  )
}

# combine all rows into a single dataframe
pielou_data_01 <- bind_rows(pielou_list_01)

# plot Pielou index as a function of pH, as a bar chart
p_01 <- ggplot(pielou_data_01, aes(x = ph, y = pielou_01)) +
  geom_col(fill = "steelblue") +
  ggtitle("Pielou Index variation with pH (tolerance 0.3 - 0.4)") +
  xlab("pH") +
  ylab("Pielou Index") +
  coord_cartesian(xlim = c(3.0, 11.0), ylim = c(0, 1)) +
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
# cat("\nColumn number of pielou_data_01 is:", ncol(pielou_data_01), "\n")
# cat("\nRow number of pielou_data_01 is:", nrow(pielou_data_01), "\n")
cat("\n#----------------------------------------------------------------------#\n")


# -------- 02 generate a dataset which differs from the previous one --------

# empty list to collect one row per pH value
pielou_list_02 <- list()

for (i in seq_along(ph_values)) {
  ph_sim <- ph_values[i]
  
  # generate the simulated dataset for this pH, with different parameters from the previous simulation
  sim_02 <- data_sim_ph_driven(
    n = 100, N = 120,
    ph = ph_sim, ph_min = 3.5, ph_max = 10.5,
    sigma_min = 0.5, sigma_max = 0.6, phi_max = 0.9,
    n_groups = 10, seed = 42
  )
  
  # compute its Pielou index and store it alongside the pH value
  pielou_list_02[[i]] <- data.frame(
    ph = ph_sim,
    pielou_02 = round(pielou_ind(sim_02$sim_counts), 2)
  )
}

# combine all rows into a single dataframe
pielou_data_02 <- bind_rows(pielou_list_02)

# plot Pielou index as a function of pH, as a bar chart
p_02 <- ggplot(pielou_data_02, aes(x = ph, y = pielou_02)) +
  geom_col(fill = "steelblue") +
  ggtitle("Pielou Index variation with pH (tolerance 0.5 - 0.6)") +
  xlab("pH") +
  ylab("Pielou Index") +
  coord_cartesian(xlim = c(3.0, 11.0), ylim = c(0, 1)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15)
  )

# print on screen the histogram, call needed because ggplots diagrams are not
# print automatically when sourced
print(p_02)

# prints on screen the calculated Pielou values along with their pH
cat("\n#----------------------------------------------------------------------#\n")
cat("\npH and Pielou index calculated are:\n")
print(pielou_data_02)
# cat("\nColumn number of pielou_data_02 is:", ncol(pielou_data_02), "\n")
# cat("\nRow number of pielou_data_02 is:", nrow(pielou_data_02), "\n")
cat("\n#----------------------------------------------------------------------#\n")


# -------- 03 - generate a dataset which differs from the previous one --------

# empty list to collect one row per pH value
pielou_list_03 <- list()

for (i in seq_along(ph_values)) {
  ph_sim <- ph_values[i]
  
  # generate the simulated dataset for this pH, with different parameters from the previous simulation
  sim_03 <- data_sim_ph_driven(
    n = 100, N = 120,
    ph = ph_sim, ph_min = 3.5, ph_max = 10.5,
    sigma_min = 0.7, sigma_max = 0.8, phi_max = 0.9,
    n_groups = 10, seed = 42
  )
  
  # compute its Pielou index and store it alongside the pH value
  pielou_list_03[[i]] <- data.frame(
    ph = ph_sim,
    pielou_03 = round(pielou_ind(sim_03$sim_counts), 2)
  )
}

# combine all rows into a single dataframe
pielou_data_03 <- bind_rows(pielou_list_03)

# plot Pielou index as a function of pH, as a bar chart
p_03 <- ggplot(pielou_data_03, aes(x = ph, y = pielou_03)) +
  geom_col(fill = "steelblue") +
  ggtitle("Pielou Index variation with pH (tolerance 0.7 - 0.8)") +
  xlab("pH") +
  ylab("Pielou Index") +
  coord_cartesian(xlim = c(3.0, 11.0), ylim = c(0, 1)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15)
  )

# print on screen the histogram, call needed because ggplots diagrams are not
# print automatically when sourced
print(p_03)

# prints on screen the calculated Pielou values along with their pH
cat("\n#----------------------------------------------------------------------#\n")
cat("\npH and Pielou index calculated are:\n")
print(pielou_data_03)
# cat("\nColumn number of pielou_data_03 is:", ncol(pielou_data_03), "\n")
# cat("\nRow number of pielou_data_03 is:", nrow(pielou_data_03), "\n")
cat("\n#----------------------------------------------------------------------#\n")


# -------- 04 - generate a dataset which differs from the previous one --------

# empty list to collect one row per pH value
pielou_list_04 <- list()

for (i in seq_along(ph_values)) {
  ph_sim <- ph_values[i]
  
  # generate the simulated dataset for this pH, with different parameters from the previous simulation
  sim_04 <- data_sim_ph_driven(
    n = 100, N = 120,
    ph = ph_sim, ph_min = 3.5, ph_max = 10.5,
    sigma_min = 0.9, sigma_max = 1.0, phi_max = 0.9,
    n_groups = 10, seed = 42
  )
  
  # compute its Pielou index and store it alongside the pH value
  pielou_list_04[[i]] <- data.frame(
    ph = ph_sim,
    pielou_04 = round(pielou_ind(sim_04$sim_counts), 2)
  )
}

# combine all rows into a single dataframe
pielou_data_04 <- bind_rows(pielou_list_04)

# plot Pielou index as a function of pH, as a bar chart
p_04 <- ggplot(pielou_data_04, aes(x = ph, y = pielou_04)) +
  geom_col(fill = "steelblue") +
  ggtitle("Pielou Index variation with pH (tolerance 0.9 - 1.0)") +
  xlab("pH") +
  ylab("Pielou Index") +
  coord_cartesian(xlim = c(3.0, 11.0), ylim = c(0, 1)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15)
  )

# print on screen the histogram, call needed because ggplots diagrams are not
# print automatically when sourced
print(p_04)

# prints on screen the calculated Pielou values along with their pH
cat("\n#----------------------------------------------------------------------#\n")
cat("\npH and Pielou index calculated are:\n")
print(pielou_data_04)
# cat("\nColumn number of pielou_data_04 is:", ncol(pielou_data_04), "\n")
# cat("\nRow number of pielou_data_04 is:", nrow(pielou_data_04), "\n")
cat("\n#----------------------------------------------------------------------#\n")

# -------- generate the final plot --------

# -------- combine all five Pielou plots into a single 3x2 figure --------


# build a simple text panel for the last (6th) cell of the grid
caption_panel <- textGrob(
  paste(
    "This shows how the Pielou index varies",
    "as a function of taxa's pH tolerance:",
    "",
    "- low tolerance means taxa live in a",
    "  narrow ecological niche, so evenness",
    "  tends to be low",
    "",
    "- high tolerance means taxa can live",
    "  across a wider pH range, so evenness",
    "  tends to be high",
    sep = "\n"
  ),
  gp = gpar(fontsize = 15),
  just = "left"
)

# generate the final plot in the "Plots" folder
png(
  filename = here("Plots", "ph_pielou_confrontation.png"),
  width = 4500, height = 3000, res = 300
)
# this is needed to arrange the 5 plots + the 6th panel with the description, in a grid 3x2
grid.arrange(
  p, p_01, p_02, p_03, p_04, caption_panel,
  ncol = 2, nrow = 3
)
# close the device
dev.off()

cat("#----------------------------------------------------------------------#")
cat("\nThe final plot has been saved in the <Plots> folder\n")
cat("#----------------------------------------------------------------------#")