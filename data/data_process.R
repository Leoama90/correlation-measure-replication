# tidyverse: data science toolkit (dplyr, ggplot2, tidyr, etc.)
# [https://tidyverse.org](https://tidyverse.org)
library(tidyverse)

# -- check for .csv data in the folders, read the data --------

# searches for .csv files in the working directory
list_csv <- list.files(pattern = "\\.csv$", recursive = TRUE) 

# clean file names and put them in a vector
clean_list_csv <- basename(list_csv)                          
                                  
# extract first three letters from the name of each file
prefixes <- substr(sub("\\.csv$", "", clean_list_csv), 1, 4)

# create list and turn the elements of the list into single variables
csv_groups <- split(clean_list_csv, prefixes)
list2env(csv_groups, envir = .GlobalEnv)

# read tables
otu_ <- read.table("data/raw/otu_HMP2_16S.csv", header = TRUE, 
                   sep = ",", row.names = 1)
meta <- read.table("data/raw/meta_HMP2.csv", header = TRUE,
                   sep = ",", row.names = 1)
taxo <- read.table("data/raw/taxonomy_HMP2_16S.csv", header = TRUE,
                   sep = ",", row.names = 1)

# -- checks --------

# checks that samples in otu_ and meta tables are the same
stopifnot(all(rownames(otu_) == rownames(meta)))
# checks that OTUs in otu_ table match the taxonomic classifications in taxo
stopifnot(all(colnames(otu_) == rownames(taxo)))
# checks that all values in otu table are numeric
stopifnot(all(apply(otu_, c(1, 2), is.numeric)))
# checks that all counts in otu table are non-negative
stopifnot(all(otu_ >= 0))

# -- generate outpusts --------

# write Data in .rds file format
saveRDS(as.matrix(otu_), "data/otu_HMP2.rds")
saveRDS(meta, "data/meta_HMP2.rds")
saveRDS(as.matrix(taxo), "data/taxonomy.rds")
