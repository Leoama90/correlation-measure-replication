# check for .csv data in the folders, read the data

# searches for .csv files in the working directory
list_csv <- list.files(pattern = "\\.csv$", recursive = TRUE) 

# cleans the file names and put them in a vector
clean_list_csv <- basename(list_csv)                          

# prints the cleaned names
# print(clean_list_csv)                                         

# extracts the first three letters from the name of each file
prefixes <- substr(sub("\\.csv$", "", clean_list_csv), 1, 4)

# creates a list of vectors and makes them single variables
csv_groups <- split(clean_list_csv, prefixes)
list2env(csv_groups, envir = .GlobalEnv)