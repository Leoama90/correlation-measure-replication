# Check for .csv data in the folders, read the data


list_csv <- list.files(pattern = "\\.csv$", recursive = TRUE) # searches for .csv files in the working directory
clean_list_csv <- basename(list_csv)                          # cleans the file names and put them in a vector
print(clean_list_csv)                                         # prints the cleaned names
