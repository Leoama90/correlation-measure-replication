library(testthat)

#expect_contains(data/raw/meta_HMP2.csv, list_csv)
#expect_contains(data/raw/otu_HMP2_16S.csv, list_csv)
#expect_contains(data/raw/taxonomy_HMP2_16S.csv, list_csv)

expect_named(clean_list_csv,
             meta_HMP2.csv,
             ignore.order = FALSE,
             ignore.case = FALSE,
             info = NULL,
             label = NULL
             )

