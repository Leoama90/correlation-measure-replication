# fil_data.R
#
# This script defines an interactive function that filters a
# metagenomic OTU count table (samples x OTUs) by prevalence, based
# on a threshold chosen by the user at runtime, and shows a
# before/after summary using datasum() (function defined in the "datasummary.R"
# script).
#
# Input:
#   - y: a data.frame, tibble, or matrix with samples on rows and
#     OTUs on columns
#
# Output:
#   - printed before/after summary (via datasum())
#   - a list with the filtered OTU table (samp_filt) and, if a taxa
#     table was found, the matching filtered taxonomy (taxa_filt)
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)
# tidyverse: data science toolkit (dplyr, ggplot2, tidyr, etc.)
# [https://tidyverse.org](https://tidyverse.org)
library(tidyverse)

# bring datasum() into scope by sourcing the file where it is defined
source(
  list.files(
    path = here(),
    pattern = "^datasummary\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)
#' Interactively filter a metagenomic OTU table by prevalence
#'
#' @param y a data.frame, tibble, or matrix with samples on rows and
#'   OTUs on columns.
#'
#' @return an invisible list with samp_filt (the filtered OTU table)
#'   and, if a taxa table is found in the environment, taxa_filt
#'   (the taxonomy table restricted to the surviving OTUs).
#' @export
#
filt_data <- function(y) {
  cat("#-----------------------------------------------------------------------# \n")
  # show the "before" summary; datasum() prints its own stats and
  # invisibly returns them, so no extra cat() is needed around it
  datasum(y)

  # keep asking until the user provides a valid number between 0 and 1
  repeat {
    question <- readline("How much filter would you like to apply? (type a number between 0 and 1) ")
    question_num <- as.numeric(question)

    # if the input is valid, exit the loop
    if (!is.na(question_num) && question_num >= 0 && question_num <= 1) {
      break
    }

    # otherwise, scold the user and loop back to ask again
    cat("I appreciate the enthusiasm, but we only need a number between 0 and 1!\n")
  }

  # filter the OTU table by the prevalence threshold chosen by the user
  samp_filt <- y[, colSums(y > 0) / nrow(y) >= question_num]

  # median of non-zero values >= 5 reads
  samp_filt <- samp_filt[, apply(samp_filt, 2, function(x) median(x[x > 0]) >= 5)]

  cat("#-----------------------------------------------------------------------# \n")

  # show the "after" summary
  datasum(samp_filt)
  cat("Used prevalence threshold was", question_num, "\n")

  # if a taxa table exists in the calling environment, align it to
  # the OTUs that survived filtering; otherwise skip this step
  result <- list(samp_filt = samp_filt)

  if (exists("taxa")) {
    taxa_filt <- taxa[colnames(samp_filt), ]
    result$taxa_filt <- taxa_filt
  }

  invisible(result)
}
