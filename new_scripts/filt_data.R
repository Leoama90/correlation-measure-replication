# filt_data.R
#
# This script defines an interactive function that filters a
# metagenomic OTU count table (samples x OTUs) by prevalence, based
# on a threshold chosen by the user at runtime, and shows a
# before/after summary using datasum() (function defined in the "datasummary.R"
# script).
#
# Input:
#   - x: a data.frame, tibble, or matrix with samples on rows and
#     OTUs on columns
#
# Output:
#   - printed before/after summary (via datasum())
#   - a list with the filtered OTU table (samp_filt) and, if a taxa
#     table was found, the matching filtered taxonomy (taxa_filt)
#
# Used functions:
#   - datasum from the script "datasummary.R"
# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)


# -------- recall datasum function from the datasummary.R script --------

# bring datasum() into scope by sourcing the file where it is defined
# (explicit relative path instead of a recursive list.files() search,
# to avoid ambiguity if more than one datasummary.R exists in the project)
source(here("new_scripts", "datasummary.R"))


#' Interactively filter a metagenomic OTU table by prevalence
#'
#' @param x a data.frame, tibble, or matrix with samples on rows and
#'   OTUs on columns.
#'
#' @return an invisible list with samp_filt (the filtered OTU table)
#'   and, if a taxa table is found in the environment, taxa_filt
#'   (the taxonomy table restricted to the surviving OTUs).
#' @export
#


# -------- body of the function --------

filt_data <- function(x) {
  cat("#-----------------------------------------------------------------------# \n")
  # show the "before" summary; datasum() prints its own stats and
  # invisibly returns them, so no extra cat() is needed around it
  cat("\n", "###-----###", "Summary of data BEFORE filtering the zeroes", "###-----###", "\n")
  datasum(x)
  # keep asking until the user provides a valid number between 0 and 1
  repeat {
    # explain and ask in a single prompt, so the user sees the full
    # context (what the threshold means, plus a concrete example)
    # together with the question itself
    question <- readline(
      paste0(
        "Prevalence threshold: an OTU is kept only if it is present ",
        "(non-zero) in at least this fraction of samples.\n",
        "Enter a number between 0 and 1 (e.g. 0.1 = keep OTUs present in ",
        "at least 10% of samples): "
      )
    )
    question_num <- as.numeric(question)
    # if the input is valid, exit the loop
    if (!is.na(question_num) && question_num >= 0 && question_num <= 1) {
      break
    }
    # otherwise, scold the user and loop back to ask again
    cat("I appreciate the enthusiasm, but we only need a number between 0 and 1!\n")
  }
  # filter the OTU table by the prevalence threshold chosen by the user
  # drop = FALSE keeps the result as a data.frame/matrix even if only
  # one OTU survives the filter
  samp_filt <- x[, colSums(x > 0) / nrow(x) >= question_num, drop = FALSE]
  # median of non-zero values >= 5 reads
  # (the anonymous function's argument is named "col", not "x", to avoid
  # shadowing the outer function's x parameter)
  samp_filt <- samp_filt[, apply(samp_filt, 2, function(col) median(col[col > 0]) >= 5), drop = FALSE]
  cat("#-----------------------------------------------------------------------# \n")
  # show the "after" summary
  cat("\n", "-----", "Summary of data AFTER filtering the zeroes", "-----", "\n")
  datasum(samp_filt)
  cat("Used prevalence threshold was", question_num, "\n")
  # if a taxa table exists in the calling environment, align it to
  # the OTUs that survived filtering; otherwise skip this step
  result <- list(samp_filt = samp_filt)
  if (exists("taxa")) {
    taxa_filt <- taxa[colnames(samp_filt), , drop = FALSE]
    result$taxa_filt <- taxa_filt
  }
  # return result as an invisible value
  invisible(result)
}
