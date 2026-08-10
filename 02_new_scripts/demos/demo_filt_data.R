# demo_filt_data.R
#
# Purpose:
#   Demonstrates how to use filt_data() (defined in filt_data.R) on a
#   small, hand-built OTU count table, showing the non-interactive call
#   used for automated/reproducible use (e.g. sensitivity analyses over
#   a range of threshold values), instead of the interactive default
#   behaviour (which would block waiting for console input).
#
# Input:
#   - no external files required; a small example OTU count table is
#     built inline
#
# Output:
#   - printed before/after summary of the filtering (via datasum(),
#     called internally by filt_data())
#
# Used scripts:
#   - filt_data.R (also sources datasummary.R internally)

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# bring filt_data() (and datasum(), sourced inside it) into scope
source(
  list.files(
    path = here(),
    pattern = "^filt_data\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)

# -------- build a small example OTU count table --------

# UI message describing what the script does
cat(
  "This script builds a small example OTU table (10 samples x 5 OTUs) ",
  "with varying prevalence and abundance, to demonstrate filt_data().\n"
)

# fix the seed for reproducible random counts
set.seed(42)

# simulate 10 samples x 5 OTUs from independent Poisson distributions,
# with different lambda per OTU so that some are dropped by prevalence
# (OTU2, very low lambda) and others by the abundance/median filter
# (OTU5, low lambda), leaving a mix of kept/dropped OTUs to inspect
otu_demo <- matrix(
  rpois(10 * 5, lambda = rep(c(20, 0.5, 8, 15, 1), each = 10)),
  nrow = 10, ncol = 5,
  dimnames = list(NULL, paste0("OTU", 1:5))
)

# show the raw example table before filtering
cat("\nExample OTU table (raw counts):\n")
print(otu_demo)


# -------- apply filt_data() non-interactively --------

# UI message describing the filtering step
cat(
  "\nApplying filt_data() with prevalence_threshold = 0.4 and ",
  "abundance_threshold = 2 (non-interactive, for reproducibility).\n"
)

# call filt_data() with both thresholds passed explicitly, so the
# interactive readline() prompt is skipped entirely; datasum() is
# called internally by filt_data(), both before and after filtering,
# so no separate summary call is needed here
filt_result <- filt_data(otu_demo, prevalence_threshold = 0.4, abundance_threshold = 2)

# show which OTUs survived the filtering
cat("\nOTUs kept after filtering:", paste(colnames(filt_result$samp_filt), collapse = ", "), "\n")


# -------- example with a taxa table present --------

# UI message describing the optional taxa_filt behaviour
cat(
  "\nRepeating the call with a 'taxa' table defined, to show that ",
  "taxa_filt is included and aligned to the surviving OTUs.\n"
)

# define a minimal taxonomy table, indexed by OTU name, to demonstrate
# the optional taxa_filt branch of filt_data()
taxa <- data.frame(
  genus = c("GenusA", "GenusB", "GenusC", "GenusD", "GenusE"),
  row.names = paste0("OTU", 1:5)
)

filt_result_with_taxa <- filt_data(otu_demo, prevalence_threshold = 0.4, abundance_threshold = 2)

# show the taxonomy restricted to the surviving OTUs
cat("\ntaxa_filt (taxonomy of surviving OTUs):\n")
print(filt_result_with_taxa$taxa_filt)

# clean up the demo taxa object, so it doesn't leak into other scripts
# run later in the same R session
rm(taxa)
