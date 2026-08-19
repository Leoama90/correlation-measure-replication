# demo_generate_matrix_factors.R
#
# Purpose:
#   Demonstrates the usage of generate_matrix_factors() (defined in
#   generate_matrix_factors.R), showing how the number of latent groups
#   controls the sparsity of the resulting correlation matrix, and
#   verifying its key properties (symmetry, unit diagonal, positive
#   semi-definiteness, exact zero pairs between different groups).
#
# Input:
#   - no external files required; all parameters are defined inline
#
# Output:
#   - printed explanation of each step, the resulting matrix, group
#     assignment, and a verification of its properties
#
# Used scripts:
#   - generate_matrix_factors.R

# here: builds file paths relative to the project root
# https://cran.r-project.org/web/packages/here/index.html
library(here)

# -------- load the script to demonstrate --------
source(
  list.files(
    path = here(),
    pattern = "^generate_matrix_factors\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)


# -------- generate a correlation matrix with 3 latent groups --------

cat("#-------------------------------------------------------#\n")
cat("This demo generates an 8 x 8 correlation matrix, splitting the 8\n")
cat("taxa into 3 latent groups. Taxa in the same group are correlated;\n")
cat("taxa in different groups have exactly zero correlation.\n")
cat("#-------------------------------------------------------#\n\n")

result <- generate_matrix_factors(n = 8, n_groups = 3, seed = 42)

cat("Group assignment (one group id per taxon):\n")
print(result$groups)

cat("\nResulting correlation matrix:\n")
print(round(result$mat, 2))

cat("\nExact number of off-diagonal zero entries (both symmetric halves):",
    result$n_zeroes, "\n")


# -------- verify the matrix is a valid correlation matrix --------

cat("\n#-------------------------------------------------------#\n")
cat("Verifying that the matrix is a valid correlation matrix:\n")
cat("symmetric, unit diagonal, and positive semi-definite (PSD).\n")
cat("#-------------------------------------------------------#\n\n")

is_symmetric <- isTRUE(all.equal(result$mat, t(result$mat)))
cat("Symmetric:", is_symmetric, "\n")

has_unit_diagonal <- isTRUE(all.equal(diag(result$mat), rep(1, 8)))
cat("Unit diagonal:", has_unit_diagonal, "\n")

eigenvalues <- eigen(result$mat, only.values = TRUE)$values
is_psd <- all(eigenvalues > -1e-8)
cat("Positive semi-definite (all eigenvalues >= 0):", is_psd, "\n")
cat("Eigenvalues:\n")
print(round(eigenvalues, 3))


# -------- show that taxa in different groups are exactly uncorrelated --------

cat("\n#-------------------------------------------------------#\n")
cat("Checking a pair of taxa known to be in different groups: their\n")
cat("correlation should be exactly zero, not just close to zero.\n")
cat("#-------------------------------------------------------#\n\n")

different_group_pair <- which(result$groups != result$groups[1])[1]
cat("Taxon 1 is in group", result$groups[1], "\n")
cat("Taxon", different_group_pair, "is in group", result$groups[different_group_pair], "\n")
cat("Their correlation is:", result$mat[1, different_group_pair], "\n")


# -------- compare a denser vs sparser matrix, by changing n_groups --------

cat("\n#-------------------------------------------------------#\n")
cat("Comparing sparsity as n_groups changes, keeping n = 8 fixed:\n")
cat("fewer groups -> fewer zeroes (more taxa share a group);\n")
cat("more groups -> more zeroes (taxa are split into smaller groups).\n")
cat("#-------------------------------------------------------#\n\n")

for (n_groups in c(1, 2, 4, 8)) {
  comparison <- generate_matrix_factors(n = 8, n_groups = n_groups, seed = 42)
  cat("n_groups =", n_groups, "-> n_zeroes =", comparison$n_zeroes, "\n")
}