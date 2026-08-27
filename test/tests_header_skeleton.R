# test-<test_name>.R
#
# Purpose:
#   <clever test here> ...,
#   checking:
#       - 
#       - 
#       - 
#       - 
#
# Inputs:
#   - 
#
# Outputs:
#   - test results printed to the console (pass/fail for each check)
#
# Note:
#
#
#
#
#
#
#
#
#


# -------- begin of libraries zone --------
#
# -------- end of libraries zone --------


# -------- load the function to test --------

source(
  list.files(
    path = here(),
    pattern = "^<script_name>\\.R$",
    full.names = TRUE,
    recursive = TRUE
  )
)

# EXPLANATION: This file contains the core header information shared across
# the project scripts. It was originally kept outside the project folder,
# as it was intended for reference only and is not meant to be executed.
# It was subsequently included in the project folder to support the revision
# and standardization of legacy scripts that lacked a proper header.