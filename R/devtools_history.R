# devtools


# place this file in .Rbuildignore
usethis::use_build_ignore("devtools_history.R")

# test files
usethis::use_test("testFile")
usethis::use_test("kernelFunction")

# packages
usethis::use_package("expm")
usethis::use_package("CompQuadForm")
usethis::use_package("DEoptim")

# data (read data before)
#usethis::use_data(csm, overwrite = TRUE)
#usethis::use_data(energy, overwrite = TRUE)

