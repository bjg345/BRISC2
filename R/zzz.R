.onAttach <- function(libname, pkgname) {
  #correct_version <- versionInfo()$version >= "1.1.28"
  print_message <- paste0("The ordering of inputs x (covariates) and y (response) in BRISC2_estimation has been changed BRISC2 1.0.0 onwards.
  Please check the new documentation with ?BRISC2_estimation.")
  packageStartupMessage(print_message)
}
