#kontaminacija nakon asimetrizacije
#source("generisanje_podataka.R") mora da se ucita pre ucitavanja ovog fajla

# Dodavanje u slucaju:
#  - "normal_spherical"
#  - "cauchy_spherical"

generate_sample_contamination <- function(setting, j, n) {
  if (!setting %in% c("normal_spherical", "cauchy_spherical")) {
    stop("Setting must be 'normal_spherical' or 'cauchy_spherical'.")
  }
  
  X <- generate_sample(setting, j, n)  # from generisanje_podataka.R (Azzalini-type)
  
  if (n < 2) stop("Need n >= 2 for contamination.")
  X[n-1, ] <- c(10, 10)
  X[n,   ] <- c(11,  1)
  
  X
}
