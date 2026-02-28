library(MASS)     # mvrnorm
library(mvtnorm)  # rmvt

set.seed(123)

n      <- 100
M      <- 3000
delta  <- c(0.15, 0.15)
j_vals <- 0:3

Sigma_sph  <- matrix(c(1,0,0,1),2,2)
Sigma_ell  <- matrix(c(2,1,1,3),2,2)

#### Normalno jezgro ####

gen_normal_kernel <- function(n, Sigma) {
  mvrnorm(n, mu = c(0,0), Sigma = Sigma)
}

#### Kosijevo jezgro ####

gen_cauchy_kernel <- function(n, Sigma) {
  rmvt(n, sigma = Sigma, df = 1)
}

#### Uslovljavanje na kone ####

in_cone <- function(z) {
  abs(atan2(z[,2], z[,1])) <= 0.5
}

gen_cone_kernel <- function(n, generator, Sigma) {
  Z <- matrix(NA, n, 2)
  count <- 0
  
  while(count < n) {
    cand <- generator(n, Sigma)
    keep <- cand[in_cone(cand), , drop=FALSE]
    
    if(nrow(keep) > 0) {
      take <- min(n - count, nrow(keep))
      Z[(count+1):(count+take), ] <- keep[1:take, ]
      count <- count + take
    }
  }
  Z
}

#### Skewing normal ####

skew_normal <- function(Z, j, delta) {
  n <- nrow(Z)
  U <- runif(n)
  
  lin <- j * as.vector(Z %*% delta)
  p   <- pnorm(lin)
  
  X <- Z
  X[U > p, ] <- -Z[U > p, ]
  
  X
}

#### Skewing kosi ####

skew_cauchy <- function(Z, j, delta) {
  n <- nrow(Z)
  U <- runif(n)
  
  quad <- rowSums(Z^2)
  lin  <- j * as.vector(Z %*% delta) * sqrt(3/(1 + quad))
  
  p <- pt(lin, df = 3)
  
  X <- Z
  X[U > p, ] <- -Z[U > p, ]
  
  X
}

#### Objedinjeno ####

generate_sample <- function(setting, j, n) {
  
  if(setting == "normal_spherical") {
    Z <- gen_normal_kernel(n, Sigma_sph)
    return(skew_normal(Z, j, delta))
  }
  
  if(setting == "normal_elliptical") {
    Z <- gen_normal_kernel(n, Sigma_ell)
    return(skew_normal(Z, j, delta))
  }
  
  if(setting == "normal_cone") {
    Z <- gen_cone_kernel(n, gen_normal_kernel, Sigma_sph)
    return(skew_normal(Z, j, delta))
  }
  
  if(setting == "cauchy_spherical") {
    Z <- gen_cauchy_kernel(n, Sigma_sph)
    return(skew_cauchy(Z, j, delta))
  }
  
  if(setting == "cauchy_elliptical") {
    Z <- gen_cauchy_kernel(n, Sigma_ell)
    return(skew_cauchy(Z, j, delta))
  }
  
  if(setting == "cauchy_cone") {
    Z <- gen_cone_kernel(n, gen_cauchy_kernel, Sigma_sph)
    return(skew_cauchy(Z, j, delta))
  }
}


