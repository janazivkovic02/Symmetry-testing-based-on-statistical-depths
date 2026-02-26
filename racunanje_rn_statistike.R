library(MASS)
library(ddalpha)

#### Sortiranje za različite dubine ####

sortiranje_sim <- function(x) {
  x.prosireno <- rbind(x, -x)
  x.dubine <- depth.simplicial(x,x.prosireno)
  return(x[order(x.dubine),])
}

sortiranje_l2 <- function(x) {
  x.prosireno <- rbind(x, -x)
  x.dubine <- depth.spatial(x,x.prosireno)
  return(x[order(x.dubine),])
}

sortiranje_mah <- function(x) {
  x.prosireno <- rbind(x, -x)
  x.dubine <- depth.Mahalanobis(x,x.prosireno)
  return(x[order(x.dubine),])
}

sortiranje_proj <- function(x) {
  x.prosireno <- rbind(x, -x)
  x.dubine <- depth.projection(x,x.prosireno)
  return(x[order(x.dubine),])
}

sortiranje_zonoid <- function(x) {
  x.prosireno <- rbind(x, -x)
  x.dubine <- depth.zonoid(x,x.prosireno)
  return(x[order(x.dubine),])
}

#### Indikator ####

indikator <- function(a,b,c) {
  A <- rbind(cbind(a,b,c),c(1,1,1))
  b <- c(0,0,1)
  
  lambdas <- solve(A, b)
  
  #svi lambda moraju biti >= 0
  unutar <- all(lambdas >= -0.00000000001)
  
  return(as.numeric(unutar))
}

#### Računanje Rn ####

rn_sim <- function(x){
  r <- 1
  x <- sortiranje_sim(x)
  for (i in 3:nrow(x)) {
    r <- r + as.numeric(indikator(x[i,],x[i-1,], x[i-2,]))
  }
  return(r)
}

rn_l2 <- function(x){
  r <- 1
  x <- sortiranje_l2(x)
  for (i in 3:nrow(x)) {
    r <- r + as.numeric(indikator(x[i,],x[i-1,], x[i-2,]))
  }
  return(r)
}

rn_mah <- function(x){
  r <- 1
  x <- sortiranje_mah(x)
  for (i in 3:nrow(x)) {
    r <- r + as.numeric(indikator(x[i,],x[i-1,], x[i-2,]))
  }
  return(r)
}

rn_proj <- function(x){
  r <- 1
  x <- sortiranje_proj(x)
  for (i in 3:nrow(x)) {
    r <- r + as.numeric(indikator(x[i,],x[i-1,], x[i-2,]))
  }
  return(r)
}

rn_zonoid <- function(x){
  r <- 1
  x <- sortiranje_zonoid(x)
  for (i in 3:nrow(x)) {
    r <- r + as.numeric(indikator(x[i,],x[i-1,], x[i-2,]))
  }
  return(r)
}