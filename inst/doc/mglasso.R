## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)

## ----eval = FALSE-------------------------------------------------------------
#  remotes::install_github("desanou/mglasso")
#  library(mglasso)

## -----------------------------------------------------------------------------
#  library(mglasso)
#  install_conesta()

## -----------------------------------------------------------------------------
#  library(Matrix)
#  n = 50
#  K = 3
#  p = 9
#  rho = 0.85
#  blocs <- list()
#  
#  for (j in 1:K) {
#    bloc <- matrix(rho, nrow = p/K, ncol = p/K)
#    for(i in 1:(p/K)) { bloc[i,i] <- 1 }
#    blocs[[j]] <- bloc
#  }
#  
#  mat.covariance <- Matrix::bdiag(blocs)
#  Matrix::image(mat.covariance)

## -----------------------------------------------------------------------------
#  rep(1:3, each = 3)

## -----------------------------------------------------------------------------
#  set.seed(11)
#  X <- mvtnorm::rmvnorm(n, mean = rep(0,p), sigma = as.matrix(mat.covariance))

## -----------------------------------------------------------------------------
#  X <- scale(X)
#  res <- mglasso(X, lambda1 = 0.1, lambda2_start = 0.1, fuse_thresh = 1e-3)

## ----fig.width=10, fig.height=10----------------------------------------------
#  plot_mglasso(res)

## -----------------------------------------------------------------------------
#  res$out$level9$clusters
#  res$out$level7$clusters
#  res$out$level4$clusters
#  res$out$level3$clusters
#  res$out$level1$clusters

