## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)

## ----eval = FALSE-------------------------------------------------------------
#  install.packages('reticulate')
#  reticulate::install_miniconda()

## ----eval = FALSE-------------------------------------------------------------
#  # install.packages('mglasso')
#  remotes::install_github("desanou/mglasso")

## -----------------------------------------------------------------------------
#  library(mglasso)
#  install_pylearn_parsimony(envname = "rmglasso", method = "conda")
#  reticulate::use_condaenv("rmglasso", required = TRUE)
#  reticulate::py_config()

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
#  mat.correlation <- Matrix::bdiag(blocs)
#  corrplot::corrplot(as.matrix(mat.correlation), method = "color", tl.col="black")

## -----------------------------------------------------------------------------
#  set.seed(11)
#  X <- mvtnorm::rmvnorm(n, mean = rep(0,p), sigma = as.matrix(mat.correlation))
#  colnames(X) <- LETTERS[1:9]

## -----------------------------------------------------------------------------
#  X <- scale(X)
#  res <- mglasso(X, lambda1 = 0.2*n, lambda2_start = 0.1, fuse_thresh = 1e-3, verbose = FALSE)

## -----------------------------------------------------------------------------
#  temp <- mglasso::conesta(X, lam1 = 0.2*n, lam2 = 0.1)

## -----------------------------------------------------------------------------
#  library(ggplot2)
#  library(ggrepel)
#  mglasso:::plot_clusterpath(as.matrix(X), res)

## -----------------------------------------------------------------------------
#  plot_mglasso(res)

