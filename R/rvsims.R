# ========================================================================
# rvsims  -  generate a vector of rvs from a simulation matrix
# ========================================================================
# Name:        
# Description: 
# Parameters:  
# Required:    none
# History:     2004-06-  : 
# If the sims array is 1-dimensional,
# it is taken to be the vector of simulations for one variable;
# if sims is 2-dimensional (L x n)
# it contains L simulations for each of n variables.
# if the sims argument is a list, assumes that each component of the
# list is a single draw from a distribution. A component may be a list,
# when the list will be recursively built into a list of rv objects.
#

.recycle_index <- function (length.from, length.to) {
  if (length.from>=1 && length.to>=1) {
    return(1+(0:(length.to-1)) %% length.from)
  }
  return(0)
}

resizeSims <- function (s, vector.length, n.sims) # NOEXPORT
{
  # here the simulations must be in columns
  if (!missing(vector.length)) {
    nr <- nrow(s) # params
    if (nr!=vector.length) {
      ix <- .recycle_index(nr, vector.length)
      s <- s[ix, , drop=FALSE]
    }
  }
  nc <- ncol(s) # sims
  if (nc!=1 && !missing(n.sims)) {
    if (nc!=n.sims) {
      ix <- .recycle_index(nc, n.sims)
      s <- s[, ix, drop=FALSE]
    }
  }
  return(s)
}

.rvsims.list <- function (x, n.sims=getnsims(), permute = FALSE) 
{
  # Assume that all elements in the list have the same dimensions.
  # This may be modified later -- filling with NAs
  dx <- dim(x[[1]])
  s <- sapply(x, length)
  if (!all(s==s[1])) {
    # some elements had different dimensions!
    stop("Simulation list was not consistent")
  }
  S <- t(matrix(unlist(x, recursive=FALSE), nrow=s[1]))
  if (!is.list(S)) {
    r <- rvsims(S, n.sims=n.sims, permute=permute)
    dim(r) <- dx
    names(r) <- names(x[[1]])
    return(r)
  }
  r <- list()
  for (i in 1:ncol(S)) {
    r[[i]] <- .rvsims.list(S[,i])
  }
  names(r) <- names(x[[1]])
  return(r)
}

rvsims <- function(sims, n.sims=getnsims(), permute=FALSE) {
  if (is.list(sims)) {
    return(.rvsims.list(sims, n.sims=n.sims, permute=permute))
  }
  is_factor <- (is.character(sims) || is.factor(sims))
  if (is.factor(sims)) {
    sims[TRUE] <- as.character(sims)
  }
  if (! length(dim(sims)) %in% c(0, 2)) {
    stop("rvsims: Argument must be a vector or matrix or a list")
  }
  if (length(sims) < 1) {
    return(rv(0))
  }
  if (is.null(d.s <- dim(sims))) {
    d.s <- dim(sims) <- c(length(sims), 1)
  }
  n.sims.max <- d.s[1]
  if (length(d.s)==2) { # A regular matrix of simulations
    .names <- dimnames(sims)[[2]]
  } else {
    stop("Simulation array has >2 dimensions. Don't know how to deal with that.")
  }
  if (n.sims.max > 1) {
    if  (nrow(sims)!=n.sims) {
      sims <- t(resizeSims(t(sims), n.sims=n.sims))
    }
    if (permute) {
      .order <- sample(n.sims)
      sims <- sims[.order, , drop=FALSE]
    }
  }
  vec <- split(sims, col(sims))
  names(vec) <- .names
  class(vec) <- class(rv())
  if (is_factor) {
    as.rvfactor(vec)
  } else {
    return(vec)
  }
}
