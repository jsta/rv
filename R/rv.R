# ================================================================================
# rv - simulation-based random variable class in R
# Version 0.9*
# Updated 2005-11-05.
# (c) 2004-2005 Jouni Kerman <kerman@stat.columbia.edu>
# ================================================================================
#
# NAME
#   rv-base.r  -  simulation-based random variables in R, basic methods
# AUTHOR
#   Jouni Kerman, kerman@stat.columbia.edu
# SEE ALSO
#  other files in rv/R in this package
# 

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# PACKAGE MANAGEMENT
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

detachrv <- function ()
{
  detach("package:rv")
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GLOBAL SETTINGS
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.RV <- list()

.rvRegisterFunctionSwitch <- function (FUN, fname, namespace)
{
  origFUN <- getFromNamespace(fname, namespace)
  .RV[[namespace]][[fname]] <<- list(R=origFUN, rv=FUN)
}

.rvFunctionSwitch <- function (attach=TRUE, verbose=TRUE)
{
  detach <- (!attach)
  for (namespace in names(.RV)) {
    funcs <- .RV[[namespace]]
    for (fun.name in names(funcs)) {
      f <- funcs[[fun.name]]
      FUN <- if (attach) f$rv else f$R
fullname <- paste(namespace, ":::", fun.name, sep="")
      if (verbose) {
        cat(paste(namespace, ":::", fun.name, " ", sep=""))
        if (attach) cat("replaced by a function in package:rv\n")
        else cat("restored.\n")
      }
      assignInNamespace(fun.name, FUN, namespace)
    }
  }
  if (verbose) {
    cat("Package rv ")
    if (attach) cat("attached") else cat("detached")
    cat(".\n")
  }
}


# ========================================================================
# rvnsims - get or set the default number of simulations (a global variable)
# ========================================================================

rvnsims <- function(nsims)
{
  if (!missing(nsims)) {
    if (nsims[1]>0) {
      options('rvnsims'=ceiling(nsims[1]))
   } else
      stop('Invalid number of simulations ',nsims[1])
  }
  getOption('rvnsims')
}


# ========================================================================
# nsims - get the number of simulations of the components of an rv
# ========================================================================

nsims <- function(x)
{
  sapply(x, length)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# AN OBJECT OF TYPE 'rv'
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  rv(n)
#  

# ========================================================================
# .attr.rvsim : set variables of rvsim objects (internal use only)
# ========================================================================
#

.attr.rvsim <- function(x, attrib) # NOEXPORT
{  
  a <- attr(x, 'rvsim')
  if (missing(attrib)) return(a)
  a[[attrib]]
}


".attr.rvsim<-" <- function(x, attrib=NULL, value) # NOEXPORT
{
  a <- attr(x, 'rvsim')
  if (is.null(attrib))
    attr(x, 'rvsim') <- value
  else {
    a[[attrib]] <- value
    attr(x, 'rvsim') <- a
  }
  x
}

# ========================================================================
# .rvsim : make a random variable component (internal use only)
# ========================================================================
# Synopsis: .rvsim(sims)
# Parameters:  vector of numbers, NA, or NaN.
# Practically this is only a class wrapper
#

.rvsim <- function(sims, attrib=NULL) # NOEXPORT
{
  sims <- as.vector(sims)
  class(sims) <- 'rvsim'
  if (!is.null(attrib))
    .attr.rvsim(sims, NULL) <- attrib
  sims
}

# ========================================================================
# rv :  make a random vector
# ========================================================================
# Name:        rv(length=0)
# Description: 
# Parameters:  
# Required:    none
# History:     2004-06-  : 
#

rv <- function(length=0)
{
  x <- as.list(numeric(length))
  class(x) <- 'rv'
  x
}

# ========================================================================
# rvsims  -  generate a vector of rvs from a simulation matrix
# ========================================================================
# Name:        
# Description: 
# Parameters:  
# Required:    none
# History:     2004-06-  : 
# New feature: attributes
#   user-defined attributes are stored in the attribute
#   'rvsim' in the attribute list of a rvsim object.
#   currently implemented are Rhat and  n.eff,
# and 
#   'order', which returns the original order of the scrambled 
#   useful for recovering the original MCMC simulations.
# However this may be useless since we have already thinned and
# discarded burn.in iterations. 
# 
# If the sims array is 1-dimensional,
# it is taken to be the vector of simulations for one variable;
# if sims is 2-dimensional (L x n)
# it contains L simulations for each of n variables.
# If sims is 3-dimensional (L x m x n)
# it contains m chains with L simulations each for each of the n variables,
# and the rv will contain then n variables with L x m simulations each.
# The L x m simulations will be scrambled.
#
# NOTE
#   now allow non-numeric simulations, e.g. logicals are OK
#   maybe complex numbers work ok too, they should (?)
#

.rvsims.list <- function (x, permute=FALSE)
{
  # Assume that all elements in the list have the same dimensions.
  dx <- dim(x[[1]])
  s <- sapply(x,I)
  if (is.list(s)) {
    # some elements had different dimensions!
    stop("Simulation list was not consistent")
  }
  r <- rvsims(t(s), permute=permute)
  dim(r) <- dx
  r
}

rvsims <- function(sims, n.sims=rvnsims(), permute=FALSE, save.order=FALSE)
{
  if (is.list(sims)) return(.rvsims.list(sims, permute=permute))
  if (!is.array(sims) && !is.vector(sims)) {
    stop("rvsims: Argument must be a numeric vector or array")
  }
  if (length(sims)<1)
    stop("rvsims: Invalid simulation vector: zero length")
  if (is.null(d.s <- dim(sims))) {
    d.s <- dim(sims) <- c(length(sims),1)
  }
  n.sims.max <- d.s[1]
  if (length(d.s)==2) { # A regular matrix of simulations
    .names <- dimnames(sims)[[2]]
    n.vars <- d.s[2]
    att <- NULL
  } else if (length(d.s)==3) { # 3-way array <=> mcmc time series
    .names <- dimnames(sims)[[3]]
    n.chains <- d.s[2]
    n.vars <- d.s[3]
    n.sims.max <- d.s[1]*n.chains
    last.values <- sims[d.s[1],,]
    dim(sims) <- c(n.sims.max,n.vars)
    permute <- TRUE
    att <- list(n.chains=n.chains, last.values=last.values)
  } else {
    stop("Simulation array has >3 dimensions. Don't know how to deal with that.")
  }
  vec <- rv(n.vars)
  if (n.sims.max>1) {
    nsims <- rvnsims()
    if (n.sims<n.sims.max) {
      omit <- (-1):(n.sims-n.sims.max)
      sims <- sims[omit,,drop=FALSE]
    } else if (n.sims>n.sims.max) {
      warning("n.sims is larger than the available number of simulations")
      include <- sample(1:n.sims.max,size=n.sims-n.sims.max,replace=TRUE)
      include <- c(1:n.sims.max, include)
      sims <- sims[include,,drop=FALSE]
    }
    if (permute) {
      .order <- sample(n.sims)
      if (save.order) att$order <- order(.order)
      sims <- sims[.order,,drop=FALSE]
    }
  }
  for (i in 1:n.vars) {
    vec[[i]] <- .rvsim(sims[,i], attrib=att)
  }
  names(vec) <- .names
  vec
}

  ## does not work well: won't be erased when variable is changed.
  ## attr(vec, 'summary') <- summary(vec)
## Is there a need to force sims to be nsims length???
##if (length(s)!=nsims)
##  s <- rep(s,length.out=nsims)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ACCESSING SIMULATIONS
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ========================================================================
# sims  -  return the matrix of simulations
# ========================================================================
#

sims <- function(x, n.sims=NULL, dimensions=FALSE, sim.matrix=FALSE, mc.array=FALSE)
{
  UseMethod('sims')
}

sims.rvsim <- function(x, ...) # NOEXPORT
{
  unclass(x) # don't drop attributes
}

sims.default <- function(x, ...) # NOEXPORT
{
  ## We should not set the dimension for constants.
  ## re: problems in max.rvsims: cbind(sims(1),sims(x)) does not work
  ## since sims(1) returns a 1x1 matrix but sims(x) may be L x n matrix.
  as.vector(x)  # drop attributes
}

.mcsims <- function(x, n.sims=NULL) # NOEXPORT
{
  if (!is.rv(x)) return(NULL)
  if (any(is.na(rvnchains(x)))) {
     stop(".mcsims: mc time series not available for ",deparse(substitute(x)))
  }
  s <- lapply(x, sims.rvsim)
  if (is.null(n.sims)) {
    n.sims <- max(sapply(s,length))
  }
  n.chains <- NA
  for (i in seq(along=s)) {
    sim <- s[[i]]
    a <- attr(sim, 'rvsim')
    if (is.na(n.chains)) n.chains <- a$n.chains
    .order <- a$order
    if (n.chains != .rvsim$n.chains || length(.order) != length(sim)) {
      stop("sims.rv: mc time series for simulations are unavailable",deparse(substitute(x)))
    }
    sim <- sim[.order]
    dim(sim) <- dim(.order)
    if (length(sim) != n.sims) sim <- rep(sim,length.out=n.sims)
    s[[i]] <- sim
  }
  m <- array(unlist(s), c(n.sims, n.chains, length(s)))
  dimnames(m) <- list(NULL, NULL, names(x))
  m
}


# ========================================================================
# sims.rv  -  return the matrix of simulations for an rv
# ========================================================================
# gives the simulations of a vector as a matrix,
# used e.g. for computing quantiles of a vector.
# length(x) columns, rvnsims rows.
# length of sims is adjusted to the longest sims length
# if dimensions=TRUE, will return an L x dim(x) matrix.
# (this is useful for simapply applied to random matrices, e.g. see determinant.rv)
#
#

sims.rv <- function(x, n.sims=NULL, dimensions=FALSE, sim.matrix=FALSE, mc.array=FALSE)
{
  if (mc.array) {
    # A 3-way matrix with original order.
    return(.mcsims(x, n.sims=n.sims))
  }
  xl <- sapply(x, length) # Lengths of the simulation vectors of each component.
  if (is.null(n.sims)) {
    n.sims <- max(xl)
  }
  s <- lapply(x, sims.rvsim)
  for (i in seq(along=xl)) {
    if (xl[i] != n.sims) s[[i]] <- rep(s[[i]],length.out=n.sims)
  }
  m <- matrix(unlist(s), nrow=n.sims)
  if (dimensions && !is.null(dim.x <- dim(x))) {
    n.s <- length(m)/prod(dim.x)
    cm <- m
    dimnames(cm) <- NULL
    m <- array(m, c(n.s,dim.x)) # multi-way array, 1st dimension is the dimension "rows"
    if (!is.null(dn <- dimnames(x))) {
      dimnames(m)[1+seq(along=dn)] <- dn ## put dimnames in correct positions
    } 
    if (sim.matrix) {
      attr(m, 'sim.matrix') <- cm # This is needed for simapply
    }
  } else {
    dimnames(m) <- list(NULL, names(x))
  }
  m
}


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# BASIC RELATIONS
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  is.rv
#  as.rv
#  is.na.rv
#  

# ========================================================================
# is.rv - is the argument an object of the class rv?
# ========================================================================
# Name:        is.rv(x)
# Description: 
# Parameters:  
# Required:    none
# History:     2004-06-  : 
#

is.rv <- function(x)
{
  UseMethod("is.rv")
}

is.rv.rv <- function(x)
{
  TRUE
}

is.rv.default <- function(x)
{
  FALSE
}

# ========================================================================
# .is.rvsim - is the argumentan rvsim object?
# ========================================================================

.is.rvsim <- function (x)
{
  (class(x)=='rvsim')
}

# ========================================================================
# is.random - is a component of a vector a random variable?
# ========================================================================
# Name:        is.random(x)
# Description: returns a vector of length(x) consisting of TRUE or FALSE
#              depending on whether the argument is random or not
# Note:        does not give false if all simulations of a component are equal
# Parameters:  x : vector, numeric or rv
# Required:    none
# History:     2004-06-  : 
#

is.random <- function(x)
{
  sapply(x, .is.rvsim)
}

# ========================================================================
# is.na.rv - is a component of a vector missing?
# ========================================================================
# Outputs a Bernoulli random vector of equal length to that of the argument..
# Parameters:  an rv
# Related   :  is.na.rvsim
# History:     2004-06-  : 
#

is.na.rv <- function(x)
{
  simapply(x, is.na)
}


# ========================================================================
# anyisrv - is any of the arguments a RV?
# ========================================================================
#

anyisrv <- function(...) # NOEXPORT
{
  'rv' %in% unlist(lapply(list(...),class))
}

# ========================================================================
# anyisrvsim - is any of the components in the argument a rvsim?
# ========================================================================
#

anyisrvsim <- function(x, na.rm=FALSE) # NOEXPORT
{
  'rvsim' %in% unlist(lapply(x,class))
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# COERCING
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ========================================================================
# as.constant  -  coerce components to constant if possible
# ========================================================================
# Name:        
# Description: returns NA if a component cannot be coerced to constant
# Parameters:  
# Required:    none
# History:     2004-06-  : 
#

as.constant <- function(x)
{
  if (is.numeric(x))
    return(x)
  else
    UseMethod('as.constant')
}

as.constant.rv <- function(x)
{
  sapply(x, as.constant)
}

as.constant.rvsim <- function(x)
{
  if (is.constant(x))
    return(x)
  else
    return(NA)
}

# ========================================================================
# as.rv  -  coerce argument to an rv
# ========================================================================
# Name:        
# Description: 
# Parameters:  
# Required:    none
# History:     2004-06-  : 
#

as.rv <- function(x)
{
  UseMethod("as.rv")
}

as.rv.rv <- function(x)
{
  x
}

as.rv.numeric <- function(x)
{
  r <- rvsims(matrix(x,nrow=1))
  cr <- class(r)
  attributes(r) <- attributes(x)
  class(r) <- cr
  r
}

as.rv.logical <- function(x)
{
  as.rv.numeric(x)
}

as.rv.integer <- function(x)
{
  as.rv.numeric(x)
}

as.rv.list <- function(x)
{
  y <- NULL
  for (i in 1:length(x))
    y <- c(y,as.rv(x[[i]]))
  y
}

as.rv.matrix <- function(x)
{
  dim.x <- dim(x)
  r <- as.rv.numeric(as.vector(x))
  dim(r) <- dim.x
  r
}

as.rv.default <- function(x, ...)
{
  stop('Cannot coerce object of class "', class(x), '" to rv')
}

# DEBUG:1 ?? The permutation of !is.null(dim(sims)) is not done yet
# DEBUG:2 (must take permutation of row indices and then permute.)


# ========================================================================
# as.logical  -  coerce to logical
# ========================================================================

#as.logical.rv <- function(x, ...)
#{
#  simapply(x, as.logical)
#}

# ========================================================================
# as.vector.rv  -  coerce to vector
# ========================================================================

as.vector.rv <- function(x, mode="any")
{
  dim(x) <- NULL # DEBUG as.vector.rv: Verify what as.vector REALLY does.
  x
}


# ========================================================================
# is.constant  -  which components are constants?
# ========================================================================
# Value:  a logical vector (not rv), TRUE for each component that has variance 0.
#

is.constant <- function(x)
{
  UseMethod('is.constant')
}

is.constant.rv <- function(x)
{
  sapply(x, is.constant)
}

is.constant.rvsim <- function(x)
{
  s <- sims.rvsim(x)
  if (length(s)==1)
    return(TRUE)
  return(var(as.vector(s)) == 0)
}

is.constant.default <- function(x)
{
  TRUE # Is it really true?
}


# ========================================================================
# rvmatrix - rv version of 'matrix'
# ========================================================================

rvmatrix <- function(...)
{
  x <- matrix(...)
  if (anyisrvsim(x)) {
     class(x) <- class(rv())
  }
  x
}

# ========================================================================
# rvarray - rv version of 'array'
# ========================================================================

rvarray <- function(...)
{
  x <- array(...)
  if (anyisrvsim(x)) {
     class(x) <- class(rv())
  }
  x
}

# ========================================================================
# length.rv - length of a random variable
# ========================================================================
# Length and dimension may be random!

#length.rv <- function(x)
#{
#  
#}


# ========================================================================
# dim.rv - dimensions of a random variable
# ========================================================================

#dim.rv <- function(x)
#{
#  
#}


# ========================================================================
# rvattach - attach rv's, based on their names.
# ========================================================================

rvattach <- function (what, name='rvattach', ...)
{
  ## Avoid duplicates
  if (any(search()==name)) eval(substitute(detach(name),list(name=as.name(name))))
  if (!is.rv(what)) stop("Argument must be an rv")
  a <- split(what)
  if (is.null(a)) return(NULL)
  for (j in 1:length(a)) { # Code from A. Gelman's attach.all.
    if (names(a)[j] %in% ls(.GlobalEnv))
      remove (list=names(a)[j], envir=.GlobalEnv)
  }
  attach (a, name=name, ...)
}


# ----------------
# end of rv-base.R
# ----------------

# rv-rbase.R

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# VECTOR INDEXING AND ASSIGNMENT OPERATIONS
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  [<-.def
#  [.rv
#  [<-.rv
#  $.rv
#  $<-.rv
#

# ========================================================================
# [<- - new replacement method for "[" assignment
# ========================================================================
# Intercepted because we want to change plain vectors into rvs if assigned a rv to
# Description: 
# Parameters:  
# Required:    none
# History:     2004-06-  : 
# 2004.05.10 !! How to allow other classes to use this function??? Should work???
# 2004.05.10 seems to work, keep testing



.rvBracketAssignment <- function(x, ..., value=NULL)
{
  v <- value # Necessary, or else 'value' will be evaluated twice:
  # Try for example a<-array(NA,c(1,5)); a[1,] <- { cat("hi!"); rnorm(5) }
  # will print hi! twice!!!
  if (is.rv(value))
    call <- substitute("[<-.rv"(x, ..., value=v))
  else
    call <- substitute(.Primitive("[<-")(x, ..., value=v))
  eval.parent(call)
}

"[<-" <- .rvBracketAssignment # EXPORT "[<-"

.rvRegisterFunctionSwitch(.rvBracketAssignment, "[<-", 'base')


# ========================================================================
# [.rv  -  rv component retrieval by index
# ========================================================================
# Name:        
# Description: 
# Parameters:  
# Required:    none
# History:     2004-06-  : 
#
# These functions must warn if l.x and l.y are not compatible lengthwise!
# 


"[.rv" <- function(x, ..., drop=TRUE)
{
  if (missing(x)) { ### ?
    if (drop) 
      return(drop(x))
     else
       return(x)
  }
  call <- substitute("["(unclass(x), ..., drop=drop))
  v <- eval.parent(call)
  v[ sapply(v, is.null) ] <- NA
  class(v) <- class(rv())
  v
}

# ========================================================================
# [<-.rv  -  rv component assignment by index
# ========================================================================
#

"[<-.rv" <- function(x, ..., value)
{
  call <- substitute(.Primitive("[<-")(unclass(x), ..., value=value))
  y <- eval.parent(call)
  class(y) <- 'rv' ##
  y
}

"[<<-.rv" <- function(x, ..., value)
{
  call <- substitute(.Primitive("[<<-")(unclass(x), ..., value=value))
  y <- eval.parent(call)
  class(y) <- 'rv' ##
  y
}

# ========================================================================
# $.rv  -  obtain rv components by a grep pattern
# ========================================================================
# It's better to implement $.rv with this kind of simple
# program rather than to use "$"(unclass(x),name),
# since "$" strips the names and de-lists x.
#

"$.rv" <- function(x, name) 
{
  if (missing(name))
    return(x)
  if (length(name)==0)
    return(NULL)
  nas <- names(x)
  if (is.null(nas))
    return(NULL)
  vec <- rv()
  class(vec) <- class(x)
  pattern <- paste("^", name, sep="")
  i <- grep(pattern, nas)
  if (length(i)==0) return(NULL)
  vec <- x[i]
  vec
} # $.rv


# ========================================================================
# $<-.rv  -  assign to rv components by name
# ========================================================================
# Name:        
# Description: 
# Parameters:  
# Required:    none
# History:     2004-06-  : 
#

"$<-.rv" <- function(x, name, value) 
{
  if (missing(name))
    return(x)
  if (length(name)==0)
    return(NULL)
  nas <- names(x)
  if (is.null(nas))
    return(NULL) # Maybe warning?
  vec <- rv()
  class(vec) <- class(x)
  pattern <- paste("^", name, sep="")
  i <- grep(pattern, nas)
  if (length(i)==0) return(NULL)
  if (is.null(value))
    x <- x[-i]
  else
    x[i] <- value
  x
} # $<-.rv



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# CONCATENATION
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  c

# ========================================================================
# c  -  concatenation of objects method
# ========================================================================
# - It was necessary to trap the 'c' function since
#   c(1,x) would not have worked; of course c(x,1) would work.
#
# Unfortunately c.rv is slow compared to the primitive 'c'.
# 

.rvConcatenate <- function(..., recursive=FALSE)
{
  x <- .Primitive("c")(..., recursive=recursive)
  if (is.list(x) && anyisrvsim(x)) # Was there a 'rvsim' in there?
    class(x) <- class(rv())
  x
}

"c" <- .rvConcatenate # EXPORT "c"

.rvRegisterFunctionSwitch(.rvConcatenate, 'c', 'base')

## DON'T : see zzz.R ## c <- .rvConcatenate


# ========================================================================
# cbind  -  column bind for rvs
# ========================================================================
# Note: It's inconvenient that we cannot call the generic (default) cbind
#       if class attributes are set.
#

cbind.rv <- function(..., deparse.level = 1)
{
  if (deparse.level != 1) 
    .NotYetUsed("deparse.level != 1")
  x <- list(...)
  if (length(x)<1) return(NULL)
  v <- NULL
  for (i in seq(along=x)) {
    v <- cbind(v, unclass(x[[i]]), deparse.level=deparse.level)
  }
  class(v) <- class(rv())
  v
}


# ========================================================================
# rvbind.rv  -  row bind for rvs
# ========================================================================
# Name:        
# Description: 
# Parameters:  
# Required:    none
# History:     2004-06-  : 
#

rbind.rv <- function(..., deparse.level = 1)
{
  if (deparse.level != 1) 
    .NotYetUsed("deparse.level != 1")
  x <- list(...)
  if (length(x)<1) return(NULL)
  v <- NULL
  for (i in seq(along=x)) {
    v <- rbind(v, unclass(x[[i]]), deparse.level=deparse.level)
  }
  class(v) <- class(rv())
  v
}

 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# SOME FUNCTIONS
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   mean.rv
#   sort.rv
#   sample.rv
#

# ========================================================================
# mean.rv - mean (rowwise)
# ========================================================================

mean.rv <- function(x, ...)
{
  rvsims(rowMeans(sims(x))) # Much faster than apply(.., mean).
}


# ========================================================================
# %*% - matrix product
# ========================================================================
# DEBUG: how to make this work with outer()?
#

.rvMatrixProductOLD <- function(x,y)
{
  if (!is.rv(x) && !is.rv(y))
    return(.Primitive("%*%")(x,y))
  n.sims <- .Internal(max(nsims(x),nsims(y), na.rm=FALSE))
  xsim <- sims(as.rv(x), dimensions=TRUE, n.sims=n.sims)
  ysim <- sims(as.rv(y), dimensions=TRUE, n.sims=n.sims)
  X <- apply(xsim, 1, list) # Make into a list
  Y <- apply(ysim, 1, list) # Make into a list
  .F <- function (x, y) .Internal(crossprod(t(x[[1]]), y[[1]]))
  XY <- mapply(.F, X, Y, SIMPLIFY=FALSE) # A list of matrices.
  m <- rvsims(XY)
  m
}

.rvMatrixProduct <- function(x,y)
{
  if (!is.rv(x) && !is.rv(y)) return(.Primitive("%*%")(x,y))
  n.sims <- .Internal(max(nsims(x),nsims(y), na.rm=FALSE))
  ysim <- sims(as.rv(y), dimensions=TRUE, n.sims=n.sims)
  d <- dim(y)
  if (!is.rv(x) && (is.null(d)) || (length(d)==2 && d[2]==1)) {
    # Typical case: constant matrix times a rv vector
    XY <- t(.Primitive("%*%")(x,t(ysim)))
  } else {
    xsim <- sims(as.rv(x), dimensions=TRUE, n.sims=n.sims)
    X <- apply(xsim, 1, list) # Make into a list
    Y <- apply(ysim, 1, list) # Make into a list
    .F <- function (x, y) .Internal(crossprod(t(x[[1]]), y[[1]]))
    XY <- mapply(.F, X, Y, SIMPLIFY=FALSE) # A list of matrices.
  }
  rvsims(XY)
}

"%*%" <- .rvMatrixProduct # EXPORT "%*%"

.rvRegisterFunctionSwitch(.rvMatrixProduct, "%*%", 'base')

# ========================================================================
# split.rv - split an rv
# ========================================================================

split.rv <- function(x, f, drop=FALSE, ...)
{
  if (missing(f))
    f <- .shortnames(x)
  if (is.null(f))
    stop('Cannot split on null vector')
  split.default(x, f)
}

# ========================================================================
# sort  -  generate order statistics
# ========================================================================
#

.sort <- getFromNamespace('sort', 'base')

sort <- function(x, ...)
{
  if (is.rv(x))
    return(simapply(x, .sort, ...))
  .sort(x, ...)
}

.rvRegisterFunctionSwitch(sort, "sort", 'base')

# ========================================================================
# is.vector  -  coerce to vector
# ========================================================================

is.vector <- function(x, mode="any")
{
  if (is.rv(x)) return(is.null(dim(x)))
 .Internal(is.vector(x, mode))
}

.rvRegisterFunctionSwitch(is.vector, "is.vector", 'base')


# ========================================================================
# is.atomic.rv  -  rv is an 'atomic' type.
# ========================================================================

is.atomic.rv <- function(x)
{
  TRUE
}

# ========================================================================
# min  -  minimum
# ========================================================================

.min <- getFromNamespace('min', 'base')

min <- function(..., na.rm=FALSE)
{
  if (anyisrv(...))
    return(simapply(cbind.rv(...), .min, na.rm=na.rm))
  .min(..., na.rm=na.rm)
}

.rvRegisterFunctionSwitch(min, "min", "base")

# ========================================================================
# max  -  maximum
# ========================================================================

.max <- getFromNamespace('max', 'base')

max <- function(..., na.rm=FALSE)
{
  if (anyisrv(...))
    return(simapply(cbind.rv(...), .max, na.rm=na.rm))
  .max(..., na.rm=na.rm)
}

.rvRegisterFunctionSwitch(max, "max", "base")


# ========================================================================
# pmin  -  parallel minimum
# ========================================================================

.pmin <- getFromNamespace('pmin', 'base')

pmin <- function(..., na.rm=FALSE)
{
  if (anyisrv(...))
    return(simapply(cbind.rv(...), .pmin, na.rm=na.rm))
  .pmin(..., na.rm=na.rm)
}

.rvRegisterFunctionSwitch(pmin, "pmin", "base")

# ========================================================================
# pmax  -  parallel maximum
# ========================================================================

.pmax <- getFromNamespace('pmax', 'base')

pmax <- function(..., na.rm=FALSE)
{
  if (anyisrv(...))
    return(simapply(cbind.rv(...), .pmax, na.rm=na.rm))
  .pmax(..., na.rm=na.rm)
}

.rvRegisterFunctionSwitch(pmax, "pmax", "base")



# ----------------
# end of rv-rbase.R
# ----------------

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# UTILITIES in the package rv
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  .matrix.ragged.array(x)
#  .expand.as.matrix(M)
#  .listByName(x)
#  .permut
#  .dim.index
#  .bracket
#  .rowmax
#  .slice

# ========================================================================
# .matrix.ragged.array - make a matrix out of a ragged array
# ========================================================================
# not used now

.matrix.ragged.array <- function(x) # NOEXPORT
{
  # DEBUG: need something more efficient here!
  ncol <- max(sapply(x,length))
  s <- array(NA, c(length(x),ncol))
  for (i in seq(length(x))) {
    y <- x[[i]]
    s[i,seq(length(y))] <- y
  }
  s
}

# ========================================================================
# .expand.as.matrix - build a large matrix out of smaller matrices or scalars, diagonally
# ========================================================================
# Useful when building a large covariance matrix from a series of small ones.
#
# Parameters:  A rv or a list.
# Required:    none
# History:     2004-06-  : 
#

.expand.as.matrix <- function(M) # NOEXPORT
{
  # M is a rv or a list.
  dims <- c(0,0)
  for (i in 1:length(M)) {
    dims <- dims + if (is.null(dim(M[[i]]))) rep(length(M[[i]]),2) else dim(M[[i]])
  }
  m <- matrix(0,nrow=dims[1],ncol=dims[2])
  col.offset <- 0
  row.offset <- 0
  for (i in 1:length(M)) {
    x <- if (is.null(dim(M[[i]]))) diag(M[[i]],nrow=length(M[[i]])) else M[[i]]
    d <- dim(x)
    m[row.offset+1:d[1],col.offset+1:d[2]] <- x
    row.offset <- row.offset + d[1]
    col.offset <- col.offset + d[2]
  }
  m
}



# .permut: index of multidimensional array.
# 
# Incidentally we can get the index of the vector by
#  z <- .permut(dims)
#  index <- 1 + (z-1) %*% cumprod(z-1)
#  so index is a vector 1:prod(dims).
#

# ========================================================================
# .permut  -  make  permutations c(1,1) ... c(m,n)
# ========================================================================
# Example: 
#   .permut(c(2,3)) gives a matrix
#        [,1] [,2]
#   [1,]    1    1
#   [2,]    2    1
#   [3,]    1    2
#   [4,]    2    2
#   [5,]    1    3
#   [6,]    2    3
#

.permut <- function(dims) # NOEXPORT
{
  tl <- prod(dims)
  co <- NULL
  p <- 1
  for (i in 1:length(dims)) {
    d <- dims[i]
    x <- rep(1:d, each=p, length.out=tl)
    co <- cbind(co, x)
    p <- p*d
  }
  dimnames(co) <- NULL
  co
}

# ========================================================================
# .dim.index  -  make a list of indices of a matrix
# ========================================================================
#

.dim.index <- function(x) # NOEXPORT
{
  if (is.null(dim(x)))
    ix <-.permut(length(x))
  else 
    ix <- .permut(dim(x))
  ixt <- paste('[',apply(ix, 1, paste, collapse=','),']',sep='')
  # pretty-print
  m.ixt <- max(nchar(ixt))
  a <- nchar(ixt)<m.ixt
  if (any(a)) {
    d <- m.ixt - nchar(ixt)
    ixt[a] <- paste(' ', ixt[a], sep='')
  }
  ixt
}


# ========================================================================
# .bracket  -  x[[ i ]] but recycle if i > length(x)
# ========================================================================
# Name:        
# Description: 
# Parameters:  
# Required:    none
# History:     2004-06-  : 
#

.bracket <- function(x,k) # NOEXPORT
{
  l <- length(x)
  if(k<=l) x[[k]] else x[[1+((k-1)%%l)]]
}

# ========================================================================
# .slice  -  apply "[ row ]" to arguments but recycle if row > length(argument)
# ========================================================================
# 
# Description: Like lapply(list(...),"[",row) but recycles.
# Parameters:  
# Required:    none
# History:     2004-06-  : 
#

.slice <- function(..., row) # NOEXPORT
{
  UseMethod('.slice')
}

.slice.default <- function(..., row) # NOEXPORT
{
  v <- list(...)
  v
}

.slice.rv <- function(..., row) # NOEXPORT
{
  v <- sapply(list(...), bracket, row)
  class(v) <- class(rv())
  v
}


# ========================================================================
# .rowmax  -  this may be obsolete... but check!
# ========================================================================
#

.rowmax <- function(..., na.rm=FALSE) # NOEXPORT
{
  # Treats each argument as a column vector and
  # returns row-wise max(...) distributions.
  # We could make row-versions of quantiles, etc. too!
  # So write rowapply --- a generic version that works rowwise.
  # SO WRITE APPLY that works with RVs!!!
  #   This apply.rv should be more flexible and allows mixing of scalars and matrices.
  return(rowapply(cbind(...), max))
  x <- list(...)
  l.x <- length(x)
  l.z <- max(unlist(lapply(x,length)))
  z <- rv(l.z)
  for (i in 1:l.z)
    z[[i]] <- max(.slice(...,row=i),na.rm=na.rm)[[1]]
  z
}


# ========================================================================
# shortnames - names of the components of a rv, with brackets removed
# ========================================================================

.shortnames <- function(x)
{
  na <- names(x)
  if (is.null(na)) return(na)
  sapply(strsplit(na,'\\['),'[',1)
}

# ========================================================================
# longnames - names of the components of a rv, with brackets removed
# ========================================================================

.longnames <- function(x)
{
  na <- names(x)
  if (is.null(na)) return(na)
  sn <- sapply(strsplit(na,'\\['),'[',1) # shortnames
  spn <- split(na, sn)
  for (i in seq(along=spn)) {
    ## .dim.index(spn[[i]])
  }
  ## NOT YET DONE
}




# end rv-util.R

# rv-alias.R - Convenience functions

#
# Convenience functions.
#
# These return always numbers, not r.v.s
#
# rvmean : componentwise means
# rvmedian : componentwise medians
# rvquantile : componentwise quantiles
# rvvar  : componentwise variances
# rvsd   : componentwise standard deviations
# E    : Expected value (alias for rvmean)
# Pr   : Probability (alias for rvmean)
# rv.sample : componentwise samples
#
# Name:        
# Description: 
# Parameters:  
# Required:    none
# History:     2004-06-  : 
#

rvmean <- function(x, ...)
{
  m <- colMeans(sims(x)) # Much faster
  dx <- dim(x)
  if (!is.null(dx) && prod(dx)==length(m)) {
    dim(m) <- dx
  }
  m
}

E <- function(x, ...)
{
  m <- colMeans(sims(x)) # Much faster
  dx <- dim(x)
  if (!is.null(dx) && prod(dx)==length(m)) {
    dim(m) <- dx
  }
  m
}

Pr <- function(x, ...)
{
  m <- colMeans(sims(x)) # Much faster
  dx <- dim(x)
  if (!is.null(dx) && prod(dx)==length(m)) {
    dim(m) <- dx
  }
  m
}

rvvar <- function(x, ...) rvsimapply(x, var, ...)

rvsd <- function(x,...) rvsimapply(x, sd, ...)

rvquantile <- function(x, ...) rvsimapply(x, quantile, ...)

rvmedian <- function(x, ...)  rvsimapply(x, median, ...)

rvmin <- function(x, ...) rvsimapply(x, min, ...)

rvmax <- function(x, ...) rvsimapply(x, max, ...)

rvsample <- function(x, size=1, replace=TRUE, prob=NULL) {
  rvsimapply(x, sample, size=size, replace=replace, prob=prob)
}

rvrange <- function (x, ...) rvsimapply(x, range)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# APPLYING
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ========================================================================
# .rvlapply  -  apply a function to components of an rv
# ========================================================================
# Name:        .rvlapply(x, FUN, ...)
# Description: the 'rv' version of lapply
#

.rvlapply <- function(x, FUN, ...)
{
  attr.x <- attributes(x)
  v <- lapply(x, FUN, ...)
  attributes(v) <- attr.x
  if (!(class(rv()) %in% class(v))) class(v) <- c(class(rv()),class(v))
  v
}

# ========================================================================
# .rvMathapply  -  apply a function to components of an rv or non-rv
# ========================================================================
# Name:        .rvlapply(x, FUN, ...)
# Description: the 'rv' version of lapply
#

.rvMathapply <- function(x, FUN, ...)
{
  v <- lapply(x, FUN, ...)
  class(v) <- class(rv())
  v
}

# ========================================================================
# simapply  -  apply a (numeric) function to the simulations, rowwise, with dimensions
# ========================================================================

simapply <- function(x, FUN, ...)
{
  s <- sims(x, dimensions=TRUE, sim.matrix=TRUE)
  ## First try the function on the first row of simulations.
  ## Obtain the dimensions of the results (apply will lose them!)
  ## This will be helpful for simapply(x, eigen), simapply(x, svd) etc.
  dims <- dim(s)
  rowlen <- prod(dims[-1])
  if (!is.null(sm <- attr(s, "sim.matrix"))) {
    first.row <- (t(sm))[seq(length=rowlen)] # First row
  } else {
    first.row <- (t(s))[seq(length=rowlen)] # First row
  }
  dim(first.row) <- dim(x)
  m0 <- (match.fun(FUN))(first.row, ...)
  ## Now, do the actual applying:
  m <- apply(s, 1, FUN, ...)
  if (is.list(m)) {
      ## For an example try simapply(x, eigen) where X is a random matrix.
      s <- list()
      for (i in seq(along=m[[1]])) {
        m1 <- sapply(m,'[[',i)
        if (is.list(m1)) stop("simapply: array expected, got list instead")
        s[[i]] <- if (is.null(dim(m1))) rvsims(m1) else rvsims(t(m1))
        dim(s[[i]]) <- dim(m0[[i]]) # Restore the dimensions from m0
      }
      names(s) <- names(m[[1]])
      return(s)
  }
  if (is.null(dim(m))) {
    if (!is.null(dim(m0))) stop("[simapply] A highly unlikely error! Contact the author...")
    return(rvsims(m))
  }
  r <- rvsims(t(m))
  if (!is.null(dim(m0))) {
    dim(r) <- dim(m0)
  }
  r
}

# ========================================================================
# rvsimapply  -  apply a (numeric) function to the simulations, columnwise
# ========================================================================
# 

rvsimapply <- function(x, FUN, ...)
{
  dx <- dim(x)
  n <- length(x)
  if (n==0)
    return(NULL)
  s <- sims(x)
  m <- apply(s, 2, FUN, ...)
  # Should we have instead an extra dimension for the vector that we receive?
  if (!is.null(dx) && prod(dx)==length(m))
    dim(m) <- dx
  m
}

# ========================================================================
# matsimapply  -  apply a (numeric) function to the simulations, pairwise/rowwise
# ========================================================================
# Description: 
# Parameters:  
# Required:    none
# History:     2004-06-  : 
#
# These routines are used for the %*%, cov, and other matrix functions requiring 2 params

# This works with functions requiring two matrix arguments, such as %*%.

matsimapply <- function(X, Y, FUN, ...)
{
  n.sims <- max(nsims(X),nsims(Y))
  x <- sims(as.rv(X), dimensions=TRUE, n.sims=n.sims)
  y <- sims(as.rv(Y), dimensions=TRUE, n.sims=n.sims)
  X. <- apply(x, 1, list)
  Y. <- apply(y, 1, list)
  MF <- match.fun(FUN)
  .F <- function (x, y, ...) list(MF(x[[1]], y[[1]], ...))
  XY <- mapply(.F, X., Y., ...)
  rvsims(XY)
}


# rv-Math.R - standard math functions for the rv class

Math.rv <- function(x, ...) {
  .rvMathapply(x, .Generic, ...)
}

Math.rvsim <- function(x, ...) {
  sim <- get(.Generic)(sims.rvsim(x), ...) # Must be vectorizable
  .rvsim(sim)
}

# cumsum, cumprod, cummax, cummin

cumsum.rv <- function (x)
{
  simapply(x, cumsum)
}

cumprod.rv <- function (x)
{
  simapply(x, cumprod)
}

cummin.rv <- function (x)
{
  simapply(x, cummin)
}

cummax.rv <- function (x)
{
  simapply(x, cummax)
}


# ----------------
# end of rv-Math.R
# ----------------

# rv-Ops.R - standard math functions for the rv class
# TODO
#  


Ops.rv <- function(e1, e2=NULL)
{
  e.attr <- attributes(e1)
  if (is.null(e.attr)) e.attr <- attributes(e2) # e1 can be a plain number
  if (is.null(e2)) {
    v <- mapply(.Generic, 0, e1, SIMPLIFY=FALSE, USE.NAMES=FALSE)
  } else {
    v <- mapply(.Generic, e1, e2, SIMPLIFY=FALSE, USE.NAMES=FALSE)
  }
  attributes(v) <- e.attr
  v
}

Ops.rvsim <- function(e1, e2)
{
  sim <- get(.Generic)(sims.rvsim(e1), sims.rvsim(e2)) # Must be vectorizable
  .rvsim(sim)
}


"!.rv" <- function(e1) 
{
  .rvlapply(e1, .Generic)
}

"!.rvsim" <- function(e1) # NOEXPORT
{
  sim <- ! sims(e1)
  .rvsim(sim)
}



# ---------------
# end of rv-Ops.R
# ---------------


# rv-Summary.R
# Summary functions for the rv class
# 
# Summary.rv Applies to
#  sum, prod
#  any, all
#  min, max
#  rv, (range); 
# We have a separate method for range().
#
# These functions take an arbitrary number of arguments and apply a function
# to each row of simulations at a time
#
# (c) Jouni Kerman 2004-2005

#
# sum, prod, any, all, min, max
#


Summary.rv <- function(..., na.rm=FALSE) {
  sims <- sims(c(...)) # an L x n matrix of numbers
  m <- apply(sims,1,.Generic,na.rm)
  rvsims(m)
}

range.rv <- function(..., na.rm=FALSE, finite=FALSE) {
  sims <- sims(c(...)) # an L x n matrix of numbers
  m <- apply(sims,1,'range',na.rm,finite) # gives a 2xL matrix!!!
  r <- rvsims(t(m)) # Transpose to get an L x 2 matrix.
  names(r) <- c('min','max')
  r
}
# rv-distr.R
#  generators of iid distributions
# 

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GENERATION OF RVS
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# utilities: .nsimmat, .parammat, .rvgen

# ========================================================================
# .nsimmat - make a matrix consisting of "random number of entries"
# ========================================================================
# Arguments:
# sim.n = lengths of random variables
# s     = values 
# and a simulation matrix of parameter simulations
# creates the corresponding matrix of parameters
# Note: used in rvparam()

.nsimmat <- function(sim.n, s) # NOEXPORT
{
  n.sims <- nrow(s)
  n.max <- max(sim.n)
  sim <- matrix(rep(s, times=n.max),nrow=n.sims)
  mc <- matrix(rep(seq(from=1,length=n.max),each=length(s)),nrow=n.sims)
  nc <- matrix(rep(sim.n,times=ncol(sim)),nrow=n.sims)
  sim[nc<mc] <- NA
  sim
}

# ========================================================================
# .parammat - make a parameter matrix (used in .rvgen)
# ========================================================================

# 

.parammat <- function(param, n) # NOEXPORT
{
  n.sims <- rvnsims()
  sim <- sims(as.rv(param), n.sims=n.sims) # guarantees an n.sims x K matrix
  if (is.rv(n)) {
    sim.n <- sims(n, n.sims=n.sims)
    sim <- .nsimmat(sim.n, sim)
  } else {
    if (length(n)>1) n <- length(n)
    if (n>1) sim <- rep(sim, n)
  }
  sim
}

.rvgen <- function(fun, n, ...) # NOEXPORT
{
  n <- n[1] ## DEBUG: (Random) dimensions not yet supported!!!
  if (!is.character(fun)) fun <- deparse(substitute(fun))
  paramlist <- list(...)
  arglist <- lapply(paramlist, .parammat, n)
  arglist$n <- max(sapply(arglist,length))
  ow <- options("warn")
  options(warn=-1)
    m <- do.call(fun, args=arglist) # will complain about NAs, but ignore.
  options(ow) # reset warnings
  n.sims <- rvnsims()
  sims <- matrix(m, nrow=n.sims)
  sims[is.na(sims)] <- NA # NaNs to NAs.
  rvsims(sims)
}


# ========================================================================
# rv.const  -  make a constant rv
# ========================================================================
#

rvconst <- function(n=1,x)
{
  if (missing(x)) { x=n; n=1 }
  l.x <- length(x)
  if (l.x==0)
    return(NULL)
  v <- rv(l.x)
  v[1:l.x] <- x # recycle
  v
}

# ========================================================================
# rvnorm  -  Generate independent normal random variables
# ========================================================================

rvnorm <- function (n=1, mean=0, sd=1, var=NULL, precision)
{
  if (!missing(precision)) {
    if (!is.null(dim(precision))) {
      var <- solve(precision) # matrix inverse
    } else {
      sd <- 1/sqrt(precision)
    }
  }
  if (!is.null(var)) {
    if (!is.null(dim(var))) {
       return(.rvmvnorm(n=n, mean=mean, Sigma=var))
    }
    sd <- sqrt(var)
  }
  .rvgen('rnorm', n=n, mean=mean, sd=sd)
}

# mvrnorm


.mvrnorm <- function (n = 1, mu, Sigma, tol = 1e-06, empirical = FALSE) # NOEXPORT
{ 
    p <- length(mu)
    if (!all(dim(Sigma) == c(p, p)))
        stop("incompatible arguments")
    eS <- eigen(Sigma, sym = TRUE, EISPACK=TRUE)
    if (any(is.na(eS$vectors))) {
      stop("Some eigenvectors are missing...")
    }
    ev <- eS$values
    if (!all(ev >= -tol * abs(ev[1])))
        stop("Sigma is not positive definite")
    X <- matrix(rnorm(p * n), n)
    if (empirical) {
        X <- scale(X, TRUE, FALSE)
        X <- X %*% svd(X, nu = 0)$v
        X <- scale(X, FALSE, TRUE)
    }
    X <- mu + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
    nm <- names(mu)
    if (is.null(nm) && !is.null(dn <- dimnames(Sigma)))
        nm <- dn[[1]]
    dimnames(X) <- list(nm, NULL)
    if (n == 1)
        drop(X)
    else t(X)
}


# ========================================================================
# .rvmvnorm  -  Generate multivariate normal random variables.
# ========================================================================
#
# TODO: 
#  - preserve dimensionality when arrays of n,mean,sd are given.
#  - dim(n) , dim(mean), dim(sd), in this order? at least dim(n) 
#    should be the most important.
#  - parameter cumulative = FALSE
#  - what to do with n rv? conditional d'n given n?
#  - what to do with the parameter 'n'? Should repeat, yes, how to implement?
#

.rvmvnorm <- function(n, mean, Sigma)
{
  n.sims <- rvnsims()
  dim.mean <- dim(mean)
  if (is.list(Sigma))
    Sigma <- .expand.as.matrix(Sigma)
  if (is.null(dim(Sigma)))
    return(rvnorm(n=n, mean=mean, var=Sigma)) ## Sigma <- diag(Sigma,nrow=length(Sigma))
  else if (!is.numeric(Sigma))
    stop('Invalid (non-numeric) covariance matrix Sigma')
  if (nrow(Sigma) != ncol(Sigma))
    stop('Invalid (nonsquare) covariance matrix Sigma')
  if (length(mean)>nrow(Sigma)) {
    # ???
    n.Sigmas <- (length(mean) %/% nrow(Sigma)) + (length(mean) %% nrow(Sigma)!=0)
    Sigma <- .expand.as.matrix(rep(rv(Sigma),n.Sigmas))
  }
  mean <- rep(mean, length.out=nrow(Sigma))
  if (missing(n)) {
    n <- length(mean)
  } else if (n!=length(mean)) {
    ## mean <- rep(mean,length.out=n)
  }
  if (length(mean)!=nrow(Sigma)) {
    stop("DEBUG: length of mean vector != nrow of Sigma. This is not yet implemented.")
  }
  X <- .mvrnorm(n=n.sims, mu=mean, Sigma=Sigma)
  r <- rvsims(X)
  nm <- names(mean)
  if ( is.null(nm) && !is.null(dn <- dimnames(sd)) )
    nm <- dn[[1]]
  if (!is.null(dim.mean) && length(r)!=prod(dim.mean))
    warning("Cannot set the dimension ",dim.mean[1],'x',dim.mean[2])
  else
    dim(r) <- dim.mean
  if (is.null(dim(r)))
    names <- nm
  else
    dimnames(r) <- list(nm, NULL)
  r
}



# ========================================================================
# rvmvt  -  multivariate t random variables 
# ========================================================================
#

.rvmvt <- function (n=1, Sigma, df=1)
{
  x <- sqrt(rvchisq(n=n, df=df)/df)
 # But will this work? x is of length n, but the returned thing is of length n x nrow(Sigma)!
  .rvmvnorm(n=n, mean=0, Sigma=Sigma)/x
}

# ========================================================================
# rvt  -  t random variables
# ========================================================================
#

rvt <- function (n=1, mu=0, scale=1, df, Sigma)
{
  if (!is.missing(Sigma)) {
    t <- .rvmvt(n=n, Sigma=Sigma, df=df)
  } else {
    t <- .rvgen('rt', n=n, df=df)
    if (scale != 1) t <- scale*t
  }
  if (mu!=0) t <- mu + t
  t
}

# ========================================================================
# rvcauchy  -  Cauchy random variables
# ========================================================================
#

rvcauchy <- function (n=1, location=0, scale=1)
{
  .rvgen('rcauchy', n=n, location=location, scale=scale)
}


# ========================================================================
# rvchisq  -  chi-square random variables
# ========================================================================
#

rvchisq <- function (n=1, df, ncp = 0) 
{
  if (missing(ncp)) 
    .rvgen('rchisq', n=n, df=df)
  else .rvgen('rchisq', n=n, df=df, ncp=ncp)
}


# ========================================================================
# rvinvchisq  -  inverse chi-square random variables
# ========================================================================
#

rvinvchisq <- function (n=1, df, scale=1)
{
  df*scale/rvchisq(n=n, df=df)
}


# ========================================================================
# rvgamma  -  gamma random variables
# ========================================================================


rvgamma <- function (n=1, shape, rate = 1, scale = 1/rate) 
{
  if (any(sims(shape) <= 0)) 
      stop("shape must be strictly positive")
  .rvgen('rgamma', n=n, shape=shape, scale=scale)
}


# ========================================================================
# rvexp  -  Exponential random variables
# ========================================================================

rvexp <- function (n=1, rate = 1)
{
  .rvgen('rexp', n=n, rate=rate)
}


# ========================================================================
# rvpois  -  Poisson random variables
# ========================================================================

rvpois <- function (n=1, lambda)
{
  .rvgen('rpois', n=n, lambda=lambda)
}

# ========================================================================
# rvunif  -  uniform random variables
# ========================================================================

rvunif <- function (n=1, min=0, max=1)
{
  .rvgen('runif', n=n, min=min, max=max)
}


# ========================================================================
# rvbern  -  Bernoulli random variables
# ========================================================================

rvbern <- function (n=1, prob)
{
  .rvgen('rbinom', n=n, size=1, prob=prob)
}


# ========================================================================
# rvbinom  -  binomial rvs
# ========================================================================

rvbinom <- function (n=1, size, prob)
{
  .rvgen('rbinom', n=n, size=size, prob=prob)
}


# ========================================================================
# rvbeta  -  beta rvs
# ========================================================================

rvbeta <- function (n=1, shape1, shape2)
{
  .rvgen('rbeta', n=n, shape1=shape1, shape2=shape2)
}


# ========================================================================
# rvdirichlet  -  dirichlet rvs
# ========================================================================

rvdirichlet <- function (n=1, alpha)
{
  # alpha is a vector
  g <- rvgamma(n=n, shape=alpha, scale=alpha )
  g/sum(g)
}

# ========================================================================
# rvboot  -  empirical (bootstrap) distribution
# ========================================================================

rvboot <- function (data)
{
  n.sims <- rvnsims()
  n <- n.sims*length(data)
  s <- matrix(sample(data, size=n, replace=TRUE), nrow=n.sims)
  r <- rvsims(s)
  dim(r) <- dim(data)
  r
}

# ========================================================================
# rvpermut  -  permutation distribution
# ========================================================================

rvpermut <- function (data)
{
  n.sims <- rvnsims()
  if (length(data)<=1) return(data)
  s <- t(sapply(1:n.sims, function(i) sample(data)))
  r <- rvsims(s)
  dim(r) <- dim(data)
  r
}






# Rewrite: simapply, sims.rv

# ========================================================================
# function  -  short description
# ========================================================================
#

#rep.rv <- function (x, times, ...) # NOEXPORT
#{
#  # We can use the default method for lists
#  x.class <- class(x)
#  if (is.rv(times)) {
#    stop('DEBUG rep.rv: is.rv(times) not yet implemented')
#  } else {
#    x.rep <- rep.default(x, times, ...)
#  }
#  class(x.rep) <- x.class
#  x.rep
#}


# ========================================================================
# rowapply  -  optimized apply for rows of rvs
# ========================================================================
#

rowapply <- function (X, FUN, ...) # 
{
  # Check whether FUN=... has been forgotten or not
  FUN <- match.fun(FUN)
  d <- dim(X)
  dl <- length(d)
  if (dl == 0) 
    stop("dim(X) must have a positive length")
  v <- rep(NA,nrow(X))
  for (i in 1:nrow(X))
    v[i] <- FUN(X[i,], ...)
  if (mode(v)=="list")
    return(rv(v))
  else
    return(v)
}



# rv-stats.R
#  cov
#  var
#  rvcov
#  quantile
#  median
#

# ========================================================================
# cov  -  short description
# ========================================================================

.cov <- getFromNamespace('cov', 'stats')

cov <- function(x, y=NULL, ...)
{
  if (anyisrv(x,y)) {
    if (!is.matrix(x)) {
      if (is.null(y)) 
        stop("supply both x and y or a matrix-like x")
        x <- as.vector(x)
    }
    return(matsimapply(x, y, .cov, ...))
  }
  .cov(x, y, ...)
}

.rvRegisterFunctionSwitch(cov, 'cov', 'stats')


## var(x) : distribution of the sample variance

# ========================================================================
# function  -  short description
# ========================================================================
#

.var <- getFromNamespace('var', 'stats')

var <- function(x, ...)
{
  if (is.rv(x)) return(simapply(x, var, ...))
  .var(x, ...)
}

.rvRegisterFunctionSwitch(var, 'var', 'stats')

# ========================================================================
# rvcov  -  covariance matrix
# ========================================================================
#

rvcov <- function(x, y=NULL, ...)
{
  if (is.null(y))
    .cov(sims(x), y, ...)
  else
    .cov(sims(x), sims(y), ...)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# QUANTILE, MEDIAN, 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ========================================================================
# quantiles of a random vector
# ========================================================================
#

quantile.rv <- function(x, ...)
{
  simapply(x, quantile, ...)
}

# ========================================================================
# median.rv - median of a random vector
# ========================================================================
#

.median <- getFromNamespace('median', 'stats')

median <- function(x, na.rm=FALSE)
{
  if (is.rv(x)) return(simapply(x, median, na.rm=na.rm))
  .median(x, na.rm=na.rm)
}

.rvRegisterFunctionSwitch(median, 'median', 'stats')



.packageName <- "rv"
#
# rv-Bayes.R
# 

# ========================================================================
# function  -  short description
# ========================================================================
# Name:        
# Description: 
# Parameters:  
# Required:    none
# History:     2004-06-  : 
#

postsim <- function(fit)
{
  UseMethod('postsim')
} 

# Old code:
#    for (sim in 1:nsim){
#      sigma[sim] <- sigma.hat*sqrt((n-k)/rchisq(1,n-k))
#      beta[sim,] <- mvrnorm (1, beta.hat, V.beta*sigma[sim]^2)
#    }
#

# ========================================================================
# postsim.lm
# ========================================================================

postsim.lm <- function(fit)
{
  # Modified version of 'sim' (from Andrew Gelman's library (gelman@stat.columbia.edu))
  #
  summ <- summary (fit)
  sigma.hat <- summ$sigma
  beta.hat <- summ$coef[,1]
  V.beta <- summ$cov.unscaled
  k <- summ$df[1]
  n <- k + summ$df[2]
  sigma <- sigma.hat*sqrt((n-k)/rvchisq(1,n-k))
  beta.0 <- as.vector(beta.hat) + sigma * rvnorm(mean=0, var=V.beta)
  names(beta.0) <- names(beta.hat)
  c(sigma=sigma, beta=beta.0)
}

# ========================================================================
# postsim.glm
# ========================================================================

postsim.glm <- function(fit)
{
  # Modified version of 'sim' (from Andrew Gelman's library (gelman@stat.columbia.edu))
  #
  summ <- summary (fit, correlation=T)
  beta.hat <- summ$coef[,1]
  sd.beta <- summ$coef[,2]
  corr.beta <- summ$corr
  k <- summ$df[1]
  n <- k + summ$df[2]
  V.beta <- corr.beta * array(sd.beta,c(k,k)) * t(array(sd.beta,c(k,k)))
  # dimnames(beta) <- list (NULL, names(beta.hat))
  rvnorm(mean=beta.hat, var=V.beta)
}

# ========================================================================
# points.rv  -  plot points and uncertainty intervals
# ========================================================================
# Name:        
# Description: 
# Parameters:  
# Required:    none
# History:     2004-06-  : 


points.rv <- function (x, y = NULL, what=NULL, rvcol=c('green','red'), ...)
{
  if (missing(y)) {
    y <- x
    x <- seq(along=y)
  }
  if (is.rv(x))
  {
    if (length(x)==1 & length(x)==length(y)) {
      points(sims(x),sims(y), pch=20, ...)
    } else {
      stop("multivariate rv-rv point plot Not yet implemented")
    }
    return(invisible(NULL))
  }
  x0 <- x
  y0 <- rvsimapply(y, quantile, c(0.025,0.975,0.25,0.75,0.5))
  col.seg <- rvcol[2]
  col.dot <- rvcol[1]
  segments(x0, y0[1,], x0, y0[2,], col=col.seg, lty='dotted', lwd=1.5, ...)
  segments(x0, y0[3,], x0, y0[4,], col=col.dot, lwd=1.73, ...)
  if (any(y.n.r <- !is.random(y))) {
    points(x0[y.n.r],E(y[y.n.r]), ...)
  }
  #segments(x0, pmin(E(y),y0[5,]), x0, pmax(E(y),y0[5,]), col='black', lwd=3, ...)
  #segments(x0, y0[5,], x0, E(y), col='black', ...)
  #points(x0, y0[5,], pch=19, ...)
  #points(x0, E(y), pch=18,  ...)
  invisible(NULL)
}

# ========================================================================
# abline.rv  -  a + bx
# ========================================================================

abline.rv <- function(x, sims=1:20, ...)
{
  rows <- as.integer(sims)
  if (any(rows<1)) stop("sims should be a vector of indices")
  s <- sims(x)[rows,,drop=FALSE]
  for (i in 1:nrow(s))
    abline(s[i,], ...)
}

# ========================================================================
# lines.rv  -  plot some random lines
# ========================================================================
# btw, "rvlines" does not make sense

lines.rv <- function(x, simrows=1:3, lcolors=c("blue", "red", "green", "orange", "yellow"), ...)
{
  rows <- as.integer(simrows)
  if (any(rows<1)) stop("bad simrows")
  s <- sims(x)[rows,,drop=FALSE]
  lc <- lcolors
  ll <- length(lcolors)
  for (i in 1:nrow(s)) {
    lines(s[i,], col=lc[1+((i-1) %% ll)], ...)
  }
}

# ========================================================================
# plot.rv  -  plot for rvs
# ========================================================================
# Not yet done: must modify.

plot.rv <- function(x, y, what=c("95%","50%","mean","median"), ylim=range(sims(y)), xlim=range(sims(x)), rvcol=c("green", "red"), ...)
{
  if (missing(y)) {
    y <- x
    x <- seq(along=y)
  }
  x. <- xlim
  y. <- ylim
  plot(x., y., type="n", ylim=ylim, xlim=xlim, ...)
  points.rv(x, y, what=what, rvcol=rvcol)
  invisible(NULL)
}

# Not yet debugged:

.NEWplot.rv <- function(x, y, what=c("95%ci","50%ci","mean","median"), ylab,
   ylim=c(min(sims(x)),max(sims(x))), ...)
{
  if (!missing(y) && !is.numeric(x))
    stop("Two-rv scatterplot Not yet implemented")
  if (missing(y)) {
    y <- x
    x <- seq(along=x)
  }
  if (missing(ylab))
    ylab <- ''
  y.means <- E(y)
  plot(x, y.means, type="n", ylab=ylab, ylim=ylim, ...)
  x0 <- x
  y0 <- rvsimapply(y, quantile, c(0.025,0.975,0.25,0.75,0.5))
  segments(x0, y0[1,], x0, y0[2,], col='red', ...)
  segments(x0, y0[3,], x0, y0[4,], col='green', ...)
  points(x0, y0[5,], pch=19, ...)
  points(x0, y.means, ...)
  ## abline(h=mean(E(x)), lty="dotted")
  invisible(NULL)
}

# ========================================================================
# rvplot  -  plot simulations of two random vectors
# ========================================================================
# Name:        
# Description: 
# Parameters:  
# Required:    none
# History:     2004-06-  : 
#

rvplot <- function(x,y,...)
{
  if (missing(y))
    plot(sims(x[1]),sims(x[2]),...)
  else
    plot(sims(x),sims(y),...)
}

# ========================================================================
# hist.rv  -  histogram, adapted for rv's
# ========================================================================
# Name:        
# Description: 
# Parameters:  
# Required:    none
# History:     2004-06-  : 
#

hist.rv <- function(x, grid=c(4,5), xlim=x.range, main=paste(xname,"simulation"), ...)
{
  par(mfrow=grid)
  #  l <- length(x)
  #  par(mfrow=c(l %% 3 + l %/% 3, 3))
  xname <- deparse(substitute(x))
  grid.n <- grid[1]*grid[2]
  s <- sims(x)
  x.range <- c(min(s),max(s))
  if (grid.n<1)
    stop("Bad grid")
  s <- s[1:grid.n,]
  for (i in 1:grid.n)
  {
    # truehist(s[i,], xlim=xlim, main, ...)
    hist(s[i,], xlim=xlim, main=main, ...)
  }
}

# ========================================================================
# rvhist  -  plot histograms of the simulations of the random components
# ========================================================================
# #

rvhist <- function(x, main=NULL, ...)
{
  if (!is.null(dim(x)))
    par(mfcol=dim(x))
  mfcol <- par('mfcol')
  n <- prod(mfcol)
  make.main <- is.null(main)
  lab <- deparse(substitute(x))
  for (i in 1:n) {
    s <- sims(x[i])
    if (make.main) {
#### Todo: if x is multidimensional, set the correct labels!!!
      x.name <- paste(lab,'[',i,']', sep='')
      main=paste('Histogram of',x.name)
      hist(s, prob=TRUE, main=main, xlab=lab, ...)
    } else {
      hist(s, prob=TRUE, xlab=lab, ...)
    }
  }
}


# ================================================================================
# hrvplot - horizontal r.v. (interval) plot
# ================================================================================

hrvplot <- function (x, labels=NULL, pos=NULL, output.file=NULL)
{
  UseMethod('hrvplot')
}

hrvplot.rv <- function(x, labels=NULL, ...)
{
  mc.array.ok <- !any(is.na(rvnchains(x)))
  probs <- c(0.025,0.25,0.5,0.75,0.975)
  debugged <- FALSE
  if (debugged && mc.array.ok) {
## Still something fishy here:
    sim <- sims(x, mc.array=TRUE)
    dim.sim <- dim(sim)
    s <- list()
    for (i in seq(length=dim.sim[2])) {
      s[[i]] <- t(apply(sim[,i,],2,quantile,probs))
      # we get an n*5 matrix
    }
    d <- dim(s[[1]])
    summ <- array(unlist(s),dim=c(d[1],dim.sim[2],d[2]))
  } else {
    summ <- t(rvquantile(x,probs=probs))
  }
  if (is.null(labels)) {
    if (is.null(labels <- names(x))) {
      labels <- paste(deparse(substitute(x)), .dim.index(x), sep='')
    }
  }
  hrvplot(summ, labels=labels, ...)
}

hrvplot.matrix <- function(x, labels, pos=NULL, output.file=NULL)
{
  summ <- x

  old.par <- par(mar=c(0,0,0,0))

  print.to.file <- ! is.null(output.file)

  if (print.to.file) {
    postscript (output.file, horizontal=TRUE, width=9)
  }

  no.labels <- missing(labels) | is.null(labels)

  cex.scale <- 0.9
  thin.style <- 0.5
  medium.style <- thin.style*2
  thick.style <- thin.style*4

  dim.sum <- dim(summ)
  sum.dim <- length(dim.sum)

  n.rows <- dim.sum[1]

  if (sum.dim==2) {
    nc <- dim.sum[2]
    lines.per.row <- 1
  } else if (sum.dim==3) {
    nc <- dim.sum[3]
    lines.per.row <- dim.sum[2]
  }

  if (is.null(pos)) {
    pos <- seq(length=n.rows)
  }

  if (no.labels) {
    labels <- as.character(seq(length=n.rows))
  }

  pos <- -pos
  bottom <- min(pos)-1

  rng <- range(summ)
  min <- rng[1]
  max <- rng[2]
  p.rng <- pretty(rng)
  width <- max-min

  left.margin <- 0.25*width

  x.left   <- min - left.margin
  x.right  <- max

  y.top    <- 2
  y.bottom <- bottom-2

  x.plotrange <- c(x.left,x.right)
  y.plotrange <- c(y.top, y.bottom)

  plot (x.plotrange, y.plotrange,
    xlab="", ylab="",
    xaxt="n", yaxt="n",
    type="n", bty="n"
  )

  y.topline <- y.top
  y.bottomline <- y.bottom+2

  text (x.left, pos, labels, adj=0, cex=cex.scale*1.1)

  lines (c(0,0), c(0,y.bottomline), lwd=thin.style) # Zero line
  lines (c(min,max), rep(0,2)) # Top line
  lines (c(min,max), rep(y.bottomline,2)) # Bottom line

  show.rng <- p.rng[ p.rng >= min & p.rng <= max ]
  tick.height <- 0.25
  x0 <- show.rng
  x1 <- x0
  y0 <- rep(0,times=length(x0))
  y1 <- rep(-tick.height,times=length(y0))
  text (x0, 1, show.rng, cex=cex.scale*1.1)
  text (x0, bottom-1, show.rng, cex=cex.scale*1.1)
  segments(x0=x0, y0=y0, x1=x1, y1=y1)
  segments(x0=x0, y0=bottom+y0, x1=x1, y1=bottom-y1)

  thin.col <- c('red','blue','yellow','orange')
  thick.col <- c('green','purple','orange', 'red')
  point.col <- c('black','dark red', 'dark red', 'dark red')
  medium.col <- thick.col

  long.line <- if (nc>=4) thin.style else medium.style
  short.line <- thick.style

  for (i in seq(length=lines.per.row)) {
    y.pos <- pos + (i-1)/8
    color <- 1 + (i-1) %% length(thin.col)
    if (sum.dim==3) {
      s <- summ[,i,]
    } else {
      s <- summ
    }
    segments(x0=s[,1], y0=y.pos, 
      x1=s[,nc], y1=y.pos,
      lwd=long.line,
      col=thin.col[color]
    )
    if (nc>=3) {
      segments(x0=s[,2], y0=y.pos, 
        x1=s[,nc-1], y1=y.pos,
        lwd=short.line,
        col=thick.col[color]
      )
    }
    if (nc %% 2 == 1) {
      midpoint <- 1 + ((nc-1) %/% 2)
      points (s[,midpoint], y.pos,
        pch=20,
        cex=cex.scale*1.25,
        col=point.col[color]
      )
    }
  }
  par(old.par)
  if (print.to.file) dev.off()
  invisible(NULL)
}


# ================================================================================
# rvplotseries  -  plot lines 
# ================================================================================
# we want to have the y range bounded by decent values

#rvplotseries <- function (x, simrows=1, lcolors=c("blue", "red", "green", "orange", "yellow"), ...)
#{
#  plot(c(1,length(x)), rvrange(x), type="n")
#  lines(x, simrows=simrows, lcolors=lcolors, ...)
#}


# ========================================================================
# plot 
# ========================================================================

# .plot <- getFromNamespace("plot", "graphics")

plot <- function (x, y, ...) 
{
    if (!missing(y) && is.rv(y) && !is.rv(x)) return(plot.rv(x, y, ...))
    if (is.null(attr(x, "class")) && is.function(x)) {
        nms <- names(list(...))
        if (missing(y)) 
            y <- {
                if (!"from" %in% nms) 
                  0
                else if (!"to" %in% nms) 
                  1
                else if (!"xlim" %in% nms) 
                  NULL
            }
        if ("ylab" %in% nms) 
            plot.function(x, y, ...)
        else plot.function(x, y, ylab = paste(deparse(substitute(x)), 
            "(x)"), ...)
    }
    else UseMethod("plot")
}


.rvRegisterFunctionSwitch(plot, "plot", "graphics")


# rv-graph.R


# ========================================================================
# mcmcarray  -  coerce rv into a LxMxN mcmc array
# ========================================================================
#  

mcmcarray <- function(x)
{
  has.chain <- !is.na(rvnchains(x))
  if (!any(has.chain))
    return(NULL)
  a <- list()
  ord <- rvsim.attr(x,'order')
  for (i in seq(along=ord)) {
    order.i <- ord[[i]]
    if (length(order.i)==1 && is.na(order.i)) {
      next
    }
    s <- sims(x[i])
    s <- s[order.i]
    dim(s) <- dim(order.i)
    ord[[i]] <- s
  }
  ord
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ATTRIBUTE ACCESS
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# ========================================================================
# rvsim.attr  -  return rvsim.attr attribute for each component
# ========================================================================
# returns a list.

rvattr <- function(x, attrib)
{
  a <- lapply(x, "attr", "rvsim")
  a <- lapply(a, "[[", attrib)
  nulls <- sapply(a, is.null)
  a[nulls] <- NA
  a
}

# rvattr(x, 'name') <- vector of values for each component or x; or,
# rvattr(x) <- list of attributes: list(attributename=vector of values, attributename2=...)
# e.g. list(Rhat=c(1.01,1.03,1.23), n.eff=c(200,100,16)).
#



"rvattr<-" <- function(x, attrib=NULL, value)
{
  if (is.null(attrib)) {
    for (a in names(value)) {
      rvattr(x, a) <- value[[a]]
    }
  } else {
    for (i in seq(along=x)) {
      a <- attr(x[[i]], "rvsim")
      if (i<=length(value)) {
        a[attrib] <- value[[i]]
      } else {
        a[attrib] <- NA
      }
      attr(x[[i]], "rvsim") <- a
    }
  }
  x
}


# ========================================================================
# rvRhat; rvneff; rvnchains - convenience functions for some attributes
# ========================================================================
#

rvRhat <- function(x)
{
  unlist(rvattr(x,'Rhat'))
}

rvneff <- function(x)
{
  unlist(rvattr(x,'n.eff'))
}

rvnchains <- function(x)
{
  unlist(rvattr(x,'n.chains'))
}


# miscellaneous functions

# ========================================================================================
# rvci - confidence/credible interval
# ========================================================================================

rvci <- function(obj, interval=0.95, one.sided=FALSE, left=TRUE)
{
  if (one.sided) {
    q <- if (left) interval else 1-interval
    ci <- t(rvquantile(obj, q))
  } else {
    lower <- (1-interval)/2
    upper <- lower + interval
    ci <- t(rvquantile(obj, c(lower,upper)))
  }
  ci
}

# ========================================================================================
# rvinci - is a given value in the confidence interval?
# ========================================================================================

rvinci <- function(obj, x, interval=0.95, one.sided=FALSE, left=FALSE)
{
  ci <- rvci(obj, interval=interval, one.sided=one.sided, left=left)
  (x>=ci[,1] && x<=ci[,2])
}

# ========================================================================================
# doc - show a pdf file from the 
# ========================================================================================

doc <- function (topic, package = NULL, lib.loc = NULL) 
{
    #
    # This is just a modified version of the function "vignette"
    # NOTHING fancy here
    # usage: doc(topic, package)
    # examples: doc('rv'); doc('R2WinBUGS')
    #
    if (is.null(package)) 
        package <- .packages(all.available = TRUE, lib.loc)
    paths <- .find.package(package, lib.loc=lib.loc)
    paths <- paths[tools::file_test("-d", file.path(paths, "doc"))]
    exts <- "pdf"
    vignettes <- lapply(paths, function(dir) {
        tools::list_files_with_exts(file.path(dir, "doc"), exts, full.names=TRUE)
    })
    if (!missing(topic)) {
        topic <- topic[1]
        vignettes <- as.character(unlist(vignettes))
        vidx <- (tools::file_path_sans_ext(basename(vignettes)) == 
            topic)
        if (any(vidx)) {
            pdf <- sub("\\.[[:alpha:]]+$", ".pdf", vignettes)
            pidx <- tools::file_test("-f", pdf)
            ok <- vidx & pidx
            if (any(ok)) {
                idx <- min(which(ok))
                if (sum(ok) > 1) {
                  warning(gettextf("vignette '%d' found more than once,\nusing the one found in '%s'", 
                    topic, dirname(pdf[idx])), call. = FALSE, 
                    domain = NA)
                }
                z <- list(file = vignettes[idx], pdf = pdf[idx])
            }
            else {
                z <- list(file = vignettes[vidx][1], pdf = character(0))
            }
            z$topic <- topic
            class(z) <- "vignette"
            return(z)
        }
        else warning(gettextf("vignette '%s' *not* found", topic), 
            call. = FALSE, domain = NA)
    }
    if (missing(topic)) {
        topic <- package
        vDB <- matrix(character(0), nr = 0, nc = 4)
        colnames(vDB) <- c("Dir", "File", "Title", "PDF")
        for (db in vignettes[sapply(vignettes, length) > 0]) {
            dir <- dirname(dirname(db[1]))
            entries <- NULL
            if (file.exists(INDEX <- file.path(dir, "Meta", "vignette.rds"))) 
                entries <- .readRDS(INDEX)
            if (NROW(entries) > 0) 
                vDB <- rbind(vDB, cbind(Dir = I(dir), entries[c("File", 
                  "Title", "PDF")]))
        }
        title <- if (NROW(vDB) > 0) {
            paste(vDB[, "Title"], paste(rep.int("(source", NROW(vDB)), 
                ifelse(vDB[, "PDF"] != "", ", pdf", ""), ")", 
                sep = ""))
        }
        else character()
        db <- cbind(Package = basename(vDB[, "Dir"]), LibPath = dirname(vDB[, 
            "Dir"]), Item = tools::file_path_sans_ext(basename(vDB[, 
            "File"])), Title = title)
        y <- list(type = "vignette", title = "Vignettes", header = NULL, 
            results = db, footer = NULL)
        class(y) <- "packageIQR"
        return(y)
    }
}


# ========================================================================
# print.rv  -  print summary of a rv on the console
# ========================================================================
# Name:        
# Description: 
# Parameters:  
# Required:    none
# History:     2004-06-  : 
#

print.rv <- function(x, ...)
{
  s <- summary(x)
  if (is.null(s))
    return(print(s))
  s <- as.data.frame(round(s,2))
  s <- cbind(s[,1:2],' '='(',s[,3:9],' '=')',s[-(1:9)])
  if (!is.null(.names <- names(x))) {
    s <- cbind(.names,s)
    names(s)[[1]] <- 'name'
  }
  print(s)
}


# ========================================================================
# function  -  short description
# ========================================================================
# Name:        
# Description: 
# Parameters:  
# Required:    none
# History:     2004-06-  : 
#

print.rvsim <- function(x, ...) # NOEXPORT
{
  print.default(x, ...)
}


# ========================================================================
# dsummary  -  summary line for summary.rv
# ========================================================================
#

dsummary <- function(v) # NOEXPORT
{
  r <- c(mean(v, na.rm=TRUE),sd(v,na.rm=TRUE),quantile(v,c(0.0,0.025,0.25,0.5,0.75,0.975,1),
    na.rm=TRUE), 100*mean(is.na(v)))
  names(r) <- c("mean","sd","Min","2.5%","25%","50%","75%","97.5%","Max","NA%")
  r
}

# ========================================================================
# summary  -  print a summary
# ========================================================================

summary.rv <- function(object, ...)
{
  x <- object
  x.len <- length(x)
  x.dim <- dim(x)
  if (is.null(x.dim))
    x.dim <- x.len
  if (x.len < 1) {
    return(NULL)
  }
  output <- t(rvsimapply(x, dsummary))
  output <- cbind(output,sims=nsims(x))
  const <- (nsims(x)==1)
  if (any(const)) output[const,'sd'] <- 0
  rnames <- dimnames(output)[[2]]
  if (all(output[,'NA%']==0)) {
    output <- output[, rnames!="NA%",drop=FALSE]
  }
  row.names <- rv:::.dim.index(x)
  dimnames(output)[[1]] <- row.names
  Rhat <- rvRhat(x)
  if (!is.null(Rhat)) {
    if (!all(is.na(Rhat)))  output <- cbind(output, Rhat, n.eff=rvneff(x))
  }
  chains <- rvnchains(x)
  if (!is.null(chains)) {
    if (!all(is.na(chains))) output <- cbind(output,chains)
  }
  output
}




as.rv.bugs <- function(x) # EXPORT
{
  # as.rv.bugs - transform a R2WinBUGS object into a random variable object
  # n.chains,
  # n.iter
  # n.burnin
  # n.thin
  # n.keep
  # n.sims
  # sims.array=sims.array[,,rank.long]
  # sims.list
  # sims.matrix=sims[,rank.long]
  # summary
  # mean=summary.mean
  # sd=summary.sd
  # median=summary.median
  # root.short
  # long.short
  # dimension.short
  # indexes.short
  # last.values
  #     sims <- c(bugs.sims(parameters.to.save, n.chains, n.iter, 
  #          n.burnin, n.thin, DIC), model.file = model.file, 
  #          is.DIC = DIC)
  #      class(sims) <- "bugs"
  #      return(sims)
  if (class(x) == 'bugs') {
    x <- x[[1]]
  } # else: plain bugs.R which does not have a class attribute
  #         and has slightly different structure
  if (is.null(x['sims.array']))
    stop("Argument does not contain a 'sims.array'")
  r <- rvsims(x$sims.array)
  rvattr(r, 'Rhat')  <- x$summary[,'Rhat']
  rvattr(r, 'n.eff') <- x$summary[,'n.eff']
  r
}


.onLoad <- function(libname,pkgname) # NOEXPORT
{
  if (is.null(rvnsims())) rvnsims(1000)
  setHook(packageEvent("rv", "detach"), 
    function (...) {
      .rvFunctionSwitch(attach=FALSE, verbose=FALSE)
    }
  )
  cat("Package rv loaded.\n")
}

.onAttach <- function(libname,pkgname) # NOEXPORT
{
  if (is.null(rvnsims())) rvnsims(200)
  .rvFunctionSwitch(attach=TRUE, verbose=FALSE)
}

