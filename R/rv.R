# ================================================================================
# rv - simulation-based random variable class in R
# Version 0.9*
# Updated 2006-08-13
# (c) 2004-2006 Jouni Kerman <jouni@kerman.com>
# ================================================================================
#
# NAME
#   rv-base0.r  -  simulation-based random variables in R, internal methods
# AUTHOR
#   Jouni Kerman, kerman@stat.columbia.edu
# SEE ALSO
#  other files in rv/R in this package
# 

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GLOBAL SETTINGS
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.RV <- list()

.rvRegisterFunctionSwitch <- function (FUN, fname, namespace, group=NULL)
{
  # 'group' is obsolete!
  #
  origFUN <- getFromNamespace(fname, namespace)
  ns.name <- paste(namespace, fname, sep=":::")
  .RV[[ns.name]] <<- list(R=origFUN, rv=FUN)
}


.rvFunctionSwitch <- function (attach=TRUE, name=NULL, verbose=FALSE)
{
  detach <- (!attach)
  if (is.null(name)) name <- names(.RV)
  RF <- .RV[name]
  full.names <- strsplit(name, ":::")
  for (i in seq(along=full.names)) {
    f <- RF[[i]]
    ns.name <- full.names[[i]]
    namespace <- ns.name[1]
    fun.name <- ns.name[2]
    FUN <- if (attach) f$rv else f$R
    if (verbose) {
       cat(name[i], " ")
       if (attach) cat("replaced by a function in package:rv\n")
         else cat("restored.\n")
     }
    assignInNamespace(fun.name, FUN, namespace)
  }
}



# ----------------
# end of rv-base0.R
# ----------------
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

# ========================================================================
# setnsims - get or set the default number of simulations (a global variable)
# ========================================================================

setnsims <- function(nsims)
{
  if (!missing(nsims)) {
    if (nsims[1]>0) {
      oldnsims <- rvpar("n.sims")
      rvpar(n.sims=ceiling(nsims[1]))
      return(oldnsims)
   } else
      stop('Invalid number of simulations ',nsims[1])
  }
  rvpar("n.sims")
}

getnsims <- function ()
{
  setnsims()
}

# ========================================================================
# rvnsims - get the number of simulations of the components of an rv
# ========================================================================

rvnsims <- function(x)
{
  sapply(x, length)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# AN OBJECT OF TYPE 'rv'
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
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
# Description: returns a random vector of given length, with NAs
# Parameters:  
# Required:    none
#

rv <- function(length=0)
{
  if (is.numeric(length)) {
    if (length>0) {
      n.sims <- setnsims()
      x <- lapply(1:length, function (x) rep(NA, n.sims))
    } else {
      x <- list()
    }
  } else {
    stop("length must be numeric")
  }
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

.rvsims.list <- function (x, permute = FALSE) 
{
  # Assume that all elements in the list have the same dimensions.
  # This may be modified later -- filling with NAs
  dx <- dim(x[[1]])
  s <- sapply(x, I)
  if (is.list(s)) {
    # some elements had different dimensions!
    stop("Simulation list was not consistent")
  }
  if (!is.null(dim(s))) s <- t(s)
  r <- rvsims(s, permute = permute)
  dim(r) <- dx
  r
}


rvsims <- function(sims, n.sims=setnsims(), permute=FALSE, save.order=FALSE)
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
    nsims <- setnsims()
    if (n.sims<n.sims.max) {
      omit <- (-1):(n.sims-n.sims.max)
      sims <- sims[omit,,drop=FALSE]
    } else if (n.sims>n.sims.max) {
      # warning("n.sims is larger than the available number of simulations")
      include <- rep(1:n.sims.max, length.out=n.sims)
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
  ##
  ## We should not set the dimension of the for constants.
  ## re: problems in max.rvsims: cbind(sims(1),sims(x)) does not work
  ## since sims(1) returns a 1x1 matrix but sims(x) may be L x n matrix.
  ## TODO: integrate better!
  ## Need to unclass first, to be consistent.
  ## DEBUG: ? why not just sims(as.rv(x, ...)) ? 
  as.vector(unclass(x))  # drop attributes
}

# ========================================================================
# .sims.as.list  -  split the simulations into a list
# ========================================================================

.sims.as.list <- function (x)
{
  # retain dimensions, and always return setnsims() simulations.
  if (is.null(x)) return(NULL)
  s <- sims(as.rv(x), n.sims=getnsims())
  s <- split(s, row(s)) ## faster than applying list to 1=rows.
  if (!is.null(d <- dim(x))) {
    dn <- dimnames(x)
    s <- lapply(s, function (x) { dim(x) <- d; dimnames(x) <- dn; x})
  }
  # The default names will be "1", "2", ... etc.
  # set names to NULL since this may interfere with mapply( ... ) in "[<-.rv"
  names(s) <- NULL
  s
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
    if (length(.order) == length(sim)) {
      sim <- sim[.order]
    }
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
 ###, as.list=FALSE)


sims.rv <- function(x, n.sims=NULL, dimensions=FALSE, sim.matrix=FALSE, mc.array=FALSE)
{
  if (length(x)<1) {
    return(NULL)
  }
  if (mc.array) {
    # A 3-way matrix
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
    m <- array(m, c(n.s, dim.x)) # multi-way array, 1st dimension is the dimension "rows"
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

# is rv numeric?

is.numeric.rv <- function (x)
{
  all(rvsimapply(x, is.numeric))
}

as.numeric.rv <- function (x, ...)
{
  simapply(x, as.numeric, ...)
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
  (rvnsims(x)>1)
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

#anyisrv <- function(...) # NOEXPORT
#{
#  'rv' %in% unlist(lapply(list(...),class))
#}

#anyisrv <- function (...) # NOEXPORT
#{
#  dots <- as.list(substitute(list(...)))[-1]
#  miss <- (sapply(dots, function (x) deparse(x)[1])=="")
#  dots[miss] <- 1
#  a <- lapply(dots, eval.parent)
#  any(sapply(a, is.rv))
#}


anyisrv <- function (...) # NOEXPORT
{
  any(sapply(list(...), is.rv))
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
  return(x)
}

as.rv.numeric <- function(x)
{
  if (is.rv(x)) return(x)
  r <- rvsims(matrix(as.vector(x), nrow=1))
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
  unlist.rv(lapply.rv(x, as.rv))
}

as.rv.matrix <- function(x)
{
  as.rv.numeric(x)
}

as.rv.default <- function(x, ...)
{
  if (is.null(x)) return(NULL)
  stop('Cannot coerce object of class "', class(x), '" to rv')
}

as.rv.xtabs <- function (x) 
{
  # NAME
  #  as.rv.xtabs
  #
  as.rv(x[])
}

# DEBUG:1 ?? The permutation of !is.null(dim(sims)) is not done yet
# DEBUG:2 (must take permutation of row indices and then permute.


# ========================================================================
# is.fuzzy  -  test whether components are logical but random
# ========================================================================

is.fuzzy <- function (x)
{
  UseMethod("is.fuzzy")
}

is.fuzzy.rv <- function (x)
{
  # NAME
  #  is.fuzzy - Is a Vector Component Logical But Random
  # 
  component.is.logical <- rvsimapply(as.rv(x), is.logical)
  component.prop <- rvmean(x)
  (component.is.logical & component.prop>0 & component.prop<1)
}

is.fuzzy.default <- function (x)
{
  FALSE
}

# ========================================================================
# as.logical, etc.
# ========================================================================



as.real.rv <- function(x, ...)
{
  simapply(x, as.real, ...)
}

as.double.rv <- function(x, ...)
{
  simapply(x, as.double, ...)
}

as.logical.rv <- function(x, ...)
{
  simapply(x, as.logical, ...)
}

as.integer.rv <- function (x, ...)
{
  simapply(x, as.integer, ...)
}

is.integer.rv <- function (x)
{
  ((is.rv(x) && all(rvsimapply(x, is.integer))) || .Primitive("is.integer")(x))
} 
.rvRegisterFunctionSwitch(is.integer.rv, "is.integer", 'base', group=2)

is.logical.rv <- function (x)
{
  ((is.rv(x) && all(rvsimapply(x, is.logical))) || .Primitive("is.logical")(x))
} 
.rvRegisterFunctionSwitch(is.logical.rv, "is.logical", 'base', group=2)

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
  (rvnsims(x)==1)
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

rvattach <- function (what, name='rvattach', overwrite=TRUE, impute=FALSE, ...)
{
  ## Avoid duplicates
  if (any(search()==name)) eval(substitute(detach(name),list(name=as.name(name))))
  if (!is.rv(what)) stop("Argument must be an rv")
  a <- splitbyname(what)
  if (is.null(a)) return(NULL)
  if (overwrite) {
    ls.GE <- ls(.GlobalEnv)
    for (j in 1:length(a)) { # Code from A. Gelman's attach.all.
      name.j <- names(a)[j]
      if (name.j %in% ls.GE) {
        if (impute) {
          v <- get(name.j, .GlobalEnv) # Value to be imputed to
          a[[name.j]] <- .impute.by.name(v, a[[name.j]])
        }
        remove (list=name.j, envir=.GlobalEnv)
      }
    }
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
#  [.rv
#  [<-.rv
#  NOT IMPLEMENTED: $.rv
#  NOT IMPLEMENTED: $<-.rv
#

# ========================================================================================
# rvcompatibility
# ========================================================================================

rvcompatibility <- function (level=NULL)
{
  old.level <- getOption("rvcompatibility")
  if (!is.null(level)) {
    if (level==0) {
      .rvFunctionSwitch(attach=TRUE)
    } else if (level==1) {
      .rvFunctionSwitch(attach=FALSE, name=c("base:::[<-", "base:::c"))
    } else {
      warning("Unknown level, ignored.")
      return(old.level)
    }
    options(rvcompatibility=level)
  }
  return(old.level)
}


# ========================================================================
# [.rv  -  rv component retrieval by index (NEW)
# ========================================================================
#


"[.rv" <- function (x, ..., drop = TRUE) 
{
    if (missing(x)) {
        if (drop) 
            return(drop(x))
        else return(x)
    }
    # Problem: cannot evaluate list(...) if we have missing values
    # (this happens often, e.g. "x[,1]")
    #
    dots <- as.list(substitute(list(...)))[-1]
    miss <- (sapply(dots, function(x) deparse(x)[1]) == "")
    # Kludge:
    if (any(miss)) {
      dots[miss] <- TRUE
      for (i in seq(length(dots))) {
          dots[[i]] <- eval.parent(dots[[i]])
      }
   } else {
     dots <- list(...)
   }
   isrv <- sapply(dots, is.rv)
   if (any(isrv)) {
      s <- .sims.as.list(x)
      dots.sims <- lapply(dots, .sims.as.list)
      list.of.sims <- do.call(mapply,
        args = c(FUN = .Primitive("["), list(s), dots.sims,
        list(MoreArgs = list(drop = drop)), SIMPLIFY = FALSE))
       v <- .rvsims.list(list.of.sims)
    }
    else {
        #call <- substitute(unclass(x)[..., drop = drop])
        #v <- eval.parent(call)
        v <- eval.parent(unclass(x)[..., drop = drop])
        if (length(v) >= 1) {
            v[sapply(v, is.null)] <- NA
        }
        else {
            v <- list()
        }
        class(v) <- class(rv())
    }
    v
}


# ========================================================================
# [<-.rv  -  rv component assignment by index (NEW)
# ========================================================================
# 2007-01-11



"[<-.rv" <- function (x, ..., value = NULL) 
{
    if (missing(x)) {
        if (drop) 
            return(drop(x))
        else return(x)
    }
    dots <- as.list(substitute(list(...)))[-1]
    miss <- (sapply(dots, function(x) deparse(x)[1]) == "")
    # Kludge:
    if (any(miss)) {
      dots[miss] <- TRUE
      for (i in seq(length(dots))) {
          ## This *might* not work in all cases, especially when called from within
          ## sapply or lapply.
          dots[[i]] <- eval.parent(dots[[i]])
      }
   } else {
     dots <- list(...)
   }
    isrv <- sapply(dots, is.rv)
    if (any(isrv)) {
        x.sims <- .sims.as.list(x)
        dots.sims <- lapply(dots, .sims.as.list)
        value.sims <- .sims.as.list(value)
        args <- c(FUN = .Primitive("[<-"), list(x.sims), dots.sims, 
            value = list(value.sims), SIMPLIFY = FALSE)
        list.of.sims <- do.call("mapply", args = args)
        v <- .rvsims.list(list.of.sims)
    }
    else {
        v <- eval.parent(.Primitive("[<-")(unclass(x), ..., value = value))
        v[sapply(v, is.null)] <- NA
        new.names <- names(v)
        dim(v) <- dim(x)
        dimnames(v) <- dimnames(x) 
        names(v) <- new.names
        class(v) <- class(rv())
    }
    v
}



# ========================================================================
# [<<-.rv  -  rv component assignment by index (NEW)
# ========================================================================
#

"[<<-.rv" <- function (x, ..., value = NULL) 
{
    if (missing(x)) {
        if (drop) 
            return(drop(x))
        else return(x)
    }
    dots <- as.list(substitute(list(...)))[-1]
    miss <- (sapply(dots, function(x) deparse(x)[1]) == "")
    # Kludge:
    if (any(miss)) {
      dots[miss] <- TRUE
      for (i in seq(length(dots))) {
          ## This *might* not work in all cases, especially when called from within
          ## sapply or lapply.
          dots[[i]] <- eval.parent(dots[[i]])
      }
   } else {
     dots <- list(...)
   }
    isrv <- sapply(dots, is.rv)
    if (any(isrv)) {
        x.sims <- .sims.as.list(x)
        dots.sims <- lapply(dots, .sims.as.list)
        value.sims <- .sims.as.list(value)
        args <- c(FUN = .Primitive("[<<-"), list(x.sims), dots.sims, 
            value = list(value.sims), SIMPLIFY = FALSE)
        list.of.sims <- do.call("mapply", args = args)
        v <- .rvsims.list(list.of.sims)
    }
    else {
        v <- eval.parent(.Primitive("[<<-")(unclass(x), ..., value = value))
        v[sapply(v, is.null)] <- NA
        new.names <- names(v)
        dim(v) <- dim(x)
        dimnames(v) <- dimnames(x) 
        names(v) <- new.names
        class(v) <- class(rv())
    }
    v
}



# ========================================================================
# [<- - new replacement method for "[" assignment
# ========================================================================


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

#"[<-.vector" <- .rvBracketAssignment # EXPORT "[<-"

.rvRegisterFunctionSwitch(.rvBracketAssignment, "[<-", 'base', group=1)



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

.rvConcatenate <- function(..., recursive=FALSE) {
  # x <- .Primitive("c")(..., recursive=recursive)
  x <- .c(..., recursive=recursive)
  if (is.list(x) && anyisrvsim(x)) # Was there a 'rvsim' in there?
    class(x) <- class(rv())
  x
}

### c.rv <- .rvConcatenate # EXPORT "c.rv"

####DEBUG#### "c" <- .rvConcatenate # EXPORT "c"

.c <- getFromNamespace("c", "base")

.rvRegisterFunctionSwitch(.rvConcatenate, 'c', 'base', group=2)

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
# updated 2007-01-11

mean.rv <- function(x, ...)
{
  rvsims(rowMeans(sims(x), ...))
}


# ========================================================================
# %*% - matrix product
# ========================================================================
# DEBUG: how to make this work with outer()?
#

.rvMatrixProduct <- function(a,b)
{
  if (!is.rv(a) && !is.rv(b)) return(.Primitive("%*%")(a,b))
  d <- dim(b)
  if (!is.rv(a) && (is.null(d)) || (length(d)==2 && d[2]==1)) {
    n.sims <- .Internal(max(rvnsims(a),rvnsims(b), na.rm=FALSE))
    bsim <- sims(as.rv(b), dimensions=TRUE, n.sims=n.sims)
    # Typical case: constant matrix times a rv vector
    AB <- t(.Primitive("%*%")(a,t(bsim)))
    rvsims(AB)
  } else {
    mapply.rv("crossprod", x=t(as.rv(a)), y=as.rv(b))
  }
}



"%*%" <- .rvMatrixProduct # EXPORT "%*%"

.rvRegisterFunctionSwitch(.rvMatrixProduct, "%*%", 'base', group=3)

# ========================================================================
# sort  -  generate order statistics
# ========================================================================
#

sort.rv <- function (x, ...) ## EXPORT sort.rv
{
  simapply(x, sort, ...)
}

# ========================================================================
# is.vector  -  coerce to vector
# ========================================================================

is.vector <- function(x, mode="any")
{
  if (is.rv(x)) return(is.null(dim(x)))
 .Internal(is.vector(x, mode))
}

.rvRegisterFunctionSwitch(is.vector, "is.vector", 'base', group=2)


# ========================================================================
# is.atomic.rv  -  rv is an 'atomic' type.
# ========================================================================

is.atomic <- function(x)
{
  (is.rv(x) || .Primitive("is.atomic")(x))
}
.rvRegisterFunctionSwitch(is.atomic, "is.atomic", 'base', group=2)



# ========================================================================
# min  -  minimum
# ========================================================================

min <- function(..., na.rm=FALSE) ## EXPORT min
{
  if (anyisrv(...)) {
    simapply(cbind.rv(...), .min, na.rm=na.rm)
  } else {
    .min(..., na.rm=na.rm)
  }
}

.min <- base:::min
.rvRegisterFunctionSwitch(min, "min", "base", group=0)

# ========================================================================
# max  -  maximum
# ========================================================================

max <- function(..., na.rm=FALSE) ## EXPORT max
{
  if (anyisrv(...)) {
    simapply(cbind.rv(...), .max, na.rm=na.rm)
  } else {
    .max(..., na.rm=na.rm)
  }
}

.max <- base:::max
.rvRegisterFunctionSwitch(max, "max", "base", group=0)


# ========================================================================
# pmin  -  parallel minimum
# ========================================================================

pmin <- function(..., na.rm=FALSE) ## EXPORT pmin
{
  if (anyisrv(...)) {
    a <- sims(cbind.rv(...), dimensions=TRUE)
    rvsims(t(apply(a, 1, function (m) apply(m, 1, .pmin))))
  } else {
    .pmin(..., na.rm=na.rm)
  }
}

.pmin <- base:::pmin
.rvRegisterFunctionSwitch(pmin, "pmin", "base", group=0)


# ========================================================================
# pmax.rv  -  parallel maximum
# ========================================================================

pmax <- function(..., na.rm=FALSE) ## EXPORT pmax.rv
{
  if (anyisrv(...)) {
    a <- sims(cbind.rv(...), dimensions=TRUE)
    rvsims(t(apply(a, 1, function (m) apply(m, 1, .pmax))))
  } else {
    .pmax(..., na.rm=na.rm)
  }
}

.pmax <- base:::pmax
.rvRegisterFunctionSwitch(pmax, "pmax", "base", group=0)

# ========================================================================
# solve.rv - solve linear systems
# ========================================================================

solve.rv <- function (a, b, ...)
{
  if (missing(b)) {
    mapply.rv("solve", a=a, ...)
  } else {
    mapply.rv("solve", a=a, b=b, ...)
  }
}


# ----------------
# end of rv-rbase.R
# ----------------

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# UTILITIES in the package rv
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  .show(x)
#  .matrix.ragged.array(x)
#  .expand.as.matrix(M)
#  .listByName(x)
#  .permut
#  .dim.index
#  .bracket
#  .impute.by.name
#  .setDimensionByName
#  .rowmax
#  .slice
#  .nodups
#
# ========================================================================
# .show - show a given variable and its value (useful for debugging purposes)
# ========================================================================

.myc <- function (x)
{
  paste("c(",paste(x,collapse=","),")", sep="")
}

.show <- function (x)
{
  cat(deparse(substitute(x)),"<-", .myc(x), "\n")
}


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
# .dimindex  -  make a list of indices of a matrix
# ========================================================================
#

.dimindex <- function(x) # NOEXPORT
{
  if (is.null(dim(x))) {
    ix <-.permut(length(x))
  } else {
    ix <- .permut(dim(x))
  }
  ixt <- paste('[',apply(ix, 1, paste, collapse=','),']',sep='')
  ixt
}

# ========================================================================
# .leftprepend  -  prepend spaces for short lines
# ========================================================================
#

.leftadjust <- function (x) # NOEXPORT
{
  m.x <- max(nchar(x))
  a <- (nchar(x)<m.x)
  if (any(a)) {
    d <- m.x - nchar(x)
    x[a] <- paste(' ', x[a], sep='') # only ONE single space
  }
  x
}

# ========================================================================
# .dim.index  -  make a list of indices of a matrix
# ========================================================================
#

.dim.index <- function(x, leftadjust=TRUE) # NOEXPORT
{
  if (is.null(dim(x)))
    ix <-.permut(length(x))
  else 
    ix <- .permut(dim(x))
  ixt <- paste('[',apply(ix, 1, paste, collapse=','),']',sep='')
  if (leftadjust) .leftadjust(ixt) else ixt
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
  v <- sapply(list(...), .bracket, row)
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
  # WRITE APPLY that works with RVs.
  #   This apply.rv should be more flexible, allowing mixing of scalars and matrices.
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
# .bracket.indices - extract indices from the brackets
# ========================================================================
# Returns NA's for missing indices

.bracket.indices <- function(x, default.bracket=FALSE) # already in rv 0.923
{
  if (length(x)<1 || !is.character(x)) return(NULL)
  no.brackets <- (regexpr("^(.*\\[(.*)\\].*)$", x)<1)
  y <- sub("^(.*\\[(.*)\\].*)$", "\\2", x)
  y <- sub(", *$", ",NA", y)
  if (any(no.brackets)) y[no.brackets] <- if (default.bracket) "1" else "0"
  y <- strsplit(y, " *, *")
  lapply(y, as.numeric)
}




# ========================================================================
# .indices - return the single-digit indices of components 
# ========================================================================
# x : a list of vectors of indices
# dim. : dimension or length of the target matrix or vector.

.indices <- function (x, dim. = NULL) 
{
    if (is.null(x)) 
        return(numeric(0))
    ld <- length(dim.)
    f <- function(x) {
        lx <- length(x)
        if (lx == 1) 
            return(x)
        if (lx != ld || any(x > dim.)) 
            return(NA)
        1 + sum(c(x - 1, 0) * c(1, cumprod(dim.)))
    }
    pos <- sapply(x, f)
    pos
}



# ========================================================================
# .impute.by.name - impute into a vector using the names attribute
# ========================================================================

.impute.by.name <- function (x, y)
{
  # impute x <- y using the names of y
  n.y <- names(y)
  name.of.y <- deparse(substitute(y))
  if (length(n.y)<1) stop("the names attribute is required for: ", name.of.y)
  bx <- .bracket.indices(n.y)
  ix <- .indices(bx, dim(x))
  if (length(ix)==1 && ix==0) {
    # Impute the whole vector.
    return(y)
  }
  na.ix <- is.na(ix)
  if (any(na.ix)) {
    warning("Couldn't figure out index ", n.y[which(na.ix)])
    ix <- ix[!na.ix]
    y  <- y[!na.ix]
    if (length(ix)<1) return(x)
  }
  x[ix] <- y
  if (is.null(dim(x))) {
    n.x <- names(x)
    n.x[ix] <- names(y)
    names(x) <- n.x
  }
  x
}



# ========================================================================
# .setDimensionByName - 
# ========================================================================

.setDimensionByName <- function (x)
{
  # 
  names.x <- names(x)
  bix <- .bracket.indices(names.x, default.bracket=TRUE)
  if (is.null(bix)) stop("Names of x MUST be set")
  b <- sapply(bix, prod)
  max.ix <- which(b==max(b))[1]
  maxdim <- bix[[max.ix]]
  if (prod(maxdim)<1) stop("Invalid dimension")
  a <- array(NA, maxdim)
  new.x <- .impute.by.name(a, x)
  names.new.x <- names(new.x)
  if (length(maxdim)>1) {
    dim(new.x) <- maxdim
  }
  names(new.x) <- names.new.x
  new.x
}



# ========================================================================
# .make.names - make names of the components
# ========================================================================

.make.names <- function(x, name=deparse(substitute(x)))
{
  if (is.null(x) || length(x)<1) return(x)
  name <- make.names(name)
  paste(name, .dim.index(x), sep="")
}

# ========================================================================
# .shortnames - names of the components, with brackets removed
# ========================================================================

.shortnames <- function(x)
{
  na <- names(x)
  if (is.null(na)) return(na)
  sapply(strsplit(na,'\\['),'[',1)
}

# ========================================================================
# make.fullnames - make eventually available!!
# ========================================================================

#.make.fullnames <- function(x)
#{
#  if (is.null(x) || length(x)<1) return(x)
#  name <- make.names(name)
#  paste(name, .dim.index(x), sep="")
#}







.list2array <- function (x, drop=FALSE)
{
  # name:
  #   list2array - coerce a list to an array, preserving dimnames
  # description:
  #   attempt to coerce the list to an array, preserving dimnames
  # arguments:
  #   x : (list) list to be coerced into one array
  # value:
  #   a vector if not all dimensions of list elements were equal;
  #   an array of dimension c(dim(x[[1]]),length(x)),
  #   if all dimensions of the elements were equal.
  # author:
  #   Jouni Kerman
  # history:
  #   2006-08-28
  #
  dm.1 <- length(x)
  if (dm.1<1) return(unlist(x))
  da <- lapply(x, function (x) if (is.null(dim(x))) length(x) else dim(x))
  dm.2 <- da[[1]]
  dn.1 <- names(x)
  dn.2 <- dimnames(x[[1]])
  # 
  a <- NULL
  for (i in seq(along=x)) {
    a <- c(a, x[[i]])
  }
  if (!all(sapply(da, function (d) all(d==dm.2)))) {
    return(a)
  }
  names(a) <- NULL
  dim(a) <- c(dm.2, dm.1)
  # If the list component does not have dimnames...
  if (is.null(dn.2)) {
    # ...and it doesn't have a dimension...
    if (is.null(dim(x[[1]]))) {
      # ...then use the names for dimensions...
      dn.2 <- names(x[[1]])
    } # ...or else just give up and use NULL for that dimension names.
    dn.a <- list(dn.2, dn.1)
  } else {
    dn.a <- c(dn.2, list(dn.1))
  }
  if (!is.null(names(dn.2))) {
    name.x <- deparse(substitute(x))
    names(dn.a) <- c(name.x, names(dn.2))
  }
  dimnames(a) <- dn.a
  if (drop) drop(a) else a
}



# .nodups

.nodups <- function (x, keep.latest=TRUE)
{
  # NAME
  # .nodups - remove list items with duplicate names (keep latest)
  # ARGUMENTS
  #   x : a list
  if (!keep.latest) stop("keep.latest=FALSE Not yet implemented")
  if (length(x)<2) return(x)
  a <- list()
  x <- rev(x) # use the fact that x[[na]] gets the *first* one
  for (na in names(x)) {
    a[[na]] <- x[[na]] # multiple assignments may occur, but that's ok.
  }
  a
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
# rvsample : componentwise samples
#
# Name:        
# Description: 
# Parameters:  
# Required:    none
# History:     2004-06-  : 
#

rvmean <- function (x, na.rm=FALSE) 
{
    m <- colMeans(sims(as.rv(x)), na.rm=na.rm)
    dx <- dim(x)
    if (!is.null(dx) && prod(dx) == length(m)) {
        dim(m) <- dx
        dimnames(m) <- dimnames(x)
    }
    names(m) <- names(x)
    m
}

E <- function (x, na.rm=FALSE) 
{
    m <- colMeans(sims(as.rv(x)), na.rm=na.rm)
    dx <- dim(x)
    if (!is.null(dx) && prod(dx) == length(m)) {
        dim(m) <- dx
        dimnames(m) <- dimnames(x)
    }
    names(m) <- names(x)
    m
}

Pr <- function (x, na.rm=FALSE) 
{
    s <- sims(as.rv(x))
    if (typeof(s) != "logical") 
        stop("Argument for Pr must be a logical statement such as 'x>0'")
    m <- colMeans(s, na.rm=na.rm)
    dx <- dim(x)
    if (!is.null(dx) && prod(dx) == length(m)) {
        dim(m) <- dx
        dimnames(m) <- dimnames(x)
    }
    names(m) <- names(x)
    m
}


rvvar <- function (x, ...)
{
    v <- rvsimapply(as.rv(x), var, na.rm=TRUE, ...)
    v[rvnsims(x)==1] <- 0
    v
}


rvsd <- function(x,...) sqrt(rvvar(x)) #rvsimapply(as.rv(x), sd, ...)

rvquantile <- function(x, ...) rvsimapply(as.rv(x), quantile, ..., na.rm=TRUE)

rvmedian <- function(x, ...)  rvsimapply(as.rv(x), median, ..., na.rm=TRUE)

rvmin <- function(x, ...) rvsimapply(as.rv(x), min, ...)

rvmax <- function(x, ...) rvsimapply(as.rv(x), max, ...)

rvrange <- function (x, ...) rvsimapply(as.rv(x), range, ...)


# rvsample is now in rv-misc.R


# end rv-alias.R

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# APPLYING
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ========================================================================
# rvlapply  -  apply a function to components of an rv, return a
# ========================================================================

.rvlapply <- function(x, FUN, ...)
{
  attr.x <- attributes(x)
  v <- lapply(x, FUN, ...)
  attributes(v) <- attr.x
  v
}

# ========================================================================
# lapply.rv  -  l-apply a function to components of an rv
# ========================================================================
# Name:        lapply.rv(x, FUN, ...)
# Description: the 'rv' version of lapply
#

lapply.rv <- function(x, FUN, ...)
{
  if (is.null(na <- names(x))) na <- paste(seq(along=x))
  lapply(split(x, na), FUN, ...)
}

# ========================================================================
# .rvMathapply  -  apply a function to components of an rv or non-rv
# ========================================================================
# Name:        .rvMathapply(x, FUN, ...)
# Description: the 'rv' version of lapply
#

.rvMathapply <- function(x, FUN, ...)
{
  v <- lapply(x, FUN, ...)
  dim(v) <- dim(x)
  dimnames(v) <- dimnames(x)
  names(v) <- names(x)
  class(v) <- class(rv())
  v
}

# ========================================================================
# simapply  -  apply a (numeric) function to the simulations, rowwise, with dimensions
# ========================================================================
# Vectorizes over the simulations of one single rv; 
# for vectorization over a group of rvs, see 'simmapply'.
#

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
    dimnames(r) <- dimnames(m0)
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
  if (!is.null(dx) && prod(dx)==length(m)) {
    dim(m) <- dx
    dimnames(m) <- dimnames(x)
  }
  names(m) <- names(x)
  m
}

# ========================================================================
# mapply.rv  -  apply any function to the simulations
# ========================================================================
# TODO: make sure that the argument names are preserved!
# Note. Won't work with functions allowing "blank" arguments
# such as "[" (e.g. x[y,,]). The functions "[" and "[<-" use
# modified versions of simmapply.
#

mapply.rv <- function (FUN, ..., MoreArgs=NULL, USE.NAMES=TRUE, SAMPLESIZE=NULL)
{
  a <- list(...)
  dim.a.names <- dimnames(a)
  a.names <- names(a)
  dimnames(a) <- dim.a.names
  names(a) <- a.names
  a <- lapply(a, .sims.as.list)
  if (!is.null(SAMPLESIZE)) {
    m <- max(sapply(a, length))
    s <- (sample(1:m, size=SAMPLESIZE, replace=TRUE)-1)
    a <- lapply(a, function (x) x[(s %% length(x))+1])
  }
  a <- c(FUN = FUN, a, SIMPLIFY = FALSE, USE.NAMES=USE.NAMES)
  a$MoreArgs <- MoreArgs
  list.of.sims <- do.call(mapply, args = a)
  list.of.sims <- lapply(list.of.sims, function (x) if (is.null(x)) NA else x)
  rvsims(list.of.sims)
}



#
# apply.rv - Apply for Random Matrices
#

apply.rv <- function (X, MARGIN, FUN, ...)
{
  FUNC <- function (x, ...) {
    class(x) <- class(rv())
    FUN(x, ...)
  }
  a <- apply(X, MARGIN, FUNC, ...)
  .list2array(a, drop=TRUE)
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
  if (is.list(v)) class(v) <- class(rv())
  v
}

Ops.rvsim <- function(e1, e2)
{
  sim <- get(.Generic)(sims.rvsim(e1), sims.rvsim(e2)) # Must be vectorizable
  .rvsim(sim)
}


"!.rv" <- function(e1) 
{
  ####.rvlapply(e1, .Generic)
  simapply(e1, .Generic)  # DEBUG ok?
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


Summary.rv <- function(..., na.rm=FALSE)
{
  sims <- sims(c(...)) # an L x n matrix of numbers
  m <- apply(sims, 1, .Generic, na.rm=na.rm)
  rvsims(m)
}

range.rv <- function(..., na.rm=FALSE, finite=FALSE)
{
  sims <- sims(c(...)) # an L x n matrix of numbers
  m <- apply(sims, 1, 'range', na.rm=na.rm, finite=finite) # gives a 2xL matrix!!!
  r <- rvsims(t(m)) # Transpose to get an L x 2 matrix.
  names(r) <- c('min','max')
  r
}
# rv-distr.R
#  generators of iid distributions
# 
# TODO
#  - truncated distributions
#  - 

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
  sim <- matrix(rep(s, times=n.max), nrow=n.sims)
  mc <- matrix(rep(seq(from=1,length=n.max), each=length(s)), nrow=n.sims)
  nc <- matrix(rep(sim.n,times=ncol(sim)), nrow=n.sims)
  sim[nc<mc] <- NA
  sim
}

# ========================================================================
# .parammat - make a parameter matrix (used in .rvgen)
# ========================================================================
# 

.parammat <- function(param, n, n.sims) # NOEXPORT
{
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
  if (length(n)>1) warning("length(n)>1? Random dimensions not yet supported")
  n <- n[1]
  if (!is.character(fun)) fun <- deparse(substitute(fun))
  paramlist <- list(...)
  n.sims <- getnsims()
  arglist <- lapply(paramlist, .parammat, n, n.sims)
  arglist$n <- max(sapply(arglist,length))
  ow <- options("warn")
  options(warn=-1)
  m <- do.call(fun, args=arglist) # will complain about NAs, but ignore.
  options(ow) # reset warnings
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
  lx <- length(x)
  if (lx==0)
    return(NULL)
  v <- rv(lx)
  v[1:lx] <- x # recycle
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
  #
  # from the MASS package
  #
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

.rvmvnorm <- function (n=1, mean, Sigma) 
{
  if (is.null(dim(Sigma))) {
    Sigma <- diag(Sigma)
  } else if (!is.numeric(Sigma)) 
    stop("Invalid (non-numeric) covariance matrix Sigma")
  else if (nrow(Sigma) != ncol(Sigma)) 
     stop("Invalid (nonsquare) covariance matrix Sigma")
  #
  if (length(mean)==1) {
    mean <- rep(mean, length.out=nrow(Sigma))
  }
  if (length(mean) != nrow(Sigma)) {
    stop("Length of mean vector != nrow of Sigma.")
  }
  #
  nm <- names(mean)
  if (is.rv(n)) stop("n cannot be random (to be implemented)")
  n <- n[1]
  if (is.rv(Sigma)) {
    r <- mapply.rv(.mvrnorm, n=n, mu=mean, Sigma=Sigma)
    if (n==1) r <- drop(r)
    if (!is.null(dim(r))) {
      dimnames(r)[[2]] <- nm
    }
    return(r)
  }
  n.sims <- setnsims()
  dim.mean <- dim(mean)
  if (is.list(Sigma)) {
    Sigma <- .expand.as.matrix(Sigma)
  }
  ##if (length(mean) > nrow(Sigma)) {
  ##  # DEBUG: ??
  ##  n.Sigmas <- (length(mean)%/%nrow(Sigma))+(length(mean)%%nrow(Sigma)!=0)
  ##  Sigma <- .expand.as.matrix(rep(rv(Sigma), n.Sigmas))
  ##}
  mean <- rep(mean, length.out = nrow(Sigma))
  mu <- rep(0, length(mean))
  n.all <- n*n.sims
  X <- .mvrnorm(n = n.all, mu = mu, Sigma = Sigma)
  if (n==1) {
    r <- rvsims(X)
  } else {
    dim(X) <- c(n.sims, length(X) %/% n.sims)
    r <- rvsims(X)
    dim(r) <- c(length(mean), n)
  }
  r <- t(r + mean) # Must be in this order to preserve dimensions
  if (n==1) r <- drop(r)
  if (is.null(dim(r))) {
    if (length(r)==length(nm)) names(r) <- nm
  } else {
    dimnames(r) <- list(NULL, nm)
  }
  r
}



# ========================================================================
# rvmvt  -  multivariate t random variables 
# ========================================================================
#

.rvmvt <- function (n=1, Sigma, df=1)
{
  x <- sqrt(rvchisq(n=n, df=df)/df)
  # But will this work? x is of length n,
  #   but the returned thing is of length n * nrow(Sigma)!
  .rvmvnorm(n=n, mean=0, Sigma=Sigma)/x
}

# ========================================================================
# rvt  -  t random variables
# ========================================================================
# FIXED 2006-04-21: is.missing --> missing

rvt <- function (n=1, mu=0, scale=1, df, Sigma)
{
  if (!missing(Sigma)) {
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

rvbern <- function (n=1, prob, logical=FALSE)
{
  r <- .rvgen('rbinom', n=n, size=1, prob=prob)
  if (logical) as.logical(r) else r
}


# ========================================================================
# rvbinom  -  binomial rvs
# ========================================================================

rvbinom <- function (n=1, size, prob)
{
  (.rvgen('rbinom', n=n, size=size, prob=prob))
}

# ========================================================================
# rvmultinom  -  multinomial rvs
# ========================================================================

rvmultinom <- function(n=1, size=1, prob)
{
  if (length(prob)<=1) {
    return(rvbinom(n=n, size=size, prob=prob))
  }
  if (anyisrv(n, size, prob)) {
    r <- mapply.rv("rmultinom", n=n, size=size, prob=prob)
  } else {
    n.sims <- getnsims()
    s <- rmultinom(n=n*n.sims, size=size, prob=prob)
    dim(s) <- c(length(s) %/% n.sims, n.sims)
    r <- rvsims(t(s))
    dim(r) <- c(length(prob), n)
    dimnames(r) <- list(names(prob), NULL)
  }
  r
}

# ========================================================================
# rvbeta  -  beta rvs
# ========================================================================

rvbeta <- function (n=1, shape1, shape2)
{
  .rvgen('rbeta', n=n, shape1=shape1, shape2=shape2)
}


# ========================================================================
# rvcat - categorical rvs
# ========================================================================

rvcat <- function (n=1, prob, levels=NULL)
{
  # NAME
  #  rvcat - Sample Categorical Random Variables
  if (anyisrv(n, prob)) {
    x <- rvmultinom(n=n, size=1, prob=prob)
    s <- sims(x, dimensions=TRUE)
    ds <- dim(s)
    s <- as.logical(s)
    dim(s) <- ds
    f <- function (m) row(m)[m]
    a <- apply(s, 1, f)
    r <- if (is.null(dim(a))) rvsims(a) else rvsims(t(a))
  } else {
    n.sims <- getnsims()
    s <- rmultinom(n=n*n.sims, size=1, prob=prob)
    s <- row(s)[as.logical(s)]
    dim(s) <- c(n.sims, length(s) %/% n.sims)
    r <- rvsims(s)
  }
  rvfactor(r, levels=levels)
}

# ========================================================================
# rvmix - mixtures of random variables
# ========================================================================

rvmix <- function (n=1, prob=NULL, index=NULL, components=list())
{
  # NAME
  #  rvmix - 
  # EXAMPLE
  #   rvmix(1, index=rvbern(1,0.05)+1, components=list(rvnorm(1), rvnorm(1, mean=10)))
  #
  if (is.null(index)) {
    if (is.null(prob)) stop("need either probabilities or an index variable")
    index <- rvcat(n=n, prob=prob)
  }
  mix <- components
  mls <- sapply(mix, length)
  maxl <- max(mls)
  if (!all(mls==maxl)) mix <- lapply(mix, function (x) x[1:maxl])
  mix <- unlist.rv(mix)
  dim(mix) <- c(maxl, length(mix) %/% maxl)
  tsmix <- t(sims(mix))
  f <- function (x) {
    s <- sims(x)
    .index <- s + seq(from=0, to=length(s)-1)*nrow(tsmix)
    mixsims <- tsmix[.index]
    rvsims(mixsims)
  }
  r <- unlist.rv(lapply.rv(index, f))
  names(r) <- NULL
  r
}



# ========================================================================
# rvdirichlet  -  dirichlet rvs
# ========================================================================
# FIXED 2006-04-21 : scale parameter was alpha, fixed to be 1.

rvdirichlet <- function (n = 1, alpha) 
{
  x <- NULL
  for (i in 1:n) {
    g <- rvgamma(n = 1, shape = alpha, scale = 1)
    x <- cbind.rv(x, g/sum(g))
  }
  x
}

# ========================================================================
# rvboot  -  empirical (bootstrap) distribution
# ========================================================================

rvboot <- function (data)
{
  n.sims <- setnsims()
  n <- n.sims*length(data)
  s <- matrix(sample(data, size=n, replace=TRUE), nrow=n.sims)
  r <- rvsims(s)
  dim(r) <- dim(data)
  r
}

# ========================================================================
# rvpermut  -  permutation distribution
# ========================================================================

rvpermut <- function (data, prob=NULL)
{
  n.sims <- getnsims()
  s <- t(sapply(rep(list(data), n.sims), sample, prob=prob))
  r <- rvsims(s)
  dim(r) <- dim(data)
  names(r) <- names(data)
  r
}

# ========================================================================
# 
# ========================================================================

rvdiscrete <- function (n=1, x, prob=NULL)
{
  n.sims <- getnsims()
  rvsims(matrix(sample(x=x, size=n*n.sims, prob=prob, replace=TRUE), nrow=n.sims))
}




rvdens <- function(n=1, FUN, range, unitprecision=10, ...)
{
  # NAME
  #   rvdensity - Sample from a given univariate density using a grid approximation
  # ARGUMENTS
  #   n : number of independent random vector components to draw
  #   FUN : density function
  #   range : range for the grid
  #   unitprecision : number of points per unit
  #   ... : other arguments passed to [FUN].
  #   
  grid <- seq(from=range[1], to=range[2], by=1/unitprecision)
  prob <- FUN(grid, ...)
  n.sims <- setnsims()
  s <- sample(grid, size=n*n.sims, prob=prob, replace=TRUE)
  noise <- runif(n*n.sims, -0.5/unitprecision, 0.5/unitprecision)
  rvsims(matrix(s+noise, nrow=n.sims))
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

cov.default <- getFromNamespace('cov', 'stats') # EXPORT cov.default

cov.rv <- function(x, y=NULL, ...)  ## EXPORT cov.rv
{
  if (!is.matrix(x)) {
    if (is.null(y)) {
      stop("supply both x and y or a matrix-like x")
    }
    x <- as.vector(x)
  }
  mapply.rv(cov.default, x=x, y=y, ...)
}

cov <- function(x, ...) {
  UseMethod("cov")
}

formals(cov.default) <- c(formals(cov.default), alist(... = ))

.rvRegisterFunctionSwitch(cov, "cov", "stats", group=3)


# ========================================================================
# cor  -  short description
# ========================================================================

cor.default <- getFromNamespace('cor', 'stats') # EXPORT cor.default

cor.rv <- function(x, y=NULL, ...)  ## EXPORT cor.rv
{
  if (!is.matrix(x)) {
    if (is.null(y)) {
      stop("supply both x and y or a matrix-like x")
    }
    x <- as.vector(x)
  }
  mapply.rv(cor.default, x=x, y=y, ...)
}

cor <- function(x, ...) {
  UseMethod("cor")
}

formals(cor.default) <- c(formals(cor.default), alist(... = ))

.rvRegisterFunctionSwitch(cor, "cor", "stats", group=3)


# ========================================================================
# function  -  short description
# ========================================================================
#

var.default <- getFromNamespace('var', 'stats') # EXPORT var.default

var.rv <- function(x, ...) ## EXPORT var.rv
{
  simapply(x, var.default, ...)
}

var <- function(x, ...) {
  UseMethod("var")
}

formals(var.default) <- c(formals(var.default), alist(... = ))

.rvRegisterFunctionSwitch(var, "var", "stats", group=3)

# ========================================================================
# sd.rv  -  short description
# ========================================================================

sd <- function (x, na.rm = FALSE) 
{
  if (is.rv(x))
    return(sd.rv(x, na.rm = na.rm))
  else if (is.matrix(x)) 
    apply(x, 2, sd, na.rm = na.rm)
  else if (is.vector(x)) 
    sqrt(var(x, na.rm = na.rm))
  else if (is.data.frame(x)) 
    sapply(x, sd, na.rm = na.rm)
  else sqrt(var(as.vector(x), na.rm = na.rm))
}

sd.rv <- function (x, na.rm = FALSE)
{
  if (is.matrix(x)) 
    apply.rv(x, 2, sd, na.rm = na.rm)
  else if (is.vector(x)) 
    sqrt(var.rv(x, na.rm = na.rm))
  else sqrt(var.rv(as.vector(x), na.rm = na.rm))
}

.rvRegisterFunctionSwitch(sd, "sd", "stats", group=3)

# ========================================================================
# rvcov  -  covariance matrix
# ========================================================================
#

rvcov <- function(x, y=NULL, ...)
{
  if (is.null(y))
    cov(sims(x), y, ...)
  else
    cov(sims(x), sims(y), ...)
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

median.rv <- function(x, na.rm=FALSE) ## EXPORT median.rv
{
  simapply(x, median, na.rm=na.rm)
}
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

#
# rvpar
#

rvpar <- function (...)
{
  args <- list(...)
  rvpar  <- getOption("rv")
  Pars <- rvpar$par
  if (length(args)==0) {
    return(Pars)
  }
  if (is.null(na <- names(args))) {
    args <- unlist(args)
    p <- Pars[args]
    if (length(args)==1) {
      return(p[[args]])
    } else {
      return(p)
    }
  }
  oldpar <- Pars
  for (a in na) { 
    if (nchar(a)>=1) {
      Pars[a] <- args[a]
    }
  }
  rvpar$par <- Pars
  options(rv=rvpar)
  oldpar
}





# ========================================================================
# plot 
# ========================================================================

plot <- function (x, y, ...) 
{
  if (!missing(y) && is.rv(y)) {
    return(plot.rv(x, y, ...))
  } else if (is.null(attr(x, "class")) && is.function(x)) {
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
            graphics:::plot.function(x, y, ...)
        else graphics:::plot.function(x, y, ylab = paste(deparse(substitute(x)), 
            "(x)"), ...)
    }
    else UseMethod("plot")
}

.rvRegisterFunctionSwitch(plot, "plot", "graphics", group=0)

plot.rv <- function (x, y=NULL, ...)
{ 
  .plot.default.rv(x, y, ...)
}

# ========================================================================
# plot.xy.rv  - 
# ========================================================================
#

.plot.xy.rv <- function (xy, type, pch = par("pch"), lty = par("lty"), col = par("col"), bg = NA, cex = 1, lwd = par("lwd"), ...)
{
  ##if (is.null(xy$rv)) {
  ##  return(.Internal(plot.xy(xy, type, pch, lty, col, bg, cex, lwd, ...)))
  ##}
  args <- c(xy, list(type=type, pch=pch, lty=lty, col=col, bg=bg, cex=cex, lwd=lwd, ...))
  typego <- list(p="points.rv", "l"="points.rv", "b"="points.rv")
  if (type=="n") return(invisible(NULL))
  if (is.null(plotroutine <- typego[[type]])) {
    stop("Plot type '", type, "' not yet implemented for rv objects!")
  }
  do.call(plotroutine, args)
  invisible(NULL)
}

# ========================================================================
# xy.coords.rv - modified version of xy.coords to accommodate rvs
# ========================================================================

.xy.coords.rv <- function (x, y = NULL, xlab = NULL, ylab = NULL, log = NULL, recycle = FALSE)
{
    if (is.null(y)) {
        ylab <- xlab
        if (is.language(x)) {
            if (inherits(x, "formula") && length(x) == 3) {
                ylab <- deparse(x[[2]])
                xlab <- deparse(x[[3]])
                y <- eval(x[[2]], environment(x), parent.frame())
                x <- eval(x[[3]], environment(x), parent.frame())
            }
            else stop("invalid first argument")
        }
        else if (is.complex(x)) {
            y <- Im(x)
            x <- Re(x)
            xlab <- paste("Re(", ylab, ")", sep = "")
            ylab <- paste("Im(", ylab, ")", sep = "")
        }
        else if (is.matrix(x) || is.data.frame(x)) {
            x <- data.matrix(x)
            if (ncol(x) == 1) {
                xlab <- "Index"
                y <- x[, 1]
                x <- 1:length(y)
            }
            else {
                colnames <- dimnames(x)[[2]]
                if (is.null(colnames)) {
                  xlab <- paste(ylab, "[,1]", sep = "")
                  ylab <- paste(ylab, "[,2]", sep = "")
                }
                else {
                  xlab <- colnames[1]
                  ylab <- colnames[2]
                }
                y <- x[, 2]
                x <- x[, 1]
            }
        }
        else if (is.list(x) && !is.rv(x)) { #### This is the only change ####
            xlab <- paste(ylab, "$x", sep = "")
            ylab <- paste(ylab, "$y", sep = "")
            y <- x[["y"]]
            x <- x[["x"]]
        }
        else {
            if (is.factor(x)) 
                x <- as.numeric(x)
            xlab <- "Index"
            y <- x
            x <- seq(along = x)
        }
    }
    if (inherits(x, "POSIXt")) 
        x <- as.POSIXct(x)
    if (length(x) != length(y)) {
        if (recycle) {
            if ((nx <- length(x)) < (ny <- length(y))) 
                x <- rep(x, length.out = ny)
            else y <- rep(y, length.out = nx)
        }
        else stop("'x' and 'y' lengths differ")
    }
    if (length(log) && log != "") {
        log <- strsplit(log, NULL)[[1]]
        f <- function (x) ((rvmean(x < 0)>0) & !rv.any.na(x))
        if ("x" %in% log && any(ii <- f(x))) {
            n <- as.integer(sum(ii))
            warning(sprintf(ngettext(n, "%d x value <= 0 omitted from logarithmic plot", 
                "%d x values <= 0 omitted from logarithmic plot"), 
                n), domain = NA)
            x[ii] <- NA
        }
        if ("y" %in% log && any(ii <- f(y))) {
            n <- as.integer(sum(ii))
            warning(sprintf(ngettext(n, "%d y value <= 0 omitted from logarithmic plot", 
                "%d y values <= 0 omitted from logarithmic plot"), 
                n), domain = NA)
            y[ii] <- NA
        }
    }
    return(list(x = as.real(x), y = as.real(y), xlab = xlab, ylab = ylab))
}

###.rvRegisterFunctionSwitch(xy.coords, "xy.coords", "grDevices")


# ========================================================================
# plot.default.rv - 
# ========================================================================


.plot.default.rv <- function (x, y = NULL, type = "p", xlim = NULL, ylim = NULL, log = "", main = NULL, sub = NULL, xlab = NULL, ylab = NULL, ann = par("ann"), axes = TRUE, frame.plot = axes, panel.first = NULL, panel.last = NULL, asp = NA, out.of.range.marker = rvpar("oorm"), rvlwd = rvpar("rvlwd"),  rvcol=rvpar("rvcol"), rvinterval=rvpar("rvinterval"), rvlex=rvpar("rvlex"), ...) 
{
    localAxis <- function(..., col, bg, pch, cex, lty, lwd) Axis(...)
    localBox <- function(..., col, bg, pch, cex, lty, lwd) box(...)
    localWindow <- function(..., col, bg, pch, cex, lty, lwd) plot.window(...)
    localTitle <- function(..., col, bg, pch, cex, lty, lwd) title(...)
    xlabel <- if (!missing(x)) 
        deparse(substitute(x))
    ylabel <- if (!missing(y)) 
        deparse(substitute(y))
    xy <- .xy.coords.rv(x, y, xlabel, ylabel, log)
    xlab <- if (is.null(xlab)) 
        xy$xlab
    else xlab
    ylab <- if (is.null(ylab)) 
        xy$ylab
    else ylab
    xlim <- if (is.null(xlim)) {
        x.sims <- sims(as.rv(xy$x))
        range(x.sims[is.finite(x.sims)])
    } else xlim
    ylim <- if (is.null(ylim)) {
        y.sims <- sims(as.rv(xy$y))
        range(y.sims[is.finite(y.sims)])
    } else ylim
    plot.new()
    localWindow(xlim, ylim, log, asp, ...)
    panel.first
    .plot.xy.rv(xy, type, out.of.range.marker = out.of.range.marker, rvlwd = rvlwd,  rvcol=rvcol, rvinterval=rvinterval, rvlex=rvlex, ...)
    panel.last
    if (axes) {
        localAxis(x, side = 1, ...)
        localAxis(y, side = 2, ...)
    }
    if (frame.plot) 
        localBox(...)
    if (ann) 
        localTitle(main = main, sub = sub, xlab = xlab, ylab = ylab, 
            ...)
    invisible()
}


###.rvRegisterFunctionSwitch(plot.default, "plot.default", "graphics")

# ========================================================================
# points.rv  -  plot points and uncertainty intervals
# ========================================================================
# 

points.rv <- function (x, y = NULL, type = "p", xlim = NULL, ylim = NULL, pch = 19, out.of.range.marker = rvpar("oorm"), rvlwd = rvpar("rvlwd"), rvcol = rvpar("rvcol"), rvinterval = rvpar("rvinterval"), rvlex = rvpar("rvlex"), ...) 
{
    xy <- .xy.coords.rv(x, y)
    x <- as.rv(xy$x)
    y <- as.rv(xy$y)
    arg <- list(...)
    draw.points <- (type == "p" || type == "b")
    draw.lines <- (type == "l" || type == "b")
    point.sample <- rvpar("point.sample")
    if (is.null(point.sample)) 
        point.sample <- NA
    line.sample <- rvpar("line.sample")
    if (is.null(line.sample)) 
        line.sample <- NA
    if (is.null(rvlex) || !is.function(rvlex)) {
        rvlex <- function(lwd) 1.5
    }
    if (is.null(rvcol) || is.na(rvcol)) {
        rvcol <- "default"
    }
    x.rv <- (is.random(x) & !rv.all.na(x))
    y.rv <- (is.random(y) & !rv.all.na(y))
    x.point <- (!x.rv)
    y.point <- (!y.rv)
    vertical.pair <- (x.point & y.rv)
    horizontal.pair <- (x.rv & y.point)
    iv.pair <- (vertical.pair | horizontal.pair)
    rv.pair <- (x.rv & y.rv)
    point.pair <- (x.point & y.point)
    segs <- NULL
    pts <- NULL
    cols <- NULL
    lwds <- NULL
    qx <- rvintervals(x, rvinterval)
    qy <- rvintervals(y, rvinterval)
    rvcol <- rep(rvcol, length.out = length(x))
    cols <- t(sapply(rvcol, rvcolortheme))
    dimnames(cols) <- list(NULL, rvinterval)
    cols <- cbind(default = if (is.null(arg$col)) 
        "black"
    else arg$col, cols)
    rvlwd <- rep(rvlwd, length.out = length(x))
    lwds <- t(sapply(rvlwd, function(wd) c(0, wd * rvlex(wd), 
        wd)))
    dimnames(lwds) <- list(NULL, rvinterval)
    lwds <- cbind(default = if (is.null(arg$lwd)) 
        "black"
    else arg$lwd, lwds)
    if (any(point.pair)) {
        x.pts <- E(x[point.pair])
        y.pts <- E(y[point.pair])
        pts <- list(x = c(pts$x, x.pts), y = c(pts$y, y.pts), 
            col = c(pts$col, cols[point.pair, 1]), lwd = c(pts$lwd, 
                lwds[point.pair, 1]))
    }
    if (any(rv.pair)) {
## THIS CODE IS NOT YET FULLY FUNCTIONAL
        xy.pts <- rvsample(c(x[rv.pair],y[rv.pair]), jointly=TRUE, size= point.sample)
        x.pts <- xy.pts[,1]
        y.pts <- xy.pts[,2]
        ## DEBUG here '19' is hard-coded!
        pchs <- rep(19, length.out = length(x.pts))
        cl <- rep(cols[rv.pair, 1], each = length(x.pts)/sum(rv.pair))
        pts <- list(x = c(pts$x, x.pts), y = c(pts$y, y.pts), 
            col = c(pts$col, cl), pch = c(pts$pch, pchs))
    }
    for (name in names(qx)) {
        if (is.na(name)) 
            next
        xiv <- qx[[name]]
        yiv <- qy[[name]]
        is.seg <- (nrow(xiv) == 2 || nrow(yiv) == 2)
        if (is.seg) {
            if (any(iv.pair)) {
                segs <- list(x0 = c(segs$x0, xiv[1, iv.pair]), 
                  y0 = c(segs$y0, yiv[1, iv.pair]), x1 = c(segs$x1, 
                    xiv[2, iv.pair]), y1 = c(segs$y1, yiv[2, 
                    iv.pair]), col = c(segs$col, cols[iv.pair, 
                    name]), lwd = c(segs$lwd, lwds[iv.pair, name]))
            }
        }
        else {
            pchs <- rep(if (is.null(arg$pch)) 
                19
            else arg$pch, length.out = length(x))
            pts <- list(x = c(pts$x, xiv[1, iv.pair]), y = c(pts$y, 
                yiv[1, iv.pair]), col = c(pts$col, cols[iv.pair, 
                name]), pch = c(pts$pch, pchs))
        }
    }
    if (draw.points) {
        if (!is.null(segs)) 
            do.call("segments", args = .nodups(c(arg, segs)))
        if (!is.null(pts)) 
            do.call("points", args = .nodups(c(arg, pts)))
        if (any(rv.pair)) {
            do.call("points", args = .nodups(c(arg, pts)))
        }
    }
    if (draw.lines) {
        lns <- .rvjointdrawxy(x, y, size = line.sample)
        lns$x <- rbind(lns$x, NA)
        lns$y <- rbind(lns$y, NA)
        do.call("lines", args = .nodups(c(arg, lns)))
    }
    invisible(NULL)
}


.rvjointdrawxy <- function (x, y, size=1, reject.na=TRUE)
{
  xy <- c(x, y)
  s <- rvsample(xy, size=size, jointly=TRUE, reject.na=reject.na)
  if (is.null(dim(s))) s <- t(s)
  xs <- t(s[,seq(along=x)])
  ys <- t(s[,-seq(along=x)])
  list(x=xs, y=ys)
}


# ========================================================================
# abline.rv  -  a + bx
# ========================================================================
# change this to
#   mapply.rv("abline", ..., MoreArgs=NULL)
#

.abline.default <- graphics:::abline

abline <- function (a = NULL, b = NULL, h = NULL, v = NULL, ...)
{
  if (anyisrv(a,b,h,v)) {
    line.sample <- rvpar("line.sample")
    if (!is.numeric(line.sample)) stop("rvpar('line.sample') is not a number")
    args <- list(FUN=.abline.default, a=a, b=b, h=h, v=v)
    nulls <- sapply(args, is.null)
    nullArgs <- names(nulls)[nulls]
    MoreArgs <- list(...)
    args[nullArgs] <- NULL
    MoreArgs[nullArgs] <- list(NULL)
    args$SAMPLESIZE <- line.sample
    args$MoreArgs <- MoreArgs
    do.call("mapply.rv", args=args)
  } else {
    .abline.default(a=a, b=b, h=h, v=v, ...)
  }
  invisible()
}

###.rvRegisterFunctionSwitch(plot.default, "plot.default", "graphics")
.rvRegisterFunctionSwitch(abline, "abline", "graphics", group=0)

# ========================================================================
# lines.rv  -  plot some random lines
# ========================================================================
# btw, "rvlines" does not make sense - that'd be something like a weird histogram

lines.rv <- function(x, y, type="l", ...)
{
  points.rv(x, y, type="l", ...)
}

# ========================================================================
# hist.rv  -  histogram, adapted for rv's
# ========================================================================

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
#

rvhist <- function(x, ...)
{
  if (!is.null(dim(x)))
    par(mfcol=dim(x))
  mfcol <- par('mfcol')
  n <- prod(mfcol)
  a <- list(...)
  make.main <- is.null(a$main)
  make.xlab <- is.null(a$xlab)
  lab <- deparse(substitute(x))
  x.names <- paste(lab, .dimindex(x), sep="")
  for (i in 1:n) {
    a$x <- sims(x[i])
    if (make.main || make.xlab) {
      this.name <- x.names[i]
      if (make.xlab) {
        a$xlab <- this.name
      }
      if (make.main) {
        a$main <- paste('Histogram of', this.name)
      }
    }
    if (is.null(a$breaks)) a$breaks <- "fd"
    do.call("hist", a)
  }
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


# ================================================================================
# rvcolortheme  -  
# ================================================================================


# Put the following in rv/par$rvcolorthemes

.rvcolorthemes <- list(
  default=c("grey20", "grey40", "grey60"),
  gray=c("grey20", "grey40", "grey60"),
  lightgray=c("grey20", "grey40", "grey60"),
  darkgray=c("black", "grey20", "grey40")
)

.makervcolortheme <- function (col)
{
  colors.with.numbers <- colors()[regexpr("[a-z]1$", colors())>0]
  colors.with.numbers <- sub("1$", "", colors.with.numbers)
  if (!col %in% colors.with.numbers) {
    if (col %in% colors())
      return(rep(col,3))
    else 
      return(NULL)
  } 
  paste(col, c("3", "2", ""), sep="") 
}

rvcolortheme <- function (theme) # NOEXPORT
{
  # NAME
  #  rvcolortheme - color themes for plotting uncertainty intervals
  # DETAILS
  # (please, later export this!)
  synonym <- list(grey="gray", lightgrey="lightgray", darkgrey="darkgray")
  theme <- theme[1]
  if (is.null(theme) || is.na(theme)) {
    theme <- "default"
  }
  (!is.null(co <- .rvcolorthemes[[theme]])) && return(co)
  (!is.null(stheme <- synonym[[theme]])) && return(.rvcolorthemes[[stheme]])
  co <- .makervcolortheme(theme)
  if (is.null(co)) {
    warning(paste("No rv color theme '",theme,"' found, using default",sep=""))
    co <- .rvcolorthemes[["default"]]
  } else {
    .rvcolorthemes[[theme]] <- co
  }
  
  names(co) <- c("dot", "thick", "thin")
  co
}

### MLPLOT ###


.mlplot3 <- function (x, ..., main=NULL, add=FALSE) 
{
  x.name <- deparse(substitute(x))
  dx <- dim(x)
  n <- dx[3]
  n.col <- ceiling(sqrt(n))
  n.row <- n %/% n.col
  if (n.col*n.row<n) n.row <- n.row+1
  grid <- c(n.row, n.col)
  #
  d.n <- dimnames(x)[[3]]
  if (length(d.n)!=n) {
    d.n[n] <- NA
    def.names <- paste(x.name, "[,,", 1:n, "]", sep="")
    d.n[is.na(d.n)] <- def.names[is.na(d.n)]
  }
  par(mfrow=grid, oma=c(0,0,5,0), mar=c(2,1,1,0.5)+0.1)
  for (k in 1:n) {
    main.title <- d.n[k]
    mlplot(x[,,k], ..., main=main.title, add=FALSE)
  }
  if (is.null(main)) {
    main <- x.name
  }
  mtext(main, line=1, outer=TRUE)
  invisible()
}

.order.rows <- function (x)
{
  row.order <- t(apply(x, 1, order))
  ord <- array(seq(along=x),dim(x))
  for (j in 1:nrow(ord)) {
    ord[j,] <- ord[j,row.order[j,]]
  }
  ord
}

mlplot <- function (X, y.center = TRUE, y.shift = 0, y.map = NULL, mar = par("mar"), left.margin = 3, top.axis = TRUE, exp.labels=FALSE, x.ticks = NULL, axes = NULL, xlim = NULL, ylim = NULL, xlab=deparse(substitute(X)), ylab=NULL, las = NULL, add = FALSE, ...) 
{
    dx <- dim(X)
    lx <- length(dx)
    if (lx > 3) {
        stop("Too many dimensions: ", lx)
    }
    else if (lx == 3) {
        return(.mlplot3(X, y.center=y.center, y.shift=y.shift, y.map = y.map,
             mar = mar, left.margin = left.margin, top.axis = top.axis, 
             exp.labels=exp.labels, x.ticks=x.ticks, axes = axes, 
             xlim=xlim, ylim = ylim,
             xlab=xlab, ylab = ylab,
             las = las,  
             add = add, ...))
    }
    else if (lx == 2) {
        labels <- dimnames(X)[[1]]
        if (any(rv.all.na(X))) {
            f <- function(x) {
                if (!any(is.na <- rv.all.na(x))) 
                  return(x)
                w <- which(is.na)
                c(x[-w], x[w])
            }
            X <- t(apply.rv(X, 1, f))
        }
    }
    else {
        labels <- names(X)
        dim(X) <- c(length(X), 1)
    }
    y.row.coords <- 1:nrow(X)
    ylim <- rev(range(y.row.coords) + c(-1, 1))
    if (is.null(y.map)) {
        y.map <- function(x) {
            if (y.center) {
                f <- function(i, j, n, nc) {
                  i + (j - 1)/n - 0.5 * (nc - 1)/n
                }
            }
            else {
                f <- function(i, j, n, nc) {
                  i + (j - 1)/n
                }
            }
            f(i = row(x), j = col(x), n = max(ncol(x),10)*1.5, nc = ncol(x))
        }
    }
    if (is.function(y.map)) {
        y <- y.map(X)
    }
    else if (is.numeric(y.map)) {
        y <- y.map
    }
    else {
        stop("Unknown type of y.map")
    }
    if (is.numeric(y.shift)) {
        y <- (y + y.shift)
    }
    else {
        stop("y.shift must be numeric")
    }
    if (length(y) != length(X)) {
        stop("y coordinates are not valid (check y.map!)")
    }
    if (is.null(xlim)) {
        x.sims <- sims(as.rv(X))
        rng <- range(x.sims[is.finite(x.sims)])
        if (length(rng) < 2) {
            rng <- c(-1, 1)
        }
        xlim <- rng
    }
    if (is.null(x.ticks)) {
      x.row.coords <- pretty(xlim)
    } else {
      x.row.coords <- x.ticks
    }
    if (exp.labels) {
      x.labels <- paste(signif(exp(x.row.coords),2))
    } else {
      x.labels <- paste(x.row.coords)
    }
    mar <- (mar + c(0, left.margin, 2, 0))
    oldpar <- par(mar = mar)
    on.exit(par(oldpar))
    las <- if (!is.null(las)) las else 1
    if (add) {
        points(X, y, xlim = xlim, ylim = ylim, ...)
    }
    else {
        plot(X, y, ..., las = las, xlim = xlim, ylim = ylim, 
            axes = FALSE, xlab=xlab, ylab = "")
        if (is.null(axes) || axes) {
            axis(1, at = x.row.coords, labels=x.labels)
            if (top.axis) 
                axis(3, at = x.row.coords, labels=x.labels)
            if (is.null(labels)) 
                labels <- paste(y.row.coords)
            axis(2, at = y.row.coords, labels = labels, tick = FALSE, 
                line = FALSE, pos = NA, outer = FALSE, font = NA, 
                las = 1)
        }
    }
    invisible(NULL)
}




rvintervals <- function (x, rvinterval=rvpar("rvinterval"), ...) # NOEXPORT
{
  which.quantiles <- list(
    "NA" = NA,
    mean = NA,
    median = 0.50,
    "50%" = c(0.25, 0.75),
    "80%" = c(0.10, 0.90),
    "95%" = c(0.025, 0.975)
  )
  .whichq <- function (iv) {
    if (is.numeric(iv)) {
      iv <- paste(100*iv, "%", sep="")
    } else {
      (!is.null(q <- which.quantiles[[iv]])) && return(q)
    }
    if (is.na(iv)) return(NA)
    n <- nchar(iv)
    if (substr(iv,n,n)=="%") {
      ivn <- as.numeric(substr(iv,1,n-1))
      c(0.5-ivn/200, 0.5+ivn/200)
    } else {
      NA
    }
  }
  .length <- function (iv) {
    lg <- .whichq(iv)
    if (length(lg)<=1) 0 else diff(lg)
  }
  probs <- as.vector(na.omit(unlist(lapply(rvinterval, .whichq))))
  if (length(probs)<=1) {
    # A trick to force probs into a named array
    # (won't otherwise return names if we have only one quantile, e.g. 0.50)
    probs <- c(probs, NA)
  }
  if (!all(is.na(probs))) {
    Q <- rvquantile(x, probs=probs, ...)
  } else {
    Q <- NULL # DEBUG: will this be ignored if we have only "mean" e.g.? #
  }
  compute.what <- list(
    "NA"   = function () NA,
    mean   = function () t(as.vector(rvmean(x, na.rm=TRUE))),
    median = function () Q["50%",,drop=FALSE],
    "50%"  = function () Q[c("25%","75%"),,drop=FALSE],
    "80%"  = function () Q[c("10%","90%"),,drop=FALSE],
    "95%"  = function () Q[c("2.5%","97.5%"),,drop=FALSE]
  )
  .lbl <- function (p) { # From 'quantile.default'
    if (is.null(p) || is.na(p)) return(NA)
    dig <- max(2, getOption("digits"))
    paste(formatC(100 * p, format = "fg", wid = 1, digits = dig), "%", sep = "")
  }
  .summaries <- function (iv) {
    if (is.null(f <- compute.what[[iv]])) {
      a <- na.omit(sapply(.whichq(iv),.lbl))
      if (all(a %in% dimnames(Q)[[1]])) {
        return(Q[a,,drop=FALSE])
      } else {
        warning("Cannot understand interval '", iv, "'")
        return(NA)
      }
    }
    a <- f()
    if (is.null(dim(a))) {
       if (length(x)==1) {
          a <- t(a)
       } else {
         na <- names(a)
         dim(a) <- c(1,length(a))
         dimnames(a) <- list(iv, na)
       }
    }
    return(a)
  } 
  lgth <- rev(order(sapply(rvinterval, .length)))
  s <- lapply(rvinterval, .summaries)
  names(s) <- rvinterval
  s[lgth]
}





# rv-graph.R



# ========================================================================
# mcmcarray  -  coerce rv into a LxMxN mcmc array
# ========================================================================
#  

mcmcarray <- function(x) # NOEXPORT : OBSOLETE
{
  has.chain <- !is.na(rvnchains(x))
  if (!any(has.chain))
    return(NULL)
  a <- list()
  ord <- rvattr(x,'order')
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
# rvattr  -  return the rvsim attribute for each component
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


# ========================================================================
# rvhyper  -  hypergeometric rvs
# ========================================================================

rvhyper <- function(nn=1, m, n, k) # NOEXPORT -- not yet ready
{
  if (anyisrv(nn, m, n, k)) {
    r <- mapply.rv("rhyper", nn=n, m=m, n=n, k=k)
  } else {
    # DEBUG: won't yet work properly as a vectorizing function!
    n.sims <- getnsims()
    s <- rhyper(nn=nn*n.sims, m=m, n=n, k=k)
    nx <- (length(s) %/% n.sims)
    dim(s) <- c(nx, n.sims)
    r <- rvsims(t(s))
  }
  r
}




rvreject <- function (x)  # NOEXPORT -- not yet ready
{
  reject <- function (x) (!(x & NA))
  simapply(x, reject)
}

# ===========================================================================

rvifelse <- function (test, yes, no)  # NOEXPORT -- not yet ready
{
  # NAME
  #   rvifelse - If-Then-Else For Random Vectors
  # DESCRIPTION
  #   ifelse for random vectors
  # DETAILS
  #   Not (yet) vectorized.
  #
  len <- max(length(yes),length(no))
  if (len==1) {
    sim <- t(sims(as.rv(c(yes[1],no[1]))))
    testsim <- as.logical(sims(test[1]))
    tst <- rbind(testsim, !testsim)
    rvsims(sim[tst])
  } else {
    # DEBUG: make this as efficient as the above code:
    if (length(yes)!=length(no)) {
      YN <- rbind.rv(no[1:len],yes[1:len])
    } else {
      YN <- rbind.rv(no,yes)
    }
    YN[as.logical(test)+1,]
  }
}




rvresample <- function (x, jointly=TRUE)  # NOEXPORT -- not yet ready
{
  # NAME
  #   rvresample - get rid of NAs by resampling non-NA simulations
  # DESCRIPTION
  #   NAs are treated as 'rejected values' and not 'missing' ones
  xa <- attributes(as.rv(x))
  ns <- max(rvnsims(x))
  xs <- rvsample(as.rv(x), size=ns, jointly=jointly, reject.na=TRUE)
  x <- rvsims(xs)
  attributes(x) <- xa
  x
}

rvtruncate <- function (x, x.intervals=NULL, middle.interval=NULL, quantile.pairs=NULL, outside=FALSE)  # NOEXPORT -- not yet ready
{
  x.attr <- attributes(as.rv(x))
  x.sims <- sims(as.rv(x))
  if (!is.null(x.intervals)) {
    if (is.null(dim(x.intervals))) {
      iv <- matrix(rep(x.intervals, length(x)), ncol=length(x))
    }
  } else if (!is.null(quantile.pairs)) {
    iv <- rvquantile(x, quantile.pairs)
  } else if (!is.null(middle.interval)) {
    iv <- rvintervals(x, middle.interval)[[1]]
  } else {
    return(rvsims(x.sims[sample(nrow(x.sims)),]))
  }
  if (!is.matrix(iv)) {
    stop("Interval not valid")
  }
  if (ncol(iv)!=length(x)) {
    stop("Invalid number of columns...") # DEBUG: not very smart error messages
  }
  if (nrow(iv) %% 2 != 0) {
    stop("Invalid number of rows...")
  }
  for (i in seq(from=1, to=nrow(iv), by=2)) {
    x.sims <- t(x.sims) # cols=simulation dimension
    chosen <- (x.sims >= iv[i,] & x.sims <=iv[i+1,])
    if (outside) chosen <- (!chosen)
    x.sims[!chosen] <- NA
  }
  x.sims <- t(x.sims)
  x <- rvsims(x.sims)
  attributes(x) <- x.attr
  x
}





# ========================================================================================
# rv.any.na - does an rv have missing values?
# ========================================================================================

rv.any.na <- function (x)
{
  # NAME
  #  rv.any.na - Which components have missing values?
  #
  (colSums(is.na(sims(as.rv(x))))>0)
}



rv.all.na <- function (x)
{
  # NAME
  #  rv.all.na - Which components are completely missing?
  (colMeans(is.na(sims(as.rv(x))))==1)
}

# ========================================================================================
# rvci - uncertainty/credible interval
# ========================================================================================

rvci <- function(obj, interval=0.95, one.sided=FALSE, left=TRUE)
{
  # NAME
  #  rvci - Uncertainty (Credible) Interval
  # 
  if (length(interval)>1) {
    ci <- t(rvquantile(obj, range(interval)))
  } else if (one.sided) {
    q <- if (left) interval else (1-interval)
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
  # NAME
  #  rvinci - Is a Given Value in the Credible Interval?
  # 
  ci <- rvci(obj, interval=interval, one.sided=one.sided, left=left)
  (x>=ci[,1] && x<=ci[,2])
}

# 
# outer.rv - 
#

outer.rv <- function (X, Y=NULL, FUN="*", ...)
{
  # NAME
  #   outer.rv - 
  # 
  if (is.null(Y)) {
    mapply.rv("outer", X, X, MoreArgs=list(FUN="*", ...))
  } else {
    mapply.rv("outer", X, Y, MoreArgs=list(FUN="*", ...))
  }
}


rvoverlap <- function(obj, iv, interval=0.95, one.sided=FALSE, left=TRUE) 
{
  # NAME
  #  rvoverlap - Which Random Scalars Fall Within Which Given Interval?
  # ARGUMENTS
  #   obj : an rv object
  #   iv  : given interval
  #   interval : interval size of obj to test
  # VALUE
  #   -Inf : if the c.i. is to the left of the given interval
  #   +Inf : if the c.i. is to the right of the given interval
  #   -1   : if the c.i. is partially covered by the interval on the left
  #   +1   : if the c.i. is partially covered by the interval on the right
  #   0    : if the c.i. is fully covered by the given range
  #
  iv <- range(iv)
  ci <- rvci(obj, interval=interval, one.sided=one.sided, left=left)
  x <- rep(0, times=nrow(ci))
  x[ ci[,2] < iv[1] ] <- (-Inf)
  x[ ci[,1] <= iv[1] & ci[,2] >= iv[1] ] <- (-1)
  x[ ci[,1] > iv[2] ] <- (+Inf)
  x[ ci[,2] >= iv[2] & ci[,1] <= iv[2] ] <- (+1)
  dim(x) <- dim(obj)
  dimnames(x) <- dimnames(obj)
  names(x) <- names(obj)
  x 
}





aperm.rv <- function (x, ...)
{
  # NAME
  #   aperm.rv - Transpose a Random Array
  #   
  a <- aperm(x, ...)
  if (!is.rv(x)) return(a)
  class(a) <- class(x)
  a
}


# =========
# unlist.rv
# =========


unlist.rv <- function (x, recursive = TRUE, use.names = TRUE) 
{
  # NAME
  #   unlist.rv - Flatten Lists That (May) Contain Random Variables 
  y <- NULL
  ix <- seq(along=x)
  xn <- names(x)
  .paste <- function (name, x) {
     nbrs <- .dim.index(x, leftadjust=FALSE)
     paste(name, nbrs, names(x), sep="")
  }
  if (recursive) {
    for (i in ix) {
      nx <- xn[i]
      if (use.names && is.null(nx)) nx <- "."
      if (!is.rv(x[[i]]) && is.list(x[[i]])) {
        new.y <- unlist.rv(x[[i]], recursive=TRUE, use.names=use.names)
      } else {
        new.y <- x[[i]]
      }
      if (is.null(names(new.y))) {
        new.names <- .paste(nx, new.y)
      } else {
        new.names <- paste(nx, ".", names(new.y), sep="")
      }
      yn <- names(y)
      y <- c(y, new.y)
      names(y) <- c(yn, new.names)
    }
  } else {
    for (i in ix) {
      nx <- xn[i]
      new.y <- x[[i]]
      new.names <- .paste(nx, new.y)
      yn <- names(y)
      y <- c(y, new.y)
      names(y) <- c(yn, new.names)
    }
  }
  y
}



rep.rv <- function (x, times, ...)
{
  if (!is.rv(x)) return(rep(x, times, ...))
  a <- rep(unclass(x), times, ...)
  class(a) <- class(x)
  a
}


rvminmax <- function (x, na.rm=FALSE, finite=FALSE, defaults=c(-Inf,Inf)) # NOEXPORT
{
  .minmax <- function (x, FUN, def) {
    if (na.rm) x <- x[!is.na(x)]
    if (finite) x <- x[is.finite(x)]
    if (length(x)==0) def else FUN(x)
  }
  m0 <- rvsimapply(x, .minmax, min, defaults[1])
  m1 <- rvsimapply(x, .minmax, max, defaults[2])
  list(min=m0, max=m1)
}

finiterange <- function (x) ## NOEXPORT
{
  x <- x[!is.na(x)]
  x <- x[is.finite(x)]
  if (length(x)==0) return(c(NA,NA))
  range(x)  
}

rvfiniterange <- function (x) ## NOEXPORT
{
  rvsimapply(x, finiterange)
}


apply.rv <- function (X, MARGIN, FUN, ...)
{
  # NAME
  #   apply.rv - Apply Functions Over Random Array Margins
  FUNC <- function (x, ...) {
    class(x) <- class(rv())
    FUN(x, ...)
  }
  a <- apply(X, MARGIN, FUNC, ...)
  .list2array(a, drop=TRUE)
}


imputeIntoTemplate <- function (x, template)
{
  # NAME
  #   imputeIntoTemplate - force a vector into specified shape
  # ARGUMENTS
  #   x : (vector or array) 
  #   template : (logical vector or array)
  # DESCRIPTION 
  #   imputes 'x' into the array shape given by 'template',
  #   which is a logical array of the desired dimension 
  #   (usually the dimnames are also important!)
  # DETAILS
  #   The values in the components of the template vector/array
  #   determine how the imputation is to be done. 
  #   A 'TRUE' value indicates a component that is replaced
  #   by a value in x;
  #   Value 'FALSE' indicates a component that is NOT imputed,
  #   but skipped and later imputed an 'NA';
  #   a component with 'NA' in the template _is imputed_ but
  #   it is later masked with an 'NA'. 
  #   Thus only components with the value 'TRUE' will appear non-NA
  #   in the return value.
  # 
  if (!is.logical(template)) {
    stop("Template must be a logical vector or array")
  }
  new.x <- template
  nas <- template
  has.nas <- any(is.na(template))
  if (has.nas) {
    template[is.na(template)] <- TRUE
  }
  new.x[template] <- x
  if (has.nas) {
    new.x[is.na(nas)] <- NA
  }
  dim(new.x) <- dim(template)
  dimnames(new.x) <- dimnames(template)
  new.x[!template] <- NA
  new.x
}


splitAndImpute <- function (x, templates)
{
  a <- splitbyname(x)
  if (!is.list(templates)) {
    return(a)
  }
  lapply(names(a), function (n) imputeIntoTemplate(a[[n]], template=templates[[n]]))
}


splitbyname <- function (x)
{
  a <- split(x, f = .shortnames(x))
  lapply(a, .setDimensionByName)
}


#
# 
#


rvsample <- function(x, size=1, jointly=TRUE, reject.na=FALSE)
{
  # NAME
  #   rvsample - Draw Samples from Random Vectors
  #
  xs <- sims(x)
  ns <- nrow(xs)
  if (is.null(size) || is.na(size)) size <- ns
  if (jointly) {
    if (reject.na) {
      f <- function (x) any(is.na(x))
      is.na.xs <- apply(xs, 1, f)
      if (all(is.na.xs)) {
        s <- sample(ns, size=size, replace=TRUE, prob=is.na.xs)
      } else {
        s <- sample(ns, size=size, replace=TRUE, prob=!is.na.xs) 
      }
    } else {
      s <- sample(ns, size=size, replace=TRUE)
    }
    s <- xs[s,]
  } else {
    s <- apply(xs, 2, function (s) {
      if (all(nas <- is.na(s))) return(s)
      sample(s[!nas], size=size, replace=TRUE)
    })
  }
  return(s)
}


# ===========================================================================




####


# ========================================================================
# print.rv  -  print summary of a rv on the console
# ========================================================================
# Name:        
# Description: 
# Parameters:  
# Required:    none
# History:     2004-06-  : 
#

print.rv <- function(x, digits=NULL, ...)
{
  s <- summary(x)
  qc <- attr(s, "quantile.cols")
  if (is.null(s))
    return(print(s))
  if (is.null(digits)) digits <- rvpar("print.digits")
  s <- as.data.frame(round(s,digits))
  ds <- dimnames(s)
  if (!is.null(qc)) {
    nqc1 <- seq(from=1, to=min(qc)-1)
    nqc2 <- seq(from=max(qc)+1, to=ncol(s))
    s <- cbind(s[,nqc1], ' '='(', s[,qc], ' '=')', s[nqc2])
  }
  if (!is.null(.names <- names(x))) {
    s <- cbind(name=.names, s)
  } else if (!is.null(unlist(.dn <- dimnames(x))) && length(dim(x))==2) {
    # 'is.null(unlist(dimnames))' since we might have e.g. list(NULL, NULL) 
    sud <- rvpar("summary.dimnames")
    if (!is.null(sud) && !is.na(sud) && is.logical(sud)) {
      da <- lapply(.dn, function (na) if (is.null(na)) rep("", nrow(s)) else na)
      rw <- da[[1]][row(x)]
      cl <- da[[2]][col(x)]
      s <- cbind(row=rw, col=cl, " "=":", s)
    }
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

.rvsummary.integer <- function (x)
{
  s <- sims(x)
  m <- colMeans(s, na.rm=TRUE)
  ns <- rvnsims(x)
  sds <- sqrt((colSums(s^2)-ns*(m^2))/(ns-1))
  sds[ns==1] <- 0
  qs <- c(0, 0.01,0.025,0.25,0.5,0.75,0.975,0.99, 1)
  Q <- t(apply(s, 2, quantile, probs=qs, na.rm=TRUE))
  dimnames(Q)[[2]][1] <- "min"
  dimnames(Q)[[2]][length(qs)] <- "max"
  NAS <- apply(s, 2, is.na)
  if (!is.null(dim(NAS))) NAS <- colMeans(NAS)*100
  nax <- if (is.null(dim(x))) NULL else names(x)
  if (all(NAS==0)) {
    S <- cbind(mean=m, sd=sds, Q, sims=ns)
    qc <- (2+seq(along=qs))
  } else {
    S <- cbind(mean=m, sd=sds, Q, "NA%"=NAS, sims=ns)
    qc <- (2+seq(along=qs))
  }
  #if (!is.null(nax)) dimnames(S)[[1]] <- nax
  attr(S, "quantile.cols") <- qc
  return(S)
}

.rvsummary.continuous <- function (x)
{
  s <- sims(x)
  m <- colMeans(s, na.rm=TRUE)
  ns <- rvnsims(x)
  sds <- sqrt((colSums(s^2)-ns*(m^2))/(ns-1))
  sds[ns==1] <- 0
  qs <- c(0.01,0.025,0.25,0.5,0.75,0.975,0.99)
  Q <- t(apply(s, 2, quantile, probs=qs, na.rm=TRUE))
  NAS <- apply(s, 2, is.na)
  if (!is.null(dim(NAS))) NAS <- colMeans(NAS)*100
  nax <- if (is.null(dim(x))) NULL else names(x)
  if (all(NAS==0)) {
    S <- cbind(mean=m, sd=sds, Q, sims=ns)
    qc <- (2+seq(along=qs))
  } else {
    S <- cbind(mean=m, sd=sds, Q, "NA%"=NAS, sims=ns)
    qc <- (2+seq(along=qs))
  }
  #if (!is.null(nax)) dimnames(S)[[1]] <- nax
  
  attr(S, "quantile.cols") <- qc
  return(S)
}

.rvsummary.logical <- function (x)
{
  s <- sims(x)
  m <- colMeans(s, na.rm=TRUE)
  ns <- rvnsims(x)
  sds <- sqrt((colSums(s^2)-ns*(m^2))/(ns-1))
  sds[ns==1] <- 0
  NAS <- apply(s, 2, is.na)
  if (!is.null(dim(NAS))) NAS <- colMeans(NAS)*100
  nax <- if (is.null(dim(x))) NULL else names(x)
  if (all(NAS==0)) {
    S <- cbind(mean=m, sd=sds, sims=ns)
  } else {
    S <- cbind(mean=m, sd=sds, "NA%"=NAS, sims=ns)
  }
  #if (!is.null(nax)) dimnames(S)[[1]] <- nax
  attr(S, "quantile.cols") <- NULL
  return(S)
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
  s <- sims(x)
  if (is.logical(x)) {
    output <- .rvsummary.logical(x)
  } else if (is.integer(x)) {
    .range <- range(s, na.rm=TRUE)
    output <- .rvsummary.integer(x)
  } else {
    output <- .rvsummary.continuous(x)
  }
  row.names <- .dim.index(x)
  dimnames(output)[[1]] <- row.names
  Rhat <- rvRhat(x)
  if (!is.null(Rhat)) {
    if (!all(is.na(Rhat)))  output <- cbind(output, Rhat, n.eff=rvneff(x))
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
  if (is.null(x['sims.array']))
    stop("Argument does not contain a 'sims.array'")
  r <- rvsims(x$sims.array)
  rvattr(r, 'Rhat')  <- x$summary[,'Rhat']
  rvattr(r, 'n.eff') <- x$summary[,'n.eff']
  r
}

as.bugs.rv <- function (x, DIC=FALSE) # What happens when DIC=TRUE?
{
  require("R2WinBUGS")
  sims.array <- sims(x, mc.array=TRUE)
  parameter.names <- dimnames(x@chains[[1]])[[2]]
  parameters.to.save <- unique(sapply(strsplit(parameter.names, "\\["), "[", 1))
  d <- dim(sims.array)
  n.burnin     <- 0
  n.keep       <- d[1]
  n.chains     <- d[2]
  n.parameters <- d[3]
  n.sims       <- n.keep*n.chains
  n.iter       <- n.keep
  n.thin       <- 1
  #
  sims <- matrix(NA, n.sims, n.parameters)
  root.long <- character(n.parameters)
  indexes.long <- vector(n.parameters, mode = "list")
  for (i in 1:n.parameters) {
    temp <- R2WinBUGS:::decode.parameter.name(parameter.names[i])
    root.long[i] <- temp$root
    indexes.long[[i]] <- temp$indexes
  }
  n.roots <- length(parameters.to.save)
  left.bracket.short <- as.vector(regexpr("[[]", parameters.to.save))
  right.bracket.short <- as.vector(regexpr("[]]", parameters.to.save))
  root.short <- ifelse(left.bracket.short == -1, parameters.to.save, 
      substring(parameters.to.save, 1, left.bracket.short - 
          1))
  dimension.short <- rep(0, n.roots)
  indexes.short <- vector(n.roots, mode = "list")
  n.indexes.short <- vector(n.roots, mode = "list")
  long.short <- vector(n.roots, mode = "list")
  length.short <- numeric(n.roots)
  for (j in 1:n.roots) {
      long.short[[j]] <- (1:n.parameters)[root.long == root.short[j]]
      length.short[j] <- length(long.short[[j]])
      if (length.short[j] == 0) 
          stop(paste("parameter", root.short[[j]], "is not in the model"))
      else if (length.short[j] > 1) {
          dimension.short[j] <- length(indexes.long[[long.short[[j]][1]]])
          n.indexes.short[[j]] <- numeric(dimension.short[j])
          for (k in 1:dimension.short[j]) n.indexes.short[[j]][k] <- length(unique(unlist(lapply(indexes.long[long.short[[j]]], 
              .subset, k))))
          length.short[j] <- prod(n.indexes.short[[j]])
          if (length(long.short[[j]]) != length.short[j]) 
              stop(paste("error in parameter", root.short[[j]], 
                "in parameters.to.save"))
          indexes.short[[j]] <- as.list(numeric(length.short[j]))
          for (k in 1:length.short[j]) indexes.short[[j]][[k]] <- indexes.long[[long.short[[j]][k]]]
      }
  }
  rank.long <- unlist(long.short)
  # -----
  # yes, it's inefficient to do this, but for now I'm just letting this be as it is:
  for (k in 1:n.parameters) {
    sims[,k] <- as.vector(sims.array[,,k])
  }
  # ----
  dimnames(sims) <- list(NULL, parameter.names)
  summary <- R2WinBUGS:::monitor(sims.array, n.chains, keep.all = TRUE)
  last.values <- as.list(numeric(n.chains))
  for (i in 1:n.chains) {
    n.roots.0 <- if (DIC) 
        n.roots - 1
    else n.roots
    last.values[[i]] <- as.list(numeric(n.roots.0))
    names(last.values[[i]]) <- root.short[1:n.roots.0]
    for (j in 1:n.roots.0) {
        if (dimension.short[j] <= 1) {
            last.values[[i]][[j]] <- sims.array[n.keep, i, 
               long.short[[j]]]
            names(last.values[[i]][[j]]) <- NULL
        }
        else last.values[[i]][[j]] <- aperm(array(sims.array[n.keep, 
           i, long.short[[j]]], rev(n.indexes.short[[j]])), 
            dimension.short[j]:1)
    }
  }
  sims <- sims[sample(n.sims), ]
  sims.list <- summary.mean <- summary.sd <- summary.median <- vector(n.roots, 
      mode = "list")
  names(sims.list) <- names(summary.mean) <- names(summary.sd) <- names(summary.median) <- root.short
  for (j in 1:n.roots) {
    if (length.short[j] == 1) {
        sims.list[[j]] <- sims[, long.short[[j]]]
        summary.mean[[j]] <- summary[long.short[[j]], "mean"]
        summary.sd[[j]] <- summary[long.short[[j]], "sd"]
        summary.median[[j]] <- summary[long.short[[j]], "50%"]
    }
    else {
        temp2 <- dimension.short[j]:1
        sims.list[[j]] <- aperm(array(sims[, long.short[[j]]], 
            c(n.sims, rev(n.indexes.short[[j]]))), c(1, (dimension.short[j] + 
            1):2))
        summary.mean[[j]] <- aperm(array(summary[long.short[[j]], 
            "mean"], rev(n.indexes.short[[j]])), temp2)
        summary.sd[[j]] <- aperm(array(summary[long.short[[j]], 
            "sd"], rev(n.indexes.short[[j]])), temp2)
        summary.median[[j]] <- aperm(array(summary[long.short[[j]], 
            "50%"], rev(n.indexes.short[[j]])), temp2)
    }
  }
  summary <- summary[rank.long, ]
  all <- list(n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, 
        n.thin = n.thin, n.keep = n.keep, n.sims = n.sims, sims.array = sims.array[, 
            , rank.long, drop = FALSE], sims.list = sims.list, 
        sims.matrix = sims[, rank.long], summary = summary, mean = summary.mean, 
        sd = summary.sd, median = summary.median, root.short = root.short, 
        long.short = long.short, dimension.short = dimension.short, 
        indexes.short = indexes.short, last.values = last.values, is.DIC=DIC)
    if (DIC) {
        deviance <- all$sims.array[, , dim(sims.array)[3], drop = FALSE]
        dim(deviance) <- dim(deviance)[1:2]
        pD <- numeric(n.chains)
        DIC <- numeric(n.chains)
        for (i in 1:n.chains) {
            pD[i] <- var(deviance[, i])/2
            DIC[i] <- mean(deviance[, i]) + pD[i]
        }
        all <- c(all, list(pD = mean(pD), DIC = mean(DIC)))
    }
  class(all) <- "bugs"
  return(all)
}



# end of rv-bugs.R


# rv factors

is.numeric.rvfactor <- function (x)
{
  FALSE
}

is.rvfactor <- function (x)
{
  UseMethod("is.rvfactor")
}

is.rvfactor.rvfactor <- function (x)
{
  TRUE
} 

is.rvfactor.rv <- function (x)
{
  all(rvsimapply(x, is.factor))
} 

is.rvfactor.default <- function (x)
{
  FALSE
} 

as.rvfactor <- function (x, ...)
{
  if (is.rvfactor(x)) x else rvfactor(x)
} 


as.rv.rvfactor <- function (x)
{
  attr(x, "levels") <- NULL
  clx <- class(x)
  clx[clx=="rvfactor"] <- NULL
  class(x) <- clx
  x
}



rvfactor <- function (x, ...)
{
  UseMethod("rvfactor")
}

rvfactor.default <- function (x, levels=NULL, ...)
{
  f <- as.factor(x)
  a <- sims(as.rv(as.integer(f)))
  rvf <- rvsims(a)
  class(rvf) <- c("rvfactor", class(rvf))
  attr(rvf, "levels") <- if (is.null(levels)) levels(f) else levels
  rvf
}


rvfactor.rv <- function (x, levels=NULL, ...)
{
  a <- sims(x)
  f <- as.factor(a)
  levs <- levels(f)
  rvf <- rvsims(array(as.integer(f), dim(a)))
  class(rvf) <- c("rvfactor", class(rvf))
  attr(rvf, "levels") <- if (is.null(levels)) levs else levels
  rvf
}

print.rvfactor <- function(x, digits=NULL, ...)
{
  s <- summary(x)
  if (is.null(s))
    return(print(s))
  if (is.null(digits)) digits <- rvpar("print.digits")
  s <- as.data.frame(round(s,digits))
  ds <- dimnames(s)
  if (!is.null(.names <- names(x))) {
    s <- cbind(name=.names, s)
  } else if (!is.null(.dn <- dimnames(x))) {
    sud <- rvpar("summary.dimnames")
    if (!is.null(sud) && !is.na(sud) && is.logical(sud)) {
      da <- lapply(.dn, function (na) if (is.null(na)) rep("", nrow(s)) else na)
      rw <- da[[1]][row(x)]
      cl <- da[[2]][col(x)]
      s <- cbind(row=rw, col=cl, " "=":", s)
    }
  }
  print(s)
}

.rvsummary.rvfactor <- function (x, all.levels=FALSE)
{
  levels <- levels(x)
  llev <- length(levels)
  num.levels <- 1:llev
  #
  s <- sims(x)
  #mns <- colMeans(s, na.rm=TRUE)
  ns <- rvnsims(x)
  #sds <- sqrt((colSums(s^2)-ns*(mns^2))/(ns-1))
  #sds[ns==1] <- 0
  #
  maxlev <- if (is.null(maxlev <- rvpar("max.levels"))) { 10 } else maxlev
  too.many.levels.to.show <- if (all.levels) FALSE else (llev>maxlev)
  last.lev.no <- llev
  choose.levels <- if (too.many.levels.to.show) c(1:(maxlev-1), last.lev.no) else num.levels
  #
  a <- apply(s, 2, function (x) table(c(x, num.levels))) # ensure that all levels are tried
  if (is.null(dim(a))) {
    dim(a) <- c(ncol(s),length(choose.levels))
  }
  a <- (a-1) # And now subtract the extra counts from the matrix that was obtained.
  a <- a[choose.levels,,drop=FALSE]
  ns <- rvnsims(x)
  NAS <- apply(s, 2, is.na)
  if (!is.null(dim(NAS))) NAS <- colSums(NAS)
  nax <- if (is.null(dim(x))) NULL else names(x)
  m <- a
  if (too.many.levels.to.show) {
    S <- rbind(a[1:(maxlev-1),,drop=FALSE], "*"=0, a[maxlev,])
    rownames(S) <- c(levels[1:(maxlev-1)], "*", levels[last.lev.no])
    remaining.col <- maxlev
  } else {
    S <- a
    rownames(S) <- levels
    remaining.col <- NA
  }
  if (any(NAS>0)) S <- rbind("NA%"=NAS*100)
  remaining <- (ns-colSums(S))
  if (!is.na(remaining.col)) {
    S[remaining.col,] <- remaining
  } else if (any(remaining>0)) {
    stop("Impossible: levels won't sum up to 0")
  } 
  S <- t(S/ns)  # compute proportions in each category and transpose
  # S <- cbind(mean=mns, sd=sds, S, sims=ns)
  S <- cbind(S, sims=ns)
  attr(S, "quantile.cols") <- NULL
  return(S)
}


summary.rvfactor <- function(object, all.levels=FALSE, ...)
{
  x <- object
  x.len <- length(x)
  x.dim <- dim(x)
  if (is.null(x.dim))
    x.dim <- x.len
  if (x.len < 1) {
    return(NULL)
  }
  output <- .rvsummary.rvfactor(x, all.levels=all.levels)
  row.names <- .dim.index(x)
  dimnames(output)[[1]] <- row.names
  output
}



# MIXED:

is.rvmixed.rv <- function (x)
{
  (!is.null(attr(x, "pointmass")))
} 

is.rvmixed <- function (x)
{
  UseMethod("is.rvmixed")
} 

is.rvmixed.rvmixed <- function (x)
{
  TRUE
}

is.rvmixed.default <- function (x)
{
  FALSE
}

is.numeric.rvmixed <- function (x)
{
  TRUE
}

rvmixed <- function (x, pointmass=NULL)
{
  if (is.null(x)) {
    return(x)
  }
  if (is.na(pointmass) || !is.numeric(pointmass)) {
    stop("Point mass must be numeric")
  }
  attr(x, "pointmass") <- pointmass
  class(x) <- c("rvmixed", class(x))
  x
}

print.rvmixed <- function(x, digits=NULL, ...)
{
  s <- summary(x)
  if (is.null(s))
    return(print(s))
  if (is.null(digits)) digits <- rvpar("print.digits")
  s <- as.data.frame(round(s,digits))
  ds <- dimnames(s)
  if (!is.null(.names <- names(x))) {
    s <- cbind(name=.names, s)
  } else if (!is.null(.dn <- dimnames(x))) {
    sud <- rvpar("summary.dimnames")
    if (!is.null(sud) && !is.na(sud) && is.logical(sud)) {
      da <- lapply(.dn, function (na) if (is.null(na)) rep("", nrow(s)) else na)
      rw <- da[[1]][row(x)]
      cl <- da[[2]][col(x)]
      s <- cbind(row=rw, col=cl, " "=":", s)
    }
  }
  print(s)
}

as.rv.rvmixed <- function (x)
{
  attr(x, "pointmass") <- NULL
  clx <- class(x)
  clx[clx=="rvmixed"] <- NULL
  class(x) <- clx
  x
}

.rvsummary.rvmixed <- function (x, all.levels=FALSE)
{
  levels <- attr(x, "pointmass")
  llev <- length(levels)
  num.levels <- 1:llev
  maxlev <- if (is.null(maxlev <- rvpar("max.mixed.points"))) { 3 } else maxlev
  too.many.levels.to.show <- if (all.levels) FALSE else (llev>maxlev)
  last.lev.no <- llev
  choose.levels <- if (too.many.levels.to.show) c(1:(maxlev-1), last.lev.no) else num.levels
  s <- sims(x)
  a <- sapply(levels[choose.levels], function (lvl) colSums(s==lvl), simplify=TRUE)
  if (is.null(dim(a))) {
    dim(a) <- c(ncol(s),length(choose.levels))
  }
  if (too.many.levels.to.show) {
    S <- rbind(a[1:(maxlev-1),,drop=FALSE], "*"=0, a[maxlev,])
    rownames(S) <- c(levels[1:(maxlev-1)], "*", levels[last.lev.no])
  } else {
    S <- a
    colnames(S) <- levels
  }
  S <- (S/nrow(s)) # columns
  X <- summary.rv(x)
  qcmax <- max(attr(X, "quantile.cols"))
  S <- cbind(X[,1:qcmax,drop=FALSE], S, X[,-(1:qcmax),drop=FALSE])
  attr(S, "quantile.cols") <- attr(X, "quantile.cols")
  attr(S, "pointmass.cols") <- (qcmax+1:length(choose.levels))
  return(S)
}


summary.rvmixed <- function(object, all.levels=FALSE, ...)
{
  x <- object
  x.len <- length(x)
  x.dim <- dim(x)
  if (is.null(x.dim))
    x.dim <- x.len
  if (x.len < 1) {
    return(NULL)
  }
  output <- .rvsummary.rvmixed(x, all.levels=all.levels, ...)
  pmc <- attr(output, "pointmass.cols")
  row.names <- .dim.index(x)
  dimnames(output)[[1]] <- row.names
  attr(output, "pointmass.cols") <- pmc
  output
}

print.rvmixed <- function(x, digits=NULL, ...)
{
  #####  return(print.rv(as.rv(x), digits=digits, ...))
  s <- summary(x)
  if (is.null(s))
    return(print(s))
  pmc <- attr(s, "pointmass.cols")
  if (is.null(digits)) digits <- rvpar("print.digits")
  s <- as.data.frame(round(s,digits))
  ds <- dimnames(s)
  if (!is.null(pmc)) {
    s <- cbind(s[,1:(min(pmc)-1),drop=FALSE], " "=":", s[,pmc], s[,-(1:max(pmc)),drop=FALSE])
  }
  if (!is.null(.names <- names(x))) {
    s <- cbind(name=.names, s)
  } else if (!is.null(.dn <- dimnames(x))) {
    sud <- rvpar("summary.dimnames")
    if (!is.null(sud) && !is.na(sud) && is.logical(sud)) {
      da <- lapply(.dn, function (na) if (is.null(na)) rep("", nrow(s)) else na)
      rw <- da[[1]][row(x)]
      cl <- da[[2]][col(x)]
      s <- cbind(row=rw, col=cl, " "=":", s)
    }
  }
  print(s)
}




.onLoad <- function(libname,pkgname) # NOEXPORT
{
  if (is.null(setnsims())) setnsims(1000)
  setHook(packageEvent("rv", "detach"), 
    function (...) {
      .rvFunctionSwitch(attach=FALSE, verbose=FALSE)
    }
  )
  rvpar(
    oorm = TRUE,
    rvlwd = 2.0,
    rvcol = "default",
    rvlex = function (lwd) 1.5,
    rvinterval = c("mean", "50%", "95%"),
    point.sample = 400,
    line.sample = 20,
    summary.dimnames = TRUE,
    print.digits = 2
  )
  cat("Package rv loaded.\n")
}

.onAttach <- function(libname,pkgname) # NOEXPORT
{
  if (is.null(setnsims())) setnsims(200)
  .rvFunctionSwitch(attach=TRUE, verbose=FALSE)
}

