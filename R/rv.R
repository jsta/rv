
#
# create and test objects of type 'rv'
#

rv <- function(length=0) {
  if (is.numeric(length)) {
    x <- as.list(rep.int(NA, length))
  } else {
    stop("length must be numeric")
  }
  class(x) <- "rv"
  return(x)
}

is.rv <- function(x)
{
  inherits(x, "rv")
}

as.rv <- function(x, ...)
{
  UseMethod("as.rv")
}

as.rv.rv <- function(x, ...)
{
  return(x)
}

as.rv.numeric <- function(x, ...)
{
  if (is.rv(x)) {
    return(x)
  }
  r <- rvsims(matrix(as.vector(x), nrow=1))
  cr <- class(r)
  attributes(r) <- attributes(x)
  class(r) <- cr
  return(r)
}

as.rv.logical <- function(x, ...)
{
  as.rv.numeric(x)
}

as.rv.integer <- function(x, ...)
{
  as.rv.numeric(x)
}

as.rv.list <- function(x, ...)
{
  stop("Cannot coerce an arbitrary list to an rv object")
}

as.rv.matrix <- function(x, ...)
{
  as.rv.numeric(x)
}

as.rv.default <- function(x, ...)
{
  if (is.null(x)) return(NULL)
  stop('Cannot coerce object of class "', paste(class(x), collapse="/"), '" to rv')
}

as.rv.xtabs <- function (x, ...) 
{
  # NAME
  #  as.rv.xtabs
  #
  as.rv(x[])
}

# DEBUG:1 ?? The permutation of !is.null(dim(sims)) is not done yet
# DEBUG:2 (must take permutation of row indices and then permute.


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

