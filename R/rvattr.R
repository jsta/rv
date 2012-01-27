# ========================================================================
# rvattr  -  return the rvsim attribute for each component
# ========================================================================
# returns a list.

rvattr <- function(x, attrib=NULL)
{
  a <- lapply(unclass(x), "attr", "rvsim")
  if (!is.null(attrib)) {
    a <- lapply(a, "[[", attrib)
    nulls <- sapply(a, is.null)
    a[nulls] <- NA
  }
  return(a)
}

# rvattr(x, 'name') <- vector of values for each component or x; or,
# rvattr(x) <- list of attributes: list(attributename=vector of values, attributename2=...)
# e.g. list(Rhat=c(1.01,1.03,1.23), n.eff=c(200,100,16)).
#



"rvattr<-" <- function(x, attrib=NULL, value)
{
  value <- as.list(value)
  if (is.null(attrib)) {
    if (length(names(value))>0) {
      for (a in names(value)) {
        rvattr(x, attrib=a) <- value[[a]]
      }
    } else {
      stop("Cannot match value to rv")
    }
  } else if (length(names(value))>0 && length(names(x))>0) {
    remaining <- rep(TRUE, length(x))
    for (i in seq_along(value)) {
      name <- names(value)[i]
      w <- which(names(x)==name)
      if (!any(remaining[w])) {
        stop("Multiple assignment")
      }
      if (length(w)==1) {
        A <- attr(x[[w]], "rvsim")
        if (is.null(A)) { A <- list() }
        A[[attrib]] <- value[[i]]
        attr(x[[w]], "rvsim") <- A
        remaining[w] <- FALSE
      } else if (length(w)>1) {
        stop("Ambiguous name: '", name, "'")
      } else {
        stop("No name ", name, " in rv.")
      }
    }
  } else if (length(x)!=length(value)) {
    stop("Value must be of the same length as x or the matching must be done by name.")
  } else {
    for (i in seq(along=x)) {
      A <- attr(x[[i]], "rvsim")
      if (is.null(A)) { A <- list() }
      A[[attrib]] <- value[[i]]
      attr(x[[i]], "rvsim") <- A
    }
  }
  return(x)
}


