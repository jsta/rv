
"[.rv" <- function (x, ..., drop = TRUE)
{
  cx <- class(x)
  X <- NextMethod()
  class(X) <- cx
  return(X)
}

"[.rvsummary" <- function (x, ..., drop = TRUE)
{
  q <- attr(x, "quantiles")
  cx <- class(x)
  x <- NextMethod()
  class(x) <- cx
  attr(x, "quantiles") <- q
  return(x)
}


"[<-.rvsummary" <- function (x, i, j, ..., value = NULL)
{
  cx <- class(x)
  q <- attr(x, "quantiles")
  value <- as.rvsummary(value, quantiles=q)
  X <- .Primitive("[<-")(unclass(x), i, j, ..., value=value)
  class(X) <- cx
  return(X)
}

"[<-.rv" <- function (x, ..., value = NULL)
{
  cx <- class(x)
  value <- as.rvobj(value)
  x <- unclass(x)
  X <- .Primitive("[<-")(x, ..., value=value)
  class(X) <- cx
  return(X)
}

