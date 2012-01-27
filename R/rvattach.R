
rvattach <- function (what, name=deparse(substitute(what)), overwrite=TRUE, impute=FALSE, ...)
{
  if (missing(what) || is.null(what)) {
    return(invisible(NULL))
  }
  if (!is.list(what)) {
    stop("Argument must be a list")
  }
  while (any(search()==name)) {
    eval(substitute(detach(name), list(name=as.name(name))))
  }
  if (!is.rv(what)) {
    stop("Argument must be an rv")
  }
  a <- splitbyname(what)
  rattach(a, name=name, overwrite=overwrite, impute=impute, ...)
}

rattach <- function (what, name=deparse(substitute(what)), overwrite=TRUE, impute=FALSE, ...)
{
  if (missing(what) || is.null(what)) {
    return(invisible(NULL))
  }
  if (!is.list(what)) {
    stop("Argument must be a list")
  }
  while (any(search()==name)) {
    eval(substitute(detach(name), list(name=as.name(name))))
  }
  if (is.null(what)) return(NULL)
  if (overwrite) {
    ls.GE <- ls(.GlobalEnv)
    conf_names <- names(what)[names(what)%in%ls.GE]
    for (obj_name in conf_names) {
      if (impute) {
        v <- get(obj_name, .GlobalEnv) # Value to be imputed to
        what[[obj_name]] <- .impute.by.name(v, what[[obj_name]])
      }
      remove (list=obj_name, envir=.GlobalEnv)
    }
  }
  attach (what, name=name, ...)
}
