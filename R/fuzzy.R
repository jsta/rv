
#' @export
is.fuzzy <- function (x) {
  UseMethod("is.fuzzy")
}

#' @method is.fuzzy rv
#' @export
is.fuzzy.rv <- function (x) {
  # NAME
  #  is.fuzzy - Is a Vector Component Logical But Random
  # 
  component.is.logical <- rvsimapply(x, is.logical)
  component.prop <- rvmean(x)
  (component.is.logical & component.prop>0 & component.prop<1)
}

#' @method is.fuzzy default
#' @export
is.fuzzy.default <- function (x)
{
  return(FALSE)
}
