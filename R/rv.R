#
# create and test objects of type 'rv'
#



#' Random Vectors
#' 
#' Creates or tests for objects of type ``\code{rv}''.
#' 
#' \code{rv} creates a random vector of the specified length.  The elements of
#' the vector are all equal to \code{NA}.
#' 
#' \code{is.rv} returns TRUE if its argument is a rv object, FALSE otherwise.
#' 
#' \code{as.rv} attempts to coerce its argument to the random vector (rv) type.
#' 
#' \code{is.random} returns \code{TRUE} or \code{FALSE} for each component of
#' the argument vector, depending on whether the component is a random variable
#' object.
#' 
#' \code{is.rvobj} tests whethe its argument object is either of class
#' \code{rv} or of class \code{rvsummary}.
#' 
#' \code{as.rvobj} coerces its argument object to \code{rv} unless the object
#' is an rv object (\code{is.rvobj(x)} is \code{TRUE}).
#' 
#' @aliases rv is.rv is.rv.rv is.rv.default as.rv as.rv.rv as.rv.numeric
#' as.rv.integer as.rv.logical as.rv.list as.rv.matrix as.rv.default is.random
#' as.rvobj is.rvobj
#' @param length desired length.
#' @param \dots further arguments passed to or from other methods.
#' @return An rv object of desired length, with the single simulation value
#' \code{NA}.
#' @note rv objects are internally lists with the class attribute set to
#' ``\code{rv}".  The number of simulations in rv objects is set by
#' \code{\link{setnsims}}.  This is by default set to 2500.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @seealso For a short version of the paper, view the vignette by
#' \code{vignette("rv")}.
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- rv(1)
#'   
#' 
#' @export rv
rv <- function(length=0) {
  if (is.numeric(length)) {
    x <- as.list(rep.int(NA, length))
  } else {
    stop("length must be numeric")
  }
  class(x) <- "rv"
  return(x)
}

#' @rdname rv
#' @param x object to be coerced or tested.
is.rv <- function(x)
{
  inherits(x, "rv")
}

# DEBUG:1 ?? The permutation of !is.null(dim(sims)) is not done yet
# DEBUG:2 (must take permutation of row indices and then permute.




#' Coercing Random Vectors to Real-valued
#' 
#' Coerces random vector objects into double-valued ones.
#' 
#' \code{as.double} coerces an rv object into double-valued one.  In effect,
#' the function \code{as.double} is applied to all simulations.
#' 
#' @param x an rv object
#' @param \dots other arguments
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- as.logical(rvbern(prob=0.5)) 
#'   print(x)
#'   print(as.double(x))
#' 
#' @export
as.double.rv <- function(x, ...)
{
  simapply(x, as.double, ...)
}



#' Logical Random vectors
#' 
#' Coerces a random variable to a logical-valued one (Bernoulli r.v.)
#' 
#' In effect, the function \code{as.logical} is applied to all simulations.
#' 
#' @param x an rv object
#' @param \dots Further arguments passed on
#' @note \code{is.logical(x)} returns \code{TRUE} if and only if \emph{each}
#' component of \code{x} is logical-valued (i.e. \code{TRUE}/\code{FALSE}).
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- rvbern(prob=0.5)   # some 0/1 valued random variable
#'   print(x)
#'   is.logical(x)           # FALSE, because by default x is 'double'
#'   x <- as.logical(x)      # coerce to logical; all zeros become FALSE, ones become TRUE
#'   is.logical(x)           # TRUE
#'   print(x)                # Shows the expectations and not the quantiles
#' 
#' @export
as.logical.rv <- function(x, ...)
{
  simapply(x, as.logical, ...)
}

#' Integer Random vectors
#' 
#' Coerces a random variable to an integer-valued (discrete) one
#' 
#' In effect, the function \code{as.integer} is applied to all simulations.
#' 
#' @param x an rv object
#' @param \dots Further arguments passed on
#' @note \code{is.integer(x)} returns \code{TRUE} if and only if \emph{each}
#' component of \code{x} is integer-valued (each simulation vector is of type
#' 'integer').
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @seealso \code{\link{as.logical.rv}}.
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- rvpois(lambda=3)   # some integer-valued random variable
#'   print(x)
#'   is.integer(x)           # FALSE, because by default x is 'double'!
#'   x <- as.integer(x)      # coerce to integer
#'   is.integer(x)           # TRUE
#'   print(x)                # Shows also the 'min' and 'max' columns
#' 
#' @export
as.integer.rv <- function (x, ...)
{
  simapply(x, as.integer, ...)
}

