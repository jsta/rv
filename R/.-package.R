

#' Concatenation of random vectors
#' 
#' Concatentates random vectors.
#' 
#' NOTE: \code{recursive} has not yet been tested.
#' 
#' \code{cc} is a function that works for both non-rv and other vectors. To
#' make code compatible for both constant vectors and rv objects, one can use
#' \code{cc} instead of \code{c}.
#' 
#' @aliases c.rv c.rvsummary cc
#' @param \dots objects to be concatenated. Can be a mixture of constants and
#' rv objects.
#' @param recursive logical. If recursive = TRUE, the function recursively
#' descends through lists (and pairlists) combining all their elements into a
#' vector.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- rvnorm(2)
#'   y <- rvbern(2, prob=0.5)
#'   z <- c(x, y)
#'   print(z)
#'   z1 <- cc(1, z)
#'   z2 <- c(as.rv(1), z)
#'   z3 <- c(as.rv(1), z)
#'   print(z1)
#'   print(z2)
#'   print(z3)
#' 
NULL





#' Constant Vectors
#' 
#' Functions to coerce or test for non-random objects.
#' 
#' \code{is.constant} returns \code{TRUE} for each component of the argument
#' object if there is only one simulation (that is, the variable is
#' ``constant'').
#' 
#' Note: rv objects that merely have variance zero are not therefore
#' necessarily ``true'' constants.
#' 
#' \code{as.constant} coerces \code{rv} or \code{rvsummary} objects into
#' constant strings; \code{NA} is returned for component that is not random.
#' 
#' @aliases is.constant as.constant as.constant.rv as.constant.rvsummary
#' @param x an object, random variable (rv) or not
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   is.constant(1)         # TRUE
#'   is.constant(as.rv(1))  # TRUE
#'   setnsims(200)
#'   x <- rvbern(prob=0.001)
#'   all(sims(x)==0)        # most probably true
#'   is.constant(x)         # always FALSE
#'   x <- rvnorm(3)
#'   x[1] <- 1
#'   as.constant(x)         # 1, NA, NA
#'   all(is.random(x) & is.na(as.constant(x))) # always TRUE
#' 
NULL





#' Distributions of Various Statistics of a Random Vector and Array
#' 
#' \code{var.rv}, \code{cov.rv} and \code{cor.rv} compute the distribution of
#' the variance statistic of x and the distribution of the covariance statistic
#' or the correlation statistic of x and y if these are vectors.  If x and y
#' are matrices then the covariances (or correlations) between the columns of x
#' and the columns of y are computed.
#' 
#' 
#' These functions are compatible with \emph{both} numeric and rv objects.  To
#' make your code compatible with \code{rv} objects, use e.g. \code{sd.rv}
#' instead of \code{sd}.
#' 
#' The functions \code{cov.rv} is implemented by applying the corresponding
#' numerical function to the rows of the simulation matrices of \code{x} and
#' \code{y} and forming a new \code{rv} object from the resulting vector of
#' simulations.  Alternatively \code{x} may be a random matrix (and \code{y}
#' \code{NULL}).  %Then the numerical function \code{cov}.
#' 
#' \code{cor.rv} works similarly, but returns the distribution of the
#' correlation statistic (i.e. function).
#' 
#' \code{var.rv} computes the distribution of the variance statistic.
#' \code{sd.rv} is the square root of the result obtained by \code{var.rv}.
#' 
#' @aliases cov.rv cor.rv var.rv sd.rv
#' @param x a numeric or random vector, matrix, or a data frame
#' @param y \code{NULL} (default) or a vector, matrix or data frame with
#' compatible dimensions to x. The default is equivalent to y = x (but more
#' efficient).
#' @param \dots further arguments passed to the corresponding numeric functions
#' @return A random vector or array.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords internal
#' @examples
#' 
#'   #
#' 
NULL





#' Fuzziness
#' 
#' Tests whether an object is ``fuzzy'', i.e.  a logical random scalar that has
#' probability strictly between zero and one (not strictly true nor strictly
#' false).
#' 
#' 
#' @aliases fuzzy is.fuzzy is.fuzzy.rv
#' @param x an object, random or constant
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- as.logical(rvbern(1,0.4)) # a logical random variable
#'   is.fuzzy(x) # TRUE, since x is logical and not constant
#'   is.fuzzy(x<2) # FALSE, since x is less than 2 with probability one
#'   is.fuzzy(rvnorm(1)) # FALSE, since it's not a probability
#'   is.fuzzy(TRUE) # FALSE, since TRUE is strictly TRUE
#'   is.fuzzy(1) # FALSE, since 1 is not a logical variable
#' 
NULL





#' Random Matrix Multiplication
#' 
#' Multiplies two random matrices, if they are conformable.  If one argument is
#' a vector, it will be coerced to either a row or column matrix to make the
#' two arguments conformable.  If both are vectors it will return the inner
#' product.
#' 
#' Optimized internally for the case of random matrix multiplied by a constant
#' column vector.
#' 
#' @aliases %*%.rv %**% matmult.rv
#' @param x,y numeric or complex matrices or vectors.
#' @return The (distribution of the) matrix product.  Use \code{\link{drop}} to
#' get rid of dimensions which have only one level.
#' @seealso \code{\link{matrix}}, \code{\link{Ops}}, \code{\link{diag}}.
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords array arith
#' @examples
#' 
#' x <- 1:4
#' (z <- x %*% x)    # scalar ("inner") product (1 x 1 matrix)
#' drop(z)             # as scalar
#' 
#' y <- diag(x)
#' z <- matrix(1:12, ncol = 3, nrow = 4)
#' y %*% z
#' y %*% x
#' x %*% z
#' 
NULL





#' Numeric Random Vectors
#' 
#' Creates or coerces rv objects of type "numeric".
#' 
#' \code{is.numeric(x)} returns \code{TRUE} if and only if \emph{each}
#' component of \code{x} is numeric-valued.
#' 
#' \code{as.numeric.rv} coerces an rv object into numeric-valued one.  In
#' effect, the function \code{as.numeric} is applied to all simulations.
#' 
#' Random factors are not numeric (just as non-random factors aren't).
#' 
#' @aliases numeric.rv as.numeric.rv is.numeric.rv is.numeric.rvfactor
#' as.numeric.rvfactor
#' @param x an rv object to be coerced or tested.
#' @param \dots further arguments passed to or from other methods.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @seealso \code{\link{numeric}}.
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- as.logical(rvbern(1,0.5)) # Bernoulli rv
#'   is.numeric(x)           # FALSE
#'   x <- as.numeric(x)      # coerce to numeric; all TRUEs become ones, FALSEs zeros
#'   is.numeric(x)           # TRUE
#' 
NULL

#' Simulation-based Random Variable Objects
#' 
#' `\code{rv}' implements a simulation-based random variable object class.
#' 
#' Please refer to the vignette: \code{vignette("rv")} for details.
#' 
#' \tabular{ll}{ Package: \tab rv \cr Version: \tab 2.3.0 \cr Date: \tab
#' 2013-05-18 \cr Namespace: \tab rv \cr Depends: \tab R(>= 2.10.0), methods,
#' utils, grDevices, graphics \cr License: \tab GPL-2 \cr }
#' 
#' @name rv-package
#' @docType package
#' @author Jouni Kerman \email{jouni@@kerman.com} Package built on Sat May 18
#' 22:47:25 CEST 2013
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
NULL

#' Random Vector Summaries
#' 
#' \code{rvsummary} is a class of objects that hold the summary information on
#' each scalar component of a random variable (quantiles, mean, sd, number of
#' simulations etc.)
#' 
#' The \code{rvsummary} class provides a means to store a concise
#' representation of the marginal posterior distributions of the vector
#' components.  By default, the 201 quantiles \preformatted{ 0, 0.005, 0.01,
#' 0.015, ..., 0.990, 0.995, 1 } are saved for each vector component in an
#' \code{rvsummary} object.
#' 
#' \code{is.rvsummary} tests whether the object is an \code{rvsummary} object;
#' \code{as.rvsummary} coerces a random vector object to a \code{rvsummary}
#' object.
#' 
#' \code{as.data.frame} is another way to obtain the data frame that is
#' produced by the \code{summary} method.
#' 
#' A data frame that has the format of an \code{rv} summary can be coerced into
#' an \code{rvsummary}; if quantiles are not specified within the data frame,
#' quantiles from the Normal distribution are filled in, if the mean and s.d.
#' are given.
#' 
#' Therefore, the following (generic) functions work with \code{rvsummary}
#' objects: \code{rvmean}, \code{rvsd}, \code{rvvar}, \code{rvquantile},
#' \code{rnsims}, \code{sims}, and consequently any `rv-only' function that
#' depends only on these functions will work; e.g. \code{is.constant}, which
#' depends only on \code{rvnsims}.
#' 
#' The method \code{is.double} is provided for compatibility reasons; this is
#' needed in a function called by \code{plot.rvsummary}
#' 
#' The arithmetic operators and mathematical functions will not work with
#' \code{rvsummary} objects.
#' 
#' The \code{sims} method returns the quantiles.
#' 
#' @aliases rvsummary is.rvsummary as.rvsummary as.rvsummary.rv
#' as.rvsummary.rvsummary as.rvsummary.data.frame as.data.frame.rvsummary
#' print.rvsummary print.rvsummary_rvfactor as.double.rvsummary
#' @param x object to be coerced or tested
#' @param quantiles quantiles to calculate and store in the object
#' @param digits integer; how many digits to round the numbers to
#' @param all.levels logical; whether to print all levels or not (see below for
#' details)
#' @param \dots further arguments passed to or from other methods.
#' @return An object of class \code{rvsummary} \emph{and} of subclass
#' \code{rvsummary_numeric}, \code{rvsummary_integer},
#' \code{rvsummary_logical}, or \code{rvsummary_rvfactor}.
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @seealso \code{\link{rvfactor}}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#'   x <- rvnorm(mean=1:12)
#'   sx <- as.rvsummary(x)
#'   print(sx)          # prints the summary of the rvsummary object
#'   length(sx)         # 12
#'   dim(sx)            # NULL
#'   dim(sx) <- c(3,4)  #   
#'   dimnames(sx) <- list(1:3, 1:4)
#'   names(sx) <- 1:12  # 
#'   print(sx)          # prints the names and dimnames as well  
#' 
NULL



