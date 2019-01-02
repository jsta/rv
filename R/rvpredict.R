#' Generate predictions from models
#' 
#' Performs predictions (in the form of \code{rv} objects) from models based on
#' given covariates.
#' 
#' The \code{lm} method generates predictions of the outcome variable.  The
#' posterior coefficient estimates (the ``intercept'' and the ``betas'') are
#' estimated Bayesianly by \code{posterior(object)}; the coefficients are
#' multiplied by \code{newdata} (if omitted, the model covariate matrix is used
#' instead) to obtain the predicted model mean; lastly, the outcomes are
#' predicted from the Normal sampling model, taking into account the sampling
#' variability along with the uncertainty in the estimation of the standard
#' deviation (`sigma').
#' 
#' The covariate matrix \code{newdata} can be an \code{rv}, representing
#' additional uncertainty in the covariates.
#' 
#' @aliases rvpredict rvpredict.lm
#' @param object An object representing a statistical model fit.
#' @param newdata A data frame with new covariates to be used in the
#' predictions. The column names of the data frame must match those in the
#' model matrix (although order may be arbitrary).  If omitted, the model
#' matrix is used instead; the resulting predictions are then the
#' \emph{replications} of the data. \emph{Note:} this can be an \code{rv}
#' object to incorporate extra uncertainty into predictions.
#' @param \dots Arguments passed to and from other methods.
#' @return For the \code{lm} method, a vector as long as there are rows in the
#' data frame \code{newdata}.
#' @author J Kerman
#' @keywords models
#' @examples
#' 
#'   ## Create some fake data
#'   n <- 10
#'   ## Some covariates
#'   set.seed(1)
#'   X <- data.frame(x1=rnorm(n, mean=0), x2=rpois(n, 10) - 10)
#'   y.mean <- (1.0 + 2.0 * X$x1 + 3.0 * X$x2)
#'   y <- rnorm(n, y.mean, sd=1.5) ## n random numbers
#'   D <- cbind(data.frame(y=y), X)
#'   ## Regression model fit
#'   obj <- lm(y ~ x1 + x2, data=D)
#'   ## Bayesian estimates
#'   posterior(obj)
#'   ## Replications
#'   y.rep <- rvpredict(obj)
#'   ## Predictions at the mean of the covariates
#'   X.pred <- data.frame(x1=mean(X$x1), x2=mean(X$x2))
#'   y.pred <- rvpredict(obj, newdata=X.pred)
#'   ## Plot predictions
#'   plot.rv(D$x1, y.rep)
#'   points(D$x1, D$y, col="red")
#'   ## `Perturb' (add uncertainty to) covariate x1
#'   X.pred2 <- X
#'   X.pred2$x1 <- rnorm(n=n, mean=X.pred2$x1, sd=sd(X.pred2$x1))
#'   y.pred2 <- rvpredict(obj, newdata=X.pred2)
#' 
#' @export rvpredict
rvpredict <- function (object, ...) {
  UseMethod("rvpredict")
}

.X.and.offset <- function (object, newdata) {
  offset <- object$offset
  tt <- terms(object)
  if (missing(newdata) || is.null(newdata)) {
    X <- model.matrix(object)
    offset <- object$offset
  } else {
    Terms <- delete.response(tt)
    m <- rvmodel.frame(Terms, newdata, na.action = na.fail, xlev = object$xlevels)
    ##
    if (!is.null(cl <- attr(Terms, "dataClasses"))) {
      stats::.checkMFClasses(cl, m)
    }
    X <- rvmodel.matrix(Terms, m, contrasts.arg = object$contrasts)
    offset <- rep(0, nrow(X))
    if (!is.null(off.num <- attr(tt, "offset"))) {
      for (i in off.num) {
        offset <- offset + eval(attr(tt, "variables")[[i + 1]], newdata)
      }
      if (!is.null(object$call$offset)) {
        offset <- offset + eval(object$call$offset, newdata)
      }
    }
  }
  list(X=X, offset=offset)
}

rvpredict.lm <- function (object, newdata, ...) {
  L <- .X.and.offset(object, newdata)
  X <- L$X
  ## end code from predict.lm
  post <- posterior(object)
  X.names <- colnames(X)
  all.coeff <- post$beta
  if (! all(ok <- (X.names %in% names(all.coeff)))) {
    stop("Some names are missing: ", paste(X.names[! ok], collapse=", "))
  }
  beta <- all.coeff[X.names]
  mu.pred <- (X %**% beta) ## rv-compatible multiplication
  mu.pred <- (mu.pred + L$offset)
  sigma <- post$sigma
  y.pred <- rvnorm(mean=mu.pred, sd=sigma)
  return(y.pred)
}

rvmodel.frame <- function (formula, data=NULL, ...) {
  if (is.null(data)) {
    return(model.frame(formula, data=data, ...))
  }
  L <- as.list(data)
  if (! any(sapply(L, is.rv))) {
    m <- match.call()
    m[[1L]] <- as.name("model.frame")
    return(eval(m, envir=parent.frame()))
  }
  terms <- NULL
  .mf <- function (data.matrix, formula, ...) {
    df <- as.data.frame(data.matrix)
    mf <- model.frame(formula, data=df, ...)
    terms <<- attr(mf, "terms")
    return(mf)
  }
  W <- lapply(L, as.rv)
  W <- do.call(cbind.rv, W)
  colnames(W) <- names(L)
  x <- simapply(W, .mf, formula=formula, ...)
  attr(x, "terms") <- terms
  return(x)
}

rvmodel.matrix <- function(object, data, ...) {
  L <- as.list(data)
  if (! any(sapply(L, is.rv))) {
    m <- match.call()
    m[[1L]] <- as.name("model.matrix")
    return(eval(m, envir=parent.frame()))
  }
  terms <- attr(data, "terms")
  mm.names <- NULL
  .mm <- function (data.matrix, object, ...) {
    df <- structure(as.data.frame(data.matrix), terms=terms)
    mm <- model.matrix(object, data=df, ...)
    mm.names <<- colnames(mm)
    return(mm)
  }
  W <- lapply(L, as.rv)
  W <- do.call(cbind.rv, W)
  colnames(W) <- names(L)
  x <- simapply(W, .mm, object=object, ...)
  colnames(x) <- mm.names
  return(x)
}

##
