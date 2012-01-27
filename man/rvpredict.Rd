\name{rvpredict}
\alias{rvpredict}
\alias{rvpredict.lm}
\title{
  Generate predictions from models
}
\description{
  Performs predictions (in the form of \code{rv} objects) from models
  based on given covariates. 
}
\usage{
  rvpredict(object, ...)
  \method{rvpredict}{lm}(object, newdata, ...)
}
\arguments{
  \item{object}{An object representing a statistical model fit.}
  \item{newdata}{A data frame with new covariates to be used in the
    predictions. The column names of the data frame must match those in
    the model matrix (although order may be arbitrary). 
    If omitted, the model matrix is used instead; the resulting
    predictions are then the \emph{replications} of the
    data. \emph{Note:} this can be an \code{rv} object to incorporate
    extra uncertainty into predictions.}
  \item{\dots}{ Arguments passed to and from other methods. }
}
\details{
  The \code{lm} method generates predictions of the outcome variable.
  The posterior coefficient estimates (the ``intercept'' and the
  ``betas'') are estimated Bayesianly by
  \code{posterior(object)}; the coefficients are multiplied by
  \code{newdata} (if omitted, the model covariate matrix is used
  instead) to obtain the predicted model mean; lastly, the outcomes are
  predicted from the Normal sampling model, taking into account the
  sampling variability along with the uncertainty in the estimation of
  the standard deviation (`sigma').

  The covariate matrix \code{newdata} can be an \code{rv}, representing
  additional uncertainty in the covariates. 
}
\value{
  For the \code{lm} method, a vector as long as there are rows in the
  data frame \code{newdata}.
}
\author{
  J Kerman
}
\examples{
  ## Create some fake data
  n <- 10
  ## Some covariates
  set.seed(1)
  X <- data.frame(x1=rnorm(n, mean=0), x2=rpois(n, 10) - 10)
  y.mean <- (1.0 + 2.0 * X$x1 + 3.0 * X$x2)
  y <- rnorm(n, y.mean, sd=1.5) ## n random numbers
  D <- cbind(data.frame(y=y), X)
  ## Regression model fit
  obj <- lm(y ~ x1 + x2, data=D)
  ## Bayesian estimates
  posterior(obj)
  ## Replications
  y.rep <- rvpredict(obj)
  ## Predictions at the mean of the covariates
  X.pred <- data.frame(x1=mean(X$x1), x2=mean(X$x2))
  y.pred <- rvpredict(obj, newdata=X.pred)
  ## Plot predictions
  plot.rv(D$x1, y.rep)
  points(D$x1, D$y, col="red")
  ## `Perturb' (add uncertainty to) covariate x1
  X.pred2 <- X
  X.pred2$x1 <- rnorm(n=n, mean=X.pred2$x1, sd=sd(X.pred2$x1))
  y.pred2 <- rvpredict(obj, newdata=X.pred2)
}  
\keyword{ models }
