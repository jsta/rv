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



#' Horizontal interval plot of components of a random vector
#' 
#' \code{mlplot} plots the scalar components as of the given random array or
#' vector as horizontal intervals, grouped by row.
#' 
#' \code{mlplot} plots the scalar components of a vector or an array (2 or
#' 3-dimensional) vertically (up to down) so that a component of a vector or a
#' row of a matrix is plotted at vertical points 1...nrow(x).
#' 
#' An 'mlplot' of a vector implements a ``forest plot.''
#' 
#' Scalars on the same row are plotted closely together.  The positioning of
#' the scalars within a row are controlled by the arguments \code{y.center},
#' \code{y.shift}, \code{y.map}.  These do not need to be set for the default
#' plot; if two arrays or vectors are plotted over on top of each other (using
#' \code{add=TRUE}) then you should probably change \code{y.shift} which
#' controls the vertical position of the array elements.
#' 
#' See \code{demo(mlplot)} for a detailed
#' 
#' To change the color of the random components of the vector, use
#' \code{rvcol}. Typically this is of the same length as \code{X}, giving the
#' color `theme' for each component.
#' 
#' If \code{X} is a 3-dimensional array, \code{mlplot} is called repeatedly for
#' each 2-dimensional array \code{X[,,k]} for each \code{k}.
#' 
#' \code{X} may also be a fixed numeric object.
#' 
#' \code{NA}s (or random scalars with 100\% NA) are not plotted.
#' 
#' \code{mlplot} is still experimental.
#' 
#' @aliases mlplot mlplot.default mlplot.rvsummary
#' @param X a random array or vector
#' @param \dots further arguments passed to plot and points
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords classes
#' @examples
#' 
#' 
#' \dontrun{
#'   # You can run this complete example by typing demo("mlplot")
#' 
#'   n.rows <- 4; n.cols <- 5; n <- (n.rows*n.cols)
#'   # Draw some fixed numbers
#'   mu.true <- rnorm(1:n.rows, mean=1:n.rows, sd=1)
#'   sigma.true <- 1
#'   theta <- rvmatrix(rvnorm(n=n.cols, mean=mu.true, sd=sigma.true), nrow=n.rows)
#'   #
#'   col.labels <- paste("Time", 1:n.cols, sep=":")
#'   row.labels <- paste("Unit", 1:n.rows, sep=":")
#'   dimnames(theta) <- list(row.labels, col.labels)
#'   #
#'   par(mfrow=c(2,2))
#'   mlplot(theta, main="theta")
#'   abline(v=0, lty="dotted")
#'   mlplot(t(theta), main="theta transposed")
#'   abline(v=0, lty="dotted")
#'   row.sd <- apply.rv(theta, 1, sd.rv)
#'   col.sd <- apply.rv(theta, 2, sd.rv)
#'   x.max <- max(rvquantile(c(row.sd, col.sd), 0.99))
#'   mlplot(row.sd, xlim=c(0, x.max), main="theta: within-row sd for each unit")
#'   abline(v=0)
#'   mlplot(col.sd, xlim=c(0, x.max), main="theta: between-row sd for each time point")
#'   abline(v=0)
#' }
#' 
#' @export mlplot
mlplot <- function (X, ...)
{
  UseMethod("mlplot")
}

#' @rdname mlplot
#' @param y.center center the intervals nicely at each y-coordinate?
#' @param y.shift add this amount to each y coordinate of an interval
#' @param y.map optional function to compute the y-coordinates, given \code{X}
#' @param mar the margins of the plot
#' @param left.margin offset to add to the left margin of the plot (to add
#' space for the labels)
#' @param vline if numeric, plot vertical lines at these (horizontal)
#' coordinates
#' @param top.axis (logical) plot the top axis?
#' @param exp.labels (logical) if the original scale is logarithmic, label
#' ticks in original (exp) scale?
#' @param x.ticks positions for the ticks of the x-axis
#' @param axes (logical) plot the axes at all?
#' @param xlim x limits
#' @param ylim y limits
#' @param las the style of axis labels, see \code{\link{par}}
#' @param add (logical) add the intervals to an existing plot?
#' @param xlab x label
#' @param ylab not used (instead of labels, the row names are shown)
#' @importFrom graphics plot points
mlplot.default <- function (X, y.center = TRUE, y.shift = 0, y.map = NULL, mar = par("mar"), left.margin = 3, vline=NULL, top.axis = TRUE, exp.labels=FALSE, x.ticks = NULL, axes = NULL, xlim = NULL, ylim = NULL, xlab=deparse(substitute(X)), ylab=NULL, las = NULL, add = FALSE, ...) 
{
    if (missing(xlab)) {
      xlab <- deparse(substitute(X)) #?
    }
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
        x.sims <- sims(as.rvobj(X))
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
    if (! is.null(names(x.ticks))) {
      x.labels <- names(x.ticks)
    } else if (exp.labels) {
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
    } else {
      plot(X, y, ..., las = las, xlim = xlim, ylim = ylim, 
           axes = FALSE, xlab=xlab, ylab = "", type="n")
      if (is.null(axes) || axes) {
        axis(1, at = x.row.coords, labels=x.labels)
        if (top.axis) 
          axis(3, at = x.row.coords, labels=x.labels)
        if (is.null(labels)) 
          labels <- paste(y.row.coords)
        axis(2, at = y.row.coords, labels = labels, tick = FALSE, 
             line = FALSE, pos = NA, outer = FALSE, font = NA, 
             las = 1)
        if (! is.null(vline)) {
          if (! is.numeric(vline)) {
            stop("'vline' must be a numeric vector (or NULL if not used)")
          }
          if (is.null(names(vline))) {
            names(vline) <- rep("dotted", length(vline))
          }
          for (i in seq_along(vline)) {
            lty <- names(vline)[i]
            abline(v=vline[i], lty=lty, col="gray")
          }
        }
      }
      points(X, y, xlim = xlim, ylim = ylim, ...)
    }
    invisible(NULL)
}


mlplot_OLD_rvsummary <- function (X, y.center = TRUE, y.shift = 0, y.map = NULL, mar = par("mar"), left.margin = 3, top.axis = TRUE, exp.labels=FALSE, x.ticks = NULL, axes = NULL, xlim = NULL, ylim = NULL, xlab=deparse(substitute(X)), ylab=NULL, las = NULL, add = FALSE, ...) # NOEXPORT
{
    if (missing(xlab)) {
      xlab <- deparse(substitute(X)) #?
    }
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
        myrow <- function (x) row(array(NA, dim(x)))
        mycol <- function (x) col(array(NA, dim(x)))
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
            f(i = myrow(x), j = mycol(x), n = max(ncol(x),10)*1.5, nc = ncol(x))
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
        x.sims <- sims(X)
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


