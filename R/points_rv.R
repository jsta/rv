#' Add Points and Intervals to a Plot
#' 
#' Draw a sequence of points or uncertainty intervals at specified (fixed)
#' x-coordinates.
#' 
#' Each point with a fixed coordinate and a random coordinate is plotted as
#' an interval.  If lines are plotted (\code{type="l"} or \code{type="b"}),
#' the result is a random draw of lines connecting the coordinates.  See
#' \code{\link{lines.rv}} for details on how to set the sample size of the
#' random draw.
#' 
#' Each interval consists of a maximum of three components.  (1) a dot (2)
#' thick interval (3) thin interval.  Typically the dot marks the mean or the
#' median; the thin and the thick intervals show a shorter and a longer middle
#' uncertainty interval.  The appearance of these intervals can be controlled
#' using the parameters \code{rvlwd}, \code{rvpoint}, \code{rvcol}, and
#' \code{rvlex}.
#' 
#' \code{rvlwd} sets the line width of the thin interval; \code{rvlex} sets the
#' factor to multiply \code{rvlwd} to get the line width of the thicker
#' interval.
#' 
#' \code{points} attempts to color the intervals and the dot using the color
#' given as \code{rvcol}. The basic name of the color should be given, e.g.
#' \code{"red"} or \code{"blue"}.  The thin line is colored using the basic
#' color, the thick line is colored using a darker hue (numbered '2', e.g.
#' \code{"red2"}) and the dot is colored using the darkest hue (numbered '3',
#' e.g. \code{"red3"}).  That is, for example. if \code{rvcol='red'}, the color
#' scheme generated for the dot, the thick line, and the thin line,
#' respectively, are \code{c('red3', 'red2', 'red')}.
#' 
#' Special color themes: the default \code{rvcol} color scheme is called
#' \code{"default"} and yields the color scheme \code{c("grey20", "grey40", "grey60")}.  Other special color themes: \code{"grey"}, \code{"lightgrey"},
#' \code{"darkgrey"}.  (The spellings 'gray' and 'grey' are interchangeable).
#' 
#' The parameter \code{rvpoint} is a character vector of length 3, with the
#' first component indicating what to plot as a dot (possible values: "mean",
#' "median"), the second component indicating what to plot as a "thick
#' interval" (possible values: "n\%" such as "50\%" or "80\%"), and the second
#' component indicating what to plot as a "thin interval".  Default:
#' \code{c("mean", "50\%", "95\%")}.  If you wish only to plot the mean and the
#' 95\% interval, use \code{rvpoint=c("mean", NA, "95\%")} or
#' \code{rvpoint=c("mean", "95\%", NA)}.
#' 
#' The color \code{col} is used for plotting fully fixed dots (both x and y
#' coordinates fixed) and lines (fixed and \emph{random lines} -- see
#' \code{\link{lines.rv}}).
#' 
#' NOTE. This parameterization is yet experimental, and may change.
#' 
#' It is possible to have both \code{x} and \code{y} random, but this code is
#' not yet fully functional.
#' 
#' @param x x-coordinates
#' @param y y-coordinates
#' @param type character indicating the type of plotting
#' @param rvcol colors for the intervals
#' @param xlim x-limits (optional)
#' @param ylim y-limits (optional)
#' @param rvlwd line width of the thin interval
#' @param rvpoint character vector of length 3, indicating intervals (points)
#' to print
#' @param rvlex factor to multiply \code{rvlwd} with, to get the thicker
#' interval
#' @param \dots further arguments passed to \code{points}
#' @author Jouni Kerman \email{jouni@@kerman.com}
#' @references Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
#' Posterior Simulations Using Random Variable Objects. Statistics and
#' Computing 17:3, 235-244.
#' 
#' See also \code{vignette("rv")}.
#' @keywords aplot
#' @examples
#' 
#'   x <- as.rv(1:10)
#'   y <- rvnorm(mean=x)
#'   par(mfrow=c(2,2))
#'   plot(x, y, main="Fixed x-coordinate")
#'   plot(y, x, main="Fixed y-coordinate")
#'   plot(x, y, lwd=4, rvcol="red", main="Color and line width changed")
#'   plot(x, y, type="b", main="Intervals and random lines", rvcol="blue", col="gray")
#'   \dontrun{
#'     # Don't use the rv-only parameters when plotting fixed vectors.
#'     plot(x, E(y), rvcol="blue", col="gray")
#'     plot(x, E(y), rvcol="blue", col="gray")
#'   }
#' 
#' @method points rv
points.rv <- function (x, y = NULL, type = "p", xlim = NULL, ylim = NULL, rvlwd = rvpar("rvlwd"), rvcol = rvpar("rvcol"), rvpoint = rvpar("rvpoint"), rvlex = rvpar("rvlex"), ...) {
  if (! (is.rvobj(x) || is.rvobj(y))) {
    return(points(x, y, type=type, xlim=xlim, ylim=ylim, ...))
  }
  xy <- .xy.coords.rv(x, y)
  x <- as.rvobj(xy$x)
  y <- as.rvobj(xy$y)
  arg <- list(...)
  draw.points <- (type == "p" || type == "b")
  draw.lines <- (type == "l" || type == "b")
  point.sample <- rvpar("point.sample")
  if (is.null(point.sample)) 
    point.sample <- NA
  line.sample <- rvpar("line.sample")
  if (is.null(line.sample)) 
    line.sample <- NA
  if (is.null(rvlex) || !is.function(rvlex)) {
    rvlex <- function(lwd) 1.5
  }
  if (is.null(rvcol) || is.na(rvcol)) {
    rvcol <- "default"
  }
  x.rv <- (is.random(x) & !rv.all.na(x))
  y.rv <- (is.random(y) & !rv.all.na(y))
  x.point <- (!x.rv)
  y.point <- (!y.rv)
  vertical.pair <- (x.point & y.rv)
  horizontal.pair <- (x.rv & y.point)
  iv.pair <- (vertical.pair | horizontal.pair)
  rv.pair <- (x.rv & y.rv)
  point.pair <- (x.point & y.point)
  segs <- NULL
  pts  <- NULL
  cols <- NULL
  lwds <- NULL
  qx <- rvintervals(x, rvpoint)
  qy <- rvintervals(y, rvpoint)
  rvcol <- rep(rvcol, length.out = length(x))
  cols <- t(sapply(rvcol, rvcolortheme))
  dimnames(cols) <- list(NULL, rvpoint)
  cols <- cbind(default = if (is.null(arg$col)) 
                "black"
  else arg$col, cols)
  rvlwd <- rep(rvlwd, length.out = length(x))
  lwds <- t(sapply(rvlwd, function(wd) c(0, wd * rvlex(wd), 
                                         wd)))
  dimnames(lwds) <- list(NULL, rvpoint)
  lwds <- cbind(default = if (is.null(arg$lwd)) 
        "black"
  else arg$lwd, lwds)
  if (any(point.pair)) {
    x.pts <- rvmean(x[point.pair])
    y.pts <- rvmean(y[point.pair])
    pts <- list(x = c(pts$x, x.pts), y = c(pts$y, y.pts), 
                col = c(pts$col, cols[point.pair, 1]), lwd = c(pts$lwd, 
                                                         lwds[point.pair, 1]))
  }
  if (any(rv.pair)) {
    ## THIS CODE IS NOT YET FULLY FUNCTIONAL
    xy.pts <- rvsample(c(x[rv.pair], y[rv.pair]), jointly = TRUE, 
                       size = point.sample)
    x.pts <- xy.pts[, 1]
    y.pts <- xy.pts[, 2]
    pchs <- rep(19, length.out = length(x.pts))
    cl <- rep(cols[rv.pair, 1], each = length(x.pts)/sum(rv.pair))
    pts <- list(x = c(pts$x, x.pts), y = c(pts$y, y.pts), 
                col = c(pts$col, cl), pch = c(pts$pch, pchs))
  }
  for (name in names(qx)) {
    if (is.na(name)) 
      next
    xiv <- qx[[name]]
    yiv <- qy[[name]]
    is.seg <- (nrow(xiv) == 2 || nrow(yiv) == 2)
    if (is.seg) {
      if (any(iv.pair)) {
        segs <- list(x0 = c(segs$x0, xiv[1, iv.pair]), 
                     y0 = c(segs$y0, yiv[1, iv.pair]), x1 = c(segs$x1, 
                                                         xiv[2, iv.pair]), y1 = c(segs$y1, yiv[2, 
                                                                             iv.pair]), col = c(segs$col, cols[iv.pair, 
                                                                                          name]), lwd = c(segs$lwd, lwds[iv.pair, name]))
      }
    } else {
      pchs <- rep(if (is.null(arg$pch)) 19 else arg$pch, 
                  length.out = length(x))
      pts <- list(x = c(pts$x, xiv[1, iv.pair]), y = c(pts$y, 
                                                   yiv[1, iv.pair]), col = c(pts$col, cols[iv.pair, 
                                                                       name]), pch = c(pts$pch, pchs))
    }
  }
  if (draw.points) {
    if (!is.null(segs)) 
      do.call("segments", args = .nodups(c(arg, segs)))
    if (!is.null(pts)) 
      do.call("points", args = .nodups(c(arg, pts)))
    if (any(rv.pair)) {
      do.call("points", args = .nodups(c(arg, pts)))
    }
  }
  if (draw.lines) {
    lns <- .rvjointdrawxy(x, y, size = line.sample)
    lns$x <- rbind(lns$x, NA)
    lns$y <- rbind(lns$y, NA)
    do.call("lines", args = .nodups(c(arg, lns)))
  }
  invisible(NULL)
}


points.rvsummary <- points.rv

.rvjointdrawxy <- function (x, y, size=1, reject.na=TRUE)
{
  xy <- c(x, y)
  s <- rvsample(xy, size=size, jointly=TRUE, reject.na=reject.na)
  if (is.null(dim(s))) s <- t(s)
  xs <- t(s[,seq(along=x)])
  ys <- t(s[,-seq(along=x)])
  list(x=xs, y=ys)
}

