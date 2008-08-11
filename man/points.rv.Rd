\name{points.rv}
\title{Add Points and Intervals to a Plot}
\alias{points.rv}
\usage{

\method{points}{rv}(x, y = NULL, type = "p", xlim = NULL, ylim = NULL, 
                 rvlwd = rvpar("rvlwd"), rvcol = rvpar("rvcol"), 
                 rvpoint = rvpar("rvpoint"), rvlex = rvpar("rvlex"), ...)
}
\description{
  Draw a sequence of points or uncertainty intervals at specified
  (fixed) x-coordinates.
}
\arguments{
  \item{x}{x-coordinates}
  \item{y}{y-coordinates}
  \item{type}{character indicating the type of plotting}
  \item{rvcol}{colors for the intervals}
  \item{xlim}{x-limits (optional)}
  \item{ylim}{y-limits (optional)}
%%% \item{pch}{Character code for plotting distribution means}
%%%  \item{out.of.range.marker}{(logical) mark intervals/points that are drawn out of the limits with an 'x' --- experimental, not yet used}
  \item{rvlwd}{line width of the thin interval}
  \item{rvpoint}{character vector of length 3, indicating intervals (points) to print}
  \item{rvlex}{factor to multiply \code{rvlwd} with, to get the thicker interval}
  \item{\dots}{further arguments passed to \code{points}}
}
\details{
  Each 'point' with a fixed coordinate and a random coordinate is
  plotted as an interval. 
  If "lines" are plotted (\code{type="l"} or \code{type="b"}),
  the result is a random draw of lines connecting the coordinates.
  See \code{\link{lines.rv}} for details on how to set the sample size of the random draw.

  Each interval consists of a maximum of three components.
  (1) a dot (2) thick interval (3) thin interval. 
  Typically the dot marks the mean or the median; the 
  thin and the thick intervals show a shorter and a longer middle uncertainty interval.
  The appearance of these intervals can be controlled using the parameters
  \code{rvlwd}, \code{rvpoint}, \code{rvcol}, and \code{rvlex}.


  \code{rvlwd} sets the line width of the thin interval;
  \code{rvlex} sets the factor to multiply \code{rvlwd} to get the line width
  of the thicker interval.

  \code{points} attempts to color the intervals and the dot using
  the color given as \code{rvcol}. The basic name of the color should be
  given, e.g. \code{"red"} or \code{"blue"}. 
  The thin line is colored using 
  the basic color, the thick line is colored using a darker hue 
  (numbered '2', e.g. \code{"red2"})
  and the dot is colored using the darkest hue (numbered '3', e.g. \code{"red3"}). 
  That is, for example. if \code{rvcol='red'}, the color scheme generated for the 
  dot, the thick line, and the thin line, respectively, are
  \code{c('red3', 'red2', 'red')}. 

  Special color themes: 
  the default \code{rvcol} color scheme is called \code{"default"} and
  yields the color scheme \code{c("grey20", "grey40", "grey60")}.
  Other special color themes: \code{"grey"}, \code{"lightgrey"}, \code{"darkgrey"}.
  (The spellings 'gray' and 'grey' are interchangeable).
  

  The parameter \code{rvpoint} is a character vector of length 3,
   with the first component
  indicating what to plot as a dot (possible values: "mean", "median"),
  the second component indicating what to plot as a "thick interval" 
  (possible values: "n%" such as "50\%" or "80\%"), and
  the second component indicating what to plot as a "thin interval".
  Default: \code{c("mean", "50\%", "95\%")}. 
  If you wish only to plot the mean and the 95\% interval, use 
  \code{rvpoint=c("mean", NA, "95\%")} or
  \code{rvpoint=c("mean", "95\%", NA)}. 

  The color \code{col} is used for plotting fully fixed dots (both x and y coordinates fixed)
  and lines (fixed and \emph{random lines} -- see \code{\link{lines.rv}}).

 NOTE. This parameterization is yet experimental, and may change.

 It is possible to have both \code{x} and \code{y} random, 
 but this code is not yet fully functional. 
}
\references{
Kerman, J. and Gelman, A. (2007). Manipulating and Summarizing
Posterior Simulations Using Random Variable Objects.
Statistics and Computing 17:3, 235-244.

See also \code{vignette("rv")}.
}
\author{
 Jouni Kerman
 \email{jouni@kerman.com}
}
\examples{
  x <- 1:10
  y <- rvnorm(mean=x)
  par(mfrow=c(2,2))
  plot(x, y, main="Fixed x-coordinate")
  plot(y, x, main="Fixed y-coordinate")
  plot(x, y, lwd=4, rvcol="red", main="Color and line width changed")
  plot(x, y, type="b", main="Intervals and random lines", rvcol="blue", col="gray")
  \dontrun{
    # Don't use the rv-only parameters when plotting fixed vectors.
    plot(x, E(y), rvcol="blue", col="gray")
    # This will work if x is an rv:
    plot(as.rv(x), E(y), rvcol="blue", col="gray")
  }
}

\keyword{aplot}


