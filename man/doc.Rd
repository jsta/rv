\name{doc}
\title{Show a document (pdf) file from the package directory}
\alias{doc}
\usage{
doc(topic, package = NULL, lib.loc=NULL)
}
\description{
  Show a document file (a pdf file) from the package "doc" directory 
}
\arguments{
  \item{topic}{a character string giving the (base) name of the pdf file to view. If omitted, all available pdf files are listed.}
  \item{package}{a character vector with the names of packages to search through, or NULL in which case all available packages in the library trees specified by lib.loc are searched.}
  \item{lib.loc}{a character vector of directory names of R libraries, or NULL. The default value of NULL corresponds to all libraries currently known.}
}
\details{
  If "package" is omitted, all packages are scanned for the file name.
}
\references{ 
Kerman, Jouni and Gelman, Andrew. Manipulating and Summarizing Posterior Simulations 
Using Random Variable Objects. Technical report, Columbia University, New York.
}
\author{
 Jouni Kerman
 \email{kerman@stat.columbia.edu}
 \url{http://www.stat.columbia.edu/~kerman}
}
\examples{
  \dontrun{doc("rv")}
  \dontrun{doc("R2WinBUGS")}
}

\keyword{classes}

