

# ========================================================================
# rvhyper  -  hypergeometric rvs
# ========================================================================

#' @importFrom stats rhyper
rvhyper <- function(nn=1, m, n, k) {
  warning("NOT YET READY")
  attr(nn, "n.name") <- "nn"
  rvvapply(rhyper, n.=nn, m=m, n=n, k=k)
}

