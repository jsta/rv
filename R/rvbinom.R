

rvbinom <- function (n=1, size, prob) {
  rvvapply(rbinom, n.=n, size=size, prob=prob)
}

