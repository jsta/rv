

rvgamma <- function (n=1, shape, rate = 1, scale = 1/rate)  {
  rvvapply(rgamma, n.=n, shape=shape, scale=scale)
}


