# ================================================================================
# rv - simulation-based random variable class in R
# Version 1.0
# Updated 2008-08-05
# (c) 2004-2008 Jouni Kerman <jouni@kerman.com>
# ================================================================================
#
# rv package (c) J Kerman 2008
#

.RV <- list()

.rvRegisterFunctionSwitch <- function (FUN, fname, namespace)
{
  origFUN <- getFromNamespace(fname, namespace)
  ns.name <- paste(namespace, fname, sep=":::")
  environment(FUN) <- environment(rv:::rv)
  .RV[[ns.name]] <<- list(R=origFUN, rv=FUN)
}


.rvFunctionSwitch <- function (attach=TRUE, name=NULL, verbose=FALSE)
{
  detach <- (!attach)
  if (is.null(name)) name <- names(.RV)
  RF <- .RV[name]
  full.names <- strsplit(name, ":::")
  for (i in seq(along=full.names)) {
    f <- RF[[i]]
    ns.name <- full.names[[i]]
    namespace <- ns.name[1]
    fun.name <- ns.name[2]
    FUN <- if (attach) f$rv else f$R
    if (verbose) {
       cat(name[i], " ")
       if (attach) cat("replaced by a function in package:rv\n")
         else cat("restored.\n")
     }
    assignInNamespace(fun.name, FUN, namespace)
  }
}



# ----------------
# end of 01FunctionSwitch.R
# ----------------


"[.rv" <- function (x, ..., drop = TRUE)
{
  cx <- class(x)
  X <- NextMethod()
  class(X) <- cx
  return(X)
}

"[<-.rv" <- function (x, ..., value = NULL)
{
  cx <- class(x)
  value <- as.rv(value)
  X <- .Primitive("[<-")(unclass(x), ..., value=value)
  class(X) <- cx
  return(X)
}




# rv-Math.R - standard math functions for the rv class

Math.rv <- function(x, ...) {
  # Componentwise operation
  X <- x # Preserve class and other attributes
  for (i in seq_along(x)) {
    x <- X[[i]] 
    X[[i]] <- NextMethod()
  }
  return(X)
}

# cumsum, cumprod, cummax, cummin

cumsum.rv <- function (x)
{
  simapply(x, cumsum)
}

cumprod.rv <- function (x)
{
  simapply(x, cumprod)
}

cummin.rv <- function (x)
{
  simapply(x, cummin)
}

cummax.rv <- function (x)
{
  simapply(x, cummax)
}


# ----------------
# end of rv-Math.R
# ----------------

# rv-Ops.R - standard math functions for the rv class
# TODO
#  


Ops.rv <- function(e1, e2=NULL)
{
  e1.attr <- attributes(e1)
  e2.attr <- attributes(e2)
  if (is.null(e2)) {
    v <- mapply(.Generic, 0, e1, SIMPLIFY=FALSE, USE.NAMES=FALSE)
  } else {
    v <- mapply(.Generic, e1, e2, SIMPLIFY=FALSE, USE.NAMES=FALSE)
  }
  if (length(v)==length(e1)) {
    attributes(v) <- e1.attr
  } else if (length(v)==length(e2)) {
    attributes(v) <- e2.attr
  }
  if (is.list(v)) class(v) <- class(rv())
  return(v)
}

"!.rv" <- function(e1) 
{
  v <- simapply(e1, .Generic)
  dim(v) <- dim(e1)
  names(v) <- names(e1)
  return(v)
}

# ---------------
# end of rv-Ops.R
# ---------------



Summary.rv <- function(..., na.rm=FALSE)
{
  S <- sims(c(...)) # an L x n matrix of numbers
  M <- apply(S, 1, .Generic, na.rm=na.rm)
  rvsims(M)
}
# ========================================================================
# abline.rv  -  a + bx
# ========================================================================
# change this to
#   mapply.rv("abline", ..., MoreArgs=NULL)
#

.abline.default <- graphics:::abline

abline.rv <- function (a = NULL, b = NULL, h = NULL, v = NULL, ...)
{
  #if (anyisrv(a,b,h,v)) {
    line.sample <- rvpar("line.sample")
    if (!is.numeric(line.sample)) stop("rvpar('line.sample') is not a number")
    args <- list(FUN=.abline.default, a=a, b=b, h=h, v=v)
    nulls <- sapply(args, is.null)
    nullArgs <- names(nulls)[nulls]
    MoreArgs <- list(...)
    args[nullArgs] <- NULL
    MoreArgs[nullArgs] <- list(NULL)
    args$SAMPLESIZE <- line.sample
    args$MoreArgs <- MoreArgs
    do.call(rvmapply, args=args)
  #} else {
  #  .abline.default(a=a, b=b, h=h, v=v, ...)
  #}
  invisible()
}



anyisrv <- function (...) # NOEXPORT
{
  any(unlist(lapply(list(...), is.rv)))
}






aperm.rv <- function (x, ...) # NOMETHOD
{
  # NAME
  #   aperm.rv - Transpose a Random Array
  #   
  a <- aperm(x, ...)
  if (!is.rv(x)) {
    return(a)
  }
  class(a) <- class(x)
  return(a)
}


apply.rv <- function (X, MARGIN, FUN, ...) # NOMETHOD
{
  # NAME
  #   apply.rv - Apply Functions Over Random Array Margins
  FUNC <- function (x, ...) {
    class(x) <- class(rv())
    FUN(x, ...)
  }
  a <- apply(X, MARGIN, FUNC, ...)
  .list2array(a, drop=TRUE)
}



as.bugs.rv <- function (x, DIC=FALSE)
{
#
# Quickly patched from code by A Gelman & others
#
# DEBUG: What happens when DIC=TRUE?
#
  require("R2WinBUGS")
  sims.array <- sims(x, mc.array=TRUE)
  parameter.names <- dimnames(x@chains[[1]])[[2]]
  parameters.to.save <- unique(sapply(strsplit(parameter.names, "\\["), "[", 1))
  d <- dim(sims.array)
  n.burnin     <- 0
  n.keep       <- d[1]
  n.chains     <- d[2]
  n.parameters <- d[3]
  n.sims       <- (n.keep*n.chains)
  n.iter       <- n.keep
  n.thin       <- 1
  #
  sims <- matrix(NA, n.sims, n.parameters)
  root.long <- character(n.parameters)
  indexes.long <- vector(n.parameters, mode = "list")
  for (i in 1:n.parameters) {
    temp <- R2WinBUGS:::decode.parameter.name(parameter.names[i])
    root.long[i] <- temp$root
    indexes.long[[i]] <- temp$indexes
  }
  n.roots <- length(parameters.to.save)
  left.bracket.short <- as.vector(regexpr("[[]", parameters.to.save))
  right.bracket.short <- as.vector(regexpr("[]]", parameters.to.save))
  root.short <- ifelse(left.bracket.short == -1, parameters.to.save, 
      substring(parameters.to.save, 1, left.bracket.short - 
          1))
  dimension.short <- rep(0, n.roots)
  indexes.short <- vector(n.roots, mode = "list")
  n.indexes.short <- vector(n.roots, mode = "list")
  long.short <- vector(n.roots, mode = "list")
  length.short <- numeric(n.roots)
  for (j in 1:n.roots) {
      long.short[[j]] <- (1:n.parameters)[root.long == root.short[j]]
      length.short[j] <- length(long.short[[j]])
      if (length.short[j] == 0) 
          stop(paste("parameter", root.short[[j]], "is not in the model"))
      else if (length.short[j] > 1) {
          dimension.short[j] <- length(indexes.long[[long.short[[j]][1]]])
          n.indexes.short[[j]] <- numeric(dimension.short[j])
          for (k in 1:dimension.short[j]) n.indexes.short[[j]][k] <- length(unique(unlist(lapply(indexes.long[long.short[[j]]], 
              .subset, k))))
          length.short[j] <- prod(n.indexes.short[[j]])
          if (length(long.short[[j]]) != length.short[j]) 
              stop(paste("error in parameter", root.short[[j]], 
                "in parameters.to.save"))
          indexes.short[[j]] <- as.list(numeric(length.short[j]))
          for (k in 1:length.short[j]) indexes.short[[j]][[k]] <- indexes.long[[long.short[[j]][k]]]
      }
  }
  rank.long <- unlist(long.short)
  # -----
  # yes, it's inefficient to do this, but for now I'm just letting this be as it is:
  for (k in 1:n.parameters) {
    sims[,k] <- as.vector(sims.array[,,k])
  }
  # ----
  dimnames(sims) <- list(NULL, parameter.names)
  summary <- R2WinBUGS:::monitor(sims.array, n.chains, keep.all = TRUE)
  last.values <- as.list(numeric(n.chains))
  for (i in 1:n.chains) {
    n.roots.0 <- if (DIC) 
        n.roots - 1
    else n.roots
    last.values[[i]] <- as.list(numeric(n.roots.0))
    names(last.values[[i]]) <- root.short[1:n.roots.0]
    for (j in 1:n.roots.0) {
        if (dimension.short[j] <= 1) {
            last.values[[i]][[j]] <- sims.array[n.keep, i, 
               long.short[[j]]]
            names(last.values[[i]][[j]]) <- NULL
        }
        else last.values[[i]][[j]] <- aperm(array(sims.array[n.keep, 
           i, long.short[[j]]], rev(n.indexes.short[[j]])), 
            dimension.short[j]:1)
    }
  }
  sims <- sims[sample(n.sims), ]
  sims.list <- summary.mean <- summary.sd <- summary.median <- vector(n.roots, 
      mode = "list")
  names(sims.list) <- names(summary.mean) <- names(summary.sd) <- names(summary.median) <- root.short
  for (j in 1:n.roots) {
    if (length.short[j] == 1) {
        sims.list[[j]] <- sims[, long.short[[j]]]
        summary.mean[[j]] <- summary[long.short[[j]], "mean"]
        summary.sd[[j]] <- summary[long.short[[j]], "sd"]
        summary.median[[j]] <- summary[long.short[[j]], "50%"]
    }
    else {
        temp2 <- dimension.short[j]:1
        sims.list[[j]] <- aperm(array(sims[, long.short[[j]]], 
            c(n.sims, rev(n.indexes.short[[j]]))), c(1, (dimension.short[j] + 
            1):2))
        summary.mean[[j]] <- aperm(array(summary[long.short[[j]], 
            "mean"], rev(n.indexes.short[[j]])), temp2)
        summary.sd[[j]] <- aperm(array(summary[long.short[[j]], 
            "sd"], rev(n.indexes.short[[j]])), temp2)
        summary.median[[j]] <- aperm(array(summary[long.short[[j]], 
            "50%"], rev(n.indexes.short[[j]])), temp2)
    }
  }
  summary <- summary[rank.long, ]
  all <- list(n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin, 
        n.thin = n.thin, n.keep = n.keep, n.sims = n.sims, sims.array = sims.array[, 
            , rank.long, drop = FALSE], sims.list = sims.list, 
        sims.matrix = sims[, rank.long], summary = summary, mean = summary.mean, 
        sd = summary.sd, median = summary.median, root.short = root.short, 
        long.short = long.short, dimension.short = dimension.short, 
        indexes.short = indexes.short, last.values = last.values, is.DIC=DIC)
    if (DIC) {
        deviance <- all$sims.array[, , dim(sims.array)[3], drop = FALSE]
        dim(deviance) <- dim(deviance)[1:2]
        pD <- numeric(n.chains)
        DIC <- numeric(n.chains)
        for (i in 1:n.chains) {
            pD[i] <- var(deviance[, i])/2
            DIC[i] <- mean(deviance[, i]) + pD[i]
        }
        all <- c(all, list(pD = mean(pD), DIC = mean(DIC)))
    }
  class(all) <- "bugs"
  return(all)
}



as.list.rv <- function (x, ...)
{
  L <- list()
  for (i in seq_along(x)) {
    L[[i]] <- x[i]
  }
  return(L)
}

as.rvsummary.bugs <- function (x, list.=TRUE, ...) # EXPORT
{
  if (list.) {
    lapply(as.rv(x, list.=TRUE, ...), as.rvsummary)
  } else {
    as.rvsummary(as.rv(x, list.=FALSE, ...))
  }
}

as.rv.bugs <- function (x, list.=TRUE, ...) # EXPORT
{
  # as.rv.bugs - transform a R2WinBUGS object into a random variable object
  # n.chains,
  # n.iter
  # n.burnin
  # n.thin
  # n.keep
  # n.sims
  # sims.array=sims.array[,,rank.long]
  # sims.list
  # sims.matrix=sims[,rank.long]
  # summary
  # mean=summary.mean
  # sd=summary.sd
  # median=summary.median
  # root.short
  # long.short
  # dimension.short
  # indexes.short
  # last.values
  #     sims <- c(bugs.sims(parameters.to.save, n.chains, n.iter, 
  #          n.burnin, n.thin, DIC), model.file = model.file, 
  #          is.DIC = DIC)
  #      class(sims) <- "bugs"
  #      return(sims)
  if (is.null(x['sims.array']))
    stop("Argument does not contain a 'sims.array'")
  A <- x$sims.array
  dnA <- dimnames(A)
  dim(A) <- c(prod(dim(A)[1:2]), dim(A)[3])
  dimnames(A) <- list(dnA[[1]], dnA[[3]])
  r <- rvsims(A)
  rvattr(r, 'Rhat')  <- x$summary[,'Rhat']
  rvattr(r, 'n.eff') <- x$summary[,'n.eff']
  if (list.) {
    return(splitbyname(r))
  } else {
    return(r)
  }
}



as.vector.rv <- function(x, mode="any")
{
  a <- attributes(x) 
  x <- lapply(unclass(x), as.vector, mode=mode)
  attributes(x) <- a
  dim(x) <- NULL
  return(x)
}




c.rv <- function(..., recursive=FALSE)
{
  # a kludge to disable dispatching rv
  x <- .c(list(NA), ..., recursive=recursive)[-1]
  class(x) <- class(rv())
  return(x)
}


# ========================================================================
# cbind  -  column bind for rvs
# ========================================================================
# Note: It's inconvenient that we cannot call the generic (default) cbind
#       if class attributes are set.
#

# DEBUG: cbind(1, rvnorm(1), rvnorm(1)) causes an error msg
#   In cbind(v, unclass(x[[i]]), deparse.level = deparse.level) :
#     number of rows of result is not a multiple of vector length (arg 2)

cbind.rv <- function(..., deparse.level = 1)
{
  if (deparse.level != 1) 
    .NotYetUsed("deparse.level != 1")
  x <- list(...)
  if (length(x)<1) return(NULL)
  v <- NULL
  for (i in seq(along=x)) {
    v <- cbind(v, unclass(x[[i]]), deparse.level=deparse.level)
  }
  class(v) <- class(rv())
  return(v)
}


# ========================================================================
# rvbind.rv  -  row bind for rvs
# ========================================================================
#

rbind.rv <- function(..., deparse.level = 1)
{
  if (deparse.level != 1) 
    .NotYetUsed("deparse.level != 1")
  x <- list(...)
  if (length(x)<1) return(NULL)
  v <- NULL
  for (i in seq(along=x)) {
    v <- rbind(v, unclass(x[[i]]), deparse.level=deparse.level)
  }
  class(v) <- class(rv())
  return(v)
}



# Value:  a logical vector (not rv), TRUE if a component is constant w.p. 1
#

is.constant <- function(x)
{
  # Note: this corresponds to "constant with probability 1", while
  # rvvar(x)==0 would correspond to "constant almost surely"
  return(rvnsims(x)==1)
}

as.constant <- function(x)
{
  UseMethod('as.constant')
}

as.constant.rv <- function(x)
{
  unlist(lapply(x, as.constant), use.names=FALSE)
}

# ========================================================================
# cor  -  short description
# ========================================================================

cor.default <- getFromNamespace('cor', 'stats') # EXPORT cor.default

formals(cor.default) <- c(formals(cor.default), alist(... = ))

cor.rv <- function(x, y=NULL, ...)  ## EXPORT cor.rv
{
  if (!is.matrix(x)) {
    if (is.null(y)) {
      stop("supply both x and y or a matrix-like x")
    }
    x <- as.vector(x)
  }
  rvmapply(cor.default, x=x, y=y, ...)
}



# ========================================================================
# cov  -  short description
# ========================================================================

cov.default <- getFromNamespace('cov', 'stats') # EXPORT cov.default

formals(cov.default) <- c(formals(cov.default), alist(... = ))

cov.rv <- function(x, y=NULL, ...)  ## EXPORT cov.rv
{
  if (!is.matrix(x)) {
    if (is.null(y)) {
      stop("supply both x and y or a matrix-like x")
    }
    x <- as.vector(x)
  }
  rvmapply(cov.default, x=x, y=y, ...)
}





detachrv <- function ()
{
  detach("package:rv")
  unloadNamespace("rv")
}



finiterange <- function (x) ## NOEXPORT
{
  x <- x[!is.na(x)]
  x <- x[is.finite(x)]
  if (length(x)==0) return(c(NA,NA))
  range(x)  
}


is.fuzzy <- function (x)
{
  UseMethod("is.fuzzy")
}

is.fuzzy.rv <- function (x)
{
  # NAME
  #  is.fuzzy - Is a Vector Component Logical But Random
  # 
  component.is.logical <- rvsimapply(x, is.logical)
  component.prop <- rvmean(x)
  (component.is.logical & component.prop>0 & component.prop<1)
}

is.fuzzy.default <- function (x)
{
  return(FALSE)
}
# ========================================================================
# hist.rv  -  histogram, adapted for rv's
# ========================================================================

hist.rv <- function(x, grid=c(4,5), xlim=x.range, main=paste(xname,"simulation"), freq=FALSE, ...)
{
  par(mfrow=grid)
  #  l <- length(x)
  #  par(mfrow=c(l %% 3 + l %/% 3, 3))
  xname <- deparse(substitute(x))
  grid.n <- grid[1]*grid[2]
  s <- sims(x)
  x.range <- c(min(s),max(s))
  if (grid.n<1)
    stop("Bad grid")
  s <- s[1:grid.n,]
  for (i in 1:grid.n)
  {
    # truehist(s[i,], xlim=xlim, main, ...)
    hist(s[i,], xlim=xlim, main=main, freq=freq, ...)
  }
}



is.na.rv <- function(x)
{
  simapply(x, is.na)
}


# ========================================================================
# lines.rv  -  plot some random lines
# ========================================================================
# btw, "rvlines" does not make sense - that'd be something like a weird histogram

lines.rv <- function(x, y, type="l", ...)
{
  points.rv(x, y, type="l", ...)
}

# ========================================================================
# %*% - matrix product
# ========================================================================
# DEBUG: how to make this work with outer()?
#

"%*%.rv" <- function(x, y)
{
  if (!is.rv(x) && !is.rv(y)) return(.Primitive("%*%")(x, y))
  d <- dim(y)
  if (!is.rv(x) && (is.null(d)) || (length(d)==2 && d[2]==1)) {
    n.sims <- .Primitive("max")(rvnsims(x), rvnsims(y), na.rm=FALSE)
    ysim <- sims(as.rv(y), dimensions=TRUE, n.sims=n.sims)
    # Typical case: constant matrix times a rv vector
    AB <- t(.Primitive("%*%")(x, t(ysim)))
    rvsims(AB)
  } else {
    rvmapply(base:::crossprod, x=t(as.rv(x)), y=as.rv(y))
  }
}




# ========================================================================
# rvRhat; rvneff; rvnchains - convenience functions for some attributes
# ========================================================================
#

rvRhat <- function(x)
{
  unlist(rvattr(x, "Rhat"))
}

rvneff <- function(x)
{
  unlist(rvattr(x, "n.eff"))
}

rvnchains <- function(x)
{
  unlist(rvattr(x, "n.chains"))
}



mean.rv <- function(x, ...)
{
  rvsims(rowMeans(sims(x), ...))
}



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

mlplot <- function (X, ...)
{
  UseMethod("mlplot")
}

mlplot.default <- function (X, y.center = TRUE, y.shift = 0, y.map = NULL, mar = par("mar"), left.margin = 3, top.axis = TRUE, exp.labels=FALSE, x.ticks = NULL, axes = NULL, xlim = NULL, ylim = NULL, xlab=deparse(substitute(X)), ylab=NULL, las = NULL, add = FALSE, ...) 
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
        x.sims <- sims(as.rv(X))
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


mlplot.rvsummary <- function (X, y.center = TRUE, y.shift = 0, y.map = NULL, mar = par("mar"), left.margin = 3, top.axis = TRUE, exp.labels=FALSE, x.ticks = NULL, axes = NULL, xlim = NULL, ylim = NULL, xlab=deparse(substitute(X)), ylab=NULL, las = NULL, add = FALSE, ...) 
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



rvnsims <- function(x)
{
  UseMethod("rvnsims")
}

rvnsims.rv <- function(x) # NOEXPORT
{
  sapply(unclass(x), length)
}

rvnsims.rvsummary <- function (x) # NOEXPORT
{
  return(x$n.sims)
}

rvnsims.default <- function (x) # NOEXPORT
{
  if (!(is.atomic(x) || is.recursive(x))) {
    stop("rvnsims: no applicable method for class '", class(x), "'")
  }
  rep.int(1, length(x))
}

setnsims <- function(n.sims)
{
# setnsims - get or set the default number of simulations (a global variable)
  if (length(n.sims)>0 && is.numeric(n.sims) && (!is.na(n.sims[1])) && n.sims[1]>=2) {
    n.sims <- as.integer(ceiling(n.sims[1]))
    oldn.sims <- rvpar("n.sims")
    rvpar(n.sims=n.sims)
  } else {
    stop('Invalid number of simulations (must be at least 2)', n.sims[1])
  }
  return(oldn.sims)
}

getnsims <- function ()
{
  n.sims <- rvpar("n.sims")
  if (!is.integer(n.sims) || n.sims<2) {
    stop("Invalid number of simulations - please set with setnsims(...)")
  }
  return(n.sims)
}




is.numeric.rv <- function (x)
{
  all(rvsimapply(x, is.numeric))
}

as.numeric.rv <- function (x, ...)
{
  simapply(x, as.numeric, ...)
}

as.numeric.rvfactor <- function (x, ...)
{
  simapply(x, as.numeric, ...)
}


outer.rv <- function (X, Y=NULL, FUN="*", ...)
{
  # NAME
  #   outer.rv - 
  # 
  if (is.null(Y)) {
    rvmapply("outer", X, X, MoreArgs=list(FUN=FUN, ...))
  } else {
    rvmapply("outer", X, Y, MoreArgs=list(FUN=FUN, ...))
  }
}

plot.rv <- function (x, y=NULL, ...)
{ 
  .plot.default.rv(x, y, ...)
}

plot.rvsummary <- function (x, y=NULL, ...)
{ 
  .plot.default.rv(x, y, ...)
}

# ========================================================================
# plot.xy.rv  - 
# ========================================================================
#

.plot.xy.rv <- function (xy, type, pch = par("pch"), lty = par("lty"), col = par("col"), bg = NA, cex = 1, lwd = par("lwd"), ...)
{
  ##if (is.null(xy$rv)) {
  ##  return(.Internal(plot.xy(xy, type, pch, lty, col, bg, cex, lwd, ...)))
  ##}
  args <- c(xy, list(type=type, pch=pch, lty=lty, col=col, bg=bg, cex=cex, lwd=lwd, ...))
  typego <- list(p="points.rv", "l"="points.rv", "b"="points.rv")
  if (type=="n") return(invisible(NULL))
  if (is.null(plotroutine <- typego[[type]])) {
    stop("Plot type '", type, "' not yet implemented for rv objects!")
  }
  do.call(plotroutine, args)
  invisible(NULL)
}

# ========================================================================
# xy.coords.rv - modified version of xy.coords to accommodate rvs
# ========================================================================

.xy.coords.rv <- function (x, y = NULL, xlab = NULL, ylab = NULL, log = NULL, recycle = FALSE)
{
    if (is.null(y)) {
        ylab <- xlab
        if (is.language(x)) {
            if (inherits(x, "formula") && length(x) == 3) {
                ylab <- deparse(x[[2]])
                xlab <- deparse(x[[3]])
                y <- eval(x[[2]], environment(x), parent.frame())
                x <- eval(x[[3]], environment(x), parent.frame())
            }
            else stop("invalid first argument")
        }
        else if (is.complex(x)) {
            y <- Im(x)
            x <- Re(x)
            xlab <- paste("Re(", ylab, ")", sep = "")
            ylab <- paste("Im(", ylab, ")", sep = "")
        }
        else if (is.matrix(x) || is.data.frame(x)) {
            x <- data.matrix(x)
            if (ncol(x) == 1) {
                xlab <- "Index"
                y <- x[, 1]
                x <- 1:length(y)
            }
            else {
                colnames <- dimnames(x)[[2]]
                if (is.null(colnames)) {
                  xlab <- paste(ylab, "[,1]", sep = "")
                  ylab <- paste(ylab, "[,2]", sep = "")
                }
                else {
                  xlab <- colnames[1]
                  ylab <- colnames[2]
                }
                y <- x[, 2]
                x <- x[, 1]
            }
        }
        else if (is.list(x) && !is.rvobj(x)) { #### This is the only change ####
            xlab <- paste(ylab, "$x", sep = "")
            ylab <- paste(ylab, "$y", sep = "")
            y <- x[["y"]]
            x <- x[["x"]]
        }
        else {
            if (is.factor(x)) 
                x <- as.numeric(x)
            xlab <- "Index"
            y <- x
            x <- seq(along = x)
        }
    }
    if (inherits(x, "POSIXt")) 
        x <- as.POSIXct(x)
    if (length(x) != length(y)) {
        if (recycle) {
            if ((nx <- length(x)) < (ny <- length(y))) 
                x <- rep(x, length.out = ny)
            else y <- rep(y, length.out = nx)
        }
        else stop("'x' and 'y' lengths differ")
    }
    if (length(log) && log != "") {
        log <- strsplit(log, NULL)[[1]]
        f <- function (x) ((Pr(x < 0)>0) & !rv.any.na(x))
        if ("x" %in% log && any(ii <- f(x))) {
            n <- as.integer(sum(ii))
            warning(sprintf(ngettext(n, "%d x value <= 0 omitted from logarithmic plot", 
                "%d x values <= 0 omitted from logarithmic plot"), 
                n), domain = NA)
            x[ii] <- NA
        }
        if ("y" %in% log && any(ii <- f(y))) {
            n <- as.integer(sum(ii))
            warning(sprintf(ngettext(n, "%d y value <= 0 omitted from logarithmic plot", 
                "%d y values <= 0 omitted from logarithmic plot"), 
                n), domain = NA)
            y[ii] <- NA
        }
    }
    return(list(x = as.real(x), y = as.real(y), xlab = xlab, ylab = ylab))
}


.plot.default.rv <- function (x, y = NULL, type = "p", xlim = NULL, ylim = NULL, log = "", main = NULL, sub = NULL, xlab = NULL, ylab = NULL, ann = par("ann"), axes = TRUE, frame.plot = axes, panel.first = NULL, panel.last = NULL, asp = NA, rvlwd = rvpar("rvlwd"),  rvcol=rvpar("rvcol"), rvpoint=rvpar("rvpoint"), rvlex=rvpar("rvlex"), ...) 
{
    localAxis <- function(..., col, bg, pch, cex, lty, lwd) Axis(...)
    localBox <- function(..., col, bg, pch, cex, lty, lwd) box(...)
    localWindow <- function(..., col, bg, pch, cex, lty, lwd) plot.window(...)
    localTitle <- function(..., col, bg, pch, cex, lty, lwd) title(...)
    xlabel <- if (!missing(x)) 
        deparse(substitute(x))
    ylabel <- if (!missing(y)) 
        deparse(substitute(y))
    xy <- .xy.coords.rv(x, y, xlabel, ylabel, log)
    xlab <- if (is.null(xlab)) 
        xy$xlab
    else xlab
    ylab <- if (is.null(ylab)) 
        xy$ylab
    else ylab
    xlim <- if (is.null(xlim)) {
        range(rvfiniterange(xy$x))
    } else xlim
    ylim <- if (is.null(ylim)) {
        range(rvfiniterange(xy$y))
    } else ylim
    plot.new()
    localWindow(xlim, ylim, log, asp, ...)
    panel.first
    .plot.xy.rv(xy, type, rvlwd = rvlwd,  rvcol=rvcol, rvpoint=rvpoint, rvlex=rvlex, ...)
    panel.last
    if (axes) {
        localAxis(x, side = 1, ...)
        localAxis(y, side = 2, ...)
    }
    if (frame.plot) 
        localBox(...)
    if (ann) 
        localTitle(main = main, sub = sub, xlab = xlab, ylab = ylab, 
            ...)
    invisible()
}

# ========================================================================
# points.rv  -  plot points and uncertainty intervals
# ========================================================================

points.rv <- function (x, y = NULL, type = "p", xlim = NULL, ylim = NULL, rvlwd = rvpar("rvlwd"), rvcol = rvpar("rvcol"), rvpoint = rvpar("rvpoint"), rvlex = rvpar("rvlex"), ...) 
{
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
        }
        else {
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



.rvjointdrawxy <- function (x, y, size=1, reject.na=TRUE)
{
  xy <- c(x, y)
  s <- rvsample(xy, size=size, jointly=TRUE, reject.na=reject.na)
  if (is.null(dim(s))) s <- t(s)
  xs <- t(s[,seq(along=x)])
  ys <- t(s[,-seq(along=x)])
  list(x=xs, y=ys)
}

.packageName <- "rv"
#
# rv-Bayes.R
# 

# ========================================================================
# function  -  short description
# ========================================================================
# Name:        
# Description: 
# Parameters:  
# Required:    none
# History:     2004-06-  : 
#

postsim <- function(fit)
{
  UseMethod('postsim')
} 

# Old code:
#    for (sim in 1:nsim){
#      sigma[sim] <- sigma.hat*sqrt((n-k)/rchisq(1,n-k))
#      beta[sim,] <- mvrnorm (1, beta.hat, V.beta*sigma[sim]^2)
#    }
#

# ========================================================================
# postsim.lm
# ========================================================================

postsim.lm <- function(fit)
{
  # Modified version of 'sim' (from Andrew Gelman's library (gelman@stat.columbia.edu))
  #
  summ <- summary (fit)
  sigma.hat <- summ$sigma
  beta.hat <- summ$coef[,1]
  V.beta <- summ$cov.unscaled
  k <- summ$df[1]
  n <- k + summ$df[2]
  sigma <- sigma.hat*sqrt((n-k)/rvchisq(1,n-k))
  beta.0 <- as.vector(beta.hat) + sigma * rvnorm(mean=0, var=V.beta)
  names(beta.0) <- names(beta.hat)
  c(sigma=sigma, beta=beta.0)
}

# ========================================================================
# postsim.glm
# ========================================================================

postsim.glm <- function(fit)
{
  # Modified version of 'sim' (from Andrew Gelman's library (gelman@stat.columbia.edu))
  #
  summ <- summary (fit, correlation=T)
  beta.hat <- summ$coef[,1]
  sd.beta <- summ$coef[,2]
  corr.beta <- summ$corr
  k <- summ$df[1]
  n <- k + summ$df[2]
  V.beta <- corr.beta * array(sd.beta,c(k,k)) * t(array(sd.beta,c(k,k)))
  # dimnames(beta) <- list (NULL, names(beta.hat))
  rvnorm(mean=beta.hat, var=V.beta)
}

# ========================================================================
# print.rv  -  print summary of a rv on the console
# ========================================================================

print.rv <- function(x, digits=rvpar("print.digits"), ...)
{
  if (length(x)==0) {
    return(cat("rv(0)\n"))
  }
  print(summary(x, ...), digits=digits, ...)
}





print.rvfactor <- function(x, all.levels=FALSE, ...)
{
  s <- summary(x, all.levels=all.levels)
  ds <- dimnames(s)
  if (!is.null(.names <- names(x))) {
    s <- cbind(name=.names, s)
  } else if (!is.null(.dn <- dimnames(x))) {
    sud <- rvpar("summary.dimnames")
    if (!is.null(sud) && !is.na(sud) && is.logical(sud)) {
      da <- lapply(.dn, function (na) if (is.null(na)) rep("", nrow(s)) else na)
      rw <- da[[1]][row(x)]
      cl <- da[[2]][col(x)]
      s <- cbind(row=rw, col=cl, " "=":", s)
    }
  }
  print(s)
}


# ========================================================================
# quantiles of a random vector
# ========================================================================
#

quantile.rv <- function(x, ...)
{
  simapply(x, quantile, ...)
}

# ========================================================================
# median.rv - median of a random vector
# ========================================================================

median.rv <- function(x, na.rm=FALSE) ## EXPORT median.rv
{
  simapply(x, median, na.rm=na.rm)
}

is.random <- function(x)
{
  return(!is.constant(x))
}



range.rv <- function(..., na.rm=FALSE, finite=FALSE)
{
  sims <- sims(c(...)) # an L x n matrix of numbers
  m <- apply(sims, 1, 'range', na.rm=na.rm, finite=finite) # gives a 2xL matrix!!!
  r <- rvsims(t(m)) # Transpose to get an L x 2 matrix.
  names(r) <- c('min', 'max')
  return(r)
}


rep.rv <- function (x, times, ...)
{
  if (!is.rv(x)) return(rep(x, times, ...))
  if (missing(times)) {
    a <- rep(unclass(x), ...)
  } else {
    a <- rep(unclass(x), times, ...)
  }
  class(a) <- class(x)
  return(a)
}

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# UTILITIES in the package rv
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  .show(x)
#  .matrix.ragged.array(x)
#  .expand.as.matrix(M)
#  .listByName(x)
#  .permut
#  .dim.index
#  .bracket
#  .impute.by.name
#  .setDimensionByName
#  .rowmax
#  .slice
#  .nodups
#
# ========================================================================
# .show - show a given variable and its value (useful for debugging purposes)
# ========================================================================

.myc <- function (x)
{
  paste("c(",paste(x,collapse=","),")", sep="")
}

.show <- function (x)
{
  cat(deparse(substitute(x)),"<-", .myc(x), "\n")
}


# ========================================================================
# .matrix.ragged.array - make a matrix out of a ragged array
# ========================================================================
# not used now

.matrix.ragged.array <- function(x) # NOEXPORT
{
  # DEBUG: need something more efficient here!
  ncol <- max(sapply(x,length))
  s <- array(NA, c(length(x),ncol))
  for (i in seq(length(x))) {
    y <- x[[i]]
    s[i,seq(length(y))] <- y
  }
  s
}

# ========================================================================
# .expand.as.matrix - build a large matrix out of smaller matrices or scalars, diagonally
# ========================================================================
# Useful when building a large covariance matrix from a series of small ones.
#
# Parameters:  A rv or a list.
# Required:    none
# History:     2004-06-  : 
#

.expand.as.matrix <- function(M) # NOEXPORT
{
  # M is a rv or a list.
  dims <- c(0,0)
  for (i in 1:length(M)) {
    dims <- dims + if (is.null(dim(M[[i]]))) rep(length(M[[i]]),2) else dim(M[[i]])
  }
  m <- matrix(0,nrow=dims[1],ncol=dims[2])
  col.offset <- 0
  row.offset <- 0
  for (i in 1:length(M)) {
    x <- if (is.null(dim(M[[i]]))) diag(M[[i]],nrow=length(M[[i]])) else M[[i]]
    d <- dim(x)
    m[row.offset+1:d[1],col.offset+1:d[2]] <- x
    row.offset <- row.offset + d[1]
    col.offset <- col.offset + d[2]
  }
  m
}



# .permut: index of multidimensional array.
# 
# Incidentally we can get the index of the vector by
#  z <- .permut(dims)
#  index <- 1 + (z-1) %*% cumprod(z-1)
#  so index is a vector 1:prod(dims).
#

# ========================================================================
# .permut  -  make  permutations c(1,1) ... c(m,n)
# ========================================================================
# Example: 
#   .permut(c(2,3)) gives a matrix
#        [,1] [,2]
#   [1,]    1    1
#   [2,]    2    1
#   [3,]    1    2
#   [4,]    2    2
#   [5,]    1    3
#   [6,]    2    3
#

.permut <- function(dims) # NOEXPORT
{
  tl <- prod(dims)
  co <- NULL
  p <- 1
  for (i in 1:length(dims)) {
    d <- dims[i]
    x <- rep(1:d, each=p, length.out=tl)
    co <- cbind(co, x)
    p <- p*d
  }
  dimnames(co) <- NULL
  co
}

# ========================================================================
# .dimindex  -  make a list of indices of a matrix
# ========================================================================
#

.dimindex <- function(x) # NOEXPORT
{
  if (is.null(dim(x))) {
    ix <- .permut(length(x))
  } else {
    ix <- .permut(dim(x))
  }
  ixt <- paste('[',apply(ix, 1, paste, collapse=','),']',sep='')
  ixt
}

# ========================================================================
# .leftprepend  -  prepend spaces for short lines
# ========================================================================
#

.leftadjust <- function (x) # NOEXPORT
{
  m.x <- max(nchar(x))
  a <- (nchar(x)<m.x)
  if (any(a)) {
    d <- m.x - nchar(x)
    x[a] <- paste(' ', x[a], sep='') # only ONE single space
  }
  x
}

# ========================================================================
# .dim.index  -  make a list of indices of a matrix
# ========================================================================
#

.dim.index <- function(x, leftadjust=TRUE) # NOEXPORT
{
  if (is.null(dim(x)))
    ix <-.permut(length(x))
  else 
    ix <- .permut(dim(x))
  ixt <- paste('[',apply(ix, 1, paste, collapse=','),']',sep='')
  if (leftadjust) .leftadjust(ixt) else ixt
}

# ========================================================================
# .bracket  -  x[[ i ]] but recycle if i > length(x)
# ========================================================================

.bracket <- function(x,k) # NOEXPORT
{
  l <- length(x)
  if(k<=l) x[[k]] else x[[1+((k-1)%%l)]]
}

# ========================================================================
# .slice  -  apply "[ row ]" to arguments but recycle if row > length(argument)
# ========================================================================
# 
# Description: Like lapply(list(...),"[",row) but recycles.
# Parameters:  
# Required:    none
# History:     2004-06-  : 
#

.slice <- function(..., row) # NOEXPORT
{
  UseMethod('.slice')
}

.slice.default <- function(..., row) # NOEXPORT
{
  list(...)
}

.slice.rv <- function(..., row) # NOEXPORT
{
  v <- sapply(list(...), .bracket, row)
  class(v) <- class(rv())
  return(v)
}

# ========================================================================
# .bracket.indices - extract indices from the brackets
# ========================================================================
# Returns NA's for missing indices

.bracket.indices <- function(x, default.bracket=FALSE) # already in rv 0.923
{
  if (length(x)<1 || !is.character(x)) return(NULL)
  rgx_bracket <- "^(.*\\[(.*)\\].*)$"
  no.brackets <- (regexpr(rgx_bracket, x)<1)
  y <- sub("^(.*\\[(.*)\\].*)$", "\\2", x)
  y <- gsub(",[[:blank:]]*$", ",NA", y) # DEBUG: won't work for x[,,]
  if (any(no.brackets)) {
    y[no.brackets] <- if (default.bracket) "1" else "0"
  }
  y <- strsplit(y, "[[:blank:]]*,[[:blank:]]*")
  lapply(y, as.numeric)
}


# ========================================================================
# .indices - return the single-digit indices of components 
# ========================================================================
# x : a list of vectors of indices
# dim. : dimension or length of the target matrix or vector.

.indices <- function (x, dim. = NULL) 
{
    if (is.null(x)) 
        return(numeric(0))
    ld <- length(dim.)
    f <- function(x) {
        lx <- length(x)
        if (lx == 1) 
            return(x)
        if (lx != ld || any(x > dim.)) 
            return(NA)
        1 + sum(c(x - 1, 0) * c(1, cumprod(dim.)))
    }
    pos <- sapply(x, f)
    pos
}



# ========================================================================
# .impute.by.name - impute into a vector using the names attribute
# ========================================================================

.impute.by.name <- function (x, y)
{
  # impute x <- y using the names of y
  n.y <- names(y)
  name.of.y <- deparse(substitute(y))
  if (length(n.y)<1) stop("the names attribute is required for: ", name.of.y)
  bx <- .bracket.indices(n.y)
  ix <- .indices(bx, dim(x))
  if (length(ix)==1 && ix==0) {
    # Impute the whole vector.
    return(y)
  }
  na.ix <- is.na(ix)
  if (any(na.ix)) {
    warning("Couldn't figure out index ", n.y[which(na.ix)])
    ix <- ix[!na.ix]
    y  <- y[!na.ix]
    if (length(ix)<1) return(x)
  }
  x[ix] <- y
  if (is.null(dim(x)) || length(dim(x))==1) {
    dim(x) <- NULL
    n.x <- names(x)
    n.x[ix] <- names(y)
    n.x[is.na(n.x)] <- ""
    names(x) <- n.x
  }
  return(x)
}



# ========================================================================
# .setDimensionByName - 
# ========================================================================

.setDimensionByName <- function (x)
{
  # 
  names.x <- names(x)
  bix <- .bracket.indices(names.x, default.bracket=TRUE)
  if (is.null(bix)) stop("Names of x MUST be set")
  b <- sapply(bix, prod)
  max.ix <- which(b==max(b))[1]
  maxdim <- bix[[max.ix]]
  if (prod(maxdim)<1) stop("Invalid dimension")
  a <- array(NA, maxdim)
  new.x <- .impute.by.name(a, x)
  names.new.x <- names(new.x)
  if (length(maxdim)>1) {
    dim(new.x) <- maxdim
  }
  names(new.x) <- names.new.x
  new.x
}



# ========================================================================
# .make.names - make names of the components
# ========================================================================

.make.names <- function(x, name=deparse(substitute(x)))
{
  if (is.null(x) || length(x)<1) return(x)
  name <- make.names(name)
  paste(name, .dim.index(x), sep="")
}

# ========================================================================
# .shortnames - names of the components, with brackets removed
# ========================================================================

.shortnames <- function(x)
{
  na <- names(x)
  if (is.null(na)) return(na)
  unlist(lapply(strsplit(na, "[", fixed=TRUE), "[", 1), use.names=FALSE)
}

# ========================================================================
# make.fullnames - make eventually available!!
# ========================================================================

#.make.fullnames <- function(x)
#{
#  if (is.null(x) || length(x)<1) return(x)
#  name <- make.names(name)
#  paste(name, .dim.index(x), sep="")
#}


.list2array <- function (x, drop=FALSE)
{
  # name:
  #   list2array - coerce a list to an array, preserving dimnames
  # description:
  #   attempt to coerce the list to an array, preserving dimnames
  # arguments:
  #   x : (list) list to be coerced into one array
  # value:
  #   a vector if not all dimensions of list elements were equal;
  #   an array of dimension c(dim(x[[1]]),length(x)),
  #   if all dimensions of the elements were equal.
  # author:
  #   Jouni Kerman
  # history:
  #   2006-08-28
  #
  dm.1 <- length(x)
  if (dm.1<1) return(unlist(x))
  da <- lapply(x, function (x) if (is.null(dim(x))) length(x) else dim(x))
  dm.2 <- da[[1]]
  dn.1 <- names(x)
  dn.2 <- dimnames(x[[1]])
  # 
  a <- NULL
  for (i in seq(along=x)) {
    a <- c(a, x[[i]])
  }
  if (!all(sapply(da, function (d) all(d==dm.2)))) {
    return(a)
  }
  names(a) <- NULL
  dim(a) <- c(dm.2, dm.1)
  # If the list component does not have dimnames...
  if (is.null(dn.2)) {
    # ...and it doesn't have a dimension...
    if (is.null(dim(x[[1]]))) {
      # ...then use the names for dimensions...
      dn.2 <- names(x[[1]])
    } # ...or else just give up and use NULL for that dimension names.
    dn.a <- list(dn.2, dn.1)
  } else {
    dn.a <- c(dn.2, list(dn.1))
  }
  if (!is.null(names(dn.2))) {
    name.x <- deparse(substitute(x))
    names(dn.a) <- c(name.x, names(dn.2))
  }
  dimnames(a) <- dn.a
  if (drop) drop(a) else a
}



# .nodups

.nodups <- function (x, keep.latest=TRUE)
{
  # NAME
  # .nodups - remove list items with duplicate names (keep latest)
  # ARGUMENTS
  #   x : a list
  if (!keep.latest) stop("keep.latest=FALSE Not yet implemented")
  if (length(x)<2) return(x)
  a <- list()
  x <- rev(x) # use the fact that x[[na]] gets the *first* one
  for (na in names(x)) {
    a[[na]] <- x[[na]] # multiple assignments may occur, but that's ok.
  }
  return(a)
}


.dimind <- function (x, dim.=NULL, MARGIN=1)
{
  # dimension indicator (generalization of row, col)
  if (is.null(dim.)) {
    d <- dim(x)
  } else {
    d <- dim.
  }
  if (is.null(d)) {
    return(rep.int(1, length(x)))
  }
  X <- array(NA, dim=d)
  if (length(d)==1) {
    X[] <- 1
  } else if (length(d)==2) {
    X <- if (MARGIN==1) row(X) else if (MARGIN==2) col(X) else stop("subscript out of bounds")
  } else if (MARGIN==1) {
    X[] <- (1:d[MARGIN])
  } else {
    X[] <- rep(1:d[MARGIN], each=prod(d[1:(MARGIN-1)]))
  }
  return(X)
}


.rvmeansd <- function (x, names.=c("mean", "sd", "NAS", "n.sims")) # for convenience
{
  a <- list(dimnames=dimnames(x), dim=dim(x), names=names(x))
  S <- sims(as.rv(x))
  m <- colMeans(S, na.rm=TRUE)
  ns <- rvnsims(x)
  v <- ((colSums(S^2, na.rm=TRUE)-ns*(m^2))/(ns-1))
  v[rvnsims(x)==1] <- 0
  if (any(naS <- is.na(S))) {
    NAS <- (colMeans(naS)*100)
  } else {
    NAS <- rep.int(0, length(x))
  }
  s <- sqrt(v)
  attributes(m) <- a
  attributes(s) <- a
  L <- list(mean=m, sd=s, NAS=NAS, n.sims=ns)
  names(L) <- names.
  return(L)
}

# end rv-util.R



#
# create and test objects of type 'rv'
#

rv <- function(length=0)
{
  if (is.numeric(length)) {
    x <- as.list(rep.int(NA, length))
  } else {
    stop("length must be numeric")
  }
  class(x) <- "rv"
  return(x)
}

is.rv <- function(x)
{
  inherits(x, "rv")
}

as.rv <- function(x, ...)
{
  UseMethod("as.rv")
}

as.rv.rv <- function(x, ...)
{
  return(x)
}

as.rv.numeric <- function(x, ...)
{
  if (is.rv(x)) return(x)
  r <- rvsims(matrix(as.vector(x), nrow=1))
  cr <- class(r)
  attributes(r) <- attributes(x)
  class(r) <- cr
  return(r)
}

as.rv.logical <- function(x, ...)
{
  as.rv.numeric(x)
}

as.rv.integer <- function(x, ...)
{
  as.rv.numeric(x)
}

as.rv.list <- function(x, ...)
{
  stop("Cannot coerce an arbitrary list to an rv object")
}

as.rv.matrix <- function(x, ...)
{
  as.rv.numeric(x)
}

as.rv.default <- function(x, ...)
{
  if (is.null(x)) return(NULL)
  stop('Cannot coerce object of class "', paste(class(x), collapse="/"), '" to rv')
}

as.rv.xtabs <- function (x, ...) 
{
  # NAME
  #  as.rv.xtabs
  #
  as.rv(x[])
}

# DEBUG:1 ?? The permutation of !is.null(dim(sims)) is not done yet
# DEBUG:2 (must take permutation of row indices and then permute.


as.real.rv <- function(x, ...)
{
  simapply(x, as.real, ...)
}

as.double.rv <- function(x, ...)
{
  simapply(x, as.double, ...)
}

as.logical.rv <- function(x, ...)
{
  simapply(x, as.logical, ...)
}

as.integer.rv <- function (x, ...)
{
  simapply(x, as.integer, ...)
}



rvVectorize <- function (FUN, vectorize.args = arg.names, SIMPLIFY = FALSE, USE.NAMES = TRUE, SAMPLESIZE=NULL) 
{
  arg.names <- as.list(formals(FUN))
  arg.names[["..."]] <- NULL
  arg.names <- names(arg.names)
  vectorize.args <- as.character(vectorize.args)
  if (!length(vectorize.args)) {
      return(FUN)
  }
  if (!all(vectorize.args %in% arg.names)) {
      stop("must specify formal argument names to vectorize")
  }
  .FUNV <- function() {
    args <- lapply(as.list(match.call())[-1], eval, envir=parent.frame())
    names <- if (is.null(names(args))) {
      character(length(args))
    } else {
      names(args)
    }
    dovec <- names %in% vectorize.args
    Args <- .Primitive("c")(FUN = FUN, args[dovec], MoreArgs = list(args[!dovec]), 
      SIMPLIFY = SIMPLIFY, USE.NAMES = USE.NAMES, SAMPLESIZE=SAMPLESIZE)
    do.call(rvmapply, Args)
  }
  formals(.FUNV) <- formals(FUN)
  return(.FUNV)
}



rv.any.na <- function (x)
{
  # NAME
  #  rv.any.na - Which components have missing values?
  #
  if (is.rvsummary(x)) {
    return(x$NAS>0)
  } else {
    return(colSums(is.na(sims(as.rv(x))))>0)
  }
}



rv.all.na <- function (x)
{
  # NAME
  #  rv.all.na - Which components are completely missing?
  if (is.rvsummary(x)) {
    return(x$NAS==1)
  } else {
    return(colSums(is.na(sims(as.rv(x))))==1)
  }
}


rvattach <- function (what, name=deparse(substitute(what)), overwrite=TRUE, impute=FALSE, ...)
{
  if (missing(what) || is.null(what)) {
    return(invisible(NULL))
  }
  if (!is.list(what)) {
    stop("Argument must be a list")
  }
  while (any(search()==name)) {
    eval(substitute(detach(name), list(name=as.name(name))))
  }
  if (!is.rv(what)) {
    stop("Argument must be an rv")
  }
  a <- splitbyname(what)
  rattach(a, name=name, overwrite=overwrite, impute=impute, ...)
}

rattach <- function (what, name=deparse(substitute(what)), overwrite=TRUE, impute=FALSE, ...)
{
  if (missing(what) || is.null(what)) {
    return(invisible(NULL))
  }
  if (!is.list(what)) {
    stop("Argument must be a list")
  }
  while (any(search()==name)) {
    eval(substitute(detach(name), list(name=as.name(name))))
  }
  if (is.null(what)) return(NULL)
  if (overwrite) {
    ls.GE <- ls(.GlobalEnv)
    conf_names <- names(what)[names(what)%in%ls.GE]
    for (obj_name in conf_names) {
      if (impute) {
        v <- get(obj_name, .GlobalEnv) # Value to be imputed to
        what[[obj_name]] <- .impute.by.name(v, what[[obj_name]])
      }
      remove (list=obj_name, envir=.GlobalEnv)
    }
  }
  attach (what, name=name, ...)
}
# ========================================================================
# rvattr  -  return the rvsim attribute for each component
# ========================================================================
# returns a list.

rvattr <- function(x, attrib)
{
  a <- lapply(unclass(x), "attr", "rvsim")
  a <- lapply(a, "[[", attrib)
  nulls <- sapply(a, is.null)
  a[nulls] <- NA
  a
}

# rvattr(x, 'name') <- vector of values for each component or x; or,
# rvattr(x) <- list of attributes: list(attributename=vector of values, attributename2=...)
# e.g. list(Rhat=c(1.01,1.03,1.23), n.eff=c(200,100,16)).
#



"rvattr<-" <- function(x, attrib=NULL, value)
{
  if (is.null(attrib)) {
    for (a in names(value)) {
      rvattr(x, a) <- value[[a]]
    }
  } else if (length(names(value))>0 && length(names(x))>0) {
    for (name in names(value)) {
      if (name %in% names(x)) {
        A <- attr(x[[name]], "rvsim")
        if (is.null(A)) { A <- list() }
        A[[attrib]] <- value[[name]]
        attr(x[[name]], "rvsim") <- A
      } else {
        stop("No name ", name, " in rv.")
      }
    }
  } else if (length(x)!=length(value)) {
    stop("Value must be of the same length as x or the matching must be done by name.")
  } else {
    for (i in seq(along=x)) {
      A <- attr(x[[i]], "rvsim")
      if (is.null(A)) { A <- list() }
      A[[attrib]] <- value[[i]]
      attr(x[[i]], "rvsim") <- A
    }
  }
  return(x)
}




rvbern <- function (n=1, prob, logical=FALSE)
{
  r <- rvvapply(stats:::rbinom, n.=n, size=1, prob=prob)
  if (logical) {
    r <- as.logical(r)
  }
  return(r)
}



rvbeta <- function (n=1, shape1, shape2)
{
  rvvapply(stats:::rbeta, n.=n, shape1=shape1, shape2=shape2)
}



rvbinom <- function (n=1, size, prob)
{
  rvvapply(stats:::rbinom, n.=n, size=size, prob=prob)
}




rvboot <- function (data)
{
#  empirical (bootstrap) distribution
#
  n.sims <- getnsims()
  n <- n.sims*length(data)
  s <- matrix(sample(data, size=n, replace=TRUE), nrow=n.sims)
  r <- rvsims(s)
  dim(r) <- dim(data)
  return(r)
}



rvcat <- function (n=1, prob, levels=NULL)
{
  # NAME
  #  rvcat - Sample Categorical Random Variables
  if (anyisrv(n, prob)) {
    x <- rvmultinom(n=n, size=1, prob=prob)
    s <- sims(x, dimensions=TRUE)
    ds <- dim(s)
    s <- as.logical(s)
    dim(s) <- ds
    f <- function (m) row(m)[m]
    a <- apply(s, 1, f)
    r <- if (is.null(dim(a))) rvsims(a) else rvsims(t(a))
  } else {
    n.sims <- getnsims()
    s <- rmultinom(n=n*n.sims, size=1, prob=prob)
    s <- row(s)[as.logical(s)]
    dim(s) <- c(n.sims, length(s) %/% n.sims)
    r <- rvsims(s)
  }
  rvfactor(r, levels=levels)
}


rvcauchy <- function (n=1, location=0, scale=1)
{
  rvvapply(stats:::rcauchy, n.=n, location=location, scale=scale)
}




rvchisq <- function (n=1, df, ncp = 0) 
{
  if (missing(ncp)) {
    rvvapply(stats:::rchisq, n.=n, df=df)
  } else {
    rvvapply(stats:::rchisq, n.=n, df=df, ncp=ncp)
  }
}




rvci <- function(obj, interval=0.95, one.sided=FALSE, left=TRUE)
{
  # NAME
  #  rvci - Uncertainty (Credible) Interval
  # 
  if (length(interval)>1) {
    ci <- rvquantile(obj, probs=range(interval))
  } else if (one.sided) {
    q <- if (left) interval else (1-interval)
    ci <- rvquantile(obj, q)
  } else {
    lower <- (1-interval)/2
    upper <- (lower + interval)
    ci <- rvquantile(obj, c(lower,upper))
  }
  return(ci)
}



# Put the following in rv/par$rvcolorthemes

.rvcolorthemes <- list(
  default=c("grey20", "grey40", "grey60"),
  gray=c("grey20", "grey40", "grey60"),
  lightgray=c("grey20", "grey40", "grey60"),
  darkgray=c("black", "grey20", "grey40")
)

.makervcolortheme <- function (col)
{
  colors.with.numbers <- colors()[regexpr("[a-z]1$", colors())>0]
  colors.with.numbers <- sub("1$", "", colors.with.numbers)
  if (!col %in% colors.with.numbers) {
    if (col %in% colors())
      return(rep(col,3))
    else 
      return(NULL)
  } 
  paste(col, c("3", "2", ""), sep="") 
}

rvcolortheme <- function (theme) # NOEXPORT
{
  # NAME
  #  rvcolortheme - color themes for plotting uncertainty intervals
  # DETAILS
  # (please, later export this!)
  synonym <- list(grey="gray", lightgrey="lightgray", darkgrey="darkgray")
  theme <- theme[1]
  if (is.null(theme) || is.na(theme)) {
    theme <- "default"
  }
  (!is.null(co <- .rvcolorthemes[[theme]])) && return(co)
  (!is.null(stheme <- synonym[[theme]])) && return(.rvcolorthemes[[stheme]])
  co <- .makervcolortheme(theme)
  if (is.null(co)) {
    warning(paste("No rv color theme '",theme,"' found, using default",sep=""))
    co <- .rvcolorthemes[["default"]]
  } else {
    .rvcolorthemes[[theme]] <- co
  }  
  names(co) <- c("dot", "thick", "thin")
  return(co)
}


rvcompatibility <- function (level=NULL)
{
  old.level <- getOption("rvcompatibility")
  if (!is.null(level)) {
    if (level==0) {
      .rvFunctionSwitch(attach=TRUE)
    } else if (level==1) {
      .rvFunctionSwitch(attach=FALSE)
    } else {
      warning("Unknown level, ignored.")
      return(old.level)
    }
    options(rvcompatibility=level)
  }
  return(old.level)
}

# =================================================================================
.min  <- base:::min
.max  <- base:::max
.pmin <- base:::pmin
.pmax <- base:::pmax
.c    <- base:::c
# =================================================================================

###.makeRvCompatible_base <- function ()
###{
  ..is.recursive <- function (x) {
    ((!is.rv(x)) && .Primitive("is.recursive")(x))
  }
  .rvRegisterFunctionSwitch(is.recursive, "is.recursive", 'base')
  #
  ..is.atomic <- function (x)
  {
    (is.rv(x) || .Primitive("is.atomic")(x))
  }
  .rvRegisterFunctionSwitch(..is.atomic, "is.atomic", 'base')
  #
  ..is.vector <- function(x, mode="any")
  {
    if (is.rv(x)) return(is.null(dim(x)))
   .Internal(is.vector(x, mode))
  }
  .rvRegisterFunctionSwitch(..is.vector, "is.vector", 'base')
  #
  ..is.integer <- function (x) {
    ((is.rv(x) && all(rvsimapply(x, is.integer))) || .Primitive("is.integer")(x))
  } 
  .rvRegisterFunctionSwitch(..is.integer, "is.integer", 'base')
  #
  ..is.logical <- function (x) {
    ((is.rv(x) && all(rvsimapply(x, is.logical))) || .Primitive("is.logical")(x))
  } 
  .rvRegisterFunctionSwitch(..is.logical, "is.logical", 'base')
  #
  #
  #
  ..min <- function(..., na.rm=FALSE) {
    if (anyisrv(...)) {
      simapply(cbind.rv(...), .min, na.rm=na.rm)
    } else {
      .min(..., na.rm=na.rm)
    }
  }
  .rvRegisterFunctionSwitch(..min, "min", "base")
  #
  ..max <- function(..., na.rm=FALSE) {
    if (anyisrv(...)) {
      simapply(cbind.rv(...), .max, na.rm=na.rm)
    } else {
      .max(..., na.rm=na.rm)
    }
  }
  .rvRegisterFunctionSwitch(..max, "max", "base")
  #
  ..pmin <- function(..., na.rm=FALSE) {
    if (anyisrv(...)) {
      a <- sims(cbind.rv(...), dimensions=TRUE)
      rvsims(t(apply(a, 1, function (m) apply(m, 1, .pmin))))
    } else {
      .pmin(..., na.rm=na.rm)
    }
  }
  .rvRegisterFunctionSwitch(..pmin, "pmin", "base")
  #
  ..pmax <- function(..., na.rm=FALSE) {
    if (anyisrv(...)) {
      a <- sims(cbind.rv(...), dimensions=TRUE)
      rvsims(t(apply(a, 1, function (m) apply(m, 1, .pmax))))
    } else {
      .pmax(..., na.rm=na.rm)
    }
  }
  .rvRegisterFunctionSwitch(..pmax, "pmax", "base")
  #
  #
  #
  ..c <- function(..., recursive=FALSE) {
    an_rv <- anyisrv(...)
    x <- .c(..., recursive=recursive)
    if (an_rv) {
      class(x) <- class(rv())
    }
    return(x)
  }
  .rvRegisterFunctionSwitch(..c, 'c', 'base')
  #
  .rvBracketAssignment <- function(x, ..., value=NULL)
  {
    if (is.rv(x) || is.rv(value)) {
      x <- as.rv(x)
      value <- as.rv(value)
      x <- .Primitive("[<-")(x, ..., value=value)
      class(x) <- class(rv())
    } else {
      x <- .Primitive("[<-")(x, ..., value=value)  
    }
    return(x)
  }
  .rvRegisterFunctionSwitch(.rvBracketAssignment, "[<-", 'base')
  #
  .rvMatrixProduct <- get("%*%.rv")
  .rvRegisterFunctionSwitch(.rvMatrixProduct, "%*%", 'base')
###}


#
#
#

plot <- function (x, y, ...) 
{
    if (!missing(y) && is.rv(y)) {
      return(plot.rv(x, y, ...))
    } else if (is.null(attr(x, "class")) && is.function(x)) {
        nms <- names(list(...))
        if (missing(y)) 
            y <- {
                if (!"from" %in% nms) 
                  0
                else if (!"to" %in% nms) 
                  1
                else if (!"xlim" %in% nms) 
                  NULL
            }
        if ("ylab" %in% nms) 
            graphics:::plot.function(x, y, ...)
        else graphics:::plot.function(x, y, ylab = paste(deparse(substitute(x)), 
            "(x)"), ...)
    }
    else UseMethod("plot")
}
.rvRegisterFunctionSwitch(plot, "plot", "graphics")

abline <- function (a = NULL, b = NULL, h = NULL, v = NULL, ...)
{
  if (anyisrv(a,b,h,v)) {
    abline.rv(a=a, b=b, h=h, v=v, ...)
  } else {
    .abline.default(a=a, b=b, h=h, v=v, ...)
  }
  invisible()
}

.rvRegisterFunctionSwitch(abline.rv, "abline", "graphics")

var <- function(x, ...)
{
  UseMethod("var")
}
.rvRegisterFunctionSwitch(var, "var", "stats")

sd <- function (x, na.rm = FALSE)
{
    if (is.rv(x))
      return(sd.rv(x, na.rm = na.rm))
    else if (is.matrix(x)) 
      apply(x, 2, sd, na.rm = na.rm)
    else if (is.vector(x)) 
      sqrt(var(x, na.rm = na.rm))
    else if (is.data.frame(x)) 
      sapply(x, sd, na.rm = na.rm)
    else sqrt(var(as.vector(x), na.rm = na.rm))
}

.rvRegisterFunctionSwitch(sd, "sd", "stats")

cov <- function(x, ...)
{
  UseMethod("cov")
}
.rvRegisterFunctionSwitch(cov, "cov", "stats")

cor <- function(x, ...)
{
  UseMethod("cor")
}
.rvRegisterFunctionSwitch(cor, "cor", "stats")



rvconst <- function(n=1, x=0)
{
  lx <- length(x)
  v <- rv(lx)
  if (lx>0) {
    v[1:lx] <- x # recycle
  }
  return(v)
}

# ========================================================================
# rvcov  -  covariance matrix
# ========================================================================
#

rvcov <- function(x, y=NULL, ...)
{
  if (is.null(y)) {
    cov(sims(x), y, ...)
  } else {
    cov(sims(x), sims(y), ...)
  }
}



rvcut <- function (x, ...)
{
  UseMethod("rvcut")
}

rvcut.default <- function (x, ...)
{
  f <- cut(x, ...)
  levs <- levels(f)
  rvf <- rvsims(as.integer(f))
  class(rvf) <- c("rvfactor", class(rvf))
  attr(rvf, "levels") <- levs
  return(rvf)
}

rvcut.rv <- function (x, ...)
{
  a <- sims(x)
  f <- cut(a, ...)
  levs <- levels(f)
  rvf <- rvsims(array(as.integer(f), dim(a)))
  class(rvf) <- c("rvfactor", class(rvf))
  attr(rvf, "levels") <- levs
  return(rvf)
}




rvdens <- function(n=1, FUN, range, unitprecision=10, ...)
{
  # NAME
  #   rvdensity - Sample from a given univariate density using a grid approximation
  # ARGUMENTS
  #   n : number of independent random vector components to draw
  #   FUN : density function, must be vectorized
  #   range : range for the grid
  #   unitprecision : number of points per unit
  #   ... : other arguments passed to [FUN].
  #   
  grid <- seq(from=range[1], to=range[2], by=1/unitprecision)
  prob <- FUN(grid, ...)
  n.sims <- getnsims()
  s <- sample(grid, size=n*n.sims, prob=prob, replace=TRUE)
  noise <- runif(n*n.sims, -0.5/unitprecision, 0.5/unitprecision)
  rvsims(matrix(s+noise, nrow=n.sims))
}




rvdirichlet <- function (n = 1, alpha) 
{
  x <- NULL
  for (i in 1:n) {
    g <- rvgamma(n = 1, shape = alpha, scale = 1)
    x <- cbind.rv(x, g/sum(g))
  }
  return(x)
}



rvdiscrete <- function (n=1, x, prob=NULL)
{
  n.sims <- getnsims()
  rvsims(matrix(sample(x=x, size=n*n.sims, prob=prob, replace=TRUE), nrow=n.sims))
}



rvexp <- function (n=1, rate=1)
{
  rvvapply(stats:::rexp, n.=n, rate=rate)
}




rvfactor <- function (x, ...)
{
  UseMethod("rvfactor")
}

rvfactor.default <- function (x, levels=NULL, ...)
{
  f <- as.factor(x)
  a <- sims(as.rv(as.integer(f)))
  rvf <- rvsims(a)
  class(rvf) <- c("rvfactor", class(rvf))
  attr(rvf, "levels") <- if (is.null(levels)) levels(f) else levels
  return(rvf)
}


rvfactor.rv <- function (x, levels=NULL, ...)
{
  a <- sims(x)
  f <- as.factor(a)
  levs <- levels(f)
  rvf <- rvsims(array(as.integer(f), dim(a)))
  class(rvf) <- c("rvfactor", class(rvf))
  attr(rvf, "levels") <- if (is.null(levels)) levs else levels
  return(rvf)
}

is.numeric.rvfactor <- function (x)
{
  FALSE
}

is.rvfactor <- function (x)
{
  UseMethod("is.rvfactor")
}

is.rvfactor.rvfactor <- function (x)
{
  TRUE
} 

is.rvfactor.rv <- function (x)
{
  all(rvsimapply(x, is.factor))
} 

is.rvfactor.default <- function (x)
{
  FALSE
} 

as.rvfactor <- function (x, ...)
{
  if (is.rvfactor(x)) x else rvfactor(x)
} 


as.rv.rvfactor <- function (x, ...)
{
  attr(x, "levels") <- NULL
  clx <- class(x)
  clx <- clx[clx!="rvfactor"]
  class(x) <- clx
  return(x)
}



rvfiniterange <- function (x) ## NOEXPORT
{
  if (is.rvsummary(x)) {
    if (is.null(x$finiterange)) {
      apply(x$quantiles, 1, range)
    } else {
      return(x$finiterange)
    }
  } else {
    rvsimapply(x, finiterange)
  }
}




rvgamma <- function (n=1, shape, rate = 1, scale = 1/rate) 
{
  rvvapply(stats:::rgamma, n.=n, shape=shape, scale=scale)
}

# ========================================================================
# rvhist  -  plot histograms of the simulations of the random components
# ========================================================================
#

rvhist <- function (x, ...) 
{
    if (!is.null(dim(x))) 
        par(mfcol = dim(x))
    mfcol <- par("mfcol")
    n <- prod(mfcol)
    n <- min(n, length(x))
    a <- list(...)
    make.main <- is.null(a$main)
    make.xlab <- is.null(a$xlab)
    lab <- deparse(substitute(x))
    x.names <- paste(lab, .dimindex(x), sep = "")
    for (i in 1:n) {
        a$x <- sims(x[i])
        if (make.main || make.xlab) {
            this.name <- x.names[i]
            if (make.xlab) {
                a$xlab <- this.name
            }
            if (make.main) {
                a$main <- paste("Histogram of", this.name)
            }
        }
        if (is.null(a$breaks)) 
            a$breaks <- "fd"
        do.call("hist", a)
    }
}




# ========================================================================
# rvhyper  -  hypergeometric rvs
# ========================================================================

rvhyper <- function(nn=1, m, n, k) # NOEXPORT -- not yet ready
{
  attr(nn, "n.name") <- "nn"
  rvvapply(stats:::rhyper, n.=nn, m=m, n=n, k=k)
}

	

rvifelse <- function (test, yes, no)
{
  #   rvifelse - If-Then-Else For Random Vectors
  rvmapply(base:::ifelse, test=test, yes=yes, no=no)
}




rvintervals <- function (x, rvpoint=rvpar("rvpoint"), ...) # NOEXPORT
{
  which.quantiles <- list(
    "NA" = NA,
    mean = NA,
    median = 0.50,
    "50%" = c(0.25, 0.75),
    "80%" = c(0.10, 0.90),
    "95%" = c(0.025, 0.975)
  )
  .whichq <- function (iv) {
    if (is.numeric(iv)) {
      iv <- paste(100*iv, "%", sep="")
    } else {
      (!is.null(q <- which.quantiles[[iv]])) && return(q)
    }
    if (is.na(iv)) return(NA)
    n <- nchar(iv)
    if (substr(iv,n,n)=="%") {
      ivn <- as.numeric(substr(iv,1,n-1))
      c(0.5-ivn/200, 0.5+ivn/200)
    } else {
      NA
    }
  }
  .length <- function (iv) {
    lg <- .whichq(iv)
    if (length(lg)<=1) 0 else diff(lg)
  }
  probs <- as.vector(na.omit(unlist(lapply(rvpoint, .whichq))))
  if (length(probs)<=1) {
    # A trick to force probs into a named array
    # (won't otherwise return names if we have only one quantile, e.g. 0.50)
    probs <- c(probs, NA)
  }
  if (!all(is.na(probs))) {
    Q <- t(rvquantile(x, probs=probs, ...))
  } else {
    Q <- NULL # DEBUG: will this be ignored if we have only "mean" e.g.? #
  }
  compute.what <- list(
    "NA"   = function () NA,
    mean   = function () t(as.vector(rvmean(x))),
    median = function () Q["50%",,drop=FALSE],
    "50%"  = function () Q[c("25%","75%"),,drop=FALSE],
    "80%"  = function () Q[c("10%","90%"),,drop=FALSE],
    "95%"  = function () Q[c("2.5%","97.5%"),,drop=FALSE]
  )
  .lbl <- function (p) { # From 'quantile.default'
    if (is.null(p) || is.na(p)) return(NA)
    dig <- max(2, getOption("digits"))
    paste(formatC(100 * p, format = "fg", wid = 1, digits = dig), "%", sep = "")
  }
  .summaries <- function (iv) {
    if (is.null(f <- compute.what[[iv]])) {
      a <- na.omit(sapply(.whichq(iv),.lbl))
      if (all(a %in% dimnames(Q)[[1]])) {
        return(Q[a,,drop=FALSE])
      } else {
        warning("Cannot understand interval '", iv, "'")
        return(NA)
      }
    }
    a <- f()
    if (is.null(dim(a))) {
       if (length(x)==1) {
          a <- t(a)
       } else {
         na <- names(a)
         dim(a) <- c(1,length(a))
         dimnames(a) <- list(iv, na)
       }
    }
    return(a)
  } 
  lgth <- rev(order(sapply(rvpoint, .length)))
  s <- lapply(rvpoint, .summaries)
  names(s) <- rvpoint
  s[lgth]
}


rvinvchisq <- function (n=1, df, scale=1)
{
  return(df*scale/rvchisq(n=n, df=df))
}


# ========================================================================
# rvmapply  -  apply a function to multiple rv objects
# ========================================================================
# Note. Won't work with functions allowing "blank" arguments
# such as "[" (e.g. x[y,,]). The functions "[" and "[<-" use
# modified versions of simmapply.
#

rvmapply <- function (FUN, ..., MoreArgs=NULL, SIMPLIFY=FALSE, USE.NAMES=TRUE, SAMPLESIZE=NULL)
{
  a <- list(...)
  dim.a.names <- dimnames(a)
  a.names <- names(a)
  dimnames(a) <- dim.a.names
  names(a) <- a.names
  a <- lapply(a, .sims.as.list)
  if (!is.null(SAMPLESIZE)) {
    m <- max(sapply(a, length))
    s <- (sample(1:m, size=SAMPLESIZE, replace=TRUE)-1)
    a <- lapply(a, function (x) x[(s %% length(x))+1])
  }
  a <- .Primitive("c")(FUN = FUN, a, SIMPLIFY = FALSE, USE.NAMES=USE.NAMES)
  a$MoreArgs <- MoreArgs
  S <- do.call(mapply, args = a)
  S <- lapply(S, function (x) if (is.null(x)) NA else x) ## DEBUG:: OK??
  r <- rvsims(S)
  ## DEBUG: match the largest-dimensional param in list(...) and 
  ## set the dimnames if they match
  ## if (isTRUE(all.equal(dim(r), dim(x)))) {
  ##  dimnames(r) <- dimnames(x)
  ##}
  return(r)
}




rvmatrix <- function (data = NA, nrow = 1, ncol = 1, byrow = FALSE, dimnames = NULL) 
{
  data <- as.vector(data)
  if (missing(nrow)) 
      nrow <- ceiling(length(data)/ncol)
  else if (missing(ncol)) 
      ncol <- ceiling(length(data)/nrow)
  if (byrow) {
    X <- t(rvarray(data, c(ncol, nrow), dimnames=dimnames))
  } else {
    X <- rvarray(data, c(nrow, ncol), dimnames=dimnames)
  }
  return(X)
}

rvarray <- function (data = NA, dim = length(data), dimnames = NULL) 
{
  as.rv(array(data = data, dim = dim, dimnames = dimnames))
}

is.matrix.rv <- function (x)
{
  dx <- dim(x)
  return((!is.null(dx)) && length(dx)==2)
}

as.matrix.rv <- function (x, ...) 
{
  if (is.matrix(x)) {
    return(x)
  }
  dn <- if (!is.null(names(x))) list(names(x), NULL) else NULL
  rvarray(x, dim=c(length(x), 1), dimnames=dn)
}


rvmean <- function (x)
{
  UseMethod("rvmean")
}

rvmean.rv <- function (x)
{
  m <- colMeans(sims(x), na.rm=TRUE)
  names(m) <- names(x)
  dim(m) <- dim(x)
  dimnames(m) <- dimnames(x)
  return(m)
}

rvmean.rvsummary <- function (x)
{
  return(x$mean)
}

rvmean.default <- function (x)
{
  if (!is.numeric(x)) {
    x[] <- as.numeric(x)
  }
  return(x)
}

E <- function (x)
{
  rvmean(x)
}

Pr <- function (X)
{
  if (!is.logical(X)) {
      stop("Argument for Pr must be a logical statement such as 'x>0'")
  }
  rvmean(X)
}


rvmin <- function(x)
{
  rvsimapply(x, min, na.rm=FALSE)
}

rvmax <- function(x)
{
  rvsimapply(x, max, na.rm=FALSE)
}

rvrange <- function (x)
{
  rvsimapply(x, range, na.rm=TRUE)
}


rvminmax <- function (x, na.rm=FALSE, finite=FALSE, defaults=c(-Inf,Inf)) # NOEXPORT
{
  .minmax <- function (x, FUN, def) {
    if (na.rm) x <- x[!is.na(x)]
    if (finite) x <- x[is.finite(x)]
    if (length(x)==0) def else FUN(x)
  }
  m0 <- rvsimapply(x, .minmax, min, defaults[1])
  m1 <- rvsimapply(x, .minmax, max, defaults[2])
  list(min=m0, max=m1)
}



rvmultinom <- function(n=1, size=1, prob)
{
  if (length(prob)<=1) {
    return(rvbinom(n=n, size=size, prob=prob))
  }
  if (anyisrv(n, size, prob)) {
    r <- rvmapply(stats:::rmultinom, n=n, size=size, prob=prob)
  } else {
    n.sims <- getnsims()
    s <- rmultinom(n=n*n.sims, size=size, prob=prob)
    dim(s) <- c(length(s) %/% n.sims, n.sims)
    r <- rvsims(t(s))
    dim(r) <- c(length(prob), n)
    dimnames(r) <- list(names(prob), NULL)
  }
  return(r)
}


# ========================================================================
# rvnorm  -  Generate variates from a normal sampling model
# ========================================================================

rvnorm <- function (n=1, mean=0, sd=1, var=NULL, precision)
{
  if (!missing(precision)) {
    if (!is.null(dim(precision))) {
      var <- solve(precision) # matrix inverse
    } else {
      sd <- 1/sqrt(precision)
    }
  }
  if (!is.null(var)) {
    if (!is.null(dim(var))) {
       return(.rvmvnorm(n=n, mean=mean, Sigma=var))
    }
    sd <- sqrt(var)
  }
  rvvapply(stats:::rnorm, n.=n, mean=mean, sd=sd)
}

.mvrnorm <- function (n = 1, mu, Sigma, tol = 1e-06, empirical = FALSE) # NOEXPORT
{ 
  #
  # from the MASS package
  #
    p <- length(mu)
    if (!all(dim(Sigma) == c(p, p)))
        stop("incompatible arguments")
    eS <- eigen(Sigma, sym = TRUE, EISPACK=TRUE)
    if (any(is.na(eS$vectors))) {
      stop("Some eigenvectors are missing...")
    }
    ev <- eS$values
    if (!all(ev >= -tol * abs(ev[1])))
        stop("Sigma is not positive definite")
    X <- matrix(rnorm(p * n), n)
    if (empirical) {
        X <- scale(X, TRUE, FALSE)
        X <- X %*% svd(X, nu = 0)$v
        X <- scale(X, FALSE, TRUE)
    }
    X <- mu + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
    nm <- names(mu)
    if (is.null(nm) && !is.null(dn <- dimnames(Sigma)))
        nm <- dn[[1]]
    dimnames(X) <- list(nm, NULL)
    if (n == 1)
        drop(X)
    else t(X)
}

# ========================================================================
# .rvmvnorm  -  Generate multivariate normal random variables.
# ========================================================================
#
# TODO: 
#  - preserve dimensionality when arrays of n,mean,sd are given.
#  - dim(n) , dim(mean), dim(sd), in this order? at least dim(n) 
#    should be the most important.
#  - parameter cumulative = FALSE
#  - what to do with n rv? conditional d'n given n?
#  - what to do with the parameter 'n'? Should repeat, yes, how to implement?
#

.rvmvnorm <- function (n=1, mean, Sigma) 
{
  if (is.null(dim(Sigma))) {
    Sigma <- diag(Sigma)
  } else if (!is.numeric(Sigma)) 
    stop("Invalid (non-numeric) covariance matrix Sigma")
  else if (nrow(Sigma) != ncol(Sigma)) 
     stop("Invalid (nonsquare) covariance matrix Sigma")
  #
  if (length(mean)==1) {
    mean <- rep(mean, length.out=nrow(Sigma))
  }
  if (length(mean) != nrow(Sigma)) {
    stop("Length of mean vector != nrow of Sigma.")
  }
  #
  nm <- names(mean)
  if (is.rv(n)) stop("n cannot be random (to be implemented)")
  n <- n[1]
  if (is.rv(Sigma)) {
    r <- rvmapply(.mvrnorm, n=n, mu=mean, Sigma=Sigma)
    if (n==1) r <- drop(r)
    if (!is.null(dim(r))) {
      dimnames(r)[[2]] <- nm
    }
    return(r)
  }
  n.sims <- getnsims()
  dim.mean <- dim(mean)
  if (is.list(Sigma)) {
    Sigma <- .expand.as.matrix(Sigma)
  }
  ##if (length(mean) > nrow(Sigma)) {
  ##  # DEBUG: ??
  ##  n.Sigmas <- (length(mean)%/%nrow(Sigma))+(length(mean)%%nrow(Sigma)!=0)
  ##  Sigma <- .expand.as.matrix(rep(rv(Sigma), n.Sigmas))
  ##}
  mean <- rep(mean, length.out = nrow(Sigma))
  mu <- rep(0, length(mean))
  n.all <- n*n.sims
  X <- .mvrnorm(n = n.all, mu = mu, Sigma = Sigma)
  if (n==1) {
    r <- rvsims(X)
  } else {
    dim(X) <- c(n.sims, length(X) %/% n.sims)
    r <- rvsims(X)
    dim(r) <- c(length(mean), n)
  }
  r <- t(r + mean) # Must be in this order to preserve dimensions
  if (n==1) r <- drop(r)
  if (is.null(dim(r))) {
    if (length(r)==length(nm)) names(r) <- nm
  } else {
    dimnames(r) <- list(NULL, nm)
  }
  return(r)
}
as.rvobj <- function (x)
{
  if (is.rvobj(x)) {
    return(x)
  }
  as.rv(x)
}

is.rvobj <- function (x)
{
  return(is.rv(x) || is.rvsummary(x))
}



rvpar <- function (...)
{
  args <- list(...)
  rvpar  <- getOption("rv")
  Pars <- rvpar$par
  if (length(args)==0) {
    return(Pars)
  }
  if (is.null(na <- names(args))) {
    args <- unlist(args)
    p <- Pars[args]
    if (length(args)==1) {
      return(p[[args]])
    } else {
      return(p)
    }
  }
  oldpar <- Pars
  for (a in na) { 
    if (nchar(a)>=1) {
      Pars[a] <- args[a]
    }
  }
  rvpar$par <- Pars
  options(rv=rvpar)
  oldpar
}




# ========================================================================
# 
# ========================================================================

rvpermut <- function (data, prob=NULL)
{
# permutation distribution
  n.sims <- getnsims()
  s <- t(sapply(rep(list(data), n.sims), sample, prob=prob))
  r <- rvsims(s)
  dim(r) <- dim(data)
  names(r) <- names(data)
  return(r)
}



rvpois <- function (n=1, lambda)
{
  rvvapply(stats:::rpois, n.=n, lambda=lambda)
}



rvquantile <- function(x, ...)
{
  UseMethod("rvquantile")
}

rvquantile.rv <- function(x, probs=c(0.025, 0.10, 0.25, 0.50, 0.75, 0.90, 0.975), ignoreInf=FALSE, ...)
{
  if (ignoreInf) {
    .f <- function (x) { quantile(x[is.finite(x)], probs=probs, ..., na.rm=TRUE) }
    t(rvsimapply(x, .f))
  } else {
    t(rvsimapply(x, quantile, probs=probs, ..., na.rm=TRUE))
  }
}

rvquantile.rvsummary <- function(x, probs=c(0.025, 0.10, 0.25, 0.50, 0.75, 0.90, 0.975), ...)
{
  Q <- x$quantiles
  all_probs <- attr(Q, "quantiles")
  M <- NULL
  name <- character(0)
  # if (all(probs %in% all_probs)) ...
  for (p in probs) {
    ix <- (all_probs==p)
    if (any(ix)) {
      M <- cbind(M, Q[,ix,drop=FALSE])
    } else {
      name <- paste(p*100, "%", sep="")
      M <- cbind(M, NA)
      colnames(M)[ncol(M)] <- name
    }
  }
  return(M)
}

rvmedian <- function(x)
{
  UseMethod("rvmedian")
}

rvmedian.rv <- function(x)
{
  rvsimapply(x, median, na.rm=TRUE)
}

rvmedian.rvsummary <- function(x)
{
  rvquantile(x, probs=0.50)
}


rvsample <- function(x, size=1, jointly=TRUE, reject.na=FALSE)
{
  # NAME
  #   rvsample - Draw Samples from Random Vectors
  #
  xs <- sims(as.rv(x))
  ns <- nrow(xs)
  if (is.null(size) || is.na(size)) size <- ns
  if (jointly) {
    if (reject.na) {
      f <- function (x) any(is.na(x))
      is.na.xs <- apply(xs, 1, f)
      if (all(is.na.xs)) {
        s <- sample(ns, size=size, replace=TRUE, prob=is.na.xs)
      } else {
        s <- sample(ns, size=size, replace=TRUE, prob=!is.na.xs) 
      }
    } else {
      s <- sample(ns, size=size, replace=TRUE)
    }
    s <- xs[s,]
  } else {
    s <- apply(xs, 2, function (s) {
      if (all(nas <- is.na(s))) return(s)
      sample(s[!nas], size=size, replace=TRUE)
    })
  }
  return(s)
}


rvvar <- function (x)
{
  UseMethod("rvvar")
}

rvvar.rv <- function (x) # NOEXPORT
{
  S <- sims(x)
  m <- colMeans(S, na.rm=TRUE)
  ns <- rvnsims(x)
  v <- ((colSums(S^2)-ns*(m^2))/(ns-1))
  v[ns==1] <- 0
  names(v) <- names(x)
  dim(v) <- dim(x)
  dimnames(v) <- dimnames(x)
  return(v)
}

rvvar.rvsummary <- function (x) # NOEXPORT
{
  return(x$sd^2)
}

rvvar.default <- function (x) # NOEXPORT
{
  rep.int(0, length(x))
}

rvsd <- function (x)
{
  UseMethod("rvsd")
}

rvsd.rv <- function (x) # NOEXPORT
{
  sqrt(rvvar(x))
}

rvsd.rvsummary <- function (x) # NOEXPORT
{
  return(x$sd)
}

rvsd.default <- function (x) # NOEXPORT
{
  rep.int(0, length(x))
}



# ========================================================================
# rvsimapply  -  apply a function to the simulations, componentwise
# ========================================================================

rvsimapply <- function(x, FUN, ...)
{
  dx <- dim(x)
  n <- length(x)
  if (n==0) {
    return(NULL)
  }
  mv <- lapply(unclass(x), FUN, ...)
  lmv <- sapply(mv, length)
  if (all(lmv==1)) {
    m <- unlist(mv, use.names=TRUE)
    dim(m) <- dx
    dimnames(m) <- dimnames(x)
    return(m)
  } else if (all(lmv==rvnsims(x))) {
    # simulation-wise function was applied - return an object of same type
    attributes(mv) <- attributes(x)
    return(mv)
  } else if (all(lmv==lmv[1])) {
    m <- unlist(mv)
    m <- matrix(m, nrow=lmv[1], ncol=n)
    dimnames(m) <- list(names(mv[[1]]), names(x))
    return(m)
  } else {
    names(mv) <- names(x)
    return(mv)
  }
}

# ========================================================================
# rvsims  -  generate a vector of rvs from a simulation matrix
# ========================================================================
# Name:        
# Description: 
# Parameters:  
# Required:    none
# History:     2004-06-  : 
# If the sims array is 1-dimensional,
# it is taken to be the vector of simulations for one variable;
# if sims is 2-dimensional (L x n)
# it contains L simulations for each of n variables.
# if the sims argument is a list, assumes that each component of the
# list is a single draw from a distribution. A component may be a list,
# when the list will be recursively built into a list of rv objects.
#

.recycle_index <- function (length.from, length.to)
{
  if (length.from>=1 && length.to>=1) {
    return(1+(0:(length.to-1)) %% length.from)
  }
  return(0)
}

resizeSims <- function (s, vector.length, n.sims) # NOEXPORT
{
  # here the simulations must be in columns
  if (!missing(vector.length)) {
    nr <- nrow(s) # params
    if (nr!=vector.length) {
      ix <- .recycle_index(nr, vector.length)
      s <- s[ix, , drop=FALSE]
    }
  }
  nc <- ncol(s) # sims
  if (nc!=1 && !missing(n.sims)) {
    if (nc!=n.sims) {
      ix <- .recycle_index(nc, n.sims)
      s <- s[, ix, drop=FALSE]
    }
  }
  return(s)
}

.rvsims.list <- function (x, n.sims=getnsims(), permute = FALSE) 
{
  # Assume that all elements in the list have the same dimensions.
  # This may be modified later -- filling with NAs
  dx <- dim(x[[1]])
  s <- sapply(x, length)
  if (!all(s==s[1])) {
    # some elements had different dimensions!
    stop("Simulation list was not consistent")
  }
  S <- t(matrix(unlist(x, recursive=FALSE), nrow=s[1]))
  if (!is.list(S)) {
    r <- rvsims(S, n.sims=n.sims, permute=permute)
    dim(r) <- dx
    names(r) <- names(x[[1]])
    return(r)
  }
  r <- list()
  for (i in 1:ncol(S)) {
    r[[i]] <- .rvsims.list(S[,i])
  }
  names(r) <- names(x[[1]])
  return(r)
}

rvsims <- function(sims, n.sims=getnsims(), permute=FALSE)
{
  if (is.list(sims)) return(.rvsims.list(sims, n.sims=n.sims, permute=permute))
  is_factor <- (is.character(sims) || is.factor(sims))
  if (is.factor(sims)) {
    sims[] <- as.character(sims)
  }
  if (!length(dim(sims))%in%c(0,2)) {
    stop("rvsims: Argument must be a vector or matrix or a list")
  }
  if (length(sims)<1) {
    return(rv(0))
  }
  if (is.null(d.s <- dim(sims))) {
    d.s <- dim(sims) <- c(length(sims),1)
  }
  n.sims.max <- d.s[1]
  if (length(d.s)==2) { # A regular matrix of simulations
    .names <- dimnames(sims)[[2]]
  } else {
    stop("Simulation array has >2 dimensions. Don't know how to deal with that.")
  }
  if (n.sims.max>1) {
    if  (nrow(sims)!=n.sims) {
      sims <- t(resizeSims(t(sims), n.sims=n.sims))
    }
    if (permute) {
      .order <- sample(n.sims)
      sims <- sims[.order,,drop=FALSE]
    }
  }
  vec <- split(sims, col(sims))
  names(vec) <- .names
  class(vec) <- class(rv())
  if (is_factor) {
    as.rvfactor(vec)
  } else {
    return(vec)
  }
}

#
# DEBUG:::
#

as.rvsummary <- function (x, ...)
{
  UseMethod("as.rvsummary")
}

is.rvsummary <- function (x)
{
  inherits(x, "rvsummary")
}


print.rvsummary <- function (x, digits=3, ...) # METHOD
{
  for (i in attr(x, "numeric")) {
    x[[i]] <- round(x[[i]], digits=digits)
  }
  print(summary(x))
}

as.rvsummary.rv <- function (x, quantiles=(0:200/200), ...)  # NOEXPORT
{
  y <- if (is.logical(x)) {
    as.rvsummary.rvlogical(x, ...)
  } else if (is.integer(x)) {
    as.rvsummary.rvinteger(x, quantiles=quantiles, ...)
  } else {
    as.rvsummary.rvnumeric(x, quantiles=quantiles, ...)
  }
  y$dimindex <- .dim.index(x)
  y$Rhat <- rvRhat(x)
  y$n.eff <- rvneff(x)
  return(y)
}

as.rvsummary.rvsummary <- function (x, ...)  # NOEXPORT
{
  return(x)
}

as.rvsummary.rvnumeric <- function (x, quantiles=(0:200/200), ...) # NOEXPORT
{
  Q <- rvquantile(x, probs=quantiles)
  fr <- if (any(is.infinite(Q))) rvfiniterange(x) else NULL
  colnames(Q) <- paste(100*quantiles, "%", sep="")
  attr(Q, "quantiles") <- quantiles
  obj <- .rvmeansd(x, names.=c("mean", "sd", "NAS", "n.sims"))
  obj <- c(obj, list(quantiles=Q, attr=attributes(x), finiterange=fr))
  nc <- c("quantiles", "mean", "sd")
  structure(obj, class=c("rvsummary_numeric", "rvsummary"), numeric=nc)
}

as.rvsummary.rvinteger <- function (x, quantiles=(0:200/200), ...) # NOEXPORT
{
  Q <-  rvquantile(x, probs=quantiles)
  fr <- if (any(is.infinite(Q))) rvfiniterange(x) else NULL
  colnames(Q) <- paste(100*quantiles, "%", sep="")
  attr(Q, "quantiles") <- quantiles
  obj <- .rvmeansd(x, names.=c("mean", "sd", "NAS", "n.sims"))
  obj <- c(obj, list(quantiles=Q, attr=attributes(x), finiterange=fr))
  nc <- c("quantiles", "mean", "sd")
  structure(obj, class=c("rvsummary_integer", "rvsummary"), numeric=nc)
}



as.rvsummary.rvlogical <- function (x, ...) # NOEXPORT
{
  obj <- .rvmeansd(x, names.=c("mean", "sd", "NAS", "n.sims"))
  obj <- c(obj, list(quantiles=NULL, attr=attributes(x)))
  structure(obj, class=c("rvsummary_logical", "rvsummary"))
}

as.rvsummary.rvfactor <- function (x, ...) # NOEXPORT
{
  levels <- levels(x)
  llev <- length(levels)
  num.levels <- 1:llev
  #
  S <- sims(x)
  a <- apply(S, 2, function (x) table(c(x, num.levels))) # ensure that all levels are tried
  if (is.null(dim(a))) {
    dim(a) <- c(ncol(S), llev)
  }
  a <- (a-1) # And now subtract the extra counts from the matrix that was obtained.
  ns <- rvnsims(x)
  if (any(naS <- is.na(S))) {
    NAS <- (colMeans(naS)*100)
  } else {
    NAS <- rep.int(0, length(x))
  }
  nax <- if (is.null(dim(x))) NULL else names(x)
  M <- a
  rownames(M) <- levels
  remaining <- (ns-colSums(M))
  if (any(remaining>0)) {
    stop("Impossible: levels won't sum up to 0")
  } 
  P <- t(M/ns)  # compute proportions in each category and transpose
  obj <- list(proportions=P, NAs=NAS, n.sims=ns, attr=attributes(x))
  structure(obj, class=c("rvsummary_rvfactor", "rvsummary"))
}

dim.rvsummary <- function (x)
{
  return(x$attr$dim)
}

"dim<-.rvsummary" <- function (x, value)
{
  d <- (x$attr$dim)
  if (is.null(value)) {
    x$attr$dim <- NULL
    x$dimindex <- .dimindex(x)
    x$attr$dimnames <- NULL
  } else if (is.numeric(value)) {
     if ((is.null(d) && length(x)==prod(value)) || (prod(value)==prod(x$attr$dim))) {
      x$attr$dim <- value
    } else {
      stop("dimensions do not match")
    }
    x$dimindex <- .dimindex(x)
    x$attr$dimnames <- NULL
  }
  return(x)
}

dimnames.rvsummary <- function (x)
{
  return(x$attr$dimnames)
}

"dimnames<-.rvsummary" <- function (x, value)
{
  d <- dim(x)
  if (is.null(value)) {
    x$attr$dimnames <- value
    return(x)
  }
  if (!(is.list(value) && length(value)==length(d))) {
    stop("value is of different length than the dimensions of the array")
  }
  for (i in seq_along(value)) {
    if (is.null(value[[i]])) next
    if (!length(value[[i]])!=d[i]) {
      stop("dimensions of the do not match")
    }
  }
  x$attr$dimnames <- value
  return(x)
}

names.rvsummary <- function (x)
{
  return(x$attr$names)
}

"names<-.rvsummary" <- function (x, value)
{
  if (!(is.null(value) || (length(value)==length(x)))) {
    stop("value is of different length than the dimensions of the array")
  }
  x$attr$names <- value
  return(x)
}

"dimnames<-.rvsummary" <- function (x, value)
{
  d <- dim(x)
  if (!(is.list(value) && length(value)==length(d))) {
    stop("value is of different length than the dimensions of the array")
  }
  for (i in seq_along(value)) {
    if (is.null(value[[i]])) next
    if (length(value[[i]])!=d[i]) {
      stop("dimensions of do not match")
    }
  }
  x$attr$dimnames <- value
  return(x)
}

length.rvsummary <- function (x)
{
  length(x$dimindex)
}

as.double.rvsummary <- function (x, ...)
{
  if (is.null(x$quantiles)) {
    stop("Cannot coerce to double.")
  }
  return(x)
}

print.rvsummary_rvfactor <- function (x, all.levels=FALSE, ...) # METHOD
{
  print(summary(x, all.levels=all.levels, ...))
}


summary.rvsummary <- function (object, ...) # NOEXPORT
{
  # assumes that the 'summary' slot exists in the rvsummary object;
  # this function adds:
  #  1. name, NA%, n.sims
  #  2. Rhat, n.eff
  #  3. dimnames columns
  #
  x <- object
  xdim <- x$attr$dim
  xdimnames <- x$attr$dimnames
  Summary <- x$summary
   if (!is.null(.names <- x$attr$names)) {
    Summary <- cbind(name=.names, Summary)
  }
  Col <- NULL
  n.sims. <- x$n.sims
  if (all(x$NAS==0)) {
    Col <- data.frame(sims=n.sims.)
  } else {
    Col <- data.frame("NA%"=x$NAS, sims=n.sims.)
  }
  if (!all(is.na(x$Rhat))) {
    Col <- cbind(Col, Rhat=x$Rhat)
  }
  if (!all(is.na(x$n.eff))) {
    Col <- cbind(Col, n.eff=x$n.eff)
  }
  Summary <- cbind(Summary, Col)
  if (!is.null(unlist(xdimnames))) {
    # 'is.null(unlist(dimnames))' since we might have e.g. list(NULL, NULL) 
    sud <- rvpar("summary.dimnames")
    if (is.null(sud) || isTRUE(sud)) {
      .f <- function (i) {
        X <- .dimind(dim.=xdim, MARGIN=i)
        na <- xdimnames[[i]]
        if (!is.null(na)) { na <- na[X] }
        return(na)
      }
      da <- lapply(seq_along(xdim), .f)
      names(da) <- names(xdimnames)
      if (is.null(names(da))) {
        names(da) <- if (length(xdim)==2) c("row", "col") else paste("d", seq_along(da), sep="")
      }
      da <- da[!sapply(da, is.null)]
      if (length(da)>0) {
        Summary <- cbind(as.data.frame(da), " "=':', Summary)
      }
    }
  }
  return(Summary)
}

summary.rvsummary_numeric <- function (object, ...) # NOEXPORT
{
  x <- object
  if (is.null(qs <- rvpar("summary.quantiles.numeric"))) {
    qs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  }
  qa <- (attr(x$quantiles, "quantiles")%in%qs)
  q <- x$quantiles[,qa,drop=FALSE]
  S <- data.frame(mean=as.vector(x$m), sd=as.vector(x$sd))
  S <- cbind(S, as.data.frame(q))
  rownames(S) <- x$dimindex
  object$summary <- S
  NextMethod()
}

summary.rvsummary_logical <- function (object, ...) # NOEXPORT
{
  x <- object
  S <- data.frame(mean=as.vector(x$m), sd=as.vector(x$sd))
  rownames(S) <- x$dimindex
  object$summary <- S
  NextMethod()
}

summary.rvsummary_integer <- function (object, ...) # NOEXPORT
{
  x <- object
  if (is.null(qs <- rvpar("summary.quantiles.integer"))) {
    qs <- c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)
  }
  qa <- (attr(x$quantiles, "quantiles")%in%qs)
  q <- x$quantiles[,qa,drop=FALSE]
  names_quantiles <- dimnames(q)[[2]]
  names_quantiles[names_quantiles=="0%"] <- "min"
  names_quantiles[names_quantiles=="100%"] <- "max"
  dimnames(q)[[2]] <- names_quantiles
  S <- data.frame(mean=as.vector(x$m), sd=as.vector(x$sd))
  S <- cbind(S, as.data.frame(q))
  rownames(S) <- x$dimindex
  object$summary <- S
  NextMethod()
}



summary.rvsummary_rvfactor <- function (object, all.levels=TRUE, ...) # NOEXPORT
{
  x <- object
  levels <- x$attr$levels
  llev <- length(levels)
  num.levels <- 1:llev
  #
  maxlev <- if (is.null(maxlev <- rvpar("max.levels"))) { 10 } else maxlev
  too.many.levels.to.show <- ((!all.levels) && (llev>maxlev))
  last.lev.no <- llev
  if (too.many.levels.to.show) {
    P1 <- x$proportions[,1:(maxlev-1),drop=FALSE]
    P2 <- x$proportions[,last.lev.no,drop=FALSE]
    omit_levels <- (!seq_along(levels) %in% c(1:(maxlev-1), last.lev.no))
    rest <- rowSums(x$proportions[,omit_levels,drop=FALSE])
    M <- cbind(P1, "*"=rest, P2)
    colnames(M) <- c(levels[1:(maxlev-1)], "*", levels[last.lev.no])
  } else {
    M <- x$proportions
    colnames(M) <- levels
  }
  S <- cbind(as.data.frame(M), sims=x$n.sims)
  rownames(S) <- x$dimindex
  object$summary <- S
  NextMethod()
}


rvt <- function (n = 1, mu = 0, scale = 1, df, Sigma) 
{
    if (!missing(Sigma)) {
        t <- .rvmvt(n = n, Sigma = Sigma, df = df)
    }
    else {
        t <- rvvapply(stats:::rt, n = n, df = df)
        if (scale != 1) 
            t <- (t * scale)
    }
    if (all(mu != 0)) {
        t <- (t + mu) # t + mu, not mu + t (preserves names)
    }
    return(t)
}



.rvmvt <- function (n=1, Sigma, df=1)
{
  x <- sqrt(rvchisq(n=n, df=df)/df)
  # DEBUG: But will this work? x is of length n,
  #   but the returned thing is of length n * nrow(Sigma)!
  return(.rvmvnorm(n=n, mean=0, Sigma=Sigma)/x)
}


rvunif <- function (n=1, min=0, max=1)
{
  rvvapply(stats:::runif, n.=n, min=min, max=max)
}


rvvapply <- function (FUN, n., ..., constantArgs=list()) # NOEXPORT
{
# 
  args <- list(...)
  n.sims <- getnsims()
  if (length(constantArgs)>0 && any(sapply(constantArgs, is.rv))) {
    stop("constantArgs should contain only constant arguments")
  }
  max_vector_length <- max(sapply(args, length))
  .resize_sims <- function (x, ...) {
    s <- t(sims(as.rv(x)))
    resizeSims(s, ...)
  }
  args <- lapply(args, .resize_sims, vector.length=max_vector_length, n.sims=n.sims)
  n.name <- attr(n., "n.name")
  if (is.null(n.name)) { n.name <- "n" }
  random_n <- (!missing(n.) && is.random(n.))
  if (!missing(n.)) {
    if (length(n.)>1) stop("length(n.)>1? Random dimensions not supported")
    n. <- n.[1]
    args[[n.name]] <- (rvmax(n.)*max_vector_length*n.sims)
  }
  args <- c(args, constantArgs)
  FUN <- match.fun(FUN)
  S <- do.call(FUN, args = args)
  S <- matrix(S, ncol=n.sims)
  if (random_n) {
    n.max <- rvmax(n.)
    ix <- rep(1:n.max, each=max_vector_length)
    ns <- as.vector(sims(n.))
    not_observed <- sapply(ns, function (n) ix>n) # A matrix
    S[not_observed] <- NA
  }
  x <- rvsims(t(S))
  n.scalars <- length(x)
  n.rows <- (n.scalars %/% max_vector_length)
  if (n.rows>1) {
    if (max_vector_length==1) {
      dim(x) <- NULL
    } else {
      dim(x) <- c(max_vector_length, n.rows)
    }
  }
  return(x)
}

# ========================================================================
# simapply  -  apply a (numeric) function to the simulations, rowwise, with dimensions
# ========================================================================
# Vectorizes over the simulations of one single rv; 
# for vectorization over a group of rvs, see 'rvmapply', 'rvvapply'.
#

simapply <- function(x, FUN, ...)
{
  # Works pretty similarly as rvmapply does
  L <- .sims.as.list(x)
  Args <- .Primitive("c")(list(FUN=FUN, SIMPLIFY=FALSE, USE.NAMES=FALSE), list(L))
  S <- do.call(base:::mapply, Args)
  r <- rvsims(S) 
  if (isTRUE(all.equal(dim(r), dim(x)))) {
    dimnames(r) <- dimnames(x)
  }
  return(r)
}


sims <- function(x, ...)
{
  UseMethod('sims')
}

sims.rvsummary <- function(x, dimensions=FALSE,  ...) # NOEXPORT
{
  S <- (x$quantiles)
  attr(S, "quantiles") <- NULL
  if (dimensions) {
    if (!is.null(dim. <- x$attr$dim)) {
      dim(S) <- c(dim., dim(S)[2])
      ndim <- length(dim(S))
      newdim <- c(ndim, 1:(ndim-1))
      S <- aperm(S, newdim)
    }
  }
  return(S)
}

sims.default <- function(x, ...) # NOEXPORT
{
  ##
  ## We should not set the dimension of constants.
  ## re: problems in max.rvsims: cbind(sims(1),sims(x)) does not work
  ## since sims(1) returns a 1x1 matrix but sims(x) may be L x n matrix.
  ## TODO: integrate better!
  ## Need to unclass first, to be consistent.
  ## DEBUG: ? why not just sims(as.rv(x, ...)) ? 
  if (is.atomic(x)) {
    as.vector(unclass(x))  # drop attributes
  } else {
    stop("No method to extract the simulations of an object of class ", paste(class(x), collapse="/"))
  }
}

# ========================================================================
# sims.rv  -  return the matrix of simulations for an rv
# ========================================================================
# gives the simulations of a vector as a matrix,
# used e.g. for computing quantiles of a vector.
# length(x) columns, rvnsims rows.
# length of sims is adjusted to the longest sims length
# if dimensions=TRUE, will return an L x dim(x) matrix.
# (this is useful for simapply applied to random matrices, e.g. see determinant.rv)
#
#
 ###, as.list=FALSE)



sims.rv <- function(x, dimensions=FALSE, n.sims=getnsims(), ...) # NOEXPORT
{
  if (length(x)<1) {
    return(matrix(nrow=n.sims,ncol=0))
  }
  xl <- rvnsims(x)
  if (missing(n.sims)) {
    n.sims <- max(xl)
  }
  class(x) <- NULL
  for (i in which(xl!=n.sims)) {
    x[[i]] <- rep(x[[i]], length.out=n.sims)
  }
  m <- matrix(unlist(x), nrow=n.sims)
  if (dimensions && !is.null(dim.x <- dim(x))) {
    n.s <- (length(m)%/%prod(dim.x))
    cm <- m
    dimnames(cm) <- NULL
    m <- array(m, c(n.s, dim.x)) # multi-way array, 1st dimension is the dimension "rows"
    if (!is.null(dn <- dimnames(x))) {
      dimnames(m)[1+seq_along(dn)] <- dn ## put dimnames in correct positions
    } 
  } else {
    dimnames(m) <- list(NULL, names(x))
  }
  return(m)
}

# ========================================================================
# .sims.as.list  -  split the simulations into a list
# ========================================================================

.sims.as.list <- function (x)
{
  # retain dimensions, and always return getnsims() simulations.
  if (is.null(x)) return(NULL)
  s <- sims(as.rv(x), n.sims=getnsims())
  s <- split(s, row(s)) ## faster than applying list to 1=rows.
  if (!is.null(d <- dim(x))) {
    dn <- dimnames(x)
    s <- lapply(s, function (x) { dim(x) <- d; dimnames(x) <- dn; x})
  }
  # The default names will be "1", "2", ... etc.
  # set names to NULL since this may interfere with mapply( ... ) in "[<-.rv"
  names(s) <- NULL
  return(s)
}


solve.rv <- function (a, b, ...)
{
  rvmapply("solve", a=a, b=b, ...)
}


sort.rv <- function (x, ...) ## EXPORT sort.rv
{
  simapply(x, sort, ...)
}


splitbyname <- function (x)
{
  a <- split(x, f = .shortnames(x))
  lapply(a, .setDimensionByName)
}


summary.rv <- function(object, ...)
{
  summary(as.rvsummary(object), ...)
}

summary.rvfactor <- function(object, all.levels=TRUE, ...)
{
  summary(as.rvsummary(object), all.levels=all.levels, ...)
}



unlist.rv <- function (x, recursive = TRUE, use.names = TRUE) 
{
  # NAME
  #   unlist.rv - Flatten Lists That (May) Contain Random Variables 
  y <- NULL
  ix <- seq(along=x)
  xn <- names(x)
  .paste <- function (name, x) {
     nbrs <- .dim.index(x, leftadjust=FALSE)
     paste(name, nbrs, names(x), sep="")
  }
  if (recursive) {
    for (i in ix) {
      nx <- xn[i]
      if (use.names && is.null(nx)) nx <- "."
      if (!is.rv(x[[i]]) && is.list(x[[i]])) {
        new.y <- unlist.rv(x[[i]], recursive=TRUE, use.names=use.names)
      } else {
        new.y <- x[[i]]
      }
      if (is.null(names(new.y))) {
        new.names <- .paste(nx, new.y)
      } else {
        new.names <- paste(nx, ".", names(new.y), sep="")
      }
      yn <- names(y)
      y <- c(y, new.y)
      names(y) <- c(yn, new.names)
    }
  } else {
    for (i in ix) {
      nx <- xn[i]
      new.y <- x[[i]]
      new.names <- .paste(nx, new.y)
      yn <- names(y)
      y <- c(y, new.y)
      names(y) <- c(yn, new.names)
    }
  }
  return(y)
}
# ========================================================================
# function  -  short description
# ========================================================================
#

# sd, sd.rv

var.default <- getFromNamespace("var", "stats") # EXPORT var.default
formals(var.default) <- c(formals(var.default), alist(... = ))

var.rv <- function(x, ...) 
{
  simapply(x, var.default, ...)
}

sd.rv <- function (x, na.rm = FALSE)
{
  if (is.matrix(x)) 
    apply.rv(x, 2, sd, na.rm = na.rm)
  else if (is.vector(x)) 
    sqrt(var.rv(x, na.rm = na.rm))
  else sqrt(var.rv(as.vector(x), na.rm = na.rm))
}



.onLoad <- function(libname, pkgname) # NOEXPORT
{
  setHook(packageEvent("rv", "detach"), 
    function (...) {
      .rvFunctionSwitch(attach=FALSE, verbose=FALSE)
    }
  )
  rvpar(
    oorm = TRUE,
    rvlwd = 2.0,
    rvcol = "default",
    rvlex = function (lwd) 1.5,
    rvpoint = c("mean", "50%", "95%"),
    point.sample = 400,
    line.sample = 20,
    summary.dimnames = TRUE,
    summary.quantiles.numeric = c(0.01, 0.025, 0.25, 0.50, 0.75, 0.975, 0.99),
    summary.quantiles.integer = c(0, 0.01, 0.025, 0.25, 0.50, 0.75, 0.975, 0.99, 1),
    print.digits = 2
  )
  if (is.null(rvpar("n.sims"))) {
    setnsims(2500)
  }
}

.onAttach <- function(libname, pkgname) # NOEXPORT
{
  .rvFunctionSwitch(attach=TRUE, verbose=FALSE)
  cat("Package rv loaded.\n")
}

