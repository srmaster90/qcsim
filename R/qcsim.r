# Westgard rules simulator for power curves
# Author: S Master
# masters@chop.edu
# code base 2016
# last update 1/22/26


# qcsim function:
#       numiter = number of empirical trials/point
#       maxs = max sigma for calculation (e.g. maxs = 5 gives results from a delta of 0s to 5s)
#       plex = plex of the assay (e.g. 1, 10, 100)
#       levels = # of levels of QC material per time point
#       reps = # of replicates of each QC level per time point
#       xsize = # of discrete points to calculate across the domain
#       sig = significance cutoff for first-level (rule1) QC.  Default is 3s.
#       num.fails = number of failures allowed per multiplex test
#       precision = are we iterating across precision (T) or accuracy/bias (F)
#       [rule flags] = Westgard rules turned on and off
#       [non-Westgard rule flags]
#           ruleT = "T-test" rule (always byLevel)
#           cutoffT = cutoff P-value for ruleT
#           ruleF = "F-test" rule (always byLevel)
#           cutoffF = cutoff P-value for ruleF
#           FTrefLength = reference population length
#           FTtestLength = test population length
#       byLevel = additional lookback via consecutive samples at one level for 2_2S, 4_1S, 10X
#       cores = # of cores for parallel processing
#       seed = seed value for reproducible random sequence
#
# Note: running the function for bias shift returns deltaSE from 0 to maxs,
#   while running for precision shift returns deltaRE from 1 to (maxs + 1)


#' Simulate quality control rule performance
#'
#' @param numiter number of empirical trials/point
#' @param maxs max sigma for calculation (e.g. maxs = 5 gives results from a delta of 0s to 5S)
#' @param plex plex of the assay (e.g. 1, 10, 100)
#' @param levels # of levels of QC material per time point
#' @param reps # of replicates of each QC level per time point
#' @param xsize # of discrete points to calculate across the domain
#' @param sig significance cutoff for first-level (rule1) QC.  Default is 3S.
#' @param num.fails number of failures allowed per multiplex test
#' @param precision are we iterating across precision (T) or accuracy/bias (F)
#' @param ruleP_2S Westgard preliminary 2s threshold
#' @param rule1_sigS Westgard threshold at "sig" SD (default 3S)
#' @param rule2_2S Westgard 2_2S rule
#' @param ruleR_4S Westgard R_4S rule
#' @param rule4_1S Westgard 4_1S rule
#' @param rule10X Westgard 10x rule
#' @param ruleT "T-test" rule (always byLevel)
#' @param cutoff.T cutoff P-value for ruleT
#' @param ruleF "F-test" rule (always byLevel)
#' @param cutoff.F cutoff P-value for ruleF
#' @param FTrefLength reference population length (T and F test)
#' @param FTtestLength test population length (T and F test)
#' @param byLevel additional lookback via consecutive samples at one level for 2_2S, 4_1S, 10X
#' @param cores # of cores for parallel processing
#' @param seed seed value for reproducible random sequence
#'
#' @return qcsim object
#'
#' import(Rcpp)
#' @importFrom doRNG %dorng%
#' @importFrom foreach foreach
#' @importFrom Rcpp evalCpp
#' @importFrom graphics lines
#' @importFrom stats qf
#' @importFrom stats qt
#' @importFrom stats rnorm
#'
#' @useDynLib qcsim, .registration=TRUE
#' 
#' @export
#'
#' @examples
#'
#' iters <- 51 # number of x points
#' maxs <- 5
#'
#' # Note: running the function for bias shift returns deltaSE from 0 to maxs,
#' #   while running for precision shift returns deltaRE from 1 to (maxs + 1)
#'
#' par(mfrow=c(1,2))
#'
#' # 3 SD rule
#' pres.bias.3s <- qcsim(maxs=maxs, xsize=iters)
#' plot(pres.bias.3s, col=1)
#' pres.prec.3s <- qcsim(maxs=maxs, xsize=iters, precision=TRUE)
#' plot(pres.prec.3s, col=1)
#'
#' # Full Westgard rules
#' pres.bias.W <- qcsim(maxs=maxs, xsize=iters,
#'                      rule1_sigS=TRUE,rule2_2S=TRUE,ruleR_4S=TRUE,rule4_1S=TRUE,rule10X=TRUE)
#' plot(pres.bias.W)
#' pres.prec.W <- qcsim(maxs=maxs, xsize=iters, precision=TRUE,
#'                      rule1_sigS=TRUE,rule2_2S=TRUE,ruleR_4S=TRUE,rule4_1S=TRUE,rule10X=TRUE)
#' plot(pres.prec.W)
#'
#'
#' # Overlay for comparison
#' plot(pres.bias.3s)
#' lines(pres.bias.W, col=2)
qcsim <- function(numiter=200000,maxs=5,plex=1,levels=1,reps=1,xsize=51,sig=3,num.fails=1,precision=FALSE,
                  ruleP_2S=TRUE, rule1_sigS=TRUE,rule2_2S=FALSE,ruleR_4S=FALSE,rule4_1S=FALSE,rule10X=FALSE,
                  ruleT=FALSE, cutoff.T=0, ruleF=FALSE, cutoff.F=0, FTrefLength=10, FTtestLength=10,
                  byLevel=TRUE,cores=1, seed=NULL) {

  doFuture::registerDoFuture(flavor=("%dofuture%"))
  options(future.globals.onReference = "error")
  future::plan(future::multisession, workers=cores)

  qcargs <- list(numiter=numiter,maxs=maxs,plex=plex,levels=levels,reps=reps,xsize=xsize,sig=sig,num.fails=num.fails,precision=precision,
                 ruleP_2S=ruleP_2S, rule1_sigS=rule1_sigS,rule2_2S=rule2_2S,ruleR_4S=ruleR_4S,rule4_1S=rule4_1S,rule10X=rule10X,
                 ruleT=ruleT, cutoff.T=cutoff.T, ruleF=ruleF, cutoff.F=cutoff.F, FTrefLength=FTrefLength, FTtestLength=FTtestLength,
                 byLevel=byLevel,cores=cores, seed=seed)

  run.size <- levels * reps
  if (rule10X) {
    buffer.W <- 9 * run.size  # enough buffer to allow lookback for 10x rule -- buffer.W covers the multiday Westgard rules
  } else if (rule4_1S) {
    buffer.W <- 3 * run.size
  } else if (rule2_2S) {
    buffer.W <- run.size
  } else {
    buffer.W <- 0
  }
  buffer.max <- buffer.W
  if (ruleT | ruleF) {
    if ((FTtestLength < 2) | (FTrefLength < 2)) {
      error("T-test and F-test rules require test lengths > 1")
    }
    buffer.TF <- FTtestLength * run.size
    if (buffer.TF > buffer.max) {
      buffer.max <- buffer.TF
    }
    t.cutoff <- abs(qt(cutoff.T / 2, df = FTrefLength + FTtestLength - 2))    # convert P value to a T-distribution threshold
    f.cutoff.low <- qf(cutoff.F / 2, df1 = FTrefLength - 1, df2 = FTtestLength - 1)  # convert P value to F-distribution threshold
    f.cutoff.high <- qf(1 - (cutoff.F / 2), df1 = FTrefLength - 1, df2 = FTtestLength - 1)  # convert P value to F-distribution threshold
  }

  # main processing

  if (!is.null(seed)) {
    set.seed(seed)
  }

  #  qcresults <- foreach (z=1:xsize,.combine='c') %do% {   # single thread
  qcresults <- foreach (z=1:xsize,.combine='c', .packages=c("qcsim")) %dorng% {  # multicore
    between.day.rules <- function(idx, rule2_2S, rule4_1S, rule10X) {   # function to more efficiently consolidate cross-day rules
      # idx is list of indices for 10X rule
      failvec <- rep(FALSE, dim(S)[2])
      if (rule10X) {
        Ssub <- S[idx,, drop = FALSE]
        failvec <- failvec | cutoffContig(Ssub, 0, 10)
      }
      idx <- idx[7:length(idx)]  # take subset for 4_1S rule
      if (rule4_1S) {
        Ssub <- S[idx,, drop = FALSE]
        failvec <- failvec | cutoffContig(Ssub, 1, 4)
      }
      idx <- idx[3:length(idx)]  # take subset for 2_2S rule
      if (rule2_2S) {
        Ssub <- S[idx,, drop = FALSE]
        failvec <- failvec | cutoffContig(Ssub, 2, 2)
      }
      return(failvec)
    }

    results <- vector(mode="numeric",numiter)

    Srows <- (numiter*run.size)+buffer.max
    S <- matrix(data=rnorm(Srows*plex),nrow=Srows,ncol=plex)

    if (ruleT | ruleF) {  # populate reference data
      ref.TF <- array(data=rnorm(FTrefLength*plex*levels),
                      dim=c(FTrefLength,plex,levels))
    }

    # only one value in the multiplex is truly out of control
    if (precision) {
      S[,1] <- S[,1] * (1 + ((z - 1) * (maxs / (xsize - 1))))
    } else {
      S[,1] <- S[,1] + (z - 1) * (maxs / (xsize - 1))
    }
    for (k in 1:numiter) {
      rulesv <- rep(TRUE,plex)
      wr.start <- buffer.max + ((k - 1) * run.size) + 1  # first row for this within-run block
      within.run.matrix <- S[wr.start:(wr.start+run.size-1),, drop = FALSE]
      wr.max <- colMax(within.run.matrix)
      wr.min <- colMin(within.run.matrix)
      wr.abs <- abs(c(wr.max, wr.min))

      pass <- TRUE
      failvec <- rep(FALSE,plex)  # vector of pass/fails that we'll carry along
      prelim.failvec <- !failvec  # for prelim 2S rule

      # Preliminary 2S rule to determine if results can be reported
      if (ruleP_2S) {
        prelim.failvec <- (wr.abs > 2)
        if (all(!prelim.failvec)) {
          pass <- TRUE
          results[k] <- pass
          next
        }
      }

      # typically 1_3S rule
      if (rule1_sigS & pass) {
        failvec <- failvec | (wr.abs > sig)
        if (sum(failvec) >= num.fails) {
          pass <- FALSE
        }
      }

      # R 4s rule
      if (ruleR_4S & pass) { # R 4s rule
        failvec <- failvec | ((wr.max - wr.min) > 4)
        if (sum(failvec) >= num.fails) {
          pass <- FALSE
        }
      }

      # remaining Westgard rules, by level
      if ((rule2_2S | rule4_1S | rule10X) & pass & byLevel) {
        pull.chunk <- ceiling(9/reps) + 1
        pull.length <- pull.chunk * reps
        pull.offset <- -(pull.chunk - 1)
        for (j in 1:levels) {
          # build set of row indices to lookback by level
          idx <- (wr.start + ((j-1) * reps) + rep((pull.offset:0)*run.size, each=reps) + rep(1:reps) - 1)[(pull.length-reps-8):pull.length] # 8 comes from 10 (10x) - 2
          failvec <- failvec | between.day.rules(idx, rule2_2S, rule4_1S, rule10X)
          if (sum(failvec) >= num.fails) {
            pass <- FALSE
            break
          }
        }
      }

      # remaining Westgard rules, all levels consecutive
      else if ((rule2_2S | rule4_1S | rule10X) & pass) {
        idx <- (wr.start - 9):(wr.start + run.size - 1)
        failvec <- failvec | between.day.rules(idx, rule2_2S, rule4_1S, rule10X)
      }

      if ((ruleT | ruleF) & pass) {
        pull.chunk <- ceiling(FTtestLength/reps)  # how many total runs we have to cross
        pull.length <- pull.chunk * reps
        pull.offset <- -(pull.chunk - 1)
        for (j in 1:levels) {
          # build set of row indices to lookback by level
          idx <- (wr.start + ((j-1) * reps) + rep((pull.offset:0)*run.size, each=reps) + rep(1:reps) - 1)
          if (length(idx) > FTtestLength) {
            idx <- idx[(FTtestLength-length(idx)+1):length(idx)]
          }
          # for T and F test, evaluate once per QC run
          ref.TF.level <- ref.TF[,,j, drop = TRUE]
          if (plex == 1) {
            ref.TF.level <- matrix(ref.TF.level, ncol=1)  # if single-plex, the last drop=T makes this a vector instead of a 1-column matrix
          }
          SSub <- rbind(ref.TF.level, S[idx,, drop = FALSE])
          if (ruleT & ruleF) {
            failvec <- failvec | tftestMat(ref.TF.level, S[idx,, drop = FALSE],
                                           useT=TRUE, useF=TRUE,
                                           t.cutoff, f.cutoff.low, f.cutoff.high)
          } else if (ruleT) {
            failvec <- failvec | tftestMat(ref.TF.level, S[idx,, drop = FALSE],
                                           useT=TRUE, useF=FALSE,
                                           t.cutoff, f.cutoff.low, f.cutoff.high)
          } else if (ruleF) {
            failvec <- failvec | tftestMat(ref.TF.level, S[idx,, drop = FALSE],
                                           useT=FALSE, useF=TRUE,
                                           t.cutoff, f.cutoff.low, f.cutoff.high)
          }
          if (sum(failvec) >= num.fails) {
            pass <- FALSE
            ref.TF <- array(data=rnorm(FTrefLength*plex*levels),
                            dim=c(FTrefLength,plex,levels))
            break
          }
        }
      }

      if (sum(failvec & prelim.failvec) >= num.fails) {
        pass <- FALSE
      }

      if ((!pass) & (buffer.max > 0)) {
        reset.start <- wr.start - buffer.max + run.size
        reset.end <- wr.start + run.size - 1
        reset.length <- reset.end - reset.start + 1
        S[reset.start:reset.end, 1:plex] <- rnorm(reset.length * plex)
        # ...and alter the "out of control" buffer
        if (precision) {
          S[reset.start:reset.end,1] <- S[reset.start:reset.end,1] * (1 + ((z - 1) * (maxs / (xsize - 1))))
        } else {
          S[reset.start:reset.end,1] <- S[reset.start:reset.end,1] + (z - 1) * (maxs / (xsize - 1))
        }
      }
      results[k] <- pass
    }
    1 - (sum(results) / numiter)
  }
  if (precision) {
    xvals <- 1+((1:xsize)-1)*(maxs / (xsize - 1))
  } else {
    xvals <- ((1:xsize)-1)*(maxs / (xsize - 1))
  }
  output <- list(x = xvals, y = as.numeric(qcresults), args = qcargs)
  class(output) <- "qcsim"
  return(output)
}

#' Plot QC power curve
#'
#' @param x object returned from qcsim
#' @param lwd line width
#' @param ... arguments to be passed to methods, such as graphical parameters.  See base::plot() for details.
#'
#' @export
#'
#' @examples
#'
#' iters <- 51 # number of x points
#' maxs <- 5
#'
#' # 3 SD rule
#' pres.bias.3s <- qcsim(maxs=maxs, xsize=iters)
#' plot(pres.bias.3s, col=1)
#'

plot.qcsim <- function(x, lwd = 2, ...) {
  if (x$args$precision) {xtext <- "RE (multiples of S)"} else {xtext <- "SE (multiples of S)"}
  plot(x$x, x$y, ylab="Probability of rejection (P)", xlab=bquote(Delta * .(xtext)), ylim=c(0,1), type="l", lwd=lwd, ...)
}

#' Add QC power curve
#'
#' @param x object returned from qcsim
#' @param lwd line width
#' @param ... arguments to be passed to methods, such as graphical parameters.  See base::plot() for details.
#'
#' @export
#'
#' @examples
#'
#' iters <- 51 # number of x points
#' maxs <- 5
#'
#' pres.bias.3s <- qcsim(maxs=maxs, xsize=iters)
#' pres.bias.W <- qcsim(maxs=maxs, xsize=iters,
#'                      rule1_sigS=TRUE,rule2_2S=TRUE,ruleR_4S=TRUE,rule4_1S=TRUE,rule10X=TRUE)
#'
#' # Overlay for comparison
#' plot(pres.bias.3s)
#' lines(pres.bias.W, col=2)

lines.qcsim <- function(x, lwd = 2, ...) {
  lines(x$x, x$y, lwd=lwd, ...)
}

######  Example use

# iters <- 51 # number of x points
# maxs <- 5
#
# par(mfrow=c(1,2))
#
# # 3s rule
# system.time(pres.bias.3s <- qcsim(maxs=maxs, xsize=iters))
# plot(pres.bias.3s, col=1)
# system.time(pres.prec.3s <- qcsim(maxs=maxs, xsize=iters, precision=T))
# plot(pres.prec.3s, col=1)
#
# # Full Westgard
# system.time(pres.bias.W <- qcsim(maxs=maxs, xsize=iters,
#                           rule1_sigS=T,rule2_2S=T,ruleR_4S=T,rule4_1S=T,rule10X=T))
# plot(pres.bias.W)
# system.time(pres.prec.W <- qcsim(maxs=maxs, xsize=iters, precision=T,
#                            rule1_sigS=T,rule2_2S=T,ruleR_4S=T,rule4_1S=T,rule10X=T))
# plot(pres.prec.W)
#
#
# # Overlay for comparison
# plot(pres.bias.3s)
# lines(pres.bias.W, col=2)
