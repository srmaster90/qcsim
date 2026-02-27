\name{qcsim-package}
\alias{qcsim-package}
\alias{qcsim}
\docType{package}
\title{
  Monte Carlo simulation for evaluating clinical laboratory QC rules
}
\description{
  This package simulates laboratory quality control data and the effects of various combinations of QC rules,
  including Westgard rules, cutoffs at arbitrary significance (such as 5 SD), as well as t and F tests.
  It calculates power curves that assess the probability of the QC rule rejecting a run as function of
  changes in the underlying bias and imprecision of the a simulated analytes.  It supports arbitrary numbers
  of levels, QC replicates, and analytes (i.e. single- or multiplex assays).
}
\details{
  The qcsim function is used to simulate the rejection power of QC rules as a function of changes in the
  underlying bias or imprecision the simulated analyte.  Plotting functions are also supplied to handle the
  results returned by qcsim().
  Note: this package does not come with any express or implied warranty of suitability for any purpose.
  You are responsible for understanding and verifying any results used for clinical implementation.
}
\author{
Stephen Master, masters@chop.edu

Maintainer: Stephen Master <masters@chop.edu>
}
\references{
  Unsihuay D and Master SR, 2026, in preparation
}
\keyword{ package }
\seealso{
  Optional links to other man pages
}
\examples{
  \donttest{
    # 3 SD rule
    pres.bias.3s <- qcsim(maxs = maxs, xsize = iters)
    plot(pres.bias.3s, col = 1)
    pres.prec.3s <- qcsim(maxs = maxs, xsize = iters, precision = TRUE)
    plot(pres.prec.3s, col = 1)
    
    # Full Westgard rules
    pres.bias.W <- qcsim(
      maxs = maxs, xsize = iters,
      rule1_sigS = TRUE, rule2_2S = TRUE, ruleR_4S = TRUE, rule4_1S = TRUE, rule10X = TRUE
    )
    plot(pres.bias.W)
    pres.prec.W <- qcsim(
      maxs = maxs, xsize = iters, precision = TRUE,
      rule1_sigS = TRUE, rule2_2S = TRUE, ruleR_4S = TRUE, rule4_1S = TRUE, rule10X = TRUE
    )
    plot(pres.prec.W)
    
    # Overlay for comparison
    plot(pres.bias.3s)
    lines(pres.bias.W, col = 2)
  }
}
