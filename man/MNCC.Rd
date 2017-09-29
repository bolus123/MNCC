\name{MNCC}
\alias{MNCC}
\title{MNCC}
\usage{
MNCC(X, FAP = 0.1, off.diag = -1/(m - 1), alternative = '2-sided',
plot.option = TRUE, maxiter = 10000)
}
\arguments{
  \item{X}{ input and it must be a matrix }
  \item{FAP}{ false Alarm Probabilty in Phase 1 }
  \item{off.diag}{ value of equal-correlation. Default setting is -1/(m - 1) }
  \item{alternative}{ '2-sided' - equal-limits applied in the calculation; '1-sided' - only upper limit involved }
  \item{plot.option}{TRUE - draw a plot for the process; FALSE - Not draw a plot for the process}
  \item{maxiter}{ maximum iteration of computation in mvtnorm }
}
\description{
Build a control chart based on multivariate normal distribution.
}
\examples{
X <- MNCC.data()
MNCC(X)
}
