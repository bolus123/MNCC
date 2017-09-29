\name{get.cc}
\alias{get.cc}
\title{get.cc}
 \usage{
get.cc(m, nu, FAP = 0.1, off.diag = -1/(m - 1), alternative = '2-sided',
maxiter = 10000)
}
\arguments{
  \item{m}{ number of subgroups }
  \item{nu}{ degrees of freedom }
  \item{FAP}{ false Alarm Probabilty in Phase 1 }
  \item{off.diag}{ value of equal-correlation. Default setting is -1/(m - 1) }
  \item{alternative}{ '2-sided' - equal-limits applied in the calculation; '1-sided' - only upper limit involved }
  \item{maxiter}{ maximum iteration of computation in mvtnorm }

}
\description{
Obtain a charting constant by multivariate normal model.
}
\examples{
get.cc(m = 10, nu = 40, FAP = 0.05)
}
