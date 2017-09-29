\name{PH1XBAR-package}
\alias{00_PH1XBAR-package}
\docType{package}
\title{
Phase 1 X-bar Control Chart
}
\description{
The utlity of this package is to build a control chart by multivariate normal distribution.
}
\details{
\tabular{ll}{
Package: \tab PH1XBAR\cr
Type: \tab Package\cr
% Version: \tab 1.0\cr
Date: \tab 2017-09-29\cr
License: \tab GPL (>= 3)\cr
}
Function 'PH1XBAR' builds a control chart by multivariate normal distribution.
Function 'get.cc' gets a charting constant by multivariate normal distribution.
}
\author{ Yuhui Yao, Subha Chakraborti.

Maintainer: Yuhui Yao <yyao17@crimson.ua.edu>
}
\references{
Champ CW, Jones LA. Designing Phase I X-bar charts with small sample sizes. Quality and Reliability Engineering International. 2004;20:497-510
}
\examples{
%\dontrun{
# Build a control chart
X <- PH1XBAR.data()
PH1XBAR(X)

# Get a charting constant
get.cc(m = 10, nu = 40, FAP = 0.05)
%}
}
