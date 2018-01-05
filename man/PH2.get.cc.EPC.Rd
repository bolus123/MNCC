\name{PH2.get.cc.EPC}
\alias{PH2.get.cc.EPC}
\title{PH2.get.cc.EPC}
 \usage{
PH2.get.cc.EPC(
	p, 
	k, 
	nu = k - 1, 
	eps = 0.1, 
	ARL0 = 370, 
	c4.option = TRUE, 
	interval = c(1, 10), 
	u = runif(100000), 
	v = runif(100000)
)
}
\arguments{
  \item{p}{ exceedance probability }
  \item{k}{ number of subgroups }
  \item{nu}{ degrees of freedom }
  \item{eps}{ a nominal allowable changes of the in-control average run length }
  \item{ARL0}{ a nominal in-control average run length }
  \item{c4.option}{ TRUE - the standard deviation estimator corrected by c4. FALSE - no c4 involved }
  \item{interval}{ the searching range of lower prediction bound of average run length }
  \item{u}{ a seqeuence between 0 and 1, and u and v are mutually independent }
  \item{v}{ a seqeuence between 0 and 1, and u and v are mutually independent }
}
\description{
Obtain a Phase II charting constant by exceedance probability criterion.
}
\examples{
PH2.get.cc.EPC(p = 0.1,	k = 30)
}
