\name{TR}
\alias{TR}
\title{Trace ratios}
\description{Trace ratios}
\usage{
TR(Fhat, F0)
}
\arguments{
\item{Fhat}{The estimated factors.}
  \item{F0}{The true factors.}
}
\value{
a numeric value of the trace ratios.
}
\references{Manuscript: Robust factor analysis with exponential squared loss}
\author{Jiaqi Hu}
\examples{
\donttest{
dat = gendata()
Y = dat$Y
F0 = dat$F0
res = REFA(dat$Y, r = 3)
Fhat = res$Fhat
TR(Fhat, F0)
}
}
