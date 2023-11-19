\name{REFA_FN}
\alias{REFA_FN}
\title{
Estimating Factor Numbers via Modified Rank Minimization
}
\description{
Estimating Factor Numbers via Modified Rank Minimization
}
\usage{
REFA_FN(Y, rmax = 8, tau = 0.75, q = 0.1, eps = 1e-04, init = TRUE)
}
\arguments{
  \item{Y}{Input matrix, of dimension \eqn{T\times N}. Each row is an observation with \eqn{N} features at time point \eqn{t}.}
  \item{rmax}{The bound of the number of factors.}
  \item{q}{Hyper parameter in modified PCA algorithm. Default is \code{0.05}.}
  \item{eps}{The stopping criterion parameter. Default is \code{1e-5}.}
  \item{tau}{Hyper parameter in selecting \eqn{\gamma} of the loss function.}
  \item{init}{Warn start by modified PCA algorithm. Default is \code{TRUE}.}
}

\value{
  \item{rhat}{The estimated factor number.}
  \item{Fhat}{The estimated factor matrix.}
  \item{Lhat}{The estimated loading matrix.}
  \item{loss}{the value of the loss function.}

}
\references{Manuscript: Robust factor analysis with exponential squared loss}
\author{Jiaqi Hu}
\examples{
\donttest{
dat = gendata()
REFA_FN(dat$Y, rmax = 8)
}
}

