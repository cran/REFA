\name{REFA}
\alias{REFA}
\title{
Robust Exponential Factor Analysis
}
\description{
Robust Exponential Factor Analysis
}
\usage{
REFA(Y, r = 3, tau = 0.75, q = 0.05, eps = 1e-05, init = TRUE)
}
\arguments{
  \item{Y}{Input matrix, of dimension \eqn{T\times N}. Each row is an observation with \eqn{N} features at time point \eqn{t}.}
  \item{r}{A positive integer indicating the factor numbers.}
  \item{q}{Hyper parameter}
  \item{eps}{The stopping criterion parameter. The default is 1e-5.}
  \item{tau}{Hyper parameter}
  \item{init}{Warn start of the algorithm. If \code{init = TRUE}, use modified PCA initialization. If \code{init} is a list contains \code{F0} and \code{L0}, we will use this initialization. Otherwise, use traditional PCA initialization.}
}
\value{
  \item{Fhat}{The estimated factor matrix.}
  \item{Lhat}{The estimated loading matrix.}
  \item{loss}{the value of the loss function.}
}
\references{Manuscript: Robust factor analysis with exponential squared loss}
\author{Jiaqi Hu}

\examples{
\donttest{
dat = gendata()
REFA(dat$Y, r = 3)
}
}

