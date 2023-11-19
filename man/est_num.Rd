\name{est_num}
\alias{est_num}
\title{Estimating Factor Numbers Corresponding  PCA}
\description{Estimating Factor Numbers Corresponding  PCA}
\usage{
est_num(X, kmax = 8, type = "BIC3")
}
\arguments{
  \item{X}{Input matrix, of dimension \eqn{T\times N}. Each row is an observation with \eqn{N} features at time point \eqn{t}.}
  \item{kmax}{The user-supplied maximum factor numbers.}
  \item{type}{the method used.}
}
\value{
the estimated factor numbers
}
\references{Manuscript: Robust factor analysis with exponential squared loss}
\author{Jiaqi Hu}
\examples{
\donttest{
dat = gendata()
est_num(dat$Y)
}
}
