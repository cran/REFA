\name{FA}
\alias{FA}
\title{
Principal Component Analysis for Factor Models
}
\description{
Principal Component Analysis for Factor Models
}
\usage{
FA(X, r)
}
\arguments{
  \item{X}{Input matrix, of dimension \eqn{T\times N}. Each row is an observation with \eqn{N} features at time point \eqn{t}.}
  \item{r}{A positive integer indicating the factor numbers.}
}
\value{
  \item{Fhat}{The estimated factor matrix.}
  \item{Lhat}{The estimated loading matrix.}
}
\references{Manuscript: Robust factor analysis with exponential squared loss}
\author{Jiaqi Hu}
\examples{
##---- Should be DIRECTLY executable !! ----


}

