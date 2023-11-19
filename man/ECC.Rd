\name{ECC}
\alias{ECC}
\title{Estimation of errors for common component}
\description{Estimation of errors for common component}
\usage{
ECC(Chat, C)
}
\arguments{
  \item{Chat}{The estimated common component}
  \item{C}{The true common component}
}
\value{
a numeric value of the ECC
}
\references{Manuscript: Robust factor analysis with exponential squared loss}
\author{Jiaqi Hu}
\examples{
\donttest{
dat = gendata()
Y = dat$Y
F0 = dat$F0
L0 = dat$L0
C0 = F0 %*% t(L0)
res = REFA(dat$Y, r = 3)
Fhat = res$Fhat
Lhat = res$Lhat
Chat = Fhat %*% t(Lhat)
ECC(Chat, C0)
}
}

