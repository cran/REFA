\name{gendata}
\alias{gendata}
\title{Data generation process
}
\description{Generate heavy-tailed data.
}
\usage{
gendata(seed = 1, T = 50, N = 50, type = "1a")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seed}{the \code{seed} used in the data generation process.
}
  \item{T}{time dimension.
}
  \item{N}{cross-sectional dimension.
}
  \item{type}{the type of the data generation process, it can be \code{"1a", "1b", "1c", "1d", "2a", "2b", "2c", "2d".}
}
}

\value{
a list consisting of \code{Y, F0, L0}.
}
\references{Manuscript: Robust factor analysis with exponential squared loss}
\author{Jiaqi Hu}
\examples{
dat = gendata()
Y = dat$Y
head(Y)
}
