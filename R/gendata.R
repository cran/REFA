#' Data generation process
#'
#' @description Generate heavy-tailed data.
#'
#' @param seed The \code{seed} used in the data generation process.
#' @param T Time dimension.
#' @param N Cross-sectional dimension.
#' @param type The type of the data generation process. Options include "1a", "1b", "1c", "1d", "2a", "2b", "2c", "2d", "3a", "3b".
#'
#' @return A list containing:
#' \item{Y}{The observed matrix.}
#' \item{F0}{The true factor matrix.}
#' \item{L0}{The true loading matrix.}
#'
#' @references Robust factor analysis with exponential squared loss. Jiaqi Hu, Tingyin Wang, Xueqin Wang. Journal of Multivariate Analysis 2026, 213 105567; doi:10.1016/j.jmva.2025.105567
#' @author Jiaqi Hu
#'
#' @examples
#' dat = gendata()
#' Y = dat$Y
#' head(Y)
#'
#' @importFrom stats rnorm rt
#' @importFrom fMultivar rmvst
#' @export
gendata <- function(seed = 1, T = 50, N = 50, type = "1a") {
  allowed_types <- c("1a", "1b", "1c", "1d", "2a", "2b", "2c", "2d", "3a", "3b")
  stopifnot(type %in% allowed_types)
  set.seed(seed)
  L0 = matrix(rnorm(N * 3), N, 3)
  F0 = matrix(0, T, 3)

  if(type == "3a" | type == "3b"){

    if(type == "3a"){
      Omega = diag(1, 3 + N, 3 + N)
    }else if(type == "3b"){
      Omega = outer(1:(N+3), 1:(N+3), function(x, y)0.5^abs(x - y))
    }
    FE = rmvst(n = T, dim = 3 + N, Omega = Omega, df = 3)
    F0 = FE[, 1:3]
    e = FE[, 4:(N+3)]

  }else{

    F0[1, ] = rnorm(3)
    for (t in 2:T) {
      F0[t, ] = F0[t - 1] * c(0.8, 0.5, 0.2) + rnorm(3)
    }
    if (type == "1a") {
      e = matrix(rnorm(T * N), T, N)
    } else if (type == "1b") {
      e = matrix(rt(T * N, 3), T, N)
    } else if (type == "1c") {
      e = matrix(rt(T * N, 2), T, N)
    } else if (type == "1d") {
      e = matrix(rt(T * N, 1), T, N)
    } else {
      J = 3
      if (type == "2a") {
        beta = 0
        rho = 0
        v = matrix(rt(T * N, 2), T, N)
      } else if (type == "2b") {
        beta = 0
        rho = 0
        v = matrix(rt(T * N, 1), T, N)
      } else if (type == "2c") {
        beta = 0.2
        rho = 0.2
        v = matrix(rt(T * N, 3), T, N)
      } else if (type == "2d") {
        beta = 0.2
        rho = 0.2
        v = matrix(rt(T * N, 2), T, N)
      }
      g = sqrt(F0[, 1]^2 + F0[, 2]^2 + F0[, 3]^2)
      e = matrix(0, T, N)
      e[1, ] = rnorm(N) * g[1]
      for (t in 2:T) {
        for (i in 1:N) {
          e[t, i] = rho * e[t - 1, i] + (1 - beta) * v[t, i] + beta * sum(v[t, (max(1, i - J)):(min(N, i + J))])
        }
      }
      e[t, ] = g[t] * e[t, ]
    }
  }

  res = normalize(F0, L0, r = 3)
  F0 = res$Fhat
  L0 = res$Lhat
  Y = F0 %*% t(L0) + e
  return(list(Y = Y, F0 = F0, L0 = L0))
}
