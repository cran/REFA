#' Robust Exponential Factor Analysis
#'
#' @description Robust Exponential Factor Analysis
#'
#' @param Y Input matrix, of dimension \eqn{T\times N}. Each row is an observation with \eqn{N} features at time point \eqn{t}.
#' @param r A positive integer indicating the factor numbers.
#' @param tau Hyper parameter in selecting \eqn{\gamma} of the loss function.
#' @param q Hyper parameter used in the initializations, truncated PCA.
#' @param eps The stopping criterion parameter. The default is 1e-5.
#' @param init Warm start of the algorithm. If \code{init = TRUE}, use truncated PCA initialization. If \code{init} is a list contains \code{F0} and \code{L0}, we will use this initialization. Otherwise, use traditional PCA initialization.
#'
#' @return A list containing:
#' \item{Fhat}{The estimated factor matrix.}
#' \item{Lhat}{The estimated loading matrix.}
#' \item{loss}{The value of the loss function.}
#'
#' @references Robust factor analysis with exponential squared loss. Jiaqi Hu, Tingyin Wang, Xueqin Wang. Journal of Multivariate Analysis 2026, 213 105567; doi:10.1016/j.jmva.2025.105567
#' @author Jiaqi Hu
#'
#' @examples
#' \donttest{
#' # Assuming gendata() is defined in your package
#' dat = gendata()
#' REFA(dat$Y, r = 3)
#' }
#' @importFrom stats median quantile rnorm nlm
#' @export
REFA <- function(Y, r = 3, tau = 0.75, q = 0.05, eps = 1e-05, init = TRUE) {
  T = nrow(Y)
  N = ncol(Y)
  # initial L
  if (is.list(init)) {
    L = init$L0
    F = init$F0
  } else if (init == TRUE) {
    med = median(Y)
    lower = quantile(Y, q)
    upper = quantile(Y, 1 - q)
    Y_normal = Y * ((Y >= lower) & (Y <= upper)) + lower * (Y < lower) + upper * (Y > upper)
    sigma = 1/T * t(Y_normal) %*% Y_normal
    L = eigen(sigma)$vectors[, 1:r] * sqrt(N)
    L = as.matrix(L)
    F = 1/N * Y_normal %*% L
    res = normalize(F, L, r)
    F = res$Fhat
    L = res$Lhat
  } else {
    L = matrix(rnorm(N * r), N, r)
    F = 1/N * Y %*% L
    res = normalize(F, L, r)
    F = res$Fhat
    L = res$Lhat
  }
  gap = 1
  gamma = quantile((Y - F %*% t(L))^2, tau)
  loss = mean(1 - exp(-1/gamma * (Y - F %*% t(L))^2))
  while (gap > eps) {
    F_new = matrix(0, T, r)
    for (t in 1:T) {
      F_new[t, ] = nlm(update_F, p = F[t, ], y = Y[t, ], L = L, gamma = gamma, gradtol = 1e-10)$estimate
    }
    L_new = matrix(0, N, r)
    for (i in 1:N) {
      L_new[i, ] = nlm(update_L, p = L[i, ], y = Y[, i], F = F_new, gamma = gamma, gradtol = 1e-10)$estimate
    }
    loss_new = mean(1 - exp(-1/gamma * (Y - F_new %*% t(L_new))^2))
    gap = loss - loss_new
    L = L_new
    F = F_new
    res = normalize(F, L, r)
    F = res$Fhat
    L = res$Lhat
    loss = loss_new
  }
  return(list(Fhat = F, Lhat = L, loss = loss))
}

#' Estimating Factor Numbers via Modified Rank Minimization
#'
#' @description Estimating Factor Numbers via Modified Rank Minimization
#'
#' @param Y Input matrix, of dimension \eqn{T\times N}. Each row is an observation with \eqn{N} features at time point \eqn{t}.
#' @param rmax The bound of the number of factors.
#' @param tau Hyper parameter in selecting \eqn{\gamma} of the loss function.
#' @param q Hyper parameter used in the initializations, truncated PCA. Default is \code{0.05}.
#' @param eps The stopping criterion parameter. Default is \code{1e-5}.
#' @param init Warm start by truncated PCA algorithm. Default is \code{TRUE}.
#'
#' @return A list containing:
#' \item{rhat}{The estimated factor number.}
#' \item{Fhat}{The estimated factor matrix.}
#' \item{Lhat}{The estimated loading matrix.}
#' \item{loss}{The value of the loss function.}
#'
#' @references Robust factor analysis with exponential squared loss. Jiaqi Hu, Tingyin Wang, Xueqin Wang. Journal of Multivariate Analysis 2026, 213 105567; doi:10.1016/j.jmva.2025.105567
#' @author Jiaqi Hu
#'
#' @examples
#' \donttest{
#' # Assuming gendata() is defined in your package
#' dat = gendata()
#' REFA_FN(dat$Y, rmax = 8)
#' }
#' @export
REFA_FN <- function(Y, rmax = 8, tau = 0.75, q = 0.1, eps = 1e-04, init = TRUE) {
  T = nrow(Y)
  N = ncol(Y)
  LNT = sqrt(min(N, T))
  PNT = LNT^(-2/3) * log(LNT)
  j = 1
  rhat = rmax
  res_list = list()
  while (j <= rmax) {
    res = REFA(Y = Y, r = j, tau = tau, q = q, eps = eps, init = init)
    res_list[[j]] = res
    sigma = diag(t(res$Lhat) %*% res$Lhat/N)
    sigma_j = sigma[j]
    if (sigma_j < PNT) {
      rhat = j - 1
      break
    } else {
      j = j + 1
    }
  }
  if (rhat == 0) {
    res = list(rhat = 0, Fhat = NULL, Lhat = NULL, loss = 1)
    return(res)
  }
  res = res_list[[rhat]]
  res = list(rhat = rhat, Fhat = res$Fhat, Lhat = res$Lhat, loss = res$loss)
  return(res)
}
