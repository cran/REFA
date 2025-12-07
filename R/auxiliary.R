# -------------------------------------------------------------------------
# Internal Helper Functions (Not exported, No Rd file generated)
# -------------------------------------------------------------------------

#' @noRd
normalize <- function(F, L, r) {
  F = as.matrix(F)
  L = as.matrix(L)
  res = svd(F %*% t(L), nu = r, nv = r)
  F1 = sqrt(nrow(F)) * res$u
  F1 = as.matrix(F1)
  L1 = res$v %*% diag(res$d[1:r], r, r)/sqrt(nrow(F))
  L1 = as.matrix(L1)
  list(Fhat = F1, Lhat = L1)
}

#' @noRd
update_F = function(ft, y, L, gamma) {
  s = y - L %*% ft
  s = s^2
  -sum(exp(-s/gamma))
}

#' @noRd
update_L = function(l, y, F, gamma) {
  s = y - F %*% l
  s = s^2
  -sum(exp(-s/gamma))
}

# -------------------------------------------------------------------------
# Exported Functions
# -------------------------------------------------------------------------

#' Estimation of errors for common component
#'
#' @description Estimation of errors for common component
#'
#' @param Chat The estimated common component.
#' @param C0 The true common component.
#'
#' @return A numeric value of the ECC.
#'
#' @references Robust factor analysis with exponential squared loss. Jiaqi Hu, Tingyin Wang, Xueqin Wang. Journal of Multivariate Analysis 2026, 213 105567; doi:10.1016/j.jmva.2025.105567
#' @author Jiaqi Hu
#'
#' @examples
#' \donttest{
#' dat = gendata()
#' Y = dat$Y
#' F0 = dat$F0
#' L0 = dat$L0
#' C0 = F0 %*% t(L0)
#' res = REFA(dat$Y, r = 3)
#' Fhat = res$Fhat
#' Lhat = res$Lhat
#' Chat = Fhat %*% t(Lhat)
#' ECC(Chat, C0)
#' }
#' @export
ECC <- function(Chat, C0){
  norm(Chat - C0, type = "F")/norm(C0, type = "F")
}

#' Trace ratios
#'
#' @description Trace ratios
#'
#' @param Fhat The estimated factors.
#' @param F0 The true factors.
#'
#' @return A numeric value of the trace ratios.
#'
#' @references Robust factor analysis with exponential squared loss. Jiaqi Hu, Tingyin Wang, Xueqin Wang. Journal of Multivariate Analysis 2026, 213 105567; doi:10.1016/j.jmva.2025.105567
#' @author Jiaqi Hu
#'
#' @examples
#' \donttest{
#' dat = gendata()
#' Y = dat$Y
#' F0 = dat$F0
#' res = REFA(dat$Y, r = 3)
#' Fhat = res$Fhat
#' TR(Fhat, F0)
#' }
#' @export
TR <- function(Fhat, F0) {
  F0 = svd(F0)$u
  F0 = as.matrix(F0)
  Fhat = svd(Fhat)$u
  Fhat = as.matrix(Fhat)
  sum(diag(t(F0) %*% Fhat %*% solve(t(Fhat) %*% Fhat) %*% t(Fhat) %*% F0))/sum(diag(t(F0) %*% F0))
}
