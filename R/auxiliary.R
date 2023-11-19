
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

update_F = function(ft, y, L, gamma) {
  s = y - L %*% ft
  s = s^2
  -sum(exp(-s/gamma))
}

update_L = function(l, y, F, gamma) {
  s = y - F %*% l
  s = s^2
  -sum(exp(-s/gamma))
}

TR <- function(Fhat, F0) {
  F0 = svd(F0)$u
  F0 = as.matrix(F0)
  Fhat = svd(Fhat)$u
  Fhat = as.matrix(Fhat)
  sum(diag(t(F0) %*% Fhat %*% solve(t(Fhat) %*% Fhat) %*% t(Fhat) %*% F0))/sum(diag(t(F0) %*% F0))
}

ECC <- function(Chat, C){
  norm(Chat - C, type = "F")/norm(C, type = "F")
}

