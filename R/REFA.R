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
