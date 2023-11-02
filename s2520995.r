###########################################
# Smoothing with Thin Plate Splines (TPS) #
###########################################

# Author: Yinjia Chen (s2520995)

# This R code

eta <- function(r) {
  res <- 0
  if (r > 0) res <- r ** 2 * log(r)
  res
}

eta_vec <- function(x, xk) {
  sapply(sqrt(colSums((t(xk) - x) ** 2)), eta)
  # apply(xk, 1, function(p) eta(dist(rbind(x_target, p))))
}

getTPS <- function(x, k=100) {
  if (k <= 3) stop("k too small!")

  # select xk
  n <- nrow(x)
  if (k >= n) {
    k <- n
    ik <- 1: n
  } else {
    ik <- sample(n, k)
  }  # concise?
  xk <- x[ik, ]

  # generate matrix E
  E <- t(apply(x, 1, function(p) eta_vec(p, xk)))  # n * k; more concise ways?
  Ek <- E[ik, ]  # k * k; subtracted from E

  # generate matrix T
  T_ <- cbind(matrix(1, n, 1), x)  # n, 3
  Tk_ <- T_[ik, ]  # k, 3

  # generate matrix Z
  qrtk <- qr(Tk_)
  Z <- qr.Q(qrtk, complete=TRUE)[, -(1:3)]  # k, k - 3

  # generate matrix X
  X <- cbind(E %*% Z, T_)  # n, k

  # generate matrix S
  S <- matrix(0, k, k)
  S[1 : (k - 3), 1 : (k - 3)] <- t(Z) %*% Ek %*% Z  # k, k

  list(xk=xk, X=X, S=S, qrtk=qrtk)
}

getV <- function(tps, y, lambda, evrsr, R, qrx) {
  n <- length(y)
  k <- nrow(R)

  X <- tps$X

  U <- evrsr$vectors
  ev <- evrsr$values

  Diag <- diag(1 + lambda * ev)

  # eigenvalues. concise  # notice the "[1: k]"!
  betah <- solve(U %*% Diag %*% t(U) %*% R, qr.qty(qrx, y)[1: k])

  EDF <- sum(1 / diag(Diag))

  muh <- X %*% betah

  list(V=sum((y - muh) ** 2) / (n - EDF) ** 2, beta=betah, mu=muh, EDF=EDF)
}

fitTPS <- function(x, y, k=100, lsp=c(-5, 5)) {

  # set up the TPS and retrieve return values
  tps <- getTPS(x, k)
  k <- nrow(tps$xk) # renew k
  X <- tps$X
  S <- tps$S

  qrx <- qr(X)  # qr of X

  R <- qr.R(qrx)  # k, k

  Rinv <- backsolve(R, diag(k))  # k, k

  RSR <- t(Rinv) %*% S %*% Rinv

  # checking symmetricy
  if (max(abs(t(RSR) - RSR)) < 1e-7) RSR <- (RSR + t(RSR)) * 0.5
  else stop("the matrix is not symmetric!")

  evrsr <- eigen(RSR)  # more concise?

  vpackb <- list(V=Inf)
  lambdab <- 0
  edfs <- vector()

  for (i in seq(lsp[1], lsp[2], length.out=100)) {
    lambda <- exp(i)
    vpack <- getV(tps, y, lambda, evrsr, R, qrx)
    edfs <- c(edfs, vpack$EDF)

    if (vpack$V < vpackb$V) {
      lambdab <- lambda
      vpackb <- vpack
    }
  }

  tps_list <- list(
    beta=vpackb$beta, mu=vpackb$mu, medf=vpackb$EDF,
    lambda=lambdab, gcv=vpackb$V, edf=edfs,
    qrtk=tps$qrtk, xk=tps$xk
  )
  class(tps_list) <- "tps"
  tps_list
}

plot.tps <- function(tps) {
  k <- length(tps$beta)

  # obtain the parameters
  alpha <- tps$beta[(k - 2): k]  # 3
  deltaz <- tps$beta[1: (k - 3)]  # k - 3
  delta <- qr.qy(tps$qrtk, c(rep(0, 3), deltaz))  # k

  # visualize test ...
  m <- 50
  x2 <- x1 <- seq(0, 1, length=m)
  xp <- cbind(rep(x1, m), rep(x2, each=m))

  # cbind is faster than rbind
  y <- cbind(matrix(1, m * m, 1), xp) %*% alpha +
    t(apply(xp, 1, function(p) eta_vec(p, tps$xk))) %*% delta

  contour(x1, x2, matrix(y, m, m))
  persp(x1, x2, matrix(y, m, m), theta=30, phi=30)
}

# add test functions?
# E's apply can be further simplified
# plot's "y <-" should be clearer