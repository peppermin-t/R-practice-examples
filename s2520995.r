eta <- function(r) {
  res <- 0
  if (r > 0) res <- r ** 2 * log(r)
  res
}

eta_vec <- function(x_target, xk) {
  apply(xk, 1, function(p) eta(dist(rbind(x_target, p))))
}

getTPS <- function(x, k=100) {
  if (k <= 3) stop("k too small!")

  # select xk
  n <- nrow(x)
  ik <- 1: n
  if (k < n) ik <- sample(n, k)
  xk <- x[ik, ]  # if k=1, xk is a vector

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
  y <- 1:n
  betah <- solve(U %*% Diag %*% t(U) %*% R, qr.qty(qrx, y)[1: k])  # eigenvalues. concise  # notice the "[1: k]"!

  EDF <- sum(1 / diag(Diag))

  muh <- X %*% betah

  list(V=sum((y - muh) ** 2) / (n - EDF) ** 2, beta=betah, mu=muh, EDF=EDF)
}

fitTPS <- function(x, y, k=100, lsp=c(-5, 5)) {

  # set up the TPS
  tps <- getTPS(x, k)

  k <- nrow(tps$xk) # renew k

  X <- tps$X
  S <- tps$S

  qrx <- qr(X)  # qr.Q, concise

  R <- qr.R(qrx)  # k, k

  R_rev <- backsolve(R, diag(k))

  RSR <- t(R_rev) %*% S %*% R_rev

  # checking symmetricy
  if (max(abs(t(RSR) - RSR)) < 1e-7) {
    RSR <- (RSR + t(RSR)) * 0.5
  } else stop("the matrix is not symmetric!")

  evrsr <- eigen(RSR)  # more concise?

  min_V <- .Machine$integer.max
  best_V_bundle <- NULL
  best_lambda <- 0
  edfs <- vector()

  for (i in seq(lsp[1], lsp[2], length.out=100)) {
    lambda <- exp(i)
    V_bundle <- getV(tps, y, lambda, evrsr, R, qrx)
    edfs <- c(edfs, V_bundle$EDF)

    if (V_bundle$V < min_V) {
      min_V <- V_bundle$V
      best_lambda <- lambda
      best_V_bundle <- V_bundle
    }
  }

  tps_list <- list(
    beta=best_V_bundle$beta, mu=best_V_bundle$mu,
    medf=best_V_bundle$EDF, lambda=best_lambda,
    gcv=min_V, edf=edfs, qrtk=tps$qrtk, xk=tps$xk
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

  y <- cbind(matrix(1, m * m, 1), xp) %*% alpha +
    t(apply(xp, 1, function(p) eta_vec(p, tps$xk))) %*% delta
  contour(x1, x2, matrix(y, m, m))
  persp(x1, x2, matrix(y, m, m), theta=30, phi=30)
}
