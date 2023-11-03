###########################################
# Smoothing with Thin Plate Splines (TPS) #
###########################################

# Author: Yinjia Chen (s2520995)

# Overview:
# This R code is a fast implementation of a data fitting method, thin plate
# spline (TPS) (an affine transformation part plus a non-affine deformation
# part), which can create smooth surfaces from the given sparse data ponits.
# The code specifically including setting up the matrices for the spline
# (getTPS), and fitting the thin plate splines to the given data (fitTPS)
# (i.e. obtaining the value of the TPS function's parameters), and then plotting
# the contour of the fitted function (plot.tps).

# Usage:
# This file only contains functions for calling from outside with
# 'source($filename$)' in R terminal or R source files to use.

# Assumption:
# each of the input point x_i has 2 values, corresponding to a single output
# value y_i.

# Define function "eta"
# Implement the TPS's kernel function: radial basis function (RBF)
# Parameter (size: scalar) used here is the input of the RBF
eta <- function(r) {

  # initialize res as 0
  res <- 0

  # if input r is greater than 0, then assign the RBF output for r to res,
  # else remain res as 0
  if (r > 0) res <- r ** 2 * log(r)

  # return the output res (size: scalar)
  res
}

# Define funciton "eta_vec"
# Aims to compute a line of the TPS kernel, i.e. applying RBF to the
# euclidian distances of several points (xk) to a specific point (x)
# Parameter:
#   x (size: 2,): the specific point with whom to compute RBF of
#     euclidian distances
#   xk (size: k, 2): the bunch of points to compute RBF of euclidian
#     distances with the single point x
eta_vec <- function(x, xk) {

  # first using the recycling rule in R to calculate the dufference
  # between x and each point of xk (need to transpose xk in advance
  # to enable the recycling rule to work correctly). (->> size 2, k)
  # then applying sqaures, colSums and a sqare root to calculate the
  # euclidian distances between x and each point of xk. (->> size k,)
  # finally apply "eta" function to each of the distances, and return
  # the obtained vector of RBF values for x and each point of xk. (->> size k,)
  sapply(sqrt(colSums((t(xk) - x) ** 2)), eta)
}

# Define function "getTPS"
# Aims to set up the matrices for the spline and prepare other variables
# needed for fitting, including the control points xk, intermediate matrix X,
# S, and Z-related decomposition object
# Parameter:
#   x (size: n, 2): the input x of the data to be fitted, containing n
#     2-dimensional points
#   k (size: scalar): number of the control points
getTPS <- function(x, k=100) {

  # ensure k is not less than or equal to 3, otherwise stop the programme,
  # since the operations afterwards will be invalid
  if (k <= 3) stop("k is too small!")

  # retrieve the input data points' number as n
  n <- nrow(x)

  # generate the select rows index of choosing control points xk from the
  # given points x, ik (->> size: k,).
  if (k >= n) {

    # if control points' number k is greater than or equal to given points'
    # number n, set k to n
    k <- n

    # and further assign ik directly as 1: n, which means the control points
    # are exactly the original given points.
    ik <- 1: n
  } else {

    # Otherwise, randomly sample k row indices from 1:n as ik corresponding to
    # the selected control points
    ik <- sample(n, k)
  }

  # generate control points matrix xk
  xk <- x[ik, ]  # xk: [k, 2]

  # according to the x and generated xk, calculate the TPS kernel matrix E*,
  # containing the internal structual information of xk; and the matrix E,
  # containing the structual information of x and xk.
  # apply eta_vec to each point inside the x through apply iterating by rows to
  # get E (needing transpose to make it the right size)
  E <- t(apply(x, 1, function(p) eta_vec(p, xk)))  # E: [n, k]
  # select the ik rows of E to get the expected Ek (reduce repeat calculations)
  Ek <- E[ik, ]  # Ek: [k, k]

  # generate the normal and control points' model matrix T_ and Tk_ bind a
  # column vector with all ones with the given points x to get T_
  T_ <- cbind(matrix(1, n, 1), x)  # T_: [n, 3]
  # select the ik rows of T_ to get the expected Tk_ (reduce repeat
  # calculations)
  Tk_ <- T_[ik, ]  # Tk_: [k, 3]

  # generate the matrix Z from QR decomposing matrix Tk_
  # get QR decomposition information of matrix Tk_ using qr()
  qrtk <- qr(Tk_)
  # first get the Q matrix of the QR decompostion of Tk_ using qr.Q(), and then
  # remove the first 1 to 3 columns to get matrix Z
  # note that qr.Q here has to have complete=TRUE to ensure getting the complete
  # square Q matrix (->> size k, k)
  Z <- qr.Q(qrtk, complete=TRUE)[, -(1:3)]  # Z: [k, k - 3]

  # generate matrix X by combining by columns of EZ and T_
  X <- cbind(E %*% Z, T_)  # X: n, k

  # generate matrix S, whose top left k-3, k-3 submatrix is (Z^T)Ek(Z), with all
  # other positions to be 0.
  # first initialize a k, k all-0 matrix
  S <- matrix(0, k, k)
  # assign the top left k-3, k-3 submatrix with (Z^T)Ek(Z) to get final S
  S[1 : (k - 3), 1 : (k - 3)] <- t(Z) %*% Ek %*% Z  # S: [k, k]

  # return the needed objects as a list, including xk, X, S, and qrtk
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
  evrsr <- eigen(RSR, symmetric=TRUE)  # more concise?

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
    qrtk=tps$qrtk, xk=tps$xk, x_range=c(min(x), max(x))
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
  x2 <- x1 <- seq(tps$x_range[1], tps$x_range[2], length=m)
  xp <- cbind(rep(x1, m), rep(x2, each=m))

  # cbind is faster than rbind
  y <- cbind(matrix(1, m * m, 1), xp) %*% alpha +
    t(apply(xp, 1, function(p) eta_vec(p, tps$xk))) %*% delta

  contour(x1, x2, matrix(y, m, m))
  persp(x1, x2, matrix(y, m, m), theta=30, phi=30)
}

# E's apply can be further simplified