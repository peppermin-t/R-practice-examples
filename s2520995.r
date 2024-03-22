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
#   k (size: scalar): number of the control points (default is 100)
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

  # return the needed objects as a list, including xk, X, S, and Z-related
  # infromation qrtk
  list(xk=xk, X=X, S=S, qrtk=qrtk)
}

# Define function "getV"
# Aims to calculate the genralized cross validation (GCV) value, the minimizing
# objective, for a single smoothing parameter, lambda, trial with needed
# matrices provided.
# Parameter:
#   X (size: n, k): the X matrix generated during getTPS in fitting
#   y (size: n,): the output y of the data to be fitted, containing n
#     single values
#   lambda (size: scalar): the lambda value of this GCV calculation trial
#   R (size: k, k): the R matrix of the QR decomposition result of X
#   evrsr: the result information of the eigen-decomposition for (R^-T)S(R^-1)
#   qrx: the result information of the QR-decomposition of X
getV <- function(X, y, lambda, R, evrsr, qrx) {
  # obtain the up-of-date k and n values from length of y and row numbers of R
  n <- length(y)
  k <- nrow(R)

  # retrieve the corresponding eigen vector matrix, U,
  # and the eigen values, ev, from evrsr
  U <- evrsr$vectors  # k, k
  ev <- evrsr$values  # k,

  # generate the diagonal matrix Diag with 1 + lambda * ev_i using diag(),
  # i.e. I + (lambda)*(Lambda)
  Diag <- diag(1 + lambda * ev)  # k, k

  # generate the parameter beta under this lambda, by solving a linear
  # transformation using solve(), with beta unknown. the rhs is obtained by
  # calculating (Q^T)y (Q is the Q matrix of QR decomposition of X) implicitly
  # with qr.qty to improve efficiency, and the lhs computes matrix
  # multiplications, regarding U is an orthogonal matrix.
  betah <- solve(U %*% Diag %*% t(U) %*% R, qr.qty(qrx, y)[1: k])  # k,
  # note that qr.qty() use the transposed complete Q matrix to calculate
  # the dot product with y, so subtracting the first k items of the result
  # will be the expected k, size result vector

  # calculate the effective degrees of freedom (EDF) of the model under the
  # given lambda, by directly summing up the inverse of each elements of Diag's
  # diagonal elements.
  EDF <- sum(1 / diag(Diag))  # scalar

  # calculate the model prediction of y's expectation, mu, by multiplication
  # between X and beta
  muh <- X %*% betah  # n,

  # calculate the generalized cross validation (GCV) value, V, under the given
  # lambda
  V <- sum((y - muh) ** 2) / (n - EDF) ** 2  # scalar

  # return the GCV value and other values needed in this lambda trial as a list,
  # including V, beta, mu, and EDF
  list(V=V, beta=betah, mu=muh, EDF=EDF)
}

# Define function "fitTPS"
# Aims to fit the spline to best minimizing the objective function and get the
# best set of parameters, i.e. beta, of the spline, according to the given x and
# y data, through finding the smoothing parameter, lambda, best minimizing the
# generalized cross validation (GCV) value.
# Parameter:
#   x (size: n, 2): the input x of the data to be fitted, containing n
#     2-dimensional points
#   y (size: n,): the output y of the data to be fitted, containing n
#     single values
#   k (size: scalar): number of the control points (default is 100)
#   lsp (size: 2,): the log scale range on which the smoothing parameter,
#     lambda's search is spaced on, with lsp[1] being the lower bound, lsp[2]
#     being the upper bound (default is c(-5, 5))
fitTPS <- function(x, y, k=100, lsp=c(-5, 5)) {

  # first call getTPS to set up the matrices for the spline and prepare
  # other variables needed for fitting, including X, S, and xk
  tps <- getTPS(x, k)

  # obtain the up-of-date k value from length of xk of tps
  k <- nrow(tps$xk)
  # retrieve the X and S matrix from the returned tps list
  X <- tps$X  # n, k
  S <- tps$S  # k, k

  # do QR decomposition for X with qr(), save the result in qrx
  qrx <- qr(X)

  # retrieve the R matrix of X's QR decomposition with qr.R
  R <- qr.R(qrx)  # k, k

  # with R being a upper triangular matrix, use backsolve to calculate
  # its inverse, faster than normal solve
  Rinv <- backsolve(R, diag(k))  # k, k

  # calculate the following term, RSR
  RSR <- t(Rinv) %*% S %*% Rinv  # k, k

  # do eigen-decomposition on RSR with eigen(), save the result in evrsr
  # note that "symmetric=TRUE" is needed because the machine error can make
  # RSR asymmetric by having tiny errors in the symmetric positions in the
  # matrix (which suppose to be symmetric)
  evrsr <- eigen(RSR, symmetric=TRUE)

  # record the best Vpack (the list returned by getV, representing the results
  # bundle in a single lambda trial), and the best lambda, initialize them with
  # valid initial values (valid for later comparisons)
  vpackb <- list(V=Inf)
  lambdab <- 0

  # initialize an empty vector to record the EDF value of each trial
  edfs <- vector()

  # looping over the possible lambda values, which is 100 values between
  # exp(lsp[1]) and exp(lsp[2])
  for (i in seq(lsp[1], lsp[2], length.out=100)) {
    # get the lambda value this trial
    lambda <- exp(i)
    # calculate the results bundle of this trial, save it in vpack
    vpack <- getV(X, y, lambda, R, evrsr, qrx)
    # append the edf value this trial to the record, edfs
    edfs <- c(edfs, vpack$EDF)

    # if the GCV value this trial is less than the current smallest GCV value,
    # then renew the best lambda and the best results bundle
    if (vpack$V < vpackb$V) {
      lambdab <- lambda
      vpackb <- vpack
    }
  }

  # ready to return the list tps_list, including best values of parameters beta,
  # mu, EDF, lambda, and GCV, along with the edf vector, and other variables
  # needed in plotting, including qrtk, xk, and x's largest and smallest
  # coordinate values
  tps_list <- list(
    beta=vpackb$beta, mu=vpackb$mu, medf=vpackb$EDF,
    lambda=lambdab, gcv=vpackb$V, edf=edfs,
    qrtk=tps$qrtk, xk=tps$xk, x_range=c(min(x), max(x))
  )

  # set the class of the return list as "tps"
  class(tps_list) <- "tps"

  # return tps_list
  tps_list
}

# Define function "plot.tps"
# Aims to plot the fitted thin plate spline with the best parameters
# obtained from fitTPS
# Parameter:
#   tps: the list of objects of class "tps", returned by fitTPS, containing
#     beta, mu, medf, lambda, gcv, edf, qrtk, xk, x_range, as mentioned above
plot.tps <- function(tps) {
  # obtain the up-of-date k value from length of beta of tps
  k <- length(tps$beta)

  # break beta vector into 2, with the first k-3 element being deltaz, and
  # the last 3 elements being alpha (the weight params of the affine
  # transformation)
  alpha <- tps$beta[(k - 2): k]  # 3
  deltaz <- tps$beta[1: (k - 3)]  # k - 3
  # do multiplication beqween Z and deltaz to get the final parameters, delta
  # (the weight params of the non-affine deformation).
  # once again not using the explicit form of Z, instead using qr.qty, since
  # Z is the last k-3 columns of the Q matrix from a QR decomposition of
  # Tk_ to improve efficiency. In this case, deltaz should be extended 3
  # locations of 0 in the front to get the desired result.
  delta <- qr.qy(tps$qrtk, c(rep(0, 3), deltaz))  # k

  # creating testing grid of input x
  # setting the side length of the testing grid, m, to 50
  m <- 50
  # generating the testing grid xp with the side length of m, lower bound of the
  # samallest of x's coordinates, and upperbound of the largest one.
  x2 <- x1 <- seq(tps$x_range[1], tps$x_range[2], length=m)  # m,
  xp <- cbind(rep(x1, m), rep(x2, each=m))  # m * m, 2

  # generate the corresponding prediction y results with vector products:
  #   alpha and delta are mentioned above; the term multipied with alpha
  #   is the model matrix of the affine transformation part, and the term
  #   multipied with delta is the matrix similar to E, containing the structual
  #   information between new xp and former xk
  # note that using cbind is faster than rbind
  y <- cbind(matrix(1, m * m, 1), xp) %*% alpha +
    t(apply(xp, 1, function(p) eta_vec(p, tps$xk))) %*% delta  # m * m,

  # plot the contour of the generated spline according to the y and xp
  contour(x1, x2, matrix(y, m, m), xlab="x1", ylab="x2")
  # add a title to the contour
  title("TPS fitting Result for the Input Set of 2-D x and 1-D y: Contour",
        font = 4)
  # plot the perspective plot of the generated spline according to the y and xp,
  # and adjust the perspective of the plot's demonstration
  persp(x1, x2, matrix(y, m, m), theta=30, phi=30, zlab="y")
  # add a title to the perspective plot
  title("TPS fitting Result for the Input Set of 2-D x and 1-D y: Perspective 
        Plot", font = 4)
}

# final score: 22/23

# Auto-mark summary

# 2/2 code loaded? (0 failed, 1 code outside functions, 2 OK)

# Basic test with n=100 ...

# 1/1 getTPS() ran?
# 1/1 used xk as x unaltered?
# 1/1 correct model matrix?
# 1/1 correct penalty?
# 1/1 EDF vector correct?
# 1/1 optimum EDF correct?
# 1/1 mu and/or beta correct?
# 2/2 efficiency of getTPS (0 >1000*t_opt, 1 >200*t_opt, otherwise 2)
# 2/2 efficiency of fitTPS (0 >100*t_opt, 1 >20*t_opt, otherwise 2)

# Plot tests with n=500...

# 2/2 plots (2 correct; 1 - minor errors; 0 - big errors/missing)?
# 2/2 plots for x not in [0,1]^2 (2 OK; 1 - minor errors; 0 - errors/missing)?

# 20/20 weighted total pre-comment marks and adjustment
# 3/3 Comments is good

# The RSR should not be computed explicitly.
# -1