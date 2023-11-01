set.seed(10)
source("s2520995.r")

ff <- function(p) {
  exp(-(p[, 1] - .3) ^ 2 / .2 ^ 2 - (p[, 2] - .3) ^ 2 / .3 ^ 2) * .5 +
    exp(-(p[, 1] - .7) ^ 2 / .25 ^ 2 - (p[, 2] - .8) ^ 2 / .3 ^ 2)
}

n <- 500
x <- matrix(runif(n * 2), n, 2)
y <- ff(x) + rnorm(n) * .1 ## generate example data to fit

res <- fitTPS(x, y)
plot(res)