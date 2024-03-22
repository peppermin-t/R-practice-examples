################################
# Monotonic Regression in JAGS #
################################

# Author: Yinjia Chen (s2520995)

# Overview:
# This R code aims at finding the relationship between the concentration of
# vinclozolin and the human androgen receptors' (AR) activity by fitting two
# proposed models with monotonic constraints using JAGS, one assumes the same
# set of decreasing sequence for each experimental run, and the other assumes
# differnet. It is based on the given data consisting of effects of five
# experimental runs, each with 9 levels of vinclozolin doses.

# Assumptions:
# 1. The data given is stored in vin.txt, and the two model files are named
#    "M1s2520995.jags" and "M2s2520995.jags", all of which are in the same
#    directory with this R code file.
# 2. The parameters monitored do not include constant values like "mu[1]" and
#    values determinable by considered parameters, like "x", "m", and "ee".

# Loading necessary libraries, including JAGS model interface, rjags, and ggplot
library(rjags)
library(ggplot2)

# Define function "plot_data_model"
# Plot the true datapoints with effect against the dose index (1...9) (in cross
#   symbols), overlayed with the mean expected effect along with the 95%
#   credible regions calculated from samples of the model (in picewise lines),
#   all colored by experimental run.
# Parameter: dt: the dataframe of original data (with column "conc", "effect",
#             "exper")
#            sam: the sampled parameters of niter iterations
#            model_name: the model index of which 'sam' is from
plot_data_model <- function(dt, sam, model_name) {
  # obtain the column index of expected effect (ee) values in the first chain
  # of sam
  ee_iid <- grep("^ee", colnames(sam))

  # retrieve the ee values in the sampling process
  EE_mcmc <- sam[, ee_iid]

  # obtain the mean value of each ee value (w.r.t. each level of each run)
  EE_mean <- colMeans(EE_mcmc)

  # obtain the 95% HPDinterval of each ee value (w.r.t. each level
  # of each run) with HPDinterval (a more robust and intepretable means
  # comparing to quantile)
  EE_interval <- HPDinterval(EE_mcmc)

  # transform the exper column (the group of the data) to a factor, for
  # later usage of coloring
  group <- factor(dt$exper)

  # generate the plot with effects against dose indices
  # first mapped the data in dt, with effects against conc, colored by
  # experimental runs
  data_plot <- ggplot(dt, aes(x=conc, y=effect, color=group)) +
    # add dt points to the plot (cross symbol)
    geom_point(shape=4, size=3) +
    # add ee mean lines
    geom_line(aes(y=EE_mean), linewidth=1.5) +
    # add ee 95% credible interval areas
    geom_ribbon(aes(ymin=EE_interval[, 1], ymax=EE_interval[, 2], fill=group),
                alpha = 0.3) +
    # set the breaks of x axis as the dose level indices
    scale_x_continuous(breaks=unique(dt$conc)) +
    # add x and y labels, the plot title, and the color tag title
    labs(title = paste("Model", model_name, ": Effect vs. Does Index"),
         x="Dose Index", y="Effect",
         color="Experimental Run", fill="Experimental Run") +
    # set minimal style theme
    theme_minimal()

  # plot the ggplot object
  print(data_plot)

}

# load data
dt <- read.table("vin.txt", header=TRUE, sep=" ")

# preprocess data by transforming concentrantion level to level index 1:9
dt$conc <- match(dt$conc, sort(unique(dt$conc)))

# Obtain the run times (J=5) and the dose level counts (N=9)
J <- length(unique(dt$exper))
N <- length(dt$effect) / J

# Obtain the effect data reshaped in a matrix with rows representing each level
# and colums representing each run
y <- matrix(dt$effect, N, J)

# record the stochastic parameters from the samples, dicarding constants and
# parameters determined by considered parameters; record ee (M * mu) for later
# plotting to avoid repeat computation
vars <- c("M", "mu", "tau", "tau0", "ee")

# setting sampling parameters, including:
# total iteration counts, burn-in steps, and chain counts.
# The values are set after checking of parameters to ensure cocnvergence
niter <- 10000
nburnin <- 2000
nchain <- 2

# create model 1, with observed data y
mod1 <- jags.model("M1s2520995.jags",
                   data=list(y=y, N=N, J=J), n.chains=nchain)

# burn-in phase
update(mod1, n.iter=nburnin)

# sample phase for model 1, rercording interested parameters
sam1 <- coda.samples(mod1, variable.names=vars, n.iter=niter)

# when considering stochastic parameters, ee and constant mu[1]
# from the list need to be removed, find the indices
discard_vars <- grep("mu\\[1\\]|^ee", colnames(sam1[[1]]))
# obtain the effective sizes for remaining parameters based on their samples
es <- effectiveSize(sam1[, -discard_vars])

# obtain the index of the parameter with the highest effective sample size,
# generally tau0, occasionally M[1]
max_es_par <- names(which(es == max(es)))
# obtain the index of the parameter with the lowest effective sample size,
# mu[9]
min_es_par <- names(which(es == min(es)))

# plot the traceplot of the highest effective size parameter, and entitle it
traceplot(sam1[[1]][, max_es_par],
          main=paste("Traceplot for the Sampled Param. with the\n
            Highest Effective Sample Size:", max_es_par))
# plot the traceplot of the lowest effective size parameter, and entitle it
traceplot(sam1[[1]][, min_es_par],
          main=paste("Traceplot for the Sampled Param. with the\n
            Smallest Effective Sample Size:", min_es_par))

# call plot_data_model and plot the target plot with data and
# expected effect means and quantiles for model 1
plot_data_model(dt, sam1[[1]], model_name=1)

# create model 2, with observed data y
mod2 <- jags.model("M2s2520995.jags",
                   data=list(y=y, N=N, J=J), n.chains=nchain)

# burn-in phase
update(mod2, n.iter=nburnin)

# sample phase for model 2, rercording interested parameters
sam2 <- coda.samples(mod2, variable.names=vars, n.iter=niter)

# call plot_data_model and plot the target plot with data and
# expected effect means and quantiles for model 2
plot_data_model(dt, sam2[[1]], model_name=2)

# compare the two models with deviance information criterion (DIC).
# DIC is the sum of two parts: the negative log likelihood term and the
# effective degrees of freedom (EDF), also the "penalty" term, the lower
# of which indicates a better model, or closer to the true model
dic.samples(mod1, n.iter=niter)
dic.samples(mod2, n.iter=niter)

# Conclusion: the penalized deviance of the two models are different between
# different runs, but are basiclly around:
#   model1: 590.6 (581.9 + 8.709)
#   model2: 576.9 (560 + 16.92)
# (with error smaller than 1), meaning that generally the second model which
# assumes different decresing sequences for each expriemntal run is better.
