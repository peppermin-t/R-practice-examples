# Yinjia Chen (S2520995), Yiwen Xing (S2530703), Kai Wen Lian (S2593019)

# https://github.com/peppermin-t/SGD.git

# Joint contribution: Go through the initial plan of the code
# Yinjia Chen: Main developer and validator of codes
# Kai Wen Lian: Sub developer of codes, comment codes
# Yiwen Xing: Sub developer of codes, verify codes, add additional comment.

##############################
##############################

# Brief Description:

# This R code contains functions to create and train a neural network for
# classification using stochastic gradient descent.

# 1. 'netup' function: Initializes a neural network with random weights and
# biases. It takes the number of nodes in each layer as input and returns a
# list of initialized parameters and empty nodes.

# 2. 'forward' function: Implements the forward propagation in the neural
# network. It updates the node values in each layer based on input data for the
# first layer and returns the updated network.

# 3. 'backward' function: Performs backward propagation. It computes the
# derivatives of the loss function with respect to the weights, biases, and
# node values, facilitating the learning process.

# 4. 'train' function: Trains the neural network using forward and backward 
# propagation. It iteratively updates the network's parameters based on the 
# gradient of the loss function calculated from a subset of the training data.

# 5. 'test' function: Tests the trained network on a separate dataset to
# evaluate its performance, specifically the misclassification rate.

# The script also includes a practical implementation where it trains a neural
# network to classify iris flowers into species based on four features. The
# dataset is split into training and testing sets, with the network's
# performance evaluated by its misclassification rate on the test set.

##############################
##############################


# netup is a function to build a list representing the network
# Input: d (a vector giving the number of nodes in each layer of a network)
# Output: An initialised network list containing a list for each h (empty nodes 
# for each layer), w (weight matrices), b (offset vectors)

netup <- function(d) {

  # Number of layers
  n <- length(d)

  # Initialise a list to store the weight matrices (W) and offset vectors (b)
  h <- W <- b <- list()

  # Loop over the layers from the 1st to (n-1)th
  for (i in 1: (n - 1)) {
    # Create a list of nodes for each layer (h)
    h[[i]] <- rep(0, d[i])

    # Dimensions for each transformation matrix and offset which is
    # given by the number of nodes in current layer (d[i+1])
    # and that in the previous layer (d[i])

    # Then the weight matrix for each link is given by:
    W[[i]] <- matrix(runif(d[i] * d[i + 1], 0, 0.2), d[i + 1], d[i])

    # And the offset vectors for each link is given by:
    b[[i]] <- runif(d[i + 1], 0, 0.2)
  }

  # Create the list of nodes for the last layer (h_L)
  h[[n]] <- rep(0, d[n])

  list(h = h, W = W, b = b)
}

# forward is a function to calculate node values in each remaining layers
# Input: nn (a list of network returned from netup()) and
# inp (a vector of input values for the first layer)
# Output: An updated network list, with node values for h

forward <- function(nn, inp) {

  # Obtain the number of layers in nn
  l <- length(nn$h)

  # Store values for the first layer
  nn$h[[1]] <- inp

  # Loop over each remaining layer
  for (i in 1: (l - 1)) {
    # Compute transformation to obtain node values and
    #   set equal to zero if element is negative (ReLU)
    nn$h[[i + 1]] <- pmax(drop(nn$W[[i]] %*% nn$h[[i]] + nn$b[[i]]), 0)
  }

  nn
}


# backward is a function to compute the derivatives of the loss for a network
# Input: nn (a list of network returned from forward()) and
#   k (the real output label k for this run)
# Output: An updated network list, containing additional lists for each
# derivative of loss

backward <- function(nn, k) {

  # Create empty lists for derivatives
  dh <- dW <- db <- list()

  l <- length(nn$h)

  # Compute the derivative of loss w.r.t. the node values (dh) at the last layer
  dh[[l]] <- exp(nn$h[[l]]) / sum(exp(nn$h[[l]]))

  # When index of node is equal to k at each iteration, minus one from output
  dh[[l]][k] <- dh[[l]][k] - 1

  ## Apply chain rule to compute derivatives of loss
  # Loop over each layer from the last layer (using back-propagation)
  for (i in (l - 1): 1) {

    # Store the dh for (i + 1)th layer in d
    d <- dh[[i + 1]]

    # Set to zero when node value is negative
    d[which(nn$h[[i + 1]] <= 0)] <- 0

    # Compute derivative w.r.t. the nodes, dh
    dh[[i]] <- t(nn$W[[i]]) %*% d

    # Compute derivative w.r.t. the offsets, db
    db[[i]] <- d

    # Compute derivative w.r.t. the weights, dW
    dW[[i]] <- d %*% t(nn$h[[i]])
  }

  # Add derivatives to the network list
  c(nn, list(dh = dh, dW = dW, db = db))
}


# train is a function to train the network given input data and labels
# Input: inp (a matrix of input data), k (a vector of labels),
# eta (learning rate), mb (batch size, randomly sampled data number to
# compute gradient for one step optimization),
# nstep (number of optimization steps to take)
# Output: A trained network list

train <- function(nn, inp, k, eta = .01, mb = 10, nstep = 10000) {

  # Obtain the number of train data
  n <- length(k)

  # Obtain the number of layers for the network
  l <- length(nn$h)

  # Loop over each step
  for (i in 1: nstep) {

    # Determine the index of data used for gradient calculation
    iid <- sample(n, mb)

    # Extract randomly sampled data
    sub_inp <- inp[iid, ]

    # Extract its corresponding labels
    sub_k <- k[iid]

    # initialize the summing-up storage of weights and offsets
    grads_W <- grads_b <- list()
    d <- sapply(nn$h, length)
    for (i in 1: (l - 1)) {
      # The weight matrix for each link is given by:
      grads_W[[i]] <- matrix(rep(0, d[i] * d[i + 1]), d[i + 1], d[i])

      # And the offset vectors for each link is given by:
      grads_b[[i]] <- rep(0, d[i + 1])
    }

    # Loop over each sampled data
    for (i in 1: mb) {

      # Forward the input data into the network and get the nn with h updated
      nn <- forward(nn, sub_inp[i, ])

      # Compute gradients with gradient descent
      grad <- backward(nn, sub_k[i])

      # Loop over each layer until the second last layer
      for (i in 1: (l - 1)) {

        # Update weight parameters with the the gradient of weights this run
        grads_W[[i]] <- grads_W[[i]] + grad$dW[[i]]

        # Update offset parameters with the the gradient of offsets this run
        grads_b[[i]] <- grads_b[[i]] + grad$db[[i]]
      }
    }

    # Loop over each layer until the second last layer
    for (i in 1: (l - 1)) {

      # Update weight parameters (W) with corresponding means
      nn$W[[i]] <- nn$W[[i]] - eta * grads_W[[i]] / mb

      # Update offset parameters (b) with corresponding means
      nn$b[[i]] <- nn$b[[i]] - eta * grads_b[[i]] / mb
    }
  }

  nn
}


# test is a function to train a network to classify items to groups based on a
# given output class k
# Input: nn (a list of network returned from forward()),
# inp (a matrix of input data), k (a vector of labels).
# Output: Rate of misclassification

test <- function(nn, inp, k) {

  # Compute number of test data
  n <- length(k)

  # Compute number of layers
  l <- length(nn$h)

  # Initialise the number of misclassified data
  mis_class <- 0

  # Loop over the test data
  for (i in 1: n) {

    # Forward the input data into the network and get the output layer values
    output_h <- forward(nn, inp[i, ])$h[[l]]

    # Find the output class with the largest output values, which is also
    # the one with the largest score (largest probability) after normalizing
    k_ <- which(output_h == max(output_h))

    # Count the misclassified data; occur when the predicted class is not
    # the real one
    mis_class <- mis_class + (k_ != k[i])
  }

  # Compute the misclassification rate
  mis_class / n
}

# set.seed(13) # training has worked
seed <- commandArgs(trailingOnly = TRUE)[1]
set.seed(seed)

# Aim: Train a network to classify irises to species based on given features
# Define the dimensions of the neural network
d <- c(4, 8, 7, 3)

## Data preparation

# Extract the features of iris to be considered
features <- colnames(iris)[1: 4]

# Extract the label based on their species
label <- colnames(iris)[5]

# Extract unique label types in the dataset
label_options <- unique(iris[, 5])

# Extract input data of the network (a matrix)
inp <- as.matrix(iris[, features])

# Extract corresponding output labels k (a integer vector)
k <- as.integer(factor(iris[, label], level = label_options))

## Divide data into test data and training data
# Test data contains every 5th row of iris dataset
test_iid <- seq(5, nrow(iris), by = 5)

# Determine which is used for training data
inp_train <- inp[-test_iid, ]
# Extract the corresponding labels for training data
k_train <- k[-test_iid]

# Determine which is used for test data
inp_test <- inp[test_iid, ]
# Extract the corresponding labels for test data
k_test <- k[test_iid]


## Create and train a network

# Create a network based on the dimensions defined
nn <- netup(d)

# Train the network nn
system.time(nn <- train(nn, inp_train, k_train, eta=.01, mb=10, nstep=10000))

# Classify test data to species according to the class predicted
# and compute misclassification rate
mis_class_rate <- test(nn, inp_test, k_test)
print(mis_class_rate)

# Grade (16/18)

# Instructor's note:
# Nicely commented, easy to follow code. Does what it is supposed to do.

# The only deficiency was that there is quite a lot of unnecessary copying and
#   list creation, which I think is probably responsible for the slightly slow
#   speed. e.g. in 'forward' you could just modify 'nn' directly and return the
#   modified list as the return object.
