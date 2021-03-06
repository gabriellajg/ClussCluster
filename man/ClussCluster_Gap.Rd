% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClussCluster.R
\name{ClussCluster_Gap}
\alias{ClussCluster_Gap}
\title{Function to select optimal tuning parameter based on Gap statistic}
\usage{
ClussCluster_Gap(x, nclust = NULL, B = 20, centers = NULL, ws = NULL,
  nepoch.max = 10, theta = NULL, seed = 1, nstart = 20, iter.max = 50,
  verbose = FALSE, progress = TRUE)
}
\description{
The tuning parameter controls the L1 bound on w, the feature weights. A permutation approach is used to select the tuning parameter.
}
