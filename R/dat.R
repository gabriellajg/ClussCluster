#' A simulated expression data set.
#'
#' This example data set containing expressing levels for 60 cells and 200 genes. The 60 cells belong to 4 cell types with 15 cells each. Each cell type is uniquely associated with 30 signature genes, i.e., the first cell type is associated with the first 30 genes, the second cell type is associated with the next 30 genes, so on and so forth. The remaining 80 genes show indistinct expression patterns among the four cell types and are considered as noise genes.
#'
#' @docType data
#'
#' @usage data(sim_dat)
#'
#' @format A data frame with 60 cells on 200 genes.
#'
#' @return A simulated dataset used to demonstrate the application of \code{ClussCluster}.
#'
#' @keywords datasets
#'
#' @examples
#' data(sim_dat)
#' head(sim_dat)
"sim_dat"
