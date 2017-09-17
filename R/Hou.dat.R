#' A real scRNA-seq expression data set from Hou et.al (2016)
#'
#' This data contains expression levels (normalized and log-transformed) for 33 cells and 27033 genes.
#'
#' @details This data contains raw expression levels (log-transformed but not centered) for 33 HCC cells and 27033 genes. The 33 cells belongs to three different subpopulations and exhibited different biological characteristics. For descriptions of how we generated this data, please refer to the paper.
#'
#' @docType data
#'
#' @usage data(Hou)
#'
#' @format A object containing the following variables:
#' \describe{
#' \item{\code{x}}{An expression data frame of 33 HCC cells on 27033 genes.}
#' \item{\code{y}}{Numerical group indicator of all cells.}
#' \item{\code{gnames}}{Gene names of all genes.}
#' \item{\code{snames}}{Cell names of all cells.}
#' \item{\code{groups}}{Cell group names.}
#' \item{\code{note}}{A simple note of the data set.}
#' }
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65364}
#' @references Hou, Yu, et al. "Single-cell triple omics sequencing reveals genetic, epigenetic, and transcriptomic heterogeneity in hepatocellular carcinomas." Cell research 26.3 (2016): 304-319.
#'
#' @keywords datasets
#'
#' @examples
#' data(Hou)
#' data <- Hou$x
"Hou"
