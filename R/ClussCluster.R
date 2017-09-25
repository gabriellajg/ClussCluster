#' Function to perform simultaneous detection of cell types and cell-type-specific signature genes
#'
#' \code{ClussCluster} takes the single-cell transcriptome data and returns an object containing cell types and type-specific signature gene sets
#'
#'
#'Takes the normalized and log transformed number of reads mapped to genes (e.g., log(RPKM+1) or log(TPM+1) where RPKM stands for Reads Per Kilobase of transcript per Million mapped reads and TPM stands for transcripts per million) but NOT centered.
#'
#' @param x An nxp data matrix. There are n cells and p genes.
#' @param nclust Number of clusters desired if the cluster centers are not provided. If both are provided, nclust must equal the number of cluster \code{centers}.
#' @param centers A set of initial (distinct) cluster centres if the number of clusters (\code{nclust}) is null. If both are provided, the number of cluster centres must equal \code{nclust}.
#' @param ws One or multiple candidate tuning parameters to be evaluated and compared. Determines the sparsity of the selected genes. Should be greater than 1.
#' @param nepoch.max The maximum number of epochs. In one epoch, each cell will be evaluated to determine if its label needs to be updated.
#' @param theta Optional argument. If provided, \code{theta} are used as the initial cluster labels of the ClussCluster algorithm; if not, K-means is performed to produce starting cluster labels.
#' @param seed This seed is used wherever K-means is used.
#' @param nstart Argument passed to \code{kmeans}. It is the number of random sets used in \code{kmeans}.
#' @param iter.max Argument passed to \code{kmeans}. The maximum number of iterations allowed.
#' @param verbose Print the updates inside every epoch? If TRUE, the updates of cluster label and the value of objective function will be printed out.
#' @param progress Print the progress? If multiple tuning parameters are evaluated, then each tuning parameter will be printed when \code{progress=TRUE}.
#' @examples
#' \dontrun{
#' data(Hou)
#' hou.dat <-Hou$x
#' run.ft <- filter_gene(hou.dat)
#' hou.test <- ClussCluster(run.ft$dat.ft, nclust=3, ws=4, verbose = F)
#' }
#' @export
#'

ClussCluster <- function(x, nclust = NULL, centers = NULL, ws = NULL, nepoch.max = 10, theta = NULL, seed = 1, nstart = 20, iter.max = 50, verbose = FALSE, progress = TRUE)
{
  # check the coexistence of arguments: nclust and centers
  if(is.null(nclust) && is.null(centers)) stop("Must provide either K or centers.")
  if(!is.null(nclust) && !is.null(centers)){
    if(nrow(centers)!=nclust) stop("If K and centers both are provided, then nrow(centers) must equal K!!!")
    if(nrow(centers)==nclust) nclust <- NULL
  }
  if(!is.null(centers) && nrow(centers)!=nrow(x)) stop("If centers are provided, then centers must have equal dimension as the original data.")

  if(!is.null(theta) && length(theta)!=ncol(x)) stop("If initial cluster labels are provided, then the length must equal the number of observations.")
  if (is.null(theta)) {
    # use regular K-means as start
    set.seed(seed)
    if (!is.null(centers)) theta <- kmeans(t(x), centers = t(centers), nstart = nstart, iter.max = iter.max)$cluster
    if (is.null(centers)) theta <- kmeans(t(x), centers = nclust, nstart = nstart, iter.max = iter.max)$cluster
  }

  if(is.null(ws)) ws <- exp(seq(log(1.3), log(sqrt(nrow(x))), length.out=20))

  output <- list()
  for (m in 1:length(ws)){
    s <- ws[m]
    if(progress){cat("s =", s, "\n")}
    sc <- list(x=x, nclust=nclust, centers=centers, s=s, theta=theta, seed=seed, nstart=nstart, iter.max=iter.max, n=ncol(x), p=nrow(x))
    sc <- initial(sc)
    if (verbose){
      cat("initial target value:", sc$wbcss, fill=TRUE)
    }
    flag <- 0
    if.conv <- FALSE
    for (epoch in 1 : nepoch.max)
    {
      if (verbose){
        cat("\nStart epoch", epoch, "...\n")
      }
      for (i in 1 : sc$n)
      {
        if (verbose){
          cat("\tTrying sample", i, "\n")
          cat(i, ',')
        }
        up <- update.theta(sc, i)
        y.new <- up$y.new
        if (y.new == sc$theta[i])
        {
          flag <- flag + 1
          if (flag >= sc$n)
          {
            if.conv <- TRUE
            break
          }
        } else
        {
          if (verbose){
            cat("target value before updating w:", sc$wbcss, fill=TRUE)
          }
          sc.new <- up$sc.update[[y.new]]
          sc <- update.w(sc.new, i, sc$theta[i], y.new)
          if (verbose){
            cat("target value after updating w:", sc$wbcss, fill=TRUE)
            cat("theta:", sc$theta, fill=TRUE)
          }
          flag <- 0
        }
      }
      if (flag >= sc$n)
      {
        if.conv <- TRUE
        break
      }
    }
    output[[m]] <- sc
  }
  if(length(ws)==1){output = output[[1]]}
  return(output)
}

#' @export
print_ClussCluster <- function(object, ...){
  if(is.null(object$s)){
    for (i in 1:length(object)){
      cat("Tuning parameter is: ", object[[i]]$s, fill=TRUE)
      cat("Number of signature genes of each cluster: ", apply(object[[i]]$w, 2, function(x) sum(x!=0)), fill=TRUE)
      cat("Clustering: ", object[[i]]$theta, fill=TRUE)
      cat("Value of objective function: ", object[[i]]$wbcss, fill=TRUE)
      cat('\n')
    }
  } else{
    cat("Tuning parameter is: ", object$s, fill=TRUE)
    cat("Number of signature genes of each cluster: ", apply(object$w, 2, function(x) sum(x!=0)), fill=TRUE)
    cat("Clustering: ", object$theta, fill=TRUE)
    cat("Value of objective function: ", object$wbcss, fill=TRUE)
    cat('\n')
  }

}

#' Plot the results of \code{ClussCluster}
#'
#' If multiple tuning parameters are evaluated in the object, the number of signature genes is plotted against the tuning parameters. If only one is included, then \code{plot_ClussCluster} returns a venn diagram and a heatmap at this particular tuning parameter.
#'
#'
#' Takes the normalized and log transformed number of reads mapped to genes (e.g., log(RPKM+1) or log(TPM+1) where RPKM stands for Reads Per Kilobase of transcript per Million mapped reads and TPM stands for transcripts per million) but NOT centered.
#'
#' If multiple tuning parameters are evaluated in the object, the number of signature genes is computed for each cluster and is plotted against the tuning parameters. Each color and line type corresponds to a cell type.
#'
#' If only one tuning parameter is evaluated, two plots will be produced. One is the venn diagram of the cell-type-specific genes, the other is the heatmap of the data with the cells and top m signature genes. See more details in the paper.
#'
#' @param object An object that is obtained by applying the ClussCluster function to the data set.
#' @param m The number of top signature genes selected to produce the heatmap.
#' @param snames The names of the cells.
#' @param gnames The names of the genes
#' @importFrom grid grid.newpage
#' @importFrom grid grid.draw
#' @importFrom VennDiagram venn.diagram
#' @importFrom scales rescale
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @export
#' @examples
#' data(Hou)
#' run.cc <- ClussCluster(Hou$x, nclust = 3, ws = c(2.4, 5, 8.8))
#' plot_ClussCluster(run.cc, m = 5, snames=Hou$snames, gnames=Hou$gnames)
#'

plot_ClussCluster <- function(object, m = 10, snames=NULL, gnames=NULL, ...){
  if (is.null(object$s)){
    ws <- sapply(object, function(x) x$s)
    num.sig <- sapply(object, function(x) apply(x$w, 2, function(x) sum(x!=0)))
    par(oma=c(0, 0, 0, 3))
    matplot(x = ws, y = t(num.sig), type = 'b', pch = 1:length(object),
            xlab = 'Tuning parameter',
            ylab = 'Number of signature genes of each group')
    legend(par('usr')[2], par('usr')[4], bty='n', xpd=NA, title="Group", legend = 1:nrow(num.sig), pch = 1:nrow(num.sig), col = 1:nrow(num.sig), lty = 1:nrow(num.sig), seg.len = 1)
  } else {
    ind.sig <- apply(object$w, 2, function(x) which(x!=0))
    #par(oma = c(0, 0, 0, 0), mfrow = c(1, 2))
    grid::grid.newpage()
    venn.plot <- VennDiagram::venn.diagram(ind.sig, fill = 1:length(ind.sig), category.names = paste("Cluster", 1:length(ind.sig)), lty =1, filename=NULL)
    grid::grid.draw(venn.plot)
    top.m.hm(object, m, snames, gnames)
  }
}


#' @export
top.m.hm <- function(object, m, snames=NULL, gnames=NULL){
  nz.gene <- which(rowSums(object$w)!=0)
  max.wt <- apply(object$w, 1, which.max)
  sig.gene <- data.frame(nz.gene, weights = object$w[cbind(nz.gene,max.wt[nz.gene])], cluster=max.wt[nz.gene])
  K <- length(unique(object$theta))
  x <- object$x

  g.i <- order.g.i <- u.k <- list()
  over.gi <- list()
  for (k in 1:K){
    u.k[[k]] <- which(object$theta==k)
    gene.k <- sig.gene[sig.gene$cluster==k,]
    g.i[[k]] <- g.k <- gene.k[order(gene.k$weights,decreasing = T),1]
  }
  for (k in 1:K){
    x.sub <- x[g.i[[k]], ]
    for (kk in 1:K){
      over.gi[[kk]] = order(apply(x.sub[,u.k[[kk]]], 1, mean)-apply(x.sub[,-u.k[[kk]]], 1, mean), decreasing = T)[1:m]
    }
    x.k <- x.sub[over.gi[[k]], u.k[[k]]]
    order.g.i[[k]] <- g.i[[k]][ over.gi[[k]] [hclust(dist(x.k))$order] ]
  }
  x.cet <- x[rev(unlist(order.g.i)), unlist(u.k)]
  coln <- snames[unlist(u.k)]
  rown <- gnames[rev(unlist(order.g.i))]
  #coln <- 1:ncol(x.cet); rown <- 1:nrow(x.cet);
  if(!is.null(snames)){coln <- snames[unlist(u.k)]}
  if(!is.null(gnames)){rown <- gnames[rev(unlist(order.g.i))]}
  data.m <- apply(x.cet,1,scales::rescale)
  data.m <- t(data.m)
  rownames(data.m) <- 1:nrow(data.m)
  colnames(data.m) <- 1:ncol(data.m)
  data.m <- reshape2::melt(data.m)
  base_size <- 12

  HM <- ggplot2::ggplot(data.m, ggplot2::aes(Var2, Var1)) +
    ggplot2::geom_tile(ggplot2::aes(fill = value), colour = "white") +
    ggplot2::scale_fill_gradient(low = "white", high = "red") +
    ggplot2::labs(x = "", y = "") +
    ggplot2::scale_y_continuous(expand = c(0, 0), labels=rown, breaks=1:length(rown)) +
    ggplot2::scale_x_continuous(expand = c(0, 0), labels=coln, breaks=1:length(coln)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position ="none",
                   axis.ticks = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 0, colour = "grey20"),
                   axis.text.y = ggplot2::element_text(hjust = 1, colour="grey20"))
  return(HM)
}



#' Function to select optimal tuning parameter based on Gap statistic
#'
#' The tuning parameter controls the L1 bound on w, the feature weights. A permutation approach is used to select the tuning parameter.
#' @export
#'
ClussCluster_Gap <- function(x, nclust = NULL, B = 20, centers = NULL, ws = NULL, nepoch.max = 10, theta = NULL, seed = 1, nstart = 20, iter.max = 50, verbose = FALSE, progress = TRUE)
{
  # check the coexistence of arguments: nclust and centers
  if(is.null(nclust) && is.null(centers)) stop("Must provide either K or centers.")
  if(!is.null(nclust) && !is.null(centers)){
    if(nrow(centers)!=nclust) stop("If K and centers both are provided, then nrow(centers) must equal K!!!")
    if(nrow(centers)==nclust) nclust <- NULL
  }
  if(!is.null(centers) && nrow(centers)!=nrow(x)) stop("If centers are provided, then centers must have equal dimension as the original data.")

  if(!is.null(theta) && length(theta)!=ncol(x)) stop("If initial cluster labels are provided, then the length must equal the number of observations.")
  if (is.null(theta)) {
    # use regular K-means as start
    set.seed(seed)
    if (!is.null(centers)) theta <- kmeans(t(x), centers = t(centers), nstart = nstart, iter.max = iter.max)$cluster
    if (is.null(centers)) theta <- kmeans(t(x), centers = nclust, nstart = nstart, iter.max = iter.max)$cluster
  }

  if(is.null(ws)) ws <- exp(seq(log(1.3), log(sqrt(nrow(x))), length.out=20))

  X_b <- run_b <- list()
  for (b in 1:B)
  {
    set.seed(7*b)
    X_b[[b]] <- t(sapply(1:nrow(x), function(t) sample(x[t,])))
  }

  run <- ClussCluster(x, nclust = nclust, centers = centers, ws = ws, nepoch.max = nepoch.max, theta=theta, seed = seed, nstart = nstart, iter.max = iter.max, verbose = verbose, progress = progress)

  for (b in 1:B){
    cat("permutation sample ", b, "\n")
    run_b[[b]] <- ClussCluster(X_b[[b]], nclust = nclust, centers = centers, ws = ws, nepoch.max = nepoch.max, theta = theta, seed = seed, nstart = nstart, iter.max = iter.max, verbose = verbose, progress = progress)
  }

  O <- rep(NA, length(ws))
  O_b <- matrix(NA, length(ws), B)
  for (t in 1:length(ws)){
    O[t] <- run[[t]]$wbcss
    for (b in 1:B){
      O_b[t,b] <- run_b[[b]][[t]]$wbcss
    }
  }

  Gap <- log(O) - apply(log(O_b), 1, mean)
  sd.Gap <- apply(log(O_b), 1, sd)
  bestw <- ws[which.max(Gap)]
  onesd.bestw <- ws[which.min(!(Gap-sd.Gap)[which.max(Gap)]<=Gap & Gap<=(Gap+sd.Gap)[which.max(Gap)])]

  return(list(ws=ws, O=O, O_b=O_b, Gap=Gap, sd.Gap=sd.Gap, run=run, bestw=bestw, onesd.bestw=onesd.bestw))
}

#' @export
print_ClussCluster_Gap <- function(object, ...){
  cat("Tuning parameter selection results for ClussCluster:", fill=TRUE)
  K <- object$run[[1]]$nclust
  nz <- round(sapply(1:length(object$ws), function(j) sum(object$run[[j]]$w!=0))/K, 2)
  mat <- round(cbind(object$ws, K, nz, object$Gap, object$sd.Gap),4)
  dimnames(mat) <- list(1:length(object$ws), c("Wbound", "# Clusters", "# Mean Non-Zero W's", "Gap Statistic", "Standard Deviation"))
  print(mat, quote=FALSE)
  cat("Tuning parameter that leads to largest Gap statistic: ", object$bestw, fill=TRUE)
  cat("Tuning parameter within one sd of the largest Gap statistic: ", object$onesd.bestw, fill=TRUE)
}




#' Plot the results of \code{ClussCluster.Gap}
#'
#' Plot the gap statistics and number of genes selected as the tuning parameter varies.
#' @param object object obtained from \code{ClussCluster_Gap}()
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @export
#'
plot_ClussCluster_Gap <- function(object){
  df <- data.frame(object$ws, object$Gap)
  limits <- ggplot2::aes(ymax = object$Gap + object$sd.Gap, ymin=object$Gap - object$sd.Gap)
  rg <- range(object$Gap - object$sd.Gap, object$Gap + object$sd.Gap)
  p1 <- ggplot2::ggplot(df, ggplot2::aes(colour="#D55E00", y=object$Gap, x=object$ws)) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::geom_errorbar(limits, width=0.8, size=0.9) +
    ggplot2::scale_y_continuous(limits = rg) +
    ggplot2::theme_bw() +
    ggplot2::labs(y="Gap",x="") +
    ggplot2::ggtitle("Gap statistic - ClussCluster") +
    ggplot2::theme(legend.position ="none", plot.title = ggplot2::element_text(hjust = 0.5))
  return(p1)
}

#' Gene Filter
#'
#' Filter out genes that are not suitable for differential expression analysis.
#'
#' This function \code{filter.genes} takes an expression data frame that has been properly normalized but NOT centered. It returns a list with the slot \code{dat.ft} being the data set that satisfies the pre-set thresholds on minumum mean, standard deviation (sd), and proportion of zeros (n0prop) for each gene.
#'
#' If the data has already been centered, one can still apply the filters of \code{mean} and \code{sd} but not \code{n0prop}.
#'
#' @param dfname name of the expression data frame
#' @param minmean minimum mean expression for each gene
#' @param n0prop minimum proportion of zero expression (count) for each gene
#' @param minsd minimum standard deviation of expression for each gene
#' @examples
#' dat <- matrix(rnbinom(300*60, mu = 2, size = 1), 300, 60)
#' dat_filtered <- filter_gene(dat, minmean=2, n0prop=0.2, minsd=1)
#' @export
#'
filter_gene = function(dfname, minmean=2, n0prop=0.2, minsd=1)
{
  keep1 <- (rowMeans(dfname) >= minmean)
  cat(sum(keep1), "out of", length(keep1), "genes have mean expression >=", minmean, fill=TRUE)

  keep2 <- (rowMeans(dfname > 1e-5) >= n0prop)
  cat(sum(keep2), "out of", length(keep2), "genes have proportion of non-zero expression >=", n0prop, fill=TRUE)

  keep3 <- (apply(dfname, 1, sd) >= minsd)
  cat(sum(keep3), "out of", length(keep3), "genes have standard deviation of expression >=", minsd, fill=TRUE)

  keep <- (keep1 & keep2 & keep3)
  cat("Overall,", sum(keep), "out of", length(keep), "genes are kept.", fill=TRUE)

  dat.ft <- dfname[keep, ]

  return(list(dfname=dfname, dat.ft=dat.ft, index=keep, minmean=minmean, n0prop=n0prop, minsd=minsd))
}
