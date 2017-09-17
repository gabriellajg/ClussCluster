# this file contains the hidden and internal functions for the ClussCluster package

initial <- function(sc){
  if (is.null(sc$nclust)) sc$nclust <- ncol(sc$centers)
  a <- w <- matrix(0, nrow=nrow(sc$x), ncol=sc$nclust)
  x <- apply(sc$x, 2, as.numeric)
  for (k in 1 : sc$nclust)
  {
    ind <- which(sc$theta == k)
    sub.x <- x[,ind]
    for (j in 1:ncol(sub.x)){
      m <- unlist(sub.x[,j])
      rowMeans((x - m) ^ 2) - rowMeans((sub.x - m) ^ 2)
    }
    a[, k] <- rowSums(
      apply(sub.x, 2, function(m) rowMeans((x - m) ^ 2) - rowMeans((sub.x - m) ^ 2))
      )
    w[, k] <- get.w(a=a[, k], s=sc$s)
  }
  sc$a <- a
  sc$w <- w
  sc$wbcss <- sum(w * a)
  return(sc)
}

## This function only update the class label of one sample: sample i.
update.theta <- function(sc, i)
{
  y.new <- y.old <- sc$theta[i] #old label of sample i
  sc.update <-list(); for (m in 1:sc$nclust){sc.update[[m]]<-sc}  # store SCs for each k
  if (length(which(sc$theta == sc$theta[i])) > 1) # we try to change it only if it is not singleton
  {
    distsq <- rep(NA, sc$nclust)
    for (k in 1 : sc$nclust)
    {
      sc.cp <- sc # a copy is created so that we can change theta
      if (k!=y.old){
        sc.cp$theta[i] <- k
        sc.cp <- get.target.value(sc.cp,i,y.old,k)
      }
      sc.update[[k]] <- sc.cp
      distsq[k] <- sc.cp$wbcss
    }
    y.new <- which.max(distsq)
    #cat("\t\tOriginal label:", sc$theta[i], "\n")
    #cat("\t\tdistsq:", distsq, "\n")
    #cat("\t\tNew label:", y.new, "\n")
    #if (y.new != sc$theta[i])
    {
      #cat("Updated!\n")
    }
  }
  return(list(y.new=y.new, sc.update=sc.update))
}

get.target.value <- function(sc, i, y.o, y.n)
{
  a <- sc$a
  a[, c(y.o,y.n)] <- 0
  x <- apply(sc$x, 2, as.numeric)
  for (k in c(y.o,y.n))
  {
    ind <- which(sc$theta == k)
    for (i in ind)
    {
      a[, k] <- a[, k] + rowMeans((x - x[, i]) ^ 2) - rowMeans((x[, ind, drop=FALSE] - x[, i]) ^ 2)
    }
  }
  sc$a <- a
  sc$wbcss <- sum(sc$w * a)
  return(sc)
}

update.w <- function(sc, i, y.o, y.n)
{
  a <- sc$a
  w <- sc$w
  for (k in c(y.o,y.n))
  {
    w[, k] <- get.w(a=a[, k], s=sc$s)
  }
  sc$w <- w
  sc$wbcss <- sum(w * a)
  return(sc)
}

get.w <- function(a, s)
{
  ap <- pmax(a, 0)
  w <- ap / l2norm(ap)

  if (sum(w) > s)
  {
    f <- function(Df)
    {
      ss <- softthred(ap, Df)
      ss2n <- l1norm(ss) / (l2norm(ss) + 1e-6) - s
    }
    D <- uniroot(f, interval=c(0, max(ap)), tol=0.0001)$root
    ss <- softthred(ap, D)
    w <- ss / l2norm(ss)
  }

  return(w)
}

softthred <- function(x, x0) {return(sign(x) * (abs(x) - x0) * (abs(x) > x0))}
l1norm <- function(x) {return(sum(abs(x)))}
l2norm <- function(x) {return(sqrt(sum(x ^ 2)))}

