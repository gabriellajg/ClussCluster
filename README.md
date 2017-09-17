# ClussCluster

This package implements a new method ClussCluster to simultaneously perform clustering analysis and signature gene selection on high-dimensional transcriptome data sets. To do so, ClussCluster incorporates a Lasso-type regularization penalty term to the objective function of K-means so that cell-type-specific signature genes can be identified while clustering the cells.

### Installing

To install this package and load it into R, do the following:

```
devtools::install_github("gabriellajg/ClussCluster")
library(ClussCluster)
```

