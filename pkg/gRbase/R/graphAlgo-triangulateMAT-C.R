## #######################################
##
## Interface to RcppEigen function for triangulation (based on sparse adj. matrix)
##
## S�ren H�jsgaard, December, 2012
## FIXME: Add tests, .Rd files and vignette
##
## #######################################

triangulateMAT_spCpp <- function(XX_, LL_=rep(2,ncol(XX_))){
  res <- .Call("C_triangulateMAT_sp", XX_, log(LL_), 0, 0
               ,package="gRbase"
               )
  ##dimnames(res) <- dimnames(XX_)
  res
}

triangulateMAT_stCpp <- function(XX_, LL_=rep(2,ncol(XX_))){
  as(triangulateMAT_spCpp(asdgCMatrix(XX_), LL_), "matrix")
  
}
