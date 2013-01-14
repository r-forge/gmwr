mcsMAT2 <-
function (amat, vn = colnames(amat), root = NULL, index = FALSE) 
{
  vn.old <- vn
  if (!is.null(root)){
    vn   <- c(root, setdiffPrim(vn, root))
    root2 <- fmatch(vn, vn.old)-1
  } else {
    root2 <- 0:(ncol(amat)-1)
  }
  #print(root2)

  if (class(amat)=="matrix"){
    a <- mcsMAT_stCpp(amat, root2)
  } else {
    a <- mcsMAT_spCpp(amat, root2)
  }

  #print(a)
  if (index){
    if (a[1]<0){
      NA
    } else {
      a+1
    }
  } else {
    if (a[1]<0){
      character(0)
    } else {
      vn.old[a+1]
    }
  }
}

