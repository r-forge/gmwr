predict.grain <- function(object, response, predictors=setdiff(names(newdata), response),
                                   newdata, type="class", ...){


  if (!inherits(object, "compgrain"))
    object <- compile(object, propagate=TRUE)
  
  type <- match.arg(type, c("class","distribution"))
  nstate <- nodeStates(object,response)
  if (missing(predictors))
    predictors  <- setdiff(names(newdata),response)
  
  p.evec <- rep(NA, nrow(newdata))
  ans <- lapply(nstate, function(a)
                {
                  v <- matrix(NA, ncol=length(a),nrow=nrow(newdata))
                  colnames(v) <- a
                  v}
                )

  for (i in 1:nrow(newdata)){
    case      <- newdata[i,predictors, drop=FALSE]
    objecttmp1    <- setFinding(object, nodes=names(case), states=case)

    p.e       <- pFinding(objecttmp1)
    #print(p.e)

    if (p.e<1e-32){
      cat(sprintf("The finding for row %i has probability 0 in then model. Consider using the 'smooth' argument when building the network. Exiting...\n", i))
      return(NULL)
    }
    
    p.evec[i] <- p.e
    for (j in 1:length(response)){
      ##pj   <- nodeMarginal(objecttmp1, response[j])[[1]]$values
      pj   <- nodeMarginal(objecttmp1, response[j])[[1]] ## BRIS
      #print(pj)
      ans[[j]][i,] <- pj
    }
  } 

  #print(ans)
  if (type=="class"){
    ns <- nodeStates(object, response)
    for (i in 1:length(ans)){
      a<-ans[[i]]
      #print(a)
      mlc <- apply(a,1,which.max)
      #print(mlc)
      
      ans[[i]] <- ns[[i]][mlc]
    }
  }
  
  value <- list(pred=ans, pFinding=p.evec)
  return(value)
}
