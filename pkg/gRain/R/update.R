"update.cpt-grain" <- function(object,  ...){
  cl <- match.call(expand.dots=TRUE)

  args <- list(...)
  arg.names <- names(args)
  if ("cptlist" %in% arg.names){
    object$cptlist <- args$cptlist
    if (object$isCompiled){
      pot.with.1   <- .defaultPotentialList(object$rip,object$universe)     
      ## FIXME lav en .insert1 funktion og g�r s�:
      ## FIXME pot.with.1   <- .insert1(object$origCQpot)
      newCQpot            <- .insertCpt(args$cptlist, pot.with.1, object$rip, details=0)
      object$origCQpot    <- newCQpot
      object$tempCQpot    <- newCQpot
      object$equilCQpot   <- .insertNA(pot.with.1)
      object$isPropagated <- FALSE		
    }
  }

  if ("origCQpot" %in% arg.names){
    if (object$isCompiled){
      object$origCQpot    <- args$origCQpot
      object$tempCQpot    <- args$origCQpot
      object$equilCQpot   <- .insertNA(object$equilCQpot)
      object$cptlist      <- .insertNA(object$cptlist)
      object$isPropagated <- FALSE		
    }
  }

  ## FIXME: pas p� n�r objektet ikke er kompileret...
  object
}  
