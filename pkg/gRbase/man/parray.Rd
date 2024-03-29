\name{parray}
\alias{parray}
\alias{as.parray}
\alias{data2parray}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Representation of and operations on multidimensional tables}
\description{General representation of multidimensional 
  tables (by \code{parray} objects).
}

\usage{
parray(varNames, levels, values = 1, normalize = "none", smooth = 0)
as.parray(values, normalize="none", smooth=0)
data2parray(data, varNames=NULL, normalize="none", smooth=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{varNames}{Names of variables defining table; can be a right hand
  sided formula.}
  \item{levels}{Either 1) a vector with number of levels of the factors in
    varNames or 2) a list with specification of the levels of the factors in
    varNames. See 'examples' below.}
  \item{values}{Values to go into the parray}
  \item{normalize}{Either "none", "first" or "all". Should result be
    normalized, see 'Details' below.} 
  \item{smooth}{Should values be smoothed, see 'Details' below.}
  \item{data}{A dataframe, a table, an xtabs (a cross classified
    contingency table) a matrix (with dimnames) or a vector (with dimnames).}
}
\details{
  
  A \code{parray} object represents a table defined by a set of variables and
  their levels, together with the values of the table. E.g. f(a,b,c) can
  be a table with a,b,c representing levels of binary variable
  
  If \code{normalize="first"} then for each configuration of all other
  variables than the first, the probabilities are normalized to sum to one. Thus f(a,b,c)
  becomes a conditional probability table of the form p(a|b,c). If
  \code{normalize="all"} then the sum over all entries of f(a,b,c) is
  one. 
  
  If \code{smooth} is positive then \code{smooth} is added to 
  \code{values}  before normalization takes place.

  \code{as.parray} can be used for coercing an array or an xtabs to a parray
  object.

}
\value{
  An object of class \code{parray}.
}
%\references{ ~put references to the literature/web site here ~ }
\author{S�ren H�jsgaard, \email{sorenh@math.aau.dk}}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}

\seealso{
  \code{\link{tableOp}},
  \code{\link{tableMargin}}
%   \code{\link{tableOp}}
%   \code{\link{tableOp}}
}

\examples{

t1 <- parray(c("gender","answer"), list(c('male','female'),c('yes','no')), values=1:4)
t1 <- parray(~gender:answer, list(c('male','female'),c('yes','no')), values=1:4)
t1 <- parray(~gender:answer, c(2,2), values=1:4)

t2 <- parray(c("answer","category"), list(c('yes','no'),c(1,2)), values=1:4+10)
t3 <- parray(c("category","foo"), c(2,2), values=1:4+100)

varNames(t1)
nLevels(t1)
valueLabels(t1)


## Create 1-dimensional vector with dim and dimnames
x1 <- 1:5
as.parray(x1)
x2 <- parray("x", levels=length(x1), values=x1)
dim(x2)
dimnames(x2)

## Matrix
x1 <- matrix(1:6, nrow=2)
as.parray(x1)
parray(~a:b, levels=dim(x1), values=x1)

## Extract parrays from data
## 1) a dataframe
data(cad1) 
data2parray(cad1, ~Sex:AngPec:AMI)
data2parray(cad1, c("Sex","AngPec","AMI"))
data2parray(cad1, c(1,2,3))
## 2) a table
data2parray(UCBAdmissions,c(1,2), normalize="first")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
