\name{querygrain}
\alias{querygrain}
\alias{querygrain.grain}
\alias{setFinding}
\alias{retractFinding}
\alias{getFinding}
\alias{pFinding}
%\alias{print.grainFinding}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Query an independence network}
\description{
  Query an independence network, i.e. obtain the conditional
  distribution of a set of variables given finding (evidence) on other variables. 
}
\usage{
querygrain(object, nodes = nodeNames(object), normalize = TRUE, type =
c("marginal", "joint", "conditional"), result="array", details = 0)

setFinding(object, nodes=NULL, states=NULL, flist=NULL, propagate=TRUE)
retractFinding(object, nodes=NULL, propagate=TRUE)
getFinding(object)
pFinding(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A "grain" object}
  \item{nodes}{A vector of nodes}
  \item{states}{A vector of states (of the nodes given by 'nodes')}
  \item{flist}{An alternative way of specifying findings (evidence), see examples below.}
  \item{propagate}{Should the network be propagated?}
  \item{normalize}{Should the results be normalized to sum to one.}
  \item{type}{Should marginals (for each node), the joint for all nodes,
  or the conditional of the first node given the rest be returned.}
  \item{result}{If "data.frame" the result is returned as a data frame
    (or possibly as a list of dataframes).}
  \item{details}{Debugging information}

}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A list of tables with potentials
}
\references{
  S�ren H�jsgaard (2012). Graphical Independence Networks with the
  gRain Package for R. Journal of Statistical Software, 46(10), 1-26.
  \url{http://www.jstatsoft.org/v46/i10/}.
}

\author{S�ren H�jsgaard, \email{sorenh@math.aau.dk}}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{
  \code{\link{cptable}}
}
\examples{
 testfile <- system.file("huginex", "chest_clinic.net", package = "gRain")
 chest <- loadHuginNet(testfile, details=0)


 qb <- querygrain(chest)
 qb

 lapply(qb, as.numeric) % Safe
 sapply(qb, as.numeric) % Risky

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{utilities}% __ONLY ONE__ keyword per line
