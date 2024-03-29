\name{extractCPT,extractPOT}
\alias{extractCPT}
\alias{extractCPT.table}
\alias{extractCPT.data.frame}
\alias{extractPOT}
\alias{extractPOT.table}
\alias{extractPOT.data.frame}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Extract conditional probabilities and clique potentials from data}
\description{
  Extract list of conditional probability tables and list of clique potentials from data
}
\usage{
extractCPT(x, graph, smooth = 0)
extractPOT(x, graph, smooth = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A table or dataframe}
  \item{graph}{A graph represented as a graphNEL object.
    For extractCPT,  graph must be a DAG while for extractPOT, graph
    must be undirected triangulated graph.}
  %%\item{V}{The set of vertices for which CPTs should be extracted}
  \item{smooth}{See 'details' below}
}
\details{
  If \code{smooth} is non--zero then \code{smooth} is added to all cell
  counts before normalization takes place.
  
}
\value{
  extractCPT: A list of conditional probability tables
  extractPOT: A list of clique potentials
}
\references{
  S�ren H�jsgaard (2012). Graphical Independence Networks with the
  gRain Package for R. Journal of Statistical Software, 46(10), 1-26.
  \url{http://www.jstatsoft.org/v46/i10/}.
}

\author{S�ren H�jsgaard, \email{sorenh@math.aau.dk}}
% \note{ ~~further notes~~ 

%  ~Make other sections like Warning with \section{Warning }{....} ~
% }
\seealso{\code{\link{compileCPT}}, \code{\link{compilePOT}}, \code{\link{grain}}}
\examples{

## Asia (chest clinique) example:

## Version 1) Specify conditional probability tables.
yn <- c("yes","no")
a    <- cptable(~asia, values=c(1,99),levels=yn)
t.a  <- cptable(~tub+asia, values=c(5,95,1,99),levels=yn)
s    <- cptable(~smoke, values=c(5,5), levels=yn)
l.s  <- cptable(~lung+smoke, values=c(1,9,1,99), levels=yn)
b.s  <- cptable(~bronc+smoke, values=c(6,4,3,7), levels=yn)
e.lt <- cptable(~either+lung+tub,values=c(1,0,1,0,1,0,0,1),levels=yn)
x.e  <- cptable(~xray+either, values=c(98,2,5,95), levels=yn)
d.be <- cptable(~dysp+bronc+either, values=c(9,1,7,3,8,2,1,9), levels=yn)
plist <- compileCPT(list(a, t.a, s, l.s, b.s, e.lt, x.e, d.be))
pn1 <- grain(plist)
q1 <- querygrain(pn1)

## Version 2) Specify DAG and data
data(chestSim100000, package="gRbase")
dgf   <- ~asia + tub * asia + smoke + lung * smoke + bronc * smoke + either * tub * lung + xray * either + dysp * bronc * either
dg    <- dag(dgf)
pp    <- extractCPT(chestSim100000, dg)
cpp2  <- compileCPT(pp)
pn2   <- grain(cpp2)
q2    <- querygrain(pn2)

## Version 2) Specify triangulated undirected graph and data
ugf <- list(c("either", "lung", "tub"), c("either", "lung", "bronc"), 
    c("either", "xray"), c("either", "dysp", "bronc"), c("smoke", 
    "lung", "bronc"), c("asia", "tub"))
gg    <- ugList(ugf)
pp    <- extractPOT(chestSim100000, gg)
cpp3  <- compilePOT(pp)
pn3   <- grain(cpp3)
q3    <- querygrain(pn3)

## Compare results:
str(q1)
str(q2[names(q1)])
str(q3[names(q1)])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
