\name{querygraph}
\alias{querygraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Query a graph}
\description{
  \code{queryg} is a general function for querying a graph object,
  specifically graphs as created with \code{newug} and \code{newdag}.  
}
\usage{
querygraph(object, op, set = NULL, set2 = NULL, set3 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A graph object; i.e. either an unidrected graph (ugsh)
    or a directed acyclic graph (dagsh)}
  \item{op}{The query operation, see 'details' below.}
  \item{set,set2,set3}{Possible arguments to a graph query of type \code{type}}
}
\details{
  The \code{op} can be:
  \itemize{
    \item {adj:} Nodes adjacant to \code{set}
    \item {an:} Ancestors of \code{set}
    \item {ancestralGraph:} Ancestral graph induced by \code{set}
    \item {ancestralSet:} Ancestral set of \code{set}
    \item {cl:} Closure of \code{set}
    \item {ch:} Children of \code{set}
    \item {maxClique:} The cliques
    \item {connectedComp} The connected components
    \item {edges:} Edges of graph
    \item {ne:} Neighbours of \code{set}
    \item {nodes:} Nodes of graph
    \item {is.complete:}
    \item {edgeList}
    \item {vpar}
    \item {is.simplicial:}
    \item {is.triangulated:}
    \item {pa:} Parents of \code{set}
    \item {separates:} Is \code{set} and \code{set2} separated by \code{set3}
    \item {simplicialNodes:} The simplicial nodes of graph
    \item {subgraph:} Subgraph induced by \code{set}
  }
}

\value{
  Depending on the \code{type}, the output will be either a new graph or
  a vector or a list.
}
%\references{ ~put references to the literature/web site here ~ }
\author{S�ren H�jsgaard,  \email{sorenh@math.aau.dk}}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{
  \code{\link{ug}},
  \code{\link{dag}}
}
\examples{

ug0 <- ug(~a:b, ~b:c:d, ~e)

querygraph(ug0, "nodes")
querygraph(ug0, "edges")

querygraph(ug0, "subgraph", c("b","c","d","e"))

querygraph(ug0, "adj", "c")
querygraph(ug0, "closure", "c") 
querygraph(ug0, "is.simplicial", "b")
querygraph(ug0, "simplicialNodes")

querygraph(ug0, "is.complete")
querygraph(ug0, "is.complete", c("b","c","d"))
querygraph(ug0, "maxClique")

querygraph(ug0, "is.triangulated")
querygraph(ug0, "is.decomposition", "a","d",c("b","c"))

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
