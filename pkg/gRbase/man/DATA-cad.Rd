\name{cad}
\alias{cad1}
\alias{cad2}
\docType{data}
\title{Coronary artery disease data}
\description{
  A cross classified table with observational data from a Danish heart
  clinic. The response variable is CAD. 
}
\usage{data(cad1)}
\format{
  A data frame with 236 observations on the following 14 variables.
  \describe{
    \item{\code{Sex}}{a factor with levels \code{Female} \code{Male}}
    \item{\code{AngPec}}{a factor with levels \code{Atypical} \code{None} \code{Typical}}
    \item{\code{AMI}}{a factor with levels \code{Definite} \code{NotCertain}}
    \item{\code{QWave}}{a factor with levels \code{No} \code{Yes}}
    \item{\code{QWavecode}}{a factor with levels \code{Nonusable} \code{Usable}}
    \item{\code{STcode}}{a factor with levels \code{Nonusable} \code{Usable}}
    \item{\code{STchange}}{a factor with levels \code{No} \code{Yes}}
    \item{\code{SuffHeartF}}{a factor with levels \code{No} \code{Yes}}
    \item{\code{Hypertrophi}}{a factor with levels \code{No} \code{Yes}}
    \item{\code{Hyperchol}}{a factor with levels \code{No} \code{Yes}}
    \item{\code{Smoker}}{a factor with levels \code{No} \code{Yes}}
    \item{\code{Inherit}}{a factor with levels \code{No} \code{Yes}}
    \item{\code{Heartfail}}{a factor with levels \code{No} \code{Yes}}
    \item{\code{CAD}}{a factor with levels \code{No} \code{Yes}}
  }
}
\details{
  cad1: Complete dataset, 236 cases.
  cad2: Incomplete dataset, 67 cases. Information on (some of) the variables
  Hyperchol, Smoker, Inherit is missing.
}

%\source{
%  ~~ reference to a publication or URL from which the data were obtained ~~
%}

\references{
H�jsgaard, S�ren and Thiesson, Bo (1995). BIFROST - Block recursive
models Induced From Relevant knowledge, Observations and Statistical
Techniques. Computational Statistics and Data Analysis,  vol. 19, p. 155-175  

Hansen, J. F. (1980). The clinical diagnoisis of ichaeme heart disease
du to coronary artery disease. Danish Medical Bulletin
}

\examples{
data(cad1)
## maybe str(cad1) ; plot(cad1) ...
}
\keyword{datasets}
