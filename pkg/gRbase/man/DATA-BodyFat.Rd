\name{BodyFat}
\alias{BodyFat}
\docType{data}
\title{
Body Fat Data
}
\description{
Estimates of the percentage of body fat determined by underwater
weighing and various body circumference measurements for 252 men.
}
\usage{data(BodyFat)}
\format{
  A data frame with 252 observations on the following 15 variables.
  \describe{
    \item{\code{Density}}{Density determined from underwater weighing, a numeric vector}
    \item{\code{BodyFat}}{Percent body fat from Siri's (1956) equation, a numeric vector}
    \item{\code{Age}}{in years, a numeric vector}
    \item{\code{Weight}}{in lbs, a numeric vector}
    \item{\code{Height}}{in inches, a numeric vector}
    \item{\code{Neck}}{circumference in cm, a numeric vector}
    \item{\code{Chest}}{circumference in cm, a numeric vector}
    \item{\code{Abdomen}}{circumference in cm, a numeric vector}
    \item{\code{Hip}}{circumference in cm, a numeric vector}
    \item{\code{Thigh}}{circumference in cm, a numeric vector}
    \item{\code{Knee}}{circumference in cm, a numeric vector}
    \item{\code{Ankle}}{circumference in cm, a numeric vector}
    \item{\code{Biceps}}{circumference in cm, a numeric vector}
    \item{\code{Forearm}}{circumference in cm, a numeric vector}
    \item{\code{Wrist}}{circumference in cm, a numeric vector}
  }
}
%\details{
%}
\source{
For more information see http://lib.stat.cmu.edu/datasets/bodyfat
}
\references{
Bailey, Covert (1994). _Smart Exercise: Burning Fat, Getting Fit_,
Houghton-Mifflin Co., Boston, pp. 179-186.

Behnke, A.R. and Wilmore, J.H. (1974). _Evaluation and Regulation of Body
Build and Composition_, Prentice-Hall, Englewood Cliffs, N.J.

Siri, W.E. (1956), "Gross composition of the body", in _Advances in 
Biological and Medical Physics_, vol. IV, edited by J.H. Lawrence and C.A.
Tobias, Academic Press, Inc., New York.

Katch, Frank and McArdle, William (1977). _Nutrition, Weight Control, and
Exercise_, Houghton Mifflin Co., Boston.

Wilmore, Jack (1976). _Athletic Training and Physical Fitness: Physiological
Principles of the Conditioning Process_, Allyn and Bacon, Inc., Boston.
}
\examples{
data(BodyFat)
head(BodyFat)
}
\keyword{datasets}
