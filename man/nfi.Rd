\name{nfi}
\alias{nfi}
\docType{data}
\title{
  Sample plots of the Norwegian National Forest Inventory (NFI) with
  a variable  of interest and an auxiliary variable  
}
\description{
  %%  ~~ A concise (1-5 lines) description of the dataset. ~~
  A total of 131 sample plots of the Norwegian National Forest Inventory (NFI) with
  timber volume interpolated and extrapolated to the year 2011 as the variable
  of interest and mean vegetation height derived from image matching  as
  the auxiliary variable. 
}
\usage{data(nfi)}
\format{
  A data frame with 131 observations on the following 2 variables.
  \describe{
    \item{\code{vol.2011}}{numeric vector of the variable of interest}
    \item{\code{Elev.Mean}}{numeric vector of the auxiliary variable}
  }
}
\source{
  %%  ~~ reference to a publication or URL from which the data were
  %%  obtained ~~
  Breidenbach, J., Ronald E. McRoberts, Astrup, R. (2015), Empirical
  coverage of model-based variance estimators for 
  remote sensing assisted estimation of stand-level timber volume. Remote
  Sensing of Environment. In press. 
}
\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
  \code{\link{JoSAE-package}} and Vignette for more examples
}
\examples{
data(nfi)
plot(vol.2011~Elev.Mean, nfi)
}
\keyword{datasets}
