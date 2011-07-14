\name{JoSAE.sample.data}
\alias{JoSAE.sample.data}
\docType{data}
\title{
  Sample plots of the Norwegian National Forest Inventory (NNFI) with
  a variable of interest and an auxiliary variable
}
\description{
  %%  ~~ A concise (1-5 lines) description of the dataset. ~~
  Above ground forest biomass over all tree species is the variable of
  interest. Mean canopy height derived from a photogrammetric canopy
  height model of 20~cm geometric and 10~cm radiometric
  resolution is the auxiliary variable.
}
\usage{data(JoSAE.sample.data)}
\format{
  A data frame with 145 observations on the following 4 variables.
  \describe{
    \item{\code{sample.ID}}{a numeric vector}
    \item{\code{domain.ID}}{a numeric vector}
    \item{\code{biomass.ha}}{a numeric vector of the variable of interest}
    \item{\code{mean.canopy.ht}}{a numeric vector of the auxiliary variable}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
  %%  ~~ reference to a publication or URL from which the data were
  %%  obtained ~~
  Breidenbach, J. and Astrup, R. (2011)
}
\references{
  %%  ~~ possibly secondary sources and usages ~~
  Breidenbach, J. and Astrup, R. (2011)
}
\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
  \code{\link{JoSAE-package}} for more examples
}
\examples{
data(JoSAE.sample.data)
## maybe str(JoSAE.sample.data) ; plot(JoSAE.sample.data) ...
plot(biomass.ha~mean.canopy.ht,JoSAE.sample.data)
}
\keyword{datasets}
