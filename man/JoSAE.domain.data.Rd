\name{JoSAE.domain.data}
\alias{JoSAE.domain.data}
\docType{data}
\title{
  %%   ~~ data name/kind ... ~~
  Dataframe containing the domains population and number of elements and
  mean of the auxiliary variable
}
\description{
  %%  ~~ A concise (1-5 lines) description of the dataset. ~~
    Auxiliary variable: Mean canopy height derived from a photogrammetric
  canopy height model
}
\usage{data(JoSAE.domain.data)}
\format{
  A data frame with 14 observations on the following 3 variables.
  \describe{
    \item{\code{domain.ID}}{a numeric vector}
    \item{\code{N.i}}{a numeric vector - number of population elements}
    \item{\code{mean.canopy.ht.bar}}{a numeric vector - mean of the
	  auxiliary variable}
  }
}
\source{
%%  ~~ reference to a publication or URL from which the data were
%%  obtained ~~
  Breidenbach, J. and Astrup, R. (2012), Small area estimation of
  forest attributes in the Norwegian National Forest Inventory. European
  Journal of Forest Research, 131:1255-1267.
}
\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
  \code{\link{JoSAE-package}} for more examples
}
\examples{
data(JoSAE.domain.data)

str(JoSAE.domain.data)
}
\keyword{datasets}
