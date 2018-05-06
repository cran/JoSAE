\name{JoSAE-package}
\alias{JoSAE-package}
\alias{JoSAE}
\docType{package}
\title{
  Provides functions for some small area estimators and their
  mean squared errors.
}

\description{ This package implements unit-level (Battese et al. 1988)
  and area-level EBLUP (Fay and Herriot 1979), and GREG (Sarndal 1984)
  estimators as well as their variance/MSE estimators. It also contains
  data and a vignette that explain its use. Heteroscedasticity can
  be considered. }

\details{

The aim in the analysis of sample surveys is frequently to derive
estimates of subpopulation characteristics. Often, the sample available
for the subpopulation is, however, too small to allow a reliable
estimate. If an auxiliary variable exists that is correlated with the
variable of interest, small area estimation (SAE) provides methods to
solve the problem (Rao 2003, Rao and Molina 2015). 

The purpose of this package is primarily to document the functions used
in the publications Breidenbach and Astrup (2012) and Breidenbach et
al. (2018). The data used in Breidenbach and Astrup (2012) are provided.
A subset of the data used in Breidenbach et al. (2018) is also provided
for testing the functions. The vignette further documents the
publication Breidenbach et al. (2015).

You might wonder why this package is called JoSAE. Well, first of all,
JoSAE sounds good (if pronounced like the female name). The other reason
was that the packages SAE and SAE2 already exist (Gomez-Rubio,
2008). They are, however, not available on CRAN and unmaintained (as of
July 2011). They also do not seem to implement the variance estimators
that we needed. So I just combined SAE with the first part of my name.

}
\author{
  Johannes Breidenbach
  
  Maintainer: Johannes Breidenbach <job@nibio.no>
}

\note{
  All the implemented functions/estimators are described in Rao (2003). This
  package merely makes the use of the estimators easier for the users
  that are not keen on programming. Especially the EBLUP variance
  estimator would require some effort.

  Today, there are several well programmed SAE packages available on
  CRAN that also provide the functions described here. This was not the
  case when the first version of the package was uploaded. This package
  is therefore mostly to document the publications mentioned above. With
  respect to heteroskedasticity, this package also implements functions
  that, as far as I can see, are not available in other packages (as of
  May 2018).

}

\references{
  Battese, G. E., Harter, R. M. & Fuller, W. A. (1988), An error-components model
  for prediction of county crop areas using survey and satellite data
  Journal of the American Statistical Association, 83, 28-36 
  
  Breidenbach, J. and Astrup, R. (2012), Small area estimation of
  forest attributes in the Norwegian National Forest Inventory. European
  Journal of Forest Research, 131:1255-1267.

  Breidenbach, J., Ronald E. McRoberts, Astrup, R. (2016), Empirical
  coverage of model-based variance estimators for 
  remote sensing assisted estimation of stand-level timber volume. Remote
  Sensing of Environment, 173,
  274-281. https://doi.org/10.1016/j.rse.2015.07.026

  Breidenbach, J., Rahlf, J., Magnussen, S., Astrup, R. (2018)
  Unit-level and area-level small area estimation under
  heteroscedasticity using digital aerial photogrammetry data. Remote
  Sensing of Environment. In press.

  Fay, R. E., Herriot, R. A., 1979. Estimates of income for small places: an
  application of James- Stein procedures to census data. Journal of the
  American Statistical Association 74 (366a), 269-277.

  Gomez-Rubio (2008), Tutorial on small area estimation, UseR
  conference 2008, August 12-14, Technische Universitat Dortmund,
  Germany.
  
  Rao, J.N.K. (2003), Small area estimation. Wiley.

  Rao, J. N., Molina, I., (2015). Small area estimation, 2nd Edition. John Wiley & Sons. 

  Sarndal, C. (1984), Design-consistent versus model-dependent estimation for
  small domains Journal of the American Statistical Association, JSTOR,
  624-631 

  Schoch, T. (2011), rsae: Robust Small Area Estimation. R package
  version 0.1-3. 
  
}

\keyword{ package }
\seealso{
  \code{\link{eblup.mse.f.wrap}},
  \code{\link{JoSAE.sample.data}}, \code{\link{JoSAE.domain.data}}, 
  \code{\link{sae.al.f}}, \code{\link{sae.ul.f}}
}
\examples{
	#mean auxiliary variables for the populations in the domains
data(JoSAE.domain.data)
	#data for the sampled elements
data(JoSAE.sample.data)
plot(biomass.ha~mean.canopy.ht,JoSAE.sample.data)

## use the original wrapper function

    #lme model
summary(fit.lme <- lme(biomass.ha ~ mean.canopy.ht, data=JoSAE.sample.data
                       , random=~1|domain.ID))

    #domain data need to have the same column names as sample data or vice versa
d.data <- JoSAE.domain.data
names(d.data)[3] <- "mean.canopy.ht"

result <- eblup.mse.f.wrap(domain.data = d.data, lme.obj = fit.lme)
result

##END: use the original wrapper function


## the same with a newer function that can consider heteroskedasticity
res <-
    sae.ul.f(samp.data=JoSAE.sample.data,
             population.data=d.data,
				#assuming homoskedasticity
             k.ij=rep(1, nrow(JoSAE.sample.data)),
             formula=biomass.ha ~ mean.canopy.ht,
             domain.col="domain.ID",
             sample.id.col="sample.ID",
             neg.sfrac=TRUE)
res$est$est

##END: the same with a newer function that can consider heteroskedasticity

}
