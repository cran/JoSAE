\name{JoSAE-package}
\alias{JoSAE-package}
\alias{JoSAE}
\docType{package}
\title{
  Provides functions for small area estimation and the estimate of
  variances (mean squared errors).
}
\description{
  This package currently implements the unit level EBLUP (Battese et
  al. 1988) and GREG (Sarndal 1984)
  estimators as well as their variance estimators. It also contains data
  and a vignette that explain its use. 
}
\details{
\tabular{ll}{
Package: \tab JoSAE\cr
Type: \tab Package\cr
Version: \tab 0.2\cr
Date: \tab 2011-07-14\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}

The aim in the analysis of sample surveys is frequently to derive
estimates of subpopulation characteristics. Often, the sample available
for the subpopulation is, however, too small to allow a reliable
estimate. If an auxiliary variable exists that is correlated with the
variable of interest, small area estimation (SAE) provides methods to
solve the problem (Rao 2003). 

The purpose of this package is primarily to document the
functions used in the publication of Breidenbach and Astrup (2011). The
data used in that study are also provided.

You might wonder why this package is called JoSAE. Well, first of all,
JoSAE sounds good. The other problem was that the packages SAE
and SAE2 already exist (Gomez-Rubio, 2008). They are, however, not
available on CRAN and unmaintained (as of July 2011). They also do not seem to implement
the variance estimators that we needed. So I just combined SAE with the
first part of my name.

}
\author{
  Johannes Breidenbach
  
  Maintainer: Johannes Breidenbach <job@skogoglandskap.no>
}
\note{
  All the implemented functions/estimators are described in Rao (2003). This
  package merely makes the use of the estimators easier for the users
  that are not keen on programming. Especially the EBLUP variance
  estimator would require some effort.

  After uploading the first JoSAE
  version (0.1), I realized there exists another SAE package that specialzes on
  robust SAE (rsae; Schoch, 2011). It implements reampling variance
  estimators as opposed to analytical variance estimators implemented
  here. The \code{\link[rsae]{rsae}} package comes with a very good
  vignette as well as the Landsat data set 
  used by Battese et al. (1988).

  There are several points where the JoSAE could be improved:

  Only univariate unit-level models with a simple block-diagonal variance structure
  are supported so far.

  The computation is based on loops on the domain level. It would be
  more elegant to use blocked matrices.

  ... many more things that will hopefully improve once I have more
  experience with programming R-packages ...
}

\references{
  Battese, G. E., Harter, R. M. & Fuller, W. A. (1988), An error-components model
  for prediction of county crop areas using survey and satellite data
  Journal of the American Statistical Association, 83, 28-36 
  
  Breidenbach, J. and Astrup, R. (submitted 2011), Small area estimation of
  forest attributes in the Norwegian National Forest Inventory. European
  Journal of Forest Research.

  Rao, J.N.K. (2003), Small area estimation. Wiley.

  Gomez-Rubio (2008), Tutorial on small area estimation, UseR
  conference 2008, August 12-14, Technische Universitat Dortmund,
  Germany.

  Sarndal, C. (1984), Design-consistent versus model-dependent estimation for
  small domains Journal of the American Statistical Association, JSTOR,
  624-631 

  Schoch, T. (2011), rsae: Robust Small Area Estimation. R package
  version 0.1-3. 
  
}

\keyword{ package }
\seealso{
  \code{\link{eblup.mse.f}},
  \code{\link{JoSAE.sample.data}}, \code{\link{JoSAE.domain.data}}
}
\examples{
	#mean auxiliary variables for the populations in the domains
data(JoSAE.domain.data)
	#data for the sampled elements
data(JoSAE.sample.data)
plot(biomass.ha~mean.canopy.ht,JoSAE.sample.data)

## the easy way: use the wrapper function to compute EBLUP and GREG estimates and variances

    #lme model
summary(fit.lme <- lme(biomass.ha ~ mean.canopy.ht, data=JoSAE.sample.data
                       , random=~1|domain.ID))

    #domain data need to have the same column names as sample data or vice versa
d.data <- JoSAE.domain.data
names(d.data)[3] <- "mean.canopy.ht"

result <- eblup.mse.f.wrap(domain.data = d.data, lme.obj = fit.lme)
result

##END: the easy way


##the hard way: compute the EBLUP MSE components yourself
    #get an overview of the domains
        #mean of the response and predictor variables from the sample. For the response this is the sample mean estimator.
tmp <- aggregate(JoSAE.sample.data[,c("biomass.ha", "mean.canopy.ht")]
                            , by=list(domain.ID=JoSAE.sample.data$domain.ID), mean)
names(tmp)[2:3] <- paste(names(tmp)[2:3], ".bar.sample", sep="")

        #number of samples within the domains
tmp1 <- aggregate(cbind(n.i=JoSAE.sample.data$biomass.ha)
                  , by=list(domain.ID=JoSAE.sample.data$domain.ID), length)

        #combine it with the population information of the domains
overview.domains <- cbind(JoSAE.domain.data, tmp[,-1], n.i=tmp1[,-1])


    #fit the models
        #lm - the auxiliary variable explains forest biomass rather good
summary(fit.lm <- lm(biomass.ha ~ mean.canopy.ht, data=JoSAE.sample.data))

        #lme
summary(fit.lme <- lme(biomass.ha ~ mean.canopy.ht, data=JoSAE.sample.data
                       , random=~1|domain.ID))

    #mean lm residual -- needed for GREG
overview.domains$mean.resid.lm <- aggregate(resid(fit.lm)
                                            , by=list(domain.ID=JoSAE.sample.data$domain.ID)
                                            , mean)[,2]

    #mean lme residual -- needed for EBLUP.var
overview.domains$mean.resid.lme <- aggregate(resid(fit.lme, level=0)
                                             , by=list(domain.ID=JoSAE.sample.data$domain.ID)
                                             , mean)[,2]

    #synthetic estimate
overview.domains$synth <- predict(fit.lm
                                  , newdata=
                                  data.frame(mean.canopy.ht=JoSAE.domain.data$mean.canopy.ht.bar))



    #GREG estimate
overview.domains$GREG <- overview.domains$synth +  overview.domains$mean.resid.lm

    #EBLUP estimate
overview.domains$EBLUP <- predict(fit.lme
                                     , newdata=data.frame(mean.canopy.ht=JoSAE.domain.data$mean.canopy.ht.bar
                                       , domain.ID=JoSAE.domain.data$domain.ID)
                                     , level=1)

    #gamma
overview.domains$gamma.i <- eblup.mse.f.gamma.i(lme.obj=fit.lme
                                          , n.i=overview.domains$n.i)

    #variance of the sample mean estimate
overview.domains$sample.var <-
    aggregate(JoSAE.sample.data$biomass.ha
          , by=list(domain.ID=JoSAE.sample.data$domain.ID), var)[,-1]/overview.domains$n.i

    #variance of the GREG estimate
overview.domains$GREG.var <-
    aggregate(resid(fit.lm)
          , by=list(domain.ID=JoSAE.sample.data$domain.ID),var)[,-1]/overview.domains$n.i

    #variance of the EBLUP
        #compute the A.i matrices for all domains (only needed once)
domain.ID <- JoSAE.domain.data$domain.ID
            #initialize the result vector
a.i.mats <- vector(mode="list", length=length(domain.ID))
for(i in 1:length(domain.ID)){
    print(i)
    a.i.mats[[i]] <- eblup.mse.f.c2.ai(gamma.i=overview.domains$gamma.i[overview.domains$domain.ID==domain.ID[i]]
                                 , n.i=overview.domains$n.i[overview.domains$domain.ID==domain.ID[i]]
                                 , lme.obj=fit.lme
                                 , X.i=as.matrix(cbind(i=1
                                   , x=JoSAE.sample.data[JoSAE.sample.data$domain.ID==domain.ID[i], "mean.canopy.ht"]
                                   ))
         )
}
                #add all the matrices
sum.A.i.mats <- Reduce("+", a.i.mats)

        #the assymptotic var-cov matrix
asy.var.cov.mat <- eblup.mse.f.c3.asyvarcovarmat(n.i=overview.domains$n.i
                                           , lme.obj=fit.lme)


        #put together the variance components
		###### Some changes are required here, if you apply it to own data!
result <- NULL
for(i in 1:length(domain.ID)){
    print(i)
    #first comp
    mse.c1.tmp <- eblup.mse.f.c1(gamma.i=overview.domains$gamma.i[overview.domains$domain.ID==domain.ID[i]]
                           , n.i=overview.domains$n.i[overview.domains$domain.ID==domain.ID[i]]
                           , lme.obj=fit.lme)
    #second comp
    mse.c2.tmp <- eblup.mse.f.c2(gamma.i=overview.domains$gamma.i[overview.domains$domain.ID==domain.ID[i]]
                           , X.i=as.matrix(cbind(i=1##cbind!!
                                   , x=JoSAE.sample.data[JoSAE.sample.data$domain.ID==domain.ID[i]
                                     , "mean.canopy.ht"]##change to other varnames if necessary
                             ))
                           , X.bar.i =as.matrix(rbind(i=1##rbind!!
                                   , x=JoSAE.domain.data[JoSAE.domain.data$domain.ID==domain.ID[i]
                                     , "mean.canopy.ht.bar"]##change to other varnames if necessary
                             ))
                           , sum.A.i = sum.A.i.mats
                           )
    #third comp
    mse.c3.tmp <- eblup.mse.f.c3(n.i=overview.domains$n.i[overview.domains$domain.ID==domain.ID[i]]
                           , lme.obj=fit.lme
                           , asympt.var.covar=asy.var.cov.mat)
    #third star comp
    mse.c3.star.tmp <- eblup.mse.f.c3.star( n.i=overview.domains$n.i[overview.domains$domain.ID==domain.ID[i]]
                                     , lme.obj=fit.lme
                                     , mean.resid.i=overview.domains$mean.resid.lme[overview.domains$domain.ID==domain.ID[i]]
                                     , asympt.var.covar=asy.var.cov.mat)
    #save result
    result <- rbind(result, data.frame(kommune=domain.ID[i]
                                       , n.i=overview.domains$n.i[overview.domains$domain.ID==domain.ID[i]]
                                       , c1=as.numeric(mse.c1.tmp), c2=as.numeric(mse.c2.tmp)
                                       , c3=as.numeric(mse.c3.tmp), c3star=as.numeric(mse.c3.star.tmp)))
}

            #derive the actual EBLUP variances
overview.domains$EBLUP.var.1 <- result$c1 + result$c2 + 2* result$c3star
overview.domains$EBLUP.var.2 <- result$c1 + result$c2 + result$c3 + result$c3star


    #display the estimates and the sampling error (sqrt(var)) in two tables
round(data.frame(overview.domains[,c("domain.ID", "n.i", "N.i")],
                 overview.domains[,c("biomass.ha.bar.sample", "GREG", "synth", "EBLUP")]),2)
        #the sampling errors of the eblup is mostly smaller than the one of the greg estimate.
        #both are always smaller than the sample mean variance.
round(data.frame(overview.domains[,c("domain.ID", "n.i", "N.i")],
                 sqrt(overview.domains[,c("sample.var", "GREG.var", "EBLUP.var.1", "EBLUP.var.2")])),2)

##END: the hard way

}
