\name{eblup.mse.f.wrap}
\alias{eblup.mse.f.wrap}
\alias{eblup.mse.f.gamma.i}
\alias{eblup.mse.f.c1}
\alias{eblup.mse.f.c2.ai}
\alias{eblup.mse.f.c2}
\alias{eblup.mse.f.c3.asyvarcovarmat}
\alias{eblup.mse.f.c3}
\alias{eblup.mse.f.c3.star}
%%\alias{eblup.mse.f.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Functions to calculate the variance of an EBLUP estimate.
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Functions to calculate the EBLUP MSE (=variance). The wrap function
  calls all the other functions and calculates EBLUP, GREG, SRS, and Synthetic
  estimates for domain means as well as the variances of the EBLUP, GREG
  and SRS estimates.
}
\usage{
eblup.mse.f.wrap(domain.data, lme.obj, debug=F, ...)
eblup.mse.f.gamma.i(lme.obj, n.i, ...)
eblup.mse.f.c1(lme.obj, n.i, gamma.i, ...)
eblup.mse.f.c2.ai(lme.obj, n.i, gamma.i, X.i, ...)
eblup.mse.f.c2(gamma.i, X.i, X.bar.i, sum.A.i, ...)
eblup.mse.f.c3.asyvarcovarmat(lme.obj, n.i, ...)
eblup.mse.f.c3(lme.obj, asympt.var.covar, n.i, ...)
eblup.mse.f.c3.star(lme.obj, asympt.var.covar, n.i, mean.resid.i, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{domain.data}{
	data set with the mean of the auxiliary variables for every
	domain including the domain ID. Names of the variables must be the same as in the
	unit-level sample data that were used to fit the lme model.
  }  
  \item{lme.obj}{
	%%     ~~Describe \code{lme.obj} here~~
	a linear mixed-effects model generated with \code{\link[nlme]{lme}}	
  }
  \item{n.i}{
	%%     ~~Describe \code{n.i} here~~
	the number of samples within domain i	
  }
  \item{gamma.i}{
	%%     ~~Describe \code{gamma.i} here~~
	the gamma_i value resulting from the gamma.i method of this function
  }
  \item{X.i}{
	the design matrix of sampled elements in domain i
  }
  \item{X.bar.i}{
	mean of the populations elements design matrix in domain i
  }
  \item{sum.A.i}{
	sum of the domains A_i matrices resulting from \code{\link{eblup.mse.f.c2.ai}}
  }
  \item{asympt.var.covar}{
	the asymptotic variance-covariance matrix of the mixed-effects model
	resulting from \code{\link{eblup.mse.f.c3.asyvarcovarmat}}
  }
  \item{mean.resid.i}{
	the mean residual of the fixed-part of the linear mixed-effects
	model in domain i (i.e., use \code{level=0} in \code{\link[nlme]{predict.lme}})
	}  
  \item{\dots}{
	%%     ~~Describe \code{\dots} here~~
	forward attributes to other functions. Not used so far.
  }
  \item{debug}{
	details are printed if true - Only used by the wrapper function
  }  
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
  Most users will probably only use the convenient wrap
  function. Nonetheless, all components of the EBLUP variance can also
  be calculated separately.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A compontent of the EBLUP variance (aka mean squared error). Which
component depends on the function used.

The wrap function returns a data
frame with many entrys for every domain: The domain-level predictor
variables obtained from the domain.data, the mean of the predictor
variables and response (aka the sample mean or SRS estimate) observed at the samples, the number of samples (n.i.sample),
the mean residuals of a lm (fitted using the fixed part of the lme) and
the lme (mean.resid.lm, mean.resid.lme), 
the synthetic (Synth), EBLUP, and GREG estimates of the mean of the
variable of interest, the gamma_i value,
the variance of the means for the sample and GREG estimates (.var.mean), the
components of the EBLUP variance (c1-c3star), the results of the first and the second
method (cf. Rao 2003) to derive the EBLUP variance
(EBLUP.var.1, EBLUP.var.2), and the standard errors derived from the
variances (.se).
}
\references{
  %% ~put references to the literature/web site here ~
  Breidenbach and Astrup (2012), Small area estimation of
  forest attributes in the Norwegian National Forest Inventory. European
  Journal of Forest Research, 131:1255-1267.
}
\author{
Johannes Breidenbach
}
\note{
  %%  ~~further notes~~
  Currently, only random intercept mixed-effects models with homogeneous
  variance structure are supported.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
  \code{\link{JoSAE-package}} for more examples
}
\examples{
library(nlme)
data(JoSAE.sample.data)
#fit a lme 
summary(fit.lme <- lme(biomass.ha ~ mean.canopy.ht, data=JoSAE.sample.data
                       , random=~1|domain.ID))
#calculate the first component of the EBLUP variance for a domain with 5 samples
eblup.mse.f.c1(fit.lme, 5, 0.2)
}

