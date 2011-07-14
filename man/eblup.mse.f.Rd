\name{eblup.mse.f}
\alias{eblup.mse.f}
\alias{eblup.mse.f.default}
\alias{eblup.mse.f.c1}
\alias{eblup.mse.f.c2.ai}
\alias{eblup.mse.f.c2}
\alias{eblup.mse.f.c3.asyvarcovarmat}
\alias{eblup.mse.f.c3}
\alias{eblup.mse.f.c3.star}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Helper functions to derive the variance of the EBLUP estimate
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  A set of helper functions. The default returns the first component.
}
\usage{
eblup.mse.f(...)
\method{eblup.mse.f}{default}(lme.obj, n.i, gamma.i, ...)
\method{eblup.mse.f}{c1}(lme.obj, n.i, gamma.i, ...)
\method{eblup.mse.f}{c2.ai}(lme.obj, n.i, gamma.i, X.i, ...)
\method{eblup.mse.f}{c2}(gamma.i, X.i, X.bar.i, sum.A.i, ...)
\method{eblup.mse.f}{c3.asyvarcovarmat}(lme.obj, n.i, ...)
\method{eblup.mse.f}{c3}(lme.obj, asympt.var.covar, n.i, ...)
\method{eblup.mse.f}{c3.star}(lme.obj, asympt.var.covar, n.i, mean.resid.i, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
	the gamma_i value resulting from \code{\link{JoSAE.gamma.i.f}}
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
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A compontent of the EBLUP variance (aka mean squared error). Which
component depends on the method used.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Johannes Breidenbach
}
\note{
%%  ~~further notes~~
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

