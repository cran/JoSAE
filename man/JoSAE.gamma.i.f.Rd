\name{JoSAE.gamma.i.f}
\alias{JoSAE.gamma.i.f}
\alias{JoSAE.gamma.i.f.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  Returns the gamma_i value for the GREG-like design of the EBLUP estimator
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  The gamma_i value is also necessary for the variance of the EBLUP estimate
}
\usage{
JoSAE.gamma.i.f(...)

\method{JoSAE.gamma.i.f}{default}(lme.obj, n.i, ...)
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
  \item{\dots}{
	%%     ~~Describe \code{\dots} here~~
	further attributes forwarded to other functions; not used so far
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
  The gamma_i value
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
#calculate gamma_i for a domain with 5 samples
JoSAE.gamma.i.f(fit.lme, 5)
}
