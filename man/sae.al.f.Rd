\name{sae.al.f}
\alias{sae.al.f}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Area-level small area estimation under heteroscedasticity or homoscedasticity.
}
\description{
  Area-level small area estimation, possibly under
  heteroscedasticity. Assumes SRS within domains. 
}
\usage{
sae.al.f(domain.id, n.i, psi.i, formula, data, b.i, type = "RE", verbose
= F, iter = 100, a.conv = 0.001,...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{domain.id}{
	Vector of unique domain IDs corresponding to data.
}
  \item{n.i}{
	Vector of number of samples within each domain i, corresponding to data.
}
  \item{psi.i}{
	Vector of variances of the direct estimator within each domain i, corresponding to data.
}
  \item{formula}{
	Formula for fixed effects.
}
  \item{data}{
	Domain-level data.
}
  \item{b.i}{
	Parameter for modeling heteroscedasticity.
}
  \item{type}{
	Method for estimating variance parameters. Currently, only option is
	"RE" for REML.
}
  \item{verbose}{
	Boolean. Should more info be printed?
}
  \item{iter}{
	Max number of iterations. Default is 100.
}
  \item{a.conv}{
	Threshold for accepting convergency. Default is 0.001
}
  \item{\dots}{
	Parameters forwarded to other functions.
}
}
\details{
  SRS (simple random sampling) is assumed within each domain. If the
  direct estimator is not based on SRS, RV components of the results do
  not make sense.
}
\value{
 \item{results }{Area-level estimates and associated standard errors
   (SE). \code{FH} Fay-Harriot estimate (the standard area-level
   estimator), \code{RV} Rivest-Vandal estimate, considering uncertainty
 in variance parameters. \code{bias.sig.sq.v} Bias in random effect
 variance. \code{synth} Synthetic estimate. \code{trans.resid.v.i}
 Transformed residuals. The other parameters are details and are named
 as in the reference.}
 \item{beta.hat }{Estimated fixed-effects parameters.}
 \item{des.mat }{Design matrix.}
 Other list elements have descriptive names.
}
\references{
Breidenbach et al. (2018) Unit-level and area-level small area
estimation under heteroscedasticity using digital aerial photogrammetry
data. Remote Sensing of Environment.
}
\author{
  Johannes Breidenbach
}
%% \note{
%% %%  ~~further notes~~
%% }
%% %% ~Make other sections like Warning with \section{Warning }{....} ~
%% \seealso{
%% %% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
library(nlme)
#sample data
data(ulal.sub.samp.dat)
#domain means of x
data(ulal.sub.dom.dat)

#The easiest way of getting the data into the right format
#is to fit a unit-level model first.
res <-
    sae.ul.f(samp.data=ulal.sub.samp.dat,
             population.data=ulal.sub.dom.dat,
             k.ij=ulal.sub.samp.dat[,"k.ij.one"],
             formula=w.VMPRHA ~ elev.mean + elev.mean.sq,
             domain.col="stand.ID",
             sample.id.col="plot.ID",
             neg.sfrac=TRUE)

#data for AL
dat.al <- merge(res$data$samp.agg.X.pop[,c("domain.id","n.i","w.VMPRHA.ybar.i","elev.mean.X.pop")],
                res$est$se[,c("domain.id","se.srs")])

#area-level SAE under homoscedasticity
res <- sae.al.f(
    domain.id=dat.al$domain.id , n.i=dat.al$n.i , psi.i=dat.al$se.srs^2
    , formula=w.VMPRHA.ybar.i ~ elev.mean.X.pop , data=dat.al
    , b.i=rep(1, nrow(dat.al))
    , type="RE")

#area-level SAE heteroskedasticity
bi.par <- 0.39#select b.i parameter
#helper function
range01 <- function(x, ...){(x - min(x, ...)) /
                                (max(x, ...) - min(x, ...))}
res <- sae.al.f(
    domain.id=dat.al$domain.id , n.i=dat.al$n.i , psi.i=dat.al$se.srs^2
    , formula=w.VMPRHA.ybar.i ~ elev.mean.X.pop, data=dat.al
    , b.i=range01(dat.al$elev.mean.X.pop)+bi.par
    , type="RE")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~survey }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~models }% __ONLY ONE__ keyword per line
