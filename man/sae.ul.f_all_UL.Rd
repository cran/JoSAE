\name{sae.ul.f}
\alias{sae.ul.f}
\alias{ul.data.prep.f}
\alias{ul.reml.f}
\alias{ul.est.f}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Unit-level small area estimation under heteroscedasticity.
}
\description{
  The only function that should be used is \code{sae.ul.f} which wraps
  the data preparation (\code{ul.data.prep.f}), parameter estimation
  (\code{ul.reml.f}), and EBLUP and MSE 
  estimation (\code{ul.est.f}) functions. The other functions
  are helpers that do not need to be called by the user directly.
}
\usage{
sae.ul.f(...)
%
ul.data.prep.f(k.ij, samp.data, population.data, formula, domain.col,
sample.id.col, neg.sfrac, pop.r, sum.i.k.ij.sq.r, N.i, ...)
%
ul.reml.f(samp.data, formula, samp.agg.X.pop, y.name, X.names, ...)
%
ul.est.f(samp.data, samp.agg.X.pop, X.names, y.name, beta.hat,
cov.beta.hat, sig.sq.e, sig.sq.v, V.bar.ee, V.bar.vv, V.bar.ve,
neg.sfrac, resid = F, ...) 
}
\arguments{
  \item{samp.data}{
	Data frame of unit-level sample data of the response and explanatory
	variables. Transformations of x have to be calculated before hand,
	not using formula.
}
  \item{k.ij}{
	Vector of the same length and corresponding to \code{samp.data}
	containing 
	the \code{k_ij} values for modeling heteroscedasticity. A vector of
	1's if heteroscedasticity is not 
	considered. \code{k_ij} is squared and used as the input to the
	\code{varFix} variance function in \code{lme}. 
}
  \item{population.data}{
	Data frame of domain-level means of the explanatory variables. Names
	have to match to sample data.
}
  \item{formula}{
	The fixed effects formula to be used in the mixed model using \code{lme}.
}
  \item{domain.col}{
	Character string identifying the column name in \code{sample.data}
	and \code{population.data}
	containing unique domain IDs.
}
  \item{sample.id.col}{
	Character string identifying the column name in \code{sample.data}
	containing unique sample IDs.
}
  \item{neg.sfrac}{
	Boolean. Are sampling fractions negligible?
}
  \item{pop.r}{
	Only required if \code{neg.sfrac} is FALSE. Same as
	\code{population.data} but calculated with the sample elements
	removed. 
}
  \item{sum.i.k.ij.sq.r}{
	Only required if \code{neg.sfrac} is FALSE. Same as the values in the
	\code{population.k.column} but calculated with the sample elements
	removed.
}
  \item{N.i}{
	Only required if \code{neg.sfrac} is FALSE. Number of elements
	within domain i.
  }
  \item{samp.agg.X.pop}{
	Data frame containing the aggregated sample data and domain-level
	explanatory variables. 
}
  \item{y.name}{
	Character string identifying the column name of the response.
}
  \item{X.names}{
	Vector of character strings identifying the column names of the
	explanatory variables.
  }
  \item{beta.hat}{
	Vector of estimated fixed effects.
}
  \item{cov.beta.hat}{
	Covariance matrix of fixed effects.
}
  \item{sig.sq.e}{
	Residual variance.
}
  \item{sig.sq.v}{
	Random effect variance.
}
  \item{V.bar.ee}{
	Assymptotic variance of the residual variance.
}
  \item{V.bar.vv}{
	Assymptotic variance of the random effect variance.
}
  \item{V.bar.ve}{
	Assymptotic covariance of the residual and random effect variance.
}
  \item{resid}{
	Boolean. Should residuals be returned?
}
  \item{\dots}{
	Parameters forwarded to other functions.
}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
  These functions can also be used for eblups
  without heteroskedasticity as the 
  the older functions around \code{eblup.mse.f.wrap}.
}
\value{
 \item{est }{List of point estimates and standard errors based on
   different estimators. \code{est}: Point estimates of EBLUP, SRS, GREG,
   Survey regression (svreg), and EBLUP synthetic estimators. Also
   contains random effects. \code{se}: Standard errors. EBLUP global was used
   in the reference publication. \code{resids}: Different forms of
   raw and transformed residuals.}
 \item{var.pars }{List of estimated variance parameters. Among others
   the lme object fitted to the data. See list of arguments above.}
 \item{data }{List of data frames in the format required by the other
   functions.}
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
%%   %%  ~~further notes~~
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

#eblup under homoskedasticity
res <-
    sae.ul.f(samp.data=ulal.sub.samp.dat,
             population.data=ulal.sub.dom.dat,
             k.ij=ulal.sub.samp.dat[,"k.ij.one"],
             formula=w.VMPRHA ~ elev.mean + elev.mean.sq,
             domain.col="stand.ID",
             sample.id.col="plot.ID",
             neg.sfrac=TRUE)

#eblup under heteroskedasticity
res <-
    sae.ul.f(samp.data=ulal.sub.samp.dat,
             population.data=ulal.sub.dom.dat,
             k.ij=ulal.sub.samp.dat[,"k.ij.em.0.48"],
             formula=w.VMPRHA ~ elev.mean + elev.mean.sq,
             domain.col="stand.ID",
             sample.id.col="plot.ID",
             neg.sfrac=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~survey }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~models }% __ONLY ONE__ keyword per line