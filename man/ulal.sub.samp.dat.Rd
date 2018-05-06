\name{ulal.sub.samp.dat}
\alias{ulal.sub.samp.dat}
\docType{data}
\title{
  %%   ~~ data name/kind ... ~~
  Unit-level sample data.
}
\description{
  %%  ~~ A concise (1-5 lines) description of the dataset. ~~
  73 sample plots clustered in 15 stands in the
  municipality of Stokke, Norway.
}
\usage{data("ulal.sub.samp.dat")}
\format{
  A data frame with 73 observations on the following 7 variables.
  \describe{
    \item{\code{plot.ID}}{a character vector with unique plot IDs.}
    \item{\code{w.VMPRHA}}{a numeric vector. Timber volume per hectare
	  (m3/ha).}
    \item{\code{stand.ID}}{a factor with unique stand IDs.} %with levels \code{H01036} \code{H03068} \code{H04063} \code{H04102} \code{H07019} \code{H07044} \code{H10082} \code{H12106} \code{H13019} \code{H16069} \code{H16078} \code{H21097} \code{H23092} \code{H25078} \code{H29055} \code{H30097} \code{L02089} \code{L04031} \code{L04046} \code{L04063} \code{L05084} \code{L05097} \code{L06047} \code{L06093} \code{L07017} \code{L07048} \code{L08034} \code{L08079} \code{L09081} \code{L15079} \code{L16017} \code{L16035} \code{L16083} \code{L17097} \code{L18099} \code{L19089} \code{L20055} \code{L21027} \code{L21061} \code{L22105} \code{L23088} \code{L24023} \code{L24104} \code{L25075} \code{L26047} \code{L26049} \code{L26066} \code{L29041} \code{L29103} \code{L29107} \code{S00059} \code{S07099} \code{S08038} \code{S30038} \code{S39023} \code{S43098} \code{S46074} \code{S48062} \code{S52099} \code{S64049} \code{S65101} \code{S65109} \code{S74075} \code{S76031} \code{S78032}}
    \item{\code{elev.mean}}{a numeric vector. Mean height
	  obtained from digital aerial photogrammetry.}
    \item{\code{elev.mean.sq}}{a numeric vector. Mean height squared
	  obtained from digital aerial photogrammetry.}
    \item{\code{k.ij.one}}{a numeric vector of ones (used in cases where
	heteroscedasticity can be ignored)}
    \item{\code{k.ij.em.0.48}}{a numeric vector. \code{elev.mean} to the
	  power of 0.48 for modeling heteroscedasticity.}
  }
}
%% \details{
%% %%  ~~ If necessary, more details than the __description__ above ~~
%% }
\source{
Breidenbach et al. (2018) Unit-level and area-level small area
estimation under heteroscedasticity using digital aerial photogrammetry
data. Remote Sensing of Environment.
}
%% \references{
%% %%  ~~ possibly secondary sources and usages ~~
%% }
\examples{
data(ulal.sub.samp.dat)
}
\keyword{datasets}
