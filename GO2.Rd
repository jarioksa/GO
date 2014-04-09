\name{GO2}
\alias{GO2}

\title{
Unconstrained Gaussian Maximum Likelihood Ordination
}
\description{
  The function fits unconstraine maximum likelihood ordination
  with unit-width Gaussian response models.
}
\usage{
GO2(comm, k = 1, tot = max(comm), freqlim = 5, 
   family = c("poisson", "binomial"), far = 10, ...)
\method{anova}{GO2}(object, ...)
spanodev(object, ...)
\method{plot}{GO2}(mod, choices = 1, label = FALSE, marginal = FALSE,
    cex = 0.7, col = 1:6, ...)
\method{predict}{GO2}(object, newdata, type = c("response", "link"), ...)
\method{calibrate}{GO2}(object, newdata, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{comm}{
   Community data frame.
}
  \item{k}{
  Number of dimensions.
}

  \item{tot}{ Total abundance used in Binomial models. This can be
  either a single value for all data, or a vector with value for each
  row of \code{comm}. The default is to use the maximum value in
  matrix.  }

  \item{freqlim}{
  Lower frequency limit for analysed species.
}
  \item{family}{
  Error distribution. Can be either \code{"poisson"} for quasi-Poisson or
  \code{"binomial"} for quasi-Binomial (and must be quoted).
}
  \item{far}{
  Threshold for distances of Gaussian optima that are regarded as alien and
  will frozen in the estimation.
}
  \item{object, mod}{\code{GO2} result object.}
  \item{choices}{The displayed dimension (axis).}
  \item{label}{Label plotted response curves.}
  \item{marginal}{Use marginal responses instead of slices through origin.}
  \item{cex, col}{Size of text and colours of lines.}
  \item{newdata}{New data.}
  \item{type}{Type of prediction, similarly as in \code{\link{predict.glm}}.}
  \item{\dots}{
  Other arguments passed to the function. In \code{GO2} these are passed to
  \code{\link{nlm}} and can include, e.g., \code{iterlim} (which often must
  be set to higher value than the default 100).
}
}
\details{

  Function is under development and unreleased. It will be released
  under different name in \pkg{vegan}. The current version is only
  provided for review purposes. The function and its support functions
  require \pkg{vegan}, although this requirements is not explicitly
  made. The optimization is based on \code{\link{nlm}} function, and passes
  arguments to this function.

  Function \code{anova} can analyse two nested models or a single model
  against null model of flat responses using parametric tests based on
  quasi-Likelihood. Function \code{spanodev} performs similar test 
  split by species. Function \code{predict} returns estimated response
  curves, and \code{newdata} can be gradient locations. Function 
  \code{calibrate} returns estimated gradient locations, and \code{newdata}
  can be community data. 

  The \code{plot} function displays fitted respose curves against one
  ordination axis. In principle, the ordination can be rotated using 
  \pkg{vegan} function \code{\link[vegan]{MDSrotate}}, but this requires
  a version that agrees to analyse \code{GO2} results. Traditional 
  ordination plots of SU scores and species optima can be displayed
  with \code{\link[vegan]{ordiplot}} (\pkg{vegan} package).

}
\value{
  Returns an object of class \code{"GO2"}.
}
\references{
  Oksanen, J, Unconstrained Gaussian maximum likelihood ordination. 
     MS under review.
}
\author{
  Jari Oksanen <jari.oksanen@oulu.fi>
}
\note{
  The ms is under review and the function is under development, and will
  change before possible release (and also have a new name). The current
  file also contains some undocumented functions.
}


\seealso{
  \code{\link[VGAM]{cgo}} in \pkg{VGAM} package. 
}
\examples{
library(vegan) ## *must* be given before using the function
data(varespec)
mod <- GO2(varespec, k=2, far=5, tot=100, family="binomial", iterlim=1000)
plot(mod, label=TRUE)
ordiplot(mod, type="t")
ordiplot(mod, dis="si", type="t")
anova(mod)
mod1 <- update(mod, k=1)
anova(mod1, mod)
spanodev(mod1)
spanodev(mod1, mod)
}

\keyword{ multivariate }
