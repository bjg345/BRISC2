\name{BRISC2_order}
\alias{BRISC2_order}
\title{Function for ordering coordinates with BRISC2}

\description{
  The function \code{BRISC2_order} outputs the ordering for a set of coordinates, which can be used as an input for
  \code{"ordering"} argument in \code{BRISC2_estimation}. This is especially useful for avoiding often computationally
  intensive location ordering scheme in case of multiple application of \code{BRISC2_estimation} on a fixed set of
  coordinates.
}

\usage{
BRISC2_order(coords, order = "Sum_coords", verbose = TRUE)
}

\arguments{

  \item{coords}{an \eqn{n \times 2}{n x 2} matrix of the observation
                coordinates in \eqn{R^2} (e.g., easting and northing). }

  \item{order}{keyword that specifies the ordering scheme to be used in ordering the observations. Supported keywords are:
               \code{"AMMD"} and \code{"Sum_coords"} for approximate Maximum Minimum Distance and sum of coordinate based ordering,
               respectively. Default value is \code{"Sum_coords"}. \eqn{n > 65} is required for \code{"AMMD"}.}

  \item{verbose}{if \code{TRUE}, progress of the algorithm is printed to the screen. Otherwise, nothing is printed to the screen.
                 Default value is \code{TRUE}. }


}


\value{
 An integer vector of ordering of the input coordinates which can be used as an input for
  \code{"ordering"} argument in \code{BRISC2_estimation}.
}

\references{

  Saha, A., & Datta, A. (2018). BRISC2: bootstrap for rapid inference on spatial
  covariances. Stat, e184, DOI: 10.1002/sta4.184.

  Datta, A., S. Banerjee, A.O. Finley, and A.E. Gelfand. (2016)
  Hierarchical Nearest-Neighbor Gaussian process models for large
  geostatistical datasets. Journal of the American Statistical
  Association, 111:800-812.

  Guinness, J. (2018) Permutation and Grouping Methods for
  Sharpening Gaussian Process Approximations, Technometrics,
  DOI: 10.1080/00401706.2018.1437476, \cr
  https://github.com/joeguinness/gp_reorder/tree/master/R .

}

\author{
  Arkajyoti Saha \email{arkajyotisaha93@gmail.com}, \cr
  Abhirup Datta \email{abhidatta@jhu.edu}
}

\examples{

set.seed(1)
n <- 1000
coords <- cbind(runif(n,0,1), runif(n,0,1))


ordering_result <- BRISC2_order(coords)
}

\keyword{model}
