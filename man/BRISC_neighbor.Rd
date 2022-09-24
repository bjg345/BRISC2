\name{BRISC2_neighbor}
\alias{BRISC2_neighbor}
\title{Function for finding set of nearest neighbors for BRISC2}

\description{
  The function \code{BRISC2_neighbor} creates the set of nearest neighbors for a given set of coordinates, which can be used as an input
  for \code{"neighbor"} argument in \code{BRISC2_estimation}. This is especially useful for avoiding often computationally
  intensive nearest neighbor finding scheme in case of multiple application of \code{BRISC2_estimation} on a fixed set of
  coordinates.
}

\usage{
BRISC2_neighbor(coords, n.neighbors = 15, n_omp = 1,
               order = "Sum_coords", search.type = "tree",
               verbose = TRUE, ordering = NULL, tol = 12
)
}

\arguments{

  \item{coords}{an \eqn{n \times 2}{n x 2} matrix of the observation
                coordinates in \eqn{R^2} (e.g., easting and northing). }

  \item{n.neighbors}{number of neighbors used in the NNGP. Default value is 15. }

  \item{n_omp}{number of threads to be used, value can be more than 1 if source code is compiled with OpenMP support.
               Default is 1. }

  \item{order}{keyword that specifies the ordering scheme to be used in ordering the observations. Supported keywords are:
               \code{"AMMD"} and \code{"Sum_coords"} for approximate Maximum Minimum Distance and sum of coordinate based ordering,
               respectively. Default value is \code{"Sum_coords"}. \eqn{n > 65} is required for \code{"AMMD"}. Ignored, if \code{ordering} in not \code{NULL}. }

  \item{search.type}{keyword that specifies type of nearest neighbor search algorithm to be used. Supported keywords are:
                     \code{"brute"}, \code{"tree"} and \code{"cb"}. \cr \code{"brute"} and \code{"tree"} provide the same result, though
                     \code{"tree"} should be faster. \code{"cb"} implements fast code book search described in Ra and Kim (1993)
                     modified for NNGP. If locations do not have identical coordinate values on the axis used for the nearest neighbor
                     ordering (see \code{order} argument) then \code{"cb"} and \code{"brute"} should produce identical neighbor sets.
                     However, if there are identical coordinate values on the axis used for nearest neighbor ordering, then \code{"cb"}
                     and \code{"brute"} might produce different, but equally valid, neighbor sets, e.g., if data are on a grid.
                     Default value is \code{"tree"}. }

  \item{verbose}{if \code{TRUE}, information regarding OpenMP support and progress of the algorithm
                 is printed to the screen. Otherwise, nothing is printed to the screen. Default value is \code{TRUE}. }

  \item{ordering}{if not \code{NULL}, denotes the \eqn{n}{n} length integer vector of ordering of the input coordinates and is used as the ordering of the coordinates for determination of the set of nearest neighbors. }

  \item{tol}{the input observation coordinates, response and the covariates are rounded to this many places after the decimal. The default value is 12. }


}


\value{
 A list containing information regarding nearest neighbors which can be used as an input
  for \code{"neighbor"} argument in \code{BRISC2_estimation}.
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
