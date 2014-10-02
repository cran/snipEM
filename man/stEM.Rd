\name{stEM}
\alias{stEM}
\title{Snipping and trimming for location and scatter estimation with
casewise and cellwise outliers}
\description{
  Computes an estimator optimizing the Gaussian likelihood over a
  snipping and trimming set. 
}
\usage{
stEM(X, V, tol = 1e-4, maxiters = 500, maxiters.S = 1000, print.it = FALSE) 
}
\arguments{
  \item{X}{Data.}
  \item{V}{Binary matrix of the same size as \code{X}. Zeros correspond to
  initial snipped entries, rows of zeros correspond to initial trimmed
  entries.}
  \item{tol}{Tolerance for convergence. Default is \code{1e-4}.}
  \item{maxiters}{Maximum number of iterations for the SM algorithm. Default is \code{500}.}
  \item{maxiters.S}{Maximum number of iterations of the inner greedy snipping algorithm. Default is \code{1000}.}
   \item{print.it}{Logical; if \code{TRUE}, partial results are print. Default is \code{FALSE}.}
}
\details{
   This function combines computes the \code{snipEM} estimator of Farcomeni
   (2014) with trimming. Optimization over a trimming set is performed
   via usual concentration steps (Rousseeuw and van Driessen, 1999). 
   It therefore provides a robust estimate of
   location and scatter in presence of entry-wise and case-wise
   outliers. The number of snipped entries and trimmed rows is kept
   fixed throughout. \code{V} must contain at least one row of zeros
   (otherwise use \code{\link{snipEM}}). 
}
\value{
  A list with the following elements: 
	\tabular{ll}{
	\code{mu} \tab Estimated location. \cr
	\code{S} \tab Estimated scatter matrix. \cr 
	\code{V} \tab Final (optimal) V matrix. \cr 
	\code{lik} \tab Gaussian log-likelihood at convergence. \cr 
	\code{iter} \tab Number of outer iterations before convergence. \cr 
	}
}
\references{
Farcomeni, A. (2014) Snipping for robust k-means clustering under component-wise contamination, \emph{Statistics and Computing}, \bold{24}, 909-917

Farcomeni, A. (2014) Robust constrained clustering in presence of entry-wise outliers, \emph{Technometrics}, \bold{56}, 102-111   

Rousseeuw, P. J. and Van Driessen, K. (1999) A fast algorithm for the
minimum covariance determinant estimator, \emph{Technometrics}, \bold{41}, 212-223.
}
\author{Alessio Farcomeni \email{alessio.farcomeni@uniroma1.it}, Andy
Leung \email{andy.leung@stat.ubc.ca}}
\seealso{
  \code{\link{sclust}}, \code{\link{snipEM}},
  \code{\link{generate.casecontam}},
  \code{\link{generate.cellcontam}}, \code{\link{emve}}, \code{\link{sumlog}},
  \code{\link{ldmvnorm}} 
}
\examples{
\dontrun{
# 10\% cellwise outliers
Xc <- generate.cellcontam(100, 10, 100, 5, 0.10)
A <- Xc$A   # true correlation matrix
Uc <- Xc$u  # (true outlier) binary matrix 
Xc <- Xc$x

# 5\% case-wise outliers
Xct <- generate.casecontam(100, 10, 100, 10, 0.05)
Xc <- rbind(Xc, Xct$x[ rowSums(Xct$u == 1) == 10,])
trim <- rep(1,105)
trim[101:105] <- 0  # (true casewise outlier) binary vector
Uc <- rbind(Uc, matrix(1,5,10))

# initial V
V <- matrix(1, 105, 10)
# identify rows to be trimmed
km <- kmeans(Xc,2)$clust
wtrim <- which(km==which.min(table(km)))
V[wtrim,] <- 0
Vtmp <- V[-wtrim,]

# identify cells to be snipped
Vtmp[!is.na(match(as.vector(Xc[-wtrim,]),boxplot(as.vector(Xc[-wtrim,]),plot=FALSE)$out))] <- 0
Xna <- Xc[-wtrim,]
Xna[ which( Vtmp == 0) ] <- NA
resEMVE <- emve(Xna)
Vtmp <- snipEM.initialV(Xc[-wtrim,], Vtmp, resEMVE@mu, resEMVE@S)
V[-wtrim,] <- Vtmp$V

# Performance in LRT distance
resSTEM <- stEM(Xc, V)
slrt(resSTEM$S, A)
slrt(cov(Xc), A)
}}
