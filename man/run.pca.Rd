% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcanalysis.R
\name{run.pca}
\alias{run.pca}
\title{Principal Component Analysis function}
\usage{
run.pca(x, center = TRUE, scale = FALSE, svd = TRUE)
}
\arguments{
\item{x}{input data matrix or data frame (numerical)}

\item{center}{zero center the data (logical value, default is TRUE)}

\item{scale}{scale the data to unit variance before analysis (logical value, default is FALSE)}

\item{svd}{indicate whether to perform single value decomposition (logical value, default is TRUE), otherwise perform eigen decomposition}
}
\value{
\code{run.pca} returns a list with the following elements:
\itemize{
  \item{Data - }{a preview of the input data matrix}
  \item{k - }{the number of observations in the dataset}
  \item{PCs - }{the matrix of loadings containing the principal components as columns (eg, eigenvectors)}
  \item{sdev - }{standard deviations of the principal components (uses singular values of the matrix if using singular value decomposition, but takes the square root of the eigenvalues of the correlation matrix if using spectral decomposition)}
}
}
\description{
This function performs a principal component analysis on a numerical data matrix.
}
\examples{
data.pca <- run.pca(data, center = TRUE, scale = TRUE, svd = TRUE)
}
\keyword{pca}
