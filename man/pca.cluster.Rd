% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcanalysis.R
\name{pca.cluster}
\alias{pca.cluster}
\title{K-means clustering and plotting function}
\usage{
pca.cluster(x, pca, k)
}
\arguments{
\item{x}{input data matrix or data frame (numerical)}

\item{pca}{principal component object returned from the run.pca() function}

\item{k}{optimal k value (numerical)}
}
\description{
This function performs a k-means clustering analysis on a numerical data matrix with a specified k value. It assigns samples into k clusters and outputs a PCA biplot in which each sample is colored according to cluster.
}
\examples{
pca.cluster(data, data_pca, 3)
}
\keyword{kmeans}
\keyword{plot}
\keyword{plotting}
