# Principal Component Analysis Package
#
# Here are the functions for the pcanalysis package
#
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' @title Principal Component Analysis function
#' @description This function performs a principal component analysis on a numerical data matrix.
#' @param x input data matrix or data frame (numerical)
#' @param center zero center the data (logical value, default is TRUE)
#' @param scale scale the data to unit variance before analysis (logical value, default is FALSE)
#' @param svd indicate whether to perform single value decomposition (logical value, default is TRUE), otherwise perform eigen decomposition
#' @return \code{run.pca} returns a list with the following elements:
#' \itemize{
#'   \item{Data - }{a preview of the input data matrix}
#'   \item{k - }{the number of observations in the dataset}
#'   \item{PCs - }{the matrix of loadings containing the principal components as columns (eg, eigenvectors)}
#'   \item{sdev - }{standard deviations of the principal components (uses singular values of the matrix if using singular value decomposition, but takes the square root of the eigenvalues of the correlation matrix if using spectral decomposition)}
#' }
#' @keywords pca
#' @export
#' @examples
#' data.pca <- run.pca(data, center = TRUE, scale = TRUE, svd = TRUE)

#Perform the PCA on a data matrix
run.pca <- function(x, center = TRUE, scale = FALSE, svd = TRUE) {
  x <- as.matrix(x) #input as matrix
  k <- ncol(x) #number of samples/variables
  x <- scale(x, center = center, scale = scale)
  if(svd) {        #SVD method (default)
    sv <- svd(x)
    pcs <- sv$v
    d <- sv$d
    colnames(pcs) <- paste('PC', c(1:ncol(pcs)), sep = '')
    rownames(pcs) <- colnames(x)
    pca <- list(Data=head(x), k=k, PCs=pcs,    #final list for output
                sdev= d / sqrt(max(1, nrow(x) - 1)))
  } else {
    R <- cov(x)   #covarance matrix
    eig <- eigen(R)
    pcs <- eig$vectors[,1:k] #get PCs
    colnames(pcs) <- paste('PC', c(1:ncol(pcs)), sep = '')
    rownames(pcs) <- colnames(x)
    eigenvals <- eig$values[1:k]
    pca <- list(Data=head(x), k=k, PCs=pcs,     #final list for output
                Eigenvalues=eigenvals, sdev = na.omit(sqrt(eigenvals)))
  }
  return(pca)
}


#' @title Principal Component Analysis Plotting function
#' @description This function outputs biplots for principal components 1 vs. 2 and 1 vs. 3, and plots the percentage of variance explained by each component.
#' @param pca principal component object returned from the run.pca() function
#' @keywords plot pca plotting
#' @export
#' @examples
#' pca.plot(data.pca)

#Plotting function
pca.plot <- function(pca){
  plot(pca$PCs[,1], pca$PCs[,2], #plots PCs 1 and 2
       main = paste("Principal Components 1 and 2"),
       xlab = paste("PC 1"),
       ylab= paste("PC 2"),
       pch=1)
  plot(pca$PCs[,1], pca$PCs[,3], #Plots PCs 1 and 3
       main = paste("Principal Components 1 and 3"),
       xlab = paste("PC 1"),
       ylab= paste("PC 3"),
       pch=1)
  p <- 100*pca$sdev^2 / sum(pca$sdev^2)  #calculates variance for each PC and take the percentage
  percent <- data.frame(p = p, PC = as.factor(1:length(p))) #creates a data frame for easy plotting
  barplot(as.vector(percent$p), ylim=c(0,100), #Plots the percent of variance explained by each PC
          main="Percent of variance explained",
          xlab="PC",
          ylab="Percent",
          xlim=c(0,11),
          names.arg=percent$PC)
}


#' @title Screeplot for k-means analysis
#' @description This function performs a k-means analysis for each k value and outputs a scree plot of the total within sum of squares for cross-validation.
#' @param data input data matrix or data frame (numerical)
#' @param kmax maximum k value (numerical)
#' @param iter.max maximum number of iterations for k-means (numerical, default is 10)
#' @keywords screeplot kmeans
#' @export
#' @examples
#' pca.scree(data, 10, iter.max = 10)


#Creates a scree plot for cross-validation to find optimal k
pca.scree <- function(data, kmax, iter.max = 10){
  WSS <- as.vector(kmax)
  for(i in 1:kmax){
    WSS[i] <- kmeans(data, centers = i, iter.max = iter.max)$tot.withinss
  }
  plot(1:kmax, WSS, type="b",
                  xlab = "Number of clusters",
                  ylab = "Within-group SS",
                  main = "K-means Scree Plot")

  return(WSS)
}


#' @title K-means clustering and plotting function
#' @description This function performs a k-means clustering analysis on a numerical data matrix with a specified k value. It assigns samples into k clusters and outputs a PCA biplot in which each sample is colored according to cluster.
#' @param x input data matrix or data frame (numerical)
#' @param pca principal component object returned from the run.pca() function
#' @param k optimal k value (numerical)
#' @keywords kmeans plotting plot
#' @export
#' @examples
#' pca.cluster(data, data_pca, 3)

#assign samples to clusters based on PCS
#input should be the raw data matrix, the output from pcanalysis.pca, and the optimal k
pca.cluster <- function(x, pca, k) {
  data <- t(as.matrix(x))
  pca <- pca$PCs
  z <- kmeans(data, k)
  princomps <- as.data.frame(pca)
  princomps$Sample <- rownames(princomps)
  clusters <- data.frame("Sample" = z[1])
  clusters$Sample <- rownames(clusters)
  rownames(clusters) <- NULL
  nd <- merge(princomps, clusters, by="Sample")
  plot(nd$PC1, nd$PC2, col=factor(nd$cluster), main = "Clustered PCA", xlab ="PC1", ylab="PC2", pch=1)
  return(clusters)
}
