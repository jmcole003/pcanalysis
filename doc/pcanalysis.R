## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
set.seed(2746)

## ----setup--------------------------------------------------------------------
library(pcanalysis)

## ----viewdata-----------------------------------------------------------------
head(gene_exp_data)[,1:7]

## ----viewdata_t---------------------------------------------------------------
gnxp_t <- t(gene_exp_data)
head(gnxp_t)[,1:7]

## ----run_pca_1----------------------------------------------------------------
gene_pca <- run.pca(gnxp_t, center = TRUE, scale = FALSE, svd = TRUE)

## ----view_pca-----------------------------------------------------------------
head(gene_pca$Data)[,1:5]
gene_pca$k
head(gene_pca$PCs)[,1:5]
head(gene_pca$sdev)[1:5]

## ----run_pca_2----------------------------------------------------------------
gene_pca_eigen <- run.pca(gnxp_t, center = TRUE, scale = FALSE, svd = FALSE)

## ----view_pca_2---------------------------------------------------------------
head(gene_pca_eigen$PCs)[,1:5]
gene_pca_eigen$Eigenvalues[1:5]
gene_pca_eigen$sdev[1:5]

## ----plot, fig.align="center", fig.height=5, fig.width=6, message=FALSE, warning=FALSE----
pca.plot(gene_pca)

## ----screeplot, fig.align="center", fig.height=5, fig.width=6, message=FALSE, warning=FALSE----
pca.scree(gnxp_t, kmax=10, iter.max=10)

## ----cluster, fig.align="center", fig.height=5, fig.width=6, message=FALSE, warning=FALSE----
#using head() to truncate the output to the first six samples
head(pca.cluster(gnxp_t, gene_pca, k=2))

