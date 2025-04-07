#' @title Import object of class 'DimReduction' to put into CYTdata object
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param importedDimReductionObject a S4 object of class 'DimReduction' to put in CYTdata object
#'
#' @return a S4 object of class 'CYTdata'
#'
#' @export
#'

addManualDimRed <- function(CYTdata,
                            DimRedName,
                            DimRedDF,
                            checkOverwrite = TRUE) {

  CYTdata = checkValidity(CYTdata, mode = "error")

  if (checkOverwrite && (DimRedName %in% colnames(CYTdata@DimRedData@cellDimRed))) {
    reply <- readline(prompt=paste("A dimensionality reduction named", DimRedName, "is already performed, do you still want to continue and overwrite (yes or no): "))
    while (!tolower(reply) %in% c("yes", "y", "no", "n")) { reply <- readline(prompt="Reply not standard, please answer yes or no: ") }
    if (tolower(reply) %in% c("no", "n")){ stop("Function stopped, CYTdata object unchanged") }
  }

  checkmate::qassert(DimRedIDs, "D2")
  CYTdata@DimRedData@cellDimRed[[DimRedName]] = DimRedDF

  CYTdata = checkValidity(CYTdata, mode = "warning")
  return(CYTdata)
}


#' @title Dimensionality Reduction of data using dimension reduction techniques
#'
#'
#' @description This function aims to generate coordinates of data in an 2D-reduced space, stored in a CYTdata object.
#'
#' The algorithm used available are :
#' - Principal component analysis (PCA), a linear dimension reduction technique
#' - Uniform Manifold Approximation and Projection (UMAP), a non-linear dimension reduction technique
#' - t-distributed stochastic neighbor embedding (t-SNE), a non-linear dimension reduction technique
#' - LargeVis-like method (lvish), a non-linear dimension reduction technique
#' The whole set of cell markers or specific cell markers can be used during the dimensionality reduction process.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param markers a character vector providing the cell markers to use to generate the reduced data
#' @param type a character value containing the type of RD method to use. Possible values are: "UMAP", "tSNE", "lvish", "PCA" (default = UMAP)
#' @param seed a numeric value providing the random seed to use during stochastic operations
#' @param ... additional arguments passed on to method from R package.
#' For PCA, please refer to prcomp method from stats package : https://rdrr.io/r/stats/stats-package.html
#' For UMAP, please refer to umap method from uwot package : https://cran.r-project.org/web/packages/uwot/uwot.pdf
#' For t-SNE, please refer to Rtsne method from Rtsne package : https://cran.r-project.org/web/packages/Rtsne/Rtsne.pdf
#' For lvish, please refer to lvish method from uwot package : https://cran.r-project.org/web/packages/uwot/uwot.pdf
#'
#' @return a S4 object of class 'CYTdata'
#'
#' @export
#'
# Compute Laplacian Eigenmaps init tsNE
# leim = LaplacianEigenmaps()
# emb <- leim@fun(as(t(log10(expr + 1)),"dimRedData"), leim@stdpars)
# plot(emb@data@data)
# Compute PCA init tsNE
# PC = prcomp(exp.markers, center=TRUE, scale=FALSE) # Compute PCA
# tsne <- Rtsne::Rtsne(exp.markers, Y_init = PC$x[,1:2], ...)

runDimRed <- function(CYTdata,
                      name,
                      markers = NULL,
                      type = c("UMAP", "tSNE", "lvish", "PCA"),
                      seed = 42,
                      checkOverwrite = TRUE,
                      ...){

  CYTdata = checkValidity(CYTdata, mode = "error")

  if (checkOverwrite && (name %in% names(CYTdata@cellData@cellDimRed))) {
    reply <- readline(prompt=paste("A dimensionality reduction named", name, "is already performed, do you still want to continue and overwrite (yes or no): "))
    while (!tolower(reply) %in% c("yes", "y", "no", "n")) { reply <- readline(prompt="Reply not standard, please answer yes or no: ") }
    if (tolower(reply) %in% c("no", "n")){ stop("Function stopped, CYTdata object unchanged") }
  }

  type = match.arg(type)
  checkmate::qassert(type, "S1")
  markers = checkorderMarkers(CYTdata, markers = markers, order = TRUE, checkDuplicates = TRUE)

  data = CYTdata@cellData@cellExprs[,markers, drop=FALSE]
  cat("\nGenerate Dimensionality reduction using", type, " :")
  cat("\n\n - Markers used to generate 2D-reduced data space : ", paste0(markers, collapse = ", "), "\n\n")

  checkmate::qassert(seed, "N1")
  set.seed(seed)

  switch(type,
         UMAP = {
           manifold <- uwot::umap(data, ...)
         },
         tSNE = {
           manifold <- Rtsne::Rtsne(data, ...)$Y
         },
         lvish = {
           manifold <- uwot::lvish(data, ...)
         },
         PCA = {
           PC <- stats::prcomp(data, ...) # Compute PCA
           manifold = PC$x[,1:2]
         })

  coordinates = data.frame(manifold)
  colnames(coordinates) = paste(type, 1:2, sep="")

  CYTdata@cellData@cellDimRed[[name]] = coordinates

  CYTdata = checkValidity(CYTdata, mode = "warning")
  return(CYTdata)
}


#' @title
#' QC for dimensionality reduction : Computes correlation between pairwise distances in the high-dimensional space and in the embedding
#' QC for dimensionality reduction : Computes proportion of nearest neighbours preservation in reduced dimension in comparison to high dimension
#'
#' @description CPD quantifies preservation of the global, or macroscropic structure.
#' A performant RD should give a boxplot with linear relationships betweens pairwise distance in the different data spaces
#'
#' @description The fraction of k-nearest neighbours in the original high-dimensional data that are preserved as k-nearest neighbours in the embedding
#' It compute the average across all n points. KNN quantifies preservation of the local, or microscopic structure.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#'
#' @param subsetSize a numeric being the sample size used for computation
#'
#' @param method a character being the method use to compute correlation. Possible values are : "pearson", "kendall", "spearman". Default value is "pearson"
#' @param knn a numeric being the number of nearest neighbors to compute
#' @param knc a numeric being the number of nearest classes to compute
#'
#' @return a list containing correlation coefficient and boxplot of pairwise distances
#'
#' @export
#'

computeDimRedQC <- function(CYTdata,
                            DimRed,
                            markers = NULL,
                            subsetProportion = 0.2,

                            correlation = TRUE,
                            KNN = TRUE,
                            KNC = TRUE,


                            correlationMethod = c("pearson", "kendall", "spearman"),
                            knn = 5,
                            knc = 5,
                            kncClusteringName = NULL){

  CYTdata = checkValidity(CYTdata, mode = "error")

  if (!DimRedName %in% names(CYTdata@cellData@cellDimRed)) {
    stop("Error : DimRedName argument (", DimRedName, ") give a dimensionality reduction data frame missing in cellDimRed list (Dimensionality redection availbale in the list :",
         paste0(names(CYTdata@cellData@cellDimRed), collapse=", "), ").")
  }
  if (subsetSize<=0 | subsetSize>=1) {
    stop("Error : 'subsetProportion' argument must be positive integer between 0 and 1(is a proportion of cells to be subsampled for QC computation.")
  }

  subIdx = sample(rownames(CYTdata@cellData@cellExprs), round(subsetProportion*nrow(CYTdata@cellData@cellExprs)), replace = F)
  dataDR = CYTdata@cellData@cellDimRed[[DimRedName]][subIdx,]
  markers = checkorderMarkers(CYTdata, markers = markers, order = TRUE, checkDuplicates = TRUE)
  dataOR = CYTdata@cellData@cellExprs[subIdx,markers,drop=FALSE]

  resQC = list()

  if (correlation) {
    correlationMethod = match.arg(correlationMethod)
    checkmate::qassert(correlationMethod, "S1")

    ## Compute pairwise distances on the DR and original space for this subsample
    distDR = dataDR %>% dist() %>% as.vector()
    distOrigin = dataOR %>% dist() %>% as.vector()

    res = stats::cor(distDR, distOrigin, method=correlationMethod) ## Report correlation coefficient

    ## Plot distance on the original space vs distance on the DR plot for each pair of point as boxplots
    plot <- data.frame("DR" = distDR, "OR" = cut(distOrigin, 50)) %>%
      ggplot2::ggplot(ggplot2::aes(x=OR, y=DR)) +
      ggplot2::geom_boxplot(fill="slateblue", alpha=0.2) +
      ggplot2::xlab("Distance cuts in original space") + ggplot2::ylab("Distance in reduced space") +
      ggplot2::ggtitle("Pairwise distances in original space vs reduced space") +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title =  ggplot2::element_text(hjust=0.5, size = 20, face = "bold"),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_text(size = 20),
                     axis.title.x = ggplot2::element_text(size = 20))
    cat("\n\n - Pairwise distance", correlationMethod, "correlation between reduced dimension and high dimension(Computed across pairwise distance from",
        round(subsetProportion*nrow(data)), "points) : ", res, "\n")
    resQC[["correlation"]] = list("boxplot" = plot, "cor" = res)
  }

  if (KNN) {
    # finding KNN for each element of data.or and data.dr -> Matrix : nrow(data.or) rows and knn columns
    neighborMatrixOR = RANN::nn2(dataOR, dataOR, knn+1, searchtype = "standard")[[1]][,-1]
    neighborMatrixDR = RANN::nn2(dataDR, dataDR, knn+1, searchtype = "standard")[[1]][,-1]
    res = sapply(1:nrow(dataOR), function(i){ intersect(neighborMatrixOR[i,], neighborMatrixDR[i,]) %>% length() %>% return() }) %>% mean()
    cat("\n\n - Average proportion of nearest neighbours preserved in reduced dimension in comparison to high dimension nearest (Across a", knn,
        "points neighbourhood and with a computation based on a subset of ", subsetSize, "data points ) : ", res/knn)
    resQC[["KNN"]] = list("score" = res/knn)
  }

  if (KNC) {

    if (!kncClusteringName %in% colnames(CYTdata@clusteringData@cellClustering)) {
      stop("Error : kncClusteringName argument (", kncClusteringName, ") is not a clustering column present in cellClustering data frame.")
    }
    dataOR = cbind.data.frame("Clustering" = CYTdata@clusteringData@cellClustering[subIdx,kncClusteringName], dataOR)
    dataDR = cbind.data.frame("Clustering" = CYTdata@clusteringData@cellClustering[subIdx,kncClusteringName], dataDR)
    CmeansOR = dataOR %>%
      group_by(Clustering) %>%
      summarise(across(everything(), mean)) %>%
      remove_rownames() %>%
      column_to_rownames("Clustering")
    CmeansDR = dataDR %>%
      group_by(Clustering) %>%
      summarise(across(everything(), mean)) %>%
      remove_rownames() %>%
      column_to_rownames("Clustering")
    neighborMatrixOR = RANN::nn2(CmeansOR, CmeansOR, knc+1, searchtype = "standard")[[1]][,-1]
    neighborMatrixDR = RANN::nn2(CmeansDR, CmeansDR, knc+1, searchtype = "standard")[[1]][,-1]
    res = sapply(1:nrow(CmeansOR), function(i){ intersect(neighborMatrixOR[i,], neighborMatrixDR[i,]) %>% length() %>% return() }) %>% mean
    cat("\n\n - Average proportion of nearest clusters preserved in reduced dimension in comparison to high dimension nearest (Across a", knc, "clusters neighbourhood and
        with a computation based on a subset of ", subsetSize, "data points ) : ", res/knc)
    resQC[["KNC"]] = list("score" = res/knc)
  }
  return(resQC)
}







