#' @title Add Manual Dimensionality Reduction
#'
#' @description
#' This function allows users to add a manually computed dimensionality reduction (DimRed) to the CYTdata object. The provided `DimRedDF` (data frame) is stored in the `cellDimRed` slot under the specified `DimRedName`. This function ensures that if the given `DimRedName` already exists, the user is prompted for confirmation before overwriting.
#'
#' @param CYTdata An object of class \code{CYTdata}. It should contain dimensionality reduction data in the \code{DimRedData} slot.
#' @param DimRedName A character string specifying the name for the new dimensionality reduction. This will be used as the key to store the dimensionality reduction in the `cellDimRed` slot.
#' @param DimRedDF A data frame containing the result of the dimensionality reduction. It should have rows corresponding to cells and columns representing reduced dimensions.
#' @param checkOverwrite A logical value indicating whether to check for the existence of the `DimRedName` and prompt the user for confirmation before overwriting an existing entry. Default is \code{TRUE}.
#'
#' @return The updated \code{CYTdata} object with the new dimensionality reduction added.
#'
#' @details
#' This function is useful for manually adding a result of a dimensionality reduction (e.g., PCA, t-SNE, UMAP) into the `CYTdata` object. If the specified `DimRedName` already exists in the `cellDimRed` slot, the user will be asked whether to overwrite it. This ensures that important data is not unintentionally lost.
#'
#' @examples
#' # Example usage:
#' CYTdata <- addManualDimRed(CYTdata = CYTdata,
#'                            DimRedName = "PCA_Results",
#'                            DimRedDF = pca_result_df)
#'
#' @import checkmate
#'
#' @export

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


#' @title Run Dimensionality Reduction
#'
#' @description
#' This function performs dimensionality reduction on the provided cytometry data using one of the selected methods (UMAP, tSNE, LVish, or PCA). The function stores the results in the `cellDimRed` slot of the `CYTdata` object. It allows for the application of different dimensionality reduction algorithms to explore and visualize high-dimensional cytometry data in a lower-dimensional space (typically 2D or 3D).
#'
#' @param CYTdata An object of class \code{CYTdata}. It contains the cytometry data and related metadata.
#' @param name A character string specifying the name of the new dimensionality reduction result. This will be used to store the result in the \code{cellDimRed} slot.
#' @param markers A character vector specifying the marker names to be used for dimensionality reduction. If \code{NULL}, all available markers in the data will be used.
#' @param type A character string specifying the type of dimensionality reduction method to be applied. It can be one of: \code{"UMAP"}, \code{"tSNE"}, \code{"lvish"}, or \code{"PCA"}.
#' @param seed An integer value to set the seed for random number generation. Default is \code{42}.
#' @param checkOverwrite A logical value indicating whether to check for the existence of the specified \code{name} in the \code{cellDimRed} slot and prompt the user for confirmation before overwriting. Default is \code{TRUE}.
#' @param \dots Additional parameters passed to the respective dimensionality reduction functions (e.g., \code{uwot::umap}, \code{Rtsne::Rtsne}, or \code{stats::prcomp}).
#'
#' @return The updated \code{CYTdata} object, with the new dimensionality reduction result stored in the \code{cellDimRed} slot.
#'
#' @details
#' This function supports the following dimensionality reduction methods:
#' - \code{UMAP}: Uniform Manifold Approximation and Projection.
#' - \code{tSNE}: t-Distributed Stochastic Neighbor Embedding.
#' - \code{lvish}: LVish manifold learning method (a non-linear method).
#' - \code{PCA}: Principal Component Analysis.
#'
#' The results of the selected dimensionality reduction method are stored in the \code{cellDimRed} slot of the \code{CYTdata} object under the specified \code{name}. This can then be used for visualization or further analysis.
#'
#' @examples
#' # Example usage:
#' CYTdata <- runDimRed(CYTdata = CYTdata,
#'                      name = "UMAP_result",
#'                      markers = c("CD3", "CD8", "CD45"),
#'                      type = "UMAP")
#'
#' # Example with PCA
#' CYTdata <- runDimRed(CYTdata = CYTdata,
#'                      name = "PCA_result",
#'                      markers = c("CD3", "CD4", "CD56"),
#'                      type = "PCA",
#'                      seed = 123)
#'
#' @import uwot
#' @import Rtsne
#' @import checkmate
#' @import stats
#'
#' @export

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
  library(uwot)

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
