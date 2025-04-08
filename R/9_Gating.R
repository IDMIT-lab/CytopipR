#' @title Perform Gating on CYTdata Object
#'
#' @description
#' This function performs gating on a given cytometry dataset (CYTdata) by filtering cells based on specified sample IDs and clustering results. The function allows users to subset the data by selecting specific samples and/or clusters. It updates the `CYTdata` object to reflect the new, gated dataset, retaining the remaining data in the object.
#'
#' @param CYTdata An object of class \code{CYTdata} containing the cytometry data, cell metadata, clustering information, and other related attributes.
#' @param samples A character vector of sample IDs that should be included in the gating process. If \code{NULL}, clustering filtering will be applied (if provided).
#' @param clustering A character string specifying the clustering column name in the \code{cellClustering} slot to filter cells by their cluster assignment. Must be specified if \code{clusters} is provided.
#' @param clusters A vector of cluster identifiers that define which clusters should be kept in the gated data. If \code{NULL}, no clustering filtering is applied.
#'
#' @return The updated \code{CYTdata} object with the gated cells, preserving the metadata, clustering, and dimensionality reduction information after the gating process.
#'
#' @details
#' The gating process subsets the cytometry data to only include the specified samples and/or clusters. If the \code{samples} argument is provided, cells belonging to the specified sample IDs are retained. If the \code{clusters} argument is provided, only cells from the specified clusters are included. The function updates the \code{CYTdata} object and ensures that clustering, dimensionality reduction, and metadata are also adjusted to reflect the gated dataset.
#'
#' If both \code{samples} and \code{clustering} are specified, the cells will be filtered based on both criteria.
#'
#' The function also ensures that the clustering levels, metadata levels, and other related elements are appropriately updated after the gating operation.
#'
#' @examples
#' # Example usage:
#' # Gating based on specific sample IDs
#' CYTdata <- CYTdataGating(CYTdata = CYTdata, samples = c("Sample1", "Sample2"))
#'
#' # Gating based on specific clusters in a clustering
#' CYTdata <- CYTdataGating(CYTdata = CYTdata,
#'                          clustering = "cluster1",
#'                          clusters = c("ClusterA", "ClusterB"))
#'
#' @import checkmate
#' @import dplyr
#' @importFrom methods droplevels
#'
#' @export

CYTdataGating <- function(CYTdata,
                          samples = NULL,
                          clustering = NULL,
                          clusters = NULL){

  CYTdata = checkValidity(CYTdata, mode = "error")

  if (is.null(samples) && is.null(clusters)) {
    stop()
  }

  checkmate::qassert(samples, c(0, "S*"))
  samples = checkorderSamples(CYTdata, samples = samples, order = TRUE, checkDuplicates = TRUE)
  gatedIdx = CYTdata@sampleData@cellSample$SampleID %in% samples

  if (!is.null(clusters)) {
    checkmate::qassert(clustering, c(0, "S1"))
    if (is.null(clustering) || (!clustering %in% colnames(CYTdata@clusteringData@cellClustering))) {
      stop("Error : 'clusters' argument is not NULL but 'clustering' argument (", clustering, ") is NULL or not a clustering column name.")
    }
    clusters = tryCatch(
      { checkorderClustering(CYTdata, clustering = clustering, clusters = clusters, order = TRUE, checkDuplicates = TRUE) },
      error = function(e) {
        message("Error : 'clustering' argument is set to", clustering, "but, for this clustering, checkorderClustering function returns error with 'clusters' argument (",
                paste0(clusters, collapse=", "), ").\n checkorderClustering's error : ", e$message)
        return(NULL)
      })
    gatedIdx = gatedIdx & (CYTdata@clusteringData@cellClustering[,clustering] %in% clusters)
  }

  #   # cat("\nGating CYTdata object and creating one object containing cells belonging to the following", level, " :",
  #   paste0(subset, collapse=", "))
  # newmatrix.expression = subset(CYTdata@matrix.expression, gatedIdx)
  # cat("\n\n - Number of cells gated :", sum(gatedIdx))
  #
  # newsamples = CYTdata@samples[gatedIdx]
  # newsamples = droplevels(newsamples)
  # removedSpls = setdiff(levels(CYTdata@samples), levels(newsamples))
  # if(length(removedSpls)>0){
  #   cat("\n\n - Gating operation remove following samples :", paste0(removedSpls, collapse = ", "))
  # }
  #  message("\n\nRemark : Clustering results are preserved during gating operation but it is recommended to the user to run a new clustering step with parameters adapted to the gated dataset")

  cat("\n\nUpdating CYTdata object..\n")
  CYTdata@cellData@cellExprs = CYTdata@cellData@cellExprs[gatedIdx,,drop=FALSE]
  CYTdata@cellData@cellAdditionalexprs = CYTdata@cellData@cellAdditionalexprs[gatedIdx,,drop=FALSE]

  CYTdata@sampleData@cellSample = CYTdata@sampleData@cellSample[gatedIdx,,drop=FALSE]
  CYTdata@sampleData@cellSample$SampleID = droplevels(CYTdata@sampleData@cellSample$SampleID)

  CYTdata@sampleData@sampleMetadata = CYTdata@sampleData@sampleMetadata[levels(CYTdata@sampleData@cellSample$SampleID),,drop=FALSE]
  CYTdata = updateMetadataLevels(CYTdata)

  CYTdata@clusteringData@cellClustering =  CYTdata@clusteringData@cellClustering[gatedIdx,,drop=FALSE]
  CYTdata = updateClusteringLevels(CYTdata)

  for (i in seq_along(CYTdata@cellData@cellDimRed)) {
    CYTdata@cellData@cellDimRed[[i]] = CYTdata@cellData@cellDimRed[[i]][gatedIdx,,drop=FALSE]
  }

  CYTdata = checkValidity(CYTdata, mode = "warning")
  return(CYTdata)
}
