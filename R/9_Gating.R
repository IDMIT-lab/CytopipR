#' @title Gate cells from selected population stored in a CYTdata object
#'
#' @description Gate cells from selected population stored in a CYTdata object and create new CYTdata object.
#' The clustering identifiers (metaclusters, supermetaclusters also if computed)
#' Dimension reduction coordinates are not kept
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param population a character vector containing the population identifiers to select and gate. By default, all the populations are gated
#' @param levels a string value being the type of population to gate. Possible values are: "metaclusters", "clusters", "supermetaclusters" (default = metaclusters)
#'
#' @return a S4 object of class 'CYTdata'
#'
#' @export
#'

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
