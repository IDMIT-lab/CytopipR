#' @title Perform Batch Correction on CYTdata Object
#'
#' @description
#' This function applies batch correction to a `CYTdata` object, using the specified batch and covariate metadata. The batch correction is performed by utilizing the `cyCombine::batch_correct` function, and it can handle multiple types of normalizations and tie-breaking methods. The function also allows for the inclusion of optional clustering and marker information.
#'
#' @param CYTdata A CYTdata object containing the data to be batch-corrected.
#' @param batchMetadata A string specifying the name of the column in `sampleMetadata` for batch information.
#' @param covarMetadata An optional string specifying the name of the column in `sampleMetadata` for covariate information.
#' @param samples A character vector of sample identifiers to subset the data for batch correction.
#' @param clustering A string specifying the name of the column containing the clustering information (optional).
#' @param clusters A vector of clusters to consider for batch correction (optional).
#' @param markers A character vector of marker names for which the batch correction is applied.
#' @param normMethod A string specifying the normalization method. Options are `"scale"`, `"rank"`, and `"qnorm"`. Default is `"scale"`.
#' @param tiesMethod A string specifying the method for handling tied ranks. Options are `"average"`, `"first"`, `"last"`, `"random"`, `"max"`, and `"min"`. Default is `"average"`.
#' @param seed A numeric value for the random seed, to ensure reproducibility of the batch correction process. Default is `92`.
#' @param xdim The number of dimensions for the output (x-axis).
#' @param ydim The number of dimensions for the output (y-axis).
#' @param rlen The number of iterations for the batch correction.
#' @param parametric A logical indicating whether to use parametric methods for batch correction. Default is `TRUE`.
#'
#' @return A CYTdata object with the batch-corrected data.
#'
#' @details
#' This function performs batch correction on the `CYTdata` object using the specified batch and optional covariate metadata. The correction method can be customized using the `normMethod` and `tiesMethod` parameters. The function also supports subsetting the data to specific samples or clusters, and the corrected expression data is returned in the CYTdata object. The `cyCombine::batch_correct` function is used for the core batch correction process.
#'
#' @examples
#' # Example 1: Perform batch correction using the "batch" metadata
#' CYTdata_corrected <- runcyCombine(CYTdata = CYTdata, batchMetadata = "batch", markers = c("marker1", "marker2"))
#'
#' # Example 2: Perform batch correction with parametric method and custom seed
#' CYTdata_corrected <- runcyCombine(CYTdata = CYTdata, batchMetadata = "batch", markers = c("marker1", "marker2"), seed = 100)
#'
#' # Example 3: Perform batch correction for specific samples and clusters
#' CYTdata_corrected <- runcyCombine(CYTdata = CYTdata, batchMetadata = "batch", samples = c("sample1", "sample2"), clustering = "cluster", clusters = c("cluster1", "cluster2"))
#'
#' @seealso
#' \code{\link{batch_correct}} from the \pkg{cyCombine} package for batch correction functionality.
#'
#' @import checkmate
#' @import cyCombine
#' @import dplyr
#' @export

runcyCombine <- function(CYTdata,
                         batchMetadata,
                         covarMetadata = NULL,
                         samples = NULL,
                         clustering = NULL, clusters = NULL,
                         markers = NULL,
                         normMethod = c("scale","rank","qnorm"),
                         tiesMethod = c("average", "first", "last", "random", "max", "min"),
                         seed = 92,
                         xdim, ydim, rlen,

                         parametric = TRUE){

  CYTdata = checkValidity(CYTdata, mode = "error")

  checkmate::qassert(samples, c(0, "S*"))
  samples = checkorderSamples(CYTdata, samples = samples, order = TRUE, checkDuplicates = TRUE)
  subIdx = CYTdata@sampleData@cellSample$SampleID %in% samples

  if (!is.null(clusters)) {
    checkmate::qassert(clustering, c(0, "S1"))
    if (is.null(clustering) || (!clustering %in% colnames(CYTdata@clusteringData@cellClustering))) {
      stop("Error : 'clusters' argument is not null but 'clustering' argument (", clustering, ") is NULL or not a clustering column name.")
    }
    clusters = tryCatch(
      { checkorderClustering(CYTdata, clustering = clustering, clusters = clusters, order = TRUE, checkDuplicates = TRUE) },
      error = function(e) {
        message("Error : 'clustering' argument is set to", clustering, "but, for this clustering, checkorderClustering function returns error with 'clusters' argument (",
                paste0(clusters, collapse=", "), ").\n checkorderClustering's error : ", e$message)
        return(NULL)
      })
    subIdx = subIdx & (CYTdata@clusteringData@cellClustering[,clustering] %in% clusters)
  }

  checkmate::qassert(markers, c(0,"S*"))
  markers = checkorderMarkers(CYTdata, markers, order=TRUE, checkDuplicates=TRUE)

  normMethod = match.arg(normMethod)
  checkmate::qassert(normMethod, "S1")
  tiesMethod = match.arg(tiesMethod)
  checkmate::qassert(tiesMethod, "S1")
  checkmate::qassert(xdim, "N1")
  checkmate::qassert(ydim, "N1")
  checkmate::qassert(rlen, "N1")
  checkmate::qassert(seed, "N1")

  checkmate::qassert(parametric, "B1")

  cat("\n\nPreparing Data for batch correction..")

  if (nrow(CYTdata@sampleData@sampleMetadata)==0){
    stop("Error in CYTdata object, sampleMetadata dataframe from sampleData slot is empty but necessary.")
  }
  checkmate::qassert(batchMetadata, "S1")
  if (!batchMetadata %in% colnames(CYTdata@sampleData@sampleMetadata)) { stop() }

  CYTdata@sampleData@sampleMetadata[CYTdata@sampleData@cellSample$SampleID,batchMetadata]
  data = cbind.data.frame(
    CYTdata@cellData@cellExprs,
    "sample" = CYTdata@sampleData@cellSample$SampleID,
    "batch" = CYTdata@sampleData@sampleMetadata[CYTdata@sampleData@cellSample$SampleID,batchMetadata],
    "id" = 1:nrow(CYTdata@cellData@cellExprs)
  )

  checkmate::qassert(covarMetadata, c("0", "S1"))
  if(!is.null(covarMetadata)){
    if (!covarMetadata %in% colnames(CYTdata@sampleData@sampleMetadata)) { stop() }
    if (batchMetadata == covarMetadata){ stop("Error : 'batchMetadata' and 'covarMetadata' arguments are the same metadata. Must be different.") }
    covarMetadata = CYTdata@sampleData@sampleMetadata[CYTdata@sampleData@cellSample$SampleID,covarMetadata]
  }

  cat("\n\nPerforming batch correction..")

  print(data[subIdx,])

  ### Correct data
  data.corrected = cyCombine::batch_correct(data[subIdx,],
                                            markers = markers,
                                            covar = covarMetadata,
                                            parametric = parametric,
                                            norm_method = normMethod,
                                            ties.method = tiesMethod,
                                            xdim = xdim, ydim = ydim, rlen = rlen,
                                            seed = seed)

  cat("\n\nCreating new CYTdata object ..")
  CYTdata@cellData@cellExprs[subIdx, markers] = data.corrected[, markers]
  CYTdata = checkValidity(CYTdata, mode = "warning")
  return(CYTdata)
}
