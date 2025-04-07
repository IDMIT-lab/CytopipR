#' @title function for the cyCombine batch correction workflow.
#'
#' @description Compute the batch correction on the data using the ComBat algorithm.
#'  Define a covariate, either as a character vector or name of tibble column.
#'  The covariate should preferable be the cell condition types but can be any column that infers heterogeneity in the data.
#'  The function assumes that the batch information is in the "batch" column and the data contains a "sample" column with sample information.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param batchMetadata a character vector containing the names of metadata associated to batch acquisition (name of a CYTdata's metadata dataframe columns)
#' @param markers a character vector containing the names of biological markers to correct. By default, all Marker are corrected
#' @param normMethod a character being the normalization method. Should be either 'rank', 'scale' or 'qnorm'. Default: 'scale'. "rank" is recommended when combining data with heavy batch effects
#' @param tiesMethod a character being the method to handle ties, when using rank. Default: 'average'
#' @param seed a numeric being the seed to use when creating the SOM.
#' @param xdim a numeric being the x-dimension size of the SOM.
#' @param ydim a numeric being the y-dimension size of the SOM.
#' @param rlen a numeric being the number of times the data is presented to the SOM network. Consider a larger value, if results are not convincing (e.g. 100)
#' @param covarMetadata a character being the covariate ComBat should use. Must be the name of a specific metadata (name of a CYTdata's metadata dataframe columns)
#' @param parametric Default: TRUE. If TRUE, the parametric version of ComBat is used. If FALSE, the non-parametric version is used.
#'
#' @return a S4 object of class 'CYTdata'
#'
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
