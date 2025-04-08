#' @title Downsample Data in CYTdata Object
#'
#' @description
#' This function performs downsampling on a CYTdata object. The downsampling can be done in various ways:
#' - **absolute**: Downsampling the entire dataset to a fixed number of cells.
#' - **percent**: Downsampling the entire dataset by a given percentage of cells.
#' - **smallSample**: Downsampling each sample to the size of the smallest sample.
#' - **percentSample**: Downsampling each sample independently by a given percentage.
#' - **absoluteSample**: Downsampling each sample to an absolute number of cells.
#'
#' @param CYTdata A CYTdata object containing the dataset to be downsampled.
#' @param type A character string specifying the downsampling type. Valid options are:
#'   - `"absolute"`: Downsample the entire dataset to a fixed number of cells.
#'   - `"percent"`: Downsample the entire dataset by a percentage.
#'   - `"smallSample"`: Downsample each sample to the size of the smallest sample.
#'   - `"percentSample"`: Downsample each sample independently by a given percentage.
#'   - `"absoluteSample"`: Downsample each sample to a fixed number of cells.
#' @param absolute A positive integer specifying the number of cells for downsampling (used in `"absolute"`, `"absoluteSample"`).
#' @param percentage A numeric value between 0 and 1 representing the proportion of the total number of cells to retain (used in `"percent"` and `"percentSample"`).
#' @param seed A numeric value for the random seed, to ensure reproducibility of the downsampling process.
#'
#' @return A CYTdata object containing the downsampled data.
#'
#' @details
#' This function provides several methods for downsampling the data in the `CYTdata` object:
#' - If `type = "absolute"`, it will select a fixed number of cells from the entire dataset.
#' - If `type = "percent"`, it will select a fixed percentage of cells from the entire dataset.
#' - If `type = "smallSample"`, each sample will be downsampled to the size of the smallest sample.
#' - If `type = "percentSample"`, each sample will be downsampled by the specified percentage.
#' - If `type = "absoluteSample"`, each sample will be downsampled to the specified absolute number of cells.
#'
#' The function updates the following components of the CYTdata object:
#' - `cellExprs`
#' - `cellAdditionalexprs`
#' - `cellSample`
#' - `sampleMetadata`
#' - `cellClustering`
#' - `cellDimRed`
#'
#' @examples
#' # Example 1: Downsample the dataset to 50% of cells
#' CYTdata_downsampled <- downsampleData(CYTdata = CYTdata, type = "percent", percentage = 0.5)
#'
#' # Example 2: Downsample each sample to 1000 cells
#' CYTdata_downsampled <- downsampleData(CYTdata = CYTdata, type = "absoluteSample", absolute = 1000)
#'
#' # Example 3: Downsample the dataset by 30% of total cells
#' CYTdata_downsampled <- downsampleData(CYTdata = CYTdata, type = "percent", percentage = 0.3)
#'
#' # Example 4: Downsample each sample to the smallest sample size
#' CYTdata_downsampled <- downsampleData(CYTdata = CYTdata, type = "smallSample")
#'
#' @seealso
#' \code{\link{checkValidity}} for validating the CYTdata object before downsampling.
#'
#' @import checkmate
#' @import dplyr
#' @import plyr
#' @export

downsampleData <- function(CYTdata,
                           type = c("absolute", "percent", "smallSample", "percentSample", "absoluteSample"),
                           absolute = 1000,
                           percentage = 1/5,
                           seed = 42){

  CYTdata = checkValidity(CYTdata, mode = "error")

  checkmate::qassert(seed, "N1")
  set.seed(seed)
  type = match.arg(type)
  checkmate::qassert(type, "S1")
  cat("\nDownsampling will be performed uniformly.")

  switch(type,
         percentSample = {
           checkmate::qassert(percentage, "N1")
           if (!(percentage>=0 && percentage<=1)) { stop("Error : 'percentage' argument is a proportion and must be a positive numeric between 0 and 1.") }
           cat("\n'type' argument set to 'percentSample' : Downsampling will be performed independently on each samples.
               The resulting subsampled object will contain a different amount of cells by samples, equal to ",
               percentage*100, "% of sample size.")
           sizes = round(table(CYTdata@sampleData@cellSample$SampleID)*percentage)
           names(sizes) = levels(CYTdata@sampleData@cellSample$SampleID)
         },
         percent = {
           checkmate::qassert(percentage, "N1")
           if (!(percentage>=0 && percentage<=1)) { stop("Error : 'percentage' argument is a proportion and must be a positive numeric between 0 and 1.") }
           cat("\n'type' argument set to 'percent' : Downsampling will be performed on the whole dataset without distinguishing each sample
               The resulting subsampled object will contain", percentage*100, "% of total amount of cells.")
           sizes = round(length(CYTdata@sampleData@cellSample$SampleID)*percentage)
         },
         smallSample = {
           cat("\n'type' argument set to 'smallSample' : Downsampling will be performed equally across the samples.
               The resulting subsampled object will contain the same amount of cells in each samples equal to the smallest sample's size")
           sizes = rep(min(table(CYTdata@sampleData@cellSample$SampleID)), nlevels(CYTdata@sampleData@cellSample$SampleID))
           names(sizes) = unique(CYTdata@sampleData@cellSample$SampleID)
         },
         absoluteSample = {
           checkmate::qassert(absolute, "N1")
           if (absolute<0) { stop("Error : 'absolute' argument must be a positive integer.") }
           cat("\n'type' argument set to 'absoluteSample' : Downsampling will be performed equally across the samples.
               The resulting subsampled object will contain the same amount of cells in each samples equal to ", absolute, " cells.")
           sizes = rep(absolute, nlevels(CYTdata@sampleData@cellSample$SampleID))
           names(sizes) = unique(CYTdata@sampleData@cellSample$SampleID) # un
         },
         absolute = {
           checkmate::qassert(absolute, "N1")
           if (absolute<0) { stop("Error : 'absolute' argument must be a positive integer.") }
           cat("\n'type' argument set to 'absolute' : Downsampling will be performed on the whole dataset without distinguishing each sample.
               The resulting subsampled object will contain ", absolute, " cells.")
           sizes = absolute
         })

  if (length(sizes)==1) { DownsamplingIndexes = sample(rownames(CYTdata@sampleData@cellSample), size = sizes) }
  else {
    DownsamplingIndexes = c()
    for (spl in names(sizes)) {
      size = sizes[spl]
      splindexes = CYTdata@sampleData@cellSample %>% filter(SampleID == spl) %>% rownames()
      if (length(splindexes) <= size){
        cat("\n - No downsampling for sample ", spl, "(enough events already or sample size smaller than downsampling target)")
        DSindexes = splindexes
      }
      else {
        cat("\n - ", size, " events to downsample, over ", length(splindexes), ", for sample ", spl)
        DSindexes = sample(splindexes, size = size)
      }
      DownsamplingIndexes = c(DownsamplingIndexes, DSindexes)
      cat("\n - Downsampling done for sample ", spl)
    }
  }

  cat("\n\nUpdating CYTdata object..\n")
  CYTdata@cellData@cellExprs = CYTdata@cellData@cellExprs[DownsamplingIndexes,]
  CYTdata@cellData@cellAdditionalexprs = CYTdata@cellData@cellAdditionalexprs[DownsamplingIndexes,,drop=FALSE]

  CYTdata@sampleData@cellSample = CYTdata@sampleData@cellSample[DownsamplingIndexes,,drop=FALSE]
  CYTdata@sampleData@cellSample$SampleID = droplevels(CYTdata@sampleData@cellSample$SampleID)
  CYTdata@sampleData@sampleMetadata = CYTdata@sampleData@sampleMetadata[levels(CYTdata@sampleData@cellSample$SampleID),,drop=FALSE]
  CYTdata = updateMetadataLevels(CYTdata)


  CYTdata@clusteringData@cellClustering =  CYTdata@clusteringData@cellClustering[DownsamplingIndexes,]
  CYTdata = updateClusteringLevels(CYTdata)

  for (i in seq_along(CYTdata@cellData@cellDimRed)) {
    CYTdata@cellData@cellDimRed[[i]] = CYTdata@cellData@cellDimRed[[i]][DownsamplingIndexes,,drop=FALSE]
  }

  CYTdata = checkValidity(CYTdata, mode = "warning")
  return(CYTdata)
}


