
# # @title downsampling.density
# # @description Downsample indexes representing sample's events with a density-based method (spade)
# #
# # @param indexes.to.downsample a numeric vector specifying the indexes of a sample's events in the dataset
# # @param target.size a numeric value specifying the number of events targeted as downsampled, is smaller than the number of indexes to downsample
# # @param whole.data.exprs a data.frame containing the entire dataset : all the marker's expression (used for density computation)
# # @param exclude.pctile a numeric value, between 0 and 1, specifying the quantile value used to determine the exclusion threshold for density values
# # @param compensate.uniform a boolean value : Density based downsampling doesn't permit us to get the exact size specified by "target.size". So, do we add a step to adjust randomly the donsampled indexes ?
# #
# # @return downsampled.indexes a numeric vector specifying the indexes downsampled
#
# downsampling.density <- function(indexes.to.downsample,
#                                  target.size,
#                                  whole.data.exprs,
#                                  exclude.pctile,
#                                  compensate.uniform,
#                                  ncores){
#   cat("- Exact size to dowsnample : ", target.size)
#
#   data.exprs = whole.data.exprs[indexes.to.downsample,]
#
#   ### SPADE density
#   df.dists = parallelDist::parDist(x = as.matrix(data.exprs), method = "euclidean", threads = ncores) %>% as.matrix() %>% as.data.frame()
#   idx.used.df.dists = sample(1:nrow(df.dists), 3*nrow(df.dists)/4)
#   med.dists = apply(df.dists[idx.used.df.dists,], 1, min, na.rm=TRUE) %>% as.vector() %>% median()
#   kernel.width = 20*med.dists
#   boof.df.dists = (df.dists < kernel.width) %>% as.data.frame()
#   density = apply(boof.df.dists, 1, sum, na.rm=TRUE) %>% as.vector()
#   #################################
#
#   data = cbind.data.frame("density" = density, "indexes" = indexes.to.downsample,  data.exprs)
#
#   print(data)
#   exclusion.boundary = stats::quantile(density, exclude.pctile, names=FALSE)
#
#   ## First step of DS : Spade exclusion
#   data = subset(data, density > exclusion.boundary) # Selection by exclude.pctile
#
#   print(data)
#
#   n.compensate = target.size - nrow(data)
#   if (n.compensate > 0){
#     cat("- Number of cells, after spade exclusion : ", nrow(data), ". So no more density downsampling. Using uniform sampling to compensate and sample the exact size..")
#     remaining.idx = indexes.to.downsample[!indexes.to.downsample %in% data$indexes]
#     compensate.idx = sample(remaining.idx, n.compensate)
#     downsampled.indexes = c(data$indexes, compensate.idx)
#   }
#   else {
#     ## Second step of DS : Target size and boundary estimation
#     density.sorted = sort(data$density)
#     cdf = rev(cumsum(1/rev(density.sorted)))
#     target.boundary = target.size/cdf[1]
#     if (target.boundary > density.sorted[1]) {
#       targets = (target.size-seq(1,length(density.sorted))) / cdf
#       target.boundary = targets[which.min(targets-density.sorted > 0)]
#     }
#     data = subset(data, target.boundary/density > stats::runif(nrow(data))) # Selection by ~target.size
#     cat("- For a sample, the number of cells, after the finished downsampling (density) process, is : ", nrow(data))
#     n.compensate = target.size - nrow(data)
#     if (n.compensate > 0){
#       cat("- Using uniform sampling to compensate and sample the exact size..")
#       remaining.idx = indexes.to.downsample[!indexes.to.downsample %in% data$indexes]
#       compensate.idx = sample(remaining.idx, n.compensate)
#       downsampled.indexes = c(data$indexes, compensate.idx)
#     }
#     else if (n.compensate < 0){
#       cat("- Using uniform sampling to compensate and sample the exact size..")
#       n.compensate = abs(n.compensate)
#       compensate.idx = sample(data$indexes, n.compensate)
#       downsampled.indexes = data$indexes[!data$indexes %in% compensate.idx]
#     }
#   }
#   return(downsampled.indexes)
# }
# @param type a string value specifying the type of downsampling selection used ("none" if no downsampling is performed, "uniform" uniformly-based, "density" density-based)
# @param parallel.ncores a numeric value specifying the number of logical processor used for parallel computation
# @param density.exclusion.pctile a numeric value, between 0 and 1, specifying, the quantile value used to determine the exclusion threshold for density values (expected if "density" selected)
# @param densityCompensate a boolean value : Density based downsampling doesn't permit us to get the exact size specified by "target.size". So, do we add a step to adjust randomly the donsampled indexes ? (expected if "density" selected)
# type = c("uniform", "density"),
# densityPercentile = 0.01,
# densityCompensate = FALSE,
# densityNcores = 1

# switch(type,
#        uniform = {
#          downsampling.function <- function(indexes, size) {
#            return(sample(indexes, size)) # replace = FALSE by default
#          }
#        },
#        density = {
#          checkmate::qassert(density.exclusion.pctile, "N1")
#          checkmate::qassert(densityCompensate, "B1")
#          checkmate::qassert(density.ncores, "N1")
#          list.density.parameters = list("density.exclusion.pctile" = density.exclusion.pctile,
#                                         "densityCompensate" = densityCompensate,
#                                         "density.ncores" = density.ncores)
#          cat("\n\n - The downsampling density parameters are :",
#              paste0(names(list.density.parameters), "=", list.density.parameters, collapse = ", "))
#          downsampling.function <- function(indexes, size) {
#            return(downsampling.density(indexes.to.downsample = indexes,
#                                        target.size = size,
#                                        whole.data.exprs = CYTdata@matrix.expression,
#                                        exclude.pctile = density.exclusion.pctile,
#                                        compensate.uniform = densityCompensate,
#                                        ncores = density.ncores) )}
#        })

#' @title Performs the downsampling of events
#' @description Perform the downsampling, using uniformly-based or density-based random selections, of events
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param type a string value specifying the way the number of downsampled events is determined ("by.sample" to take the smallest size of sample, "by.number" to take an absolute size, "by.percent" to take a part of each sample)
#' @param absolute a numeric value specifying the exact number of downsampled events (expected if "by.number" selected)
#' @param percentage a numeric value specifying the percent of each sample size taken as number of downsampled events (expected if "by.percent" selected)
#' @param seed a numeric value : the seed for RNG
#'
#' @return a S4 object of class 'CYTdata'
#'
#' @export
#'

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
