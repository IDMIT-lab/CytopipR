

#' @title Performs the upsampling of downsampled events
#'
#' @description This function aims to perform the upsampling of downsampled events events based on an existing CYTdata object and existing cell events stored in tab-separated or FCS files.
#'
#' Importantly, the identification of cell clusters must have been performed prior to this operation.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param newCYTdata a S4 object of class 'CYTdata'
#' @param markers a character vector providing the markers to use to perform upsampling
#' @param type a character value containing the type of the method to apply. Possible values are: "KNN", "RF", "Logistic.regression" (default = KNN)
#'
#' @return a S4 object of class 'CYTdata'
#'
#' @export
#'

performUpsampling <- function(CYTdata,
                              newCYTdata,
                              markers = NULL,
                              type = c("KNN", "RF", "Logistic.regression")){

  # Check args
  checkmate::qassert(markers, c(0,"S*"))
  type = match.arg(type)
  checkmate::qassert(type, "S1")

  if (is.null(markers)){ markers = CYTdata@markers }
  markers = as.vector(markers)
  markers_interest = intersect(markers, newCYTdata@markers)
  if (length(markers_interest)<length(markers)) {
    stop("Error in performUpsampling : Some of the markers to use for Upsampling are not present in newCYTdata (",
         paste0(markers[!markers %in% markers_interest], collapse = ", "),")\n")
  }
  cat("Common markers used for Upsampling :", paste0(markers_interest, collapse=", "))

  exprs = CYTdata@matrix.expression
  new.exprs = newCYTdata@matrix.expression[, newCYTdata@markers %in% CYTdata@markers]

  raw.exprs = CYTdata@raw.matrix.expression
  new.raw.exprs = newCYTdata@raw.matrix.expression[, newCYTdata@markers %in% CYTdata@markers]

  overall.exprs = rbind.data.frame(exprs, new.exprs)
  overall.raw.exprs = rbind.data.frame(raw.exprs, new.raw.exprs)
  overall.clusters = c(CYTdata@Clustering@clusters, rep(NA, nrow(new.exprs)))
  overall.samples = c(CYTdata@samples, newCYTdata@samples)

  indexes.duplicated = duplicated(overall.exprs)

  overall.exprs = overall.exprs[indexes.duplicated,]
  overall.raw.exprs = overall.raw.exprs[indexes.duplicated,]
  overall.clusters = overall.clusters[indexes.duplicated]
  overall.samples = overall.samples[indexes.duplicated]

  data.clusters = cbind.data.frame(overall.exprs, "clusters" = overall.clusters)
  training.dataset = subset(data.clusters, !is.na(clusters))
  test.dataset = subset(data.clusters, is.na(clusters))
  test.dataset$clusters = NULL

  switch(type,
         KNN = {
           centroids = plyr::ddply(training.dataset, "clusters", function(x) {
             x$clusters <- NULL
             centers <- apply(x, 2, stats::median, rm.na = TRUE)
             return(centers)
           })
           upsampled.clusters <- FNN::knnx.index(centroids, test.dataset, k = 1, algorithm = "kd_tree")
         },
         RF = { stop("Not coded yet") },
         Logistic.regression = { stop("Not coded yet") })


  overall.clusters = c(CYTdata@Clustering@clusters, upsampled.clusters)

  # Put into "Clustering" object
  cat("\n\nCreating new Clustering object :")
  cat("\n - Computing cell cluster count matrix...")
  cellcount = compute.cell.count(clusters = overall.clusters, samples = overall.samples)
  cat("\n - Computing cell cluster abundance matrix...")
  abundance = compute.cell.abundance(count.matrix = cellcount)

  Clustering.object <- methods::new("Clustering",
                                    clusters = overall.clusters,
                                    cellcount = cellcount,
                                    abundance = abundance,
                                    parameters = append(CYTdata@Clustering@parameters, list("type.upsampling" = type,
                                                                                            "markers.upsampling" = markers)))
  cat("Creating new CYTdata object.. \n")
  CYTdata <- methods::new("CYTdata",
                          matrix.expression = overall.exprs,
                          markers = colnames(exprs),
                          samples = samples,
                          files = files,
                          raw.matrix.expression = overall.raw.exprs,
                          Clustering = Clustering.object)
  CYTdata@Clustering = Clustering.object
  validObject(CYTdata)
  return(CYTdata)

  CYTdata@matrix.expression.r <- matrix.expression.r[, colnames(matrix.expression.r) %in% colnames(downsampled.exp)]
  CYTdata@samples <- c(CYTdata@samples, upsampled.samples)
  CYTdata@identify.clusters <- c(CYTdata@identify.clusters, knn)

  message("computing cell cluster count matrix...")
  CYTdata@matrix.cell.count <- computeCellCounts(proj = CYTdata@matrix.expression,
                                                 clusters = CYTdata@identify.clusters,
                                                 samples = CYTdata@samples)

  message("computing cell cluster abundance matrix...")
  CYTdata@matrix.abundance <- computeClusterAbundances(count = CYTdata@matrix.cell.count)

  validObject(CYTdata)
  return(CYTdata)
}

