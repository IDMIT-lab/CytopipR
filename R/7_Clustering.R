#' @title Rename Clustering Columns or Cluster IDs in CYTdata
#'
#' @description
#' This function allows you to rename clustering column names or the cluster IDs in a `CYTdata` object. Depending on the input, it can either rename the column names of clustering data or rename the individual cluster IDs (levels) within a specified clustering column. Optionally, duplicate cluster levels can be merged if necessary.
#'
#' @param CYTdata A CYTdata object containing the clustering data to be renamed.
#' @param clustering A string specifying the name of the clustering column to rename cluster IDs (optional). If `NULL`, the function will rename the clustering column names instead.
#' @param merge A logical indicating whether duplicated clustering levels should be merged if they share the same name after renaming. Default is `TRUE`.
#' @param from A character vector of clustering names to be renamed (either column names or cluster IDs).
#' @param to A character vector of new names for the clustering columns or cluster IDs.
#'
#' @return A CYTdata object with renamed clustering columns or cluster IDs.
#'
#' @details
#' This function performs two operations depending on the provided parameters:
#' 1. **Renaming clustering columns**: If the `clustering` argument is `NULL`, the function will rename the columns of the `cellClustering` dataframe. The `from` argument should contain the existing column names, and `to` should contain the new column names. The lengths of `from` and `to` must be equal.
#' 2. **Renaming cluster IDs**: If the `clustering` argument is specified, the function renames the cluster IDs within that column. The `from` argument contains the old cluster IDs, and the `to` argument contains the new names for those clusters. If the new names result in duplicate clusters, the `merge` argument controls whether the duplicates should be merged.
#'
#' The function checks for various errors, such as mismatched lengths between `from` and `to`, missing or duplicated names, and invalid clustering names or IDs.
#'
#' @examples
#' # Example 1: Rename clustering columns
#' CYTdata_corrected <- renameClustering(CYTdata, from = c("oldClust1", "oldClust2"), to = c("newClust1", "newClust2"))
#'
#' # Example 2: Rename cluster IDs within a specific clustering column
#' CYTdata_corrected <- renameClustering(CYTdata, clustering = "clusterColumn", from = c("clusterA", "clusterB"), to = c("clusterX", "clusterY"))
#'
#' # Example 3: Merge duplicated cluster IDs after renaming
#' CYTdata_corrected <- renameClustering(CYTdata, clustering = "clusterColumn", from = c("clusterA", "clusterB"), to = c("clusterA", "clusterA"), merge = TRUE)
#'
#' @seealso
#' \code{\link{plyr::mapvalues}} for mapping old values to new values.
#'
#' @import checkmate
#' @import plyr
#' @export

renameClustering <- function(CYTdata, clustering = NULL, merge = TRUE, from = NULL, to){

  CYTdata = checkValidity(CYTdata, mode = "error")
  if (nrow(CYTdata@clusteringData@cellClustering)==0 && ncol(CYTdata@clusteringData@cellClustering)==0) { stop("Error : clusteringData slot is empty") }

  checkmate::qassert(clustering, c(0, "S*"))
  checkmate::qassert(to, "S*")
  checkmate::qassert(merge, "B1")

  if (is.null(clustering)) { # Rename clustering columns
    if (is.null(from)) { from = colnames(CYTdata@clusteringData@cellClustering) }
    else {
      clErr = setdiff(from, colnames(object@clusteringData@cellClustering))
      if (length(clErr)>0) { stop("Error : 'from' argument contains names not present among cellClustering's clustering column names (", paste0(clErr, collapse=", "), ".") }
    }
    if (length(from)!=length(to)) { stop("Error : Length of argument 'from' (", length(from), ") and argument 'to' (", length(to), ") must be equal") }
    dupfrom = unique(from[duplicated(from)])
    if (length(dupfrom)>0) { stop("Error : Argument 'from' contains duplicated clustering column names :", paste0(dupfrom, collapse=", ")) }
    dupto = unique(to[duplicated(to)])
    if (length(dupto)>0) { stop("Error : Argument 'to' contains duplicated names :", paste0(dupto, collapse=", ")) }

    cat("\n\nCurrent clustering column names of cellClustering dataframe are (in the order) :", paste0(colnames(CYTdata@clusteringData@cellClustering), collapse = ", "), "\n")
    cat("\n\nThe following values :")
    cat("\n - ", paste0(from, collapse=", "))
    cat("\n\nwill be renamed, in the order, by :")
    cat("\n - ", paste0(to, collapse=", "))
    colnames(CYTdata@clusteringData@cellClustering)[match(colnames(CYTdata@clusteringData@cellClustering), from)] = to
  }
  else { # Rename clusters names
    if (!clustering %in% colnames(CYTdata@clusteringData@cellClustering)) {
      stop("Error : 'clustering' argument (", clustering, ") is not NULL or a clustering column names")
    }
    oldclulev = levels(CYTdata@clusteringData@cellClustering[,clustering])
    if (is.null(from)) { from = oldclulev }
    else {
      clErr = setdiff(from, oldclulev)
      if (length(clErr)>0) { stop("Error : 'from' argument contains clustering values not present among", clustering, "clustering values (", paste0(clErr, collapse=", "), ".") }
      dupfrom = unique(from[duplicated(from)])
      if (length(dupfrom)>0) { stop("Error : Argument 'from' contains duplicated", clustering, "clustering values :", paste0(dupfrom, collapse=", ")) }
    }
    if (length(from)!=length(to)) { stop("Error : Length of argument 'from' (", length(from), ") and argument 'to' (", length(to), ") must be equal") }
    if (length(to) == 1) { to = rep(to, length(from)) }
    newclulev = oldclulev
    newclulev[match(from, newclulev)] = to

    levdup = unique(newclulev[duplicated(newclulev)])
    if (length(levdup)>0){
      message("After renaming, several ", clustering, "'s cluster IDs have the same name (",  paste0(levdup, collapse=", "), "), either by renaming to an already existing and unchanged clustering value, or by duplicate in the 'to' argument, or both.")
      if (merge) { message("\nWarning : 'merge' argument is set to TRUE. After renaming,", clustering, " column of cellClustering dataframe get its duplicated levels merged.") }
      else { stop("\nError : 'merge' argument is set to FALSE. CYTdata unchanged.") }
    }
    cat("\n\nCurrent", clustering, "'s cluster IDs (levels) are (in the order) :", paste0(oldclulev, collapse = ", "), "\n")
    cat("\n\nThe following values :")
    cat("\n - ", paste0(from, collapse=", "))
    cat("\n\nwill be renamed, in the order, by :")
    cat("\n - ", paste0(to, collapse=", "))
    CYTdata@clusteringData@cellClustering[,clustering] = plyr::mapvalues(CYTdata@clusteringData@cellClustering[,clustering], from = from, to = to)
  }
  CYTdata = checkValidity(CYTdata, mode = "warning")
  return(CYTdata)
}

#' @title Reorder Clustering Levels in CYTdata object
#'
#' @description
#' This function allows you to reorder the levels of a specified clustering column in the `CYTdata` object. The clustering levels can be reordered either alphabetically or based on a custom order provided by the user.
#'
#' @param CYTdata A `CYTdata` object containing the clustering data that will be reordered.
#' @param clustering A string specifying the name of the clustering column whose levels are to be reordered.
#' @param alphabetic A logical indicating whether the clustering levels should be sorted alphabetically. Default is `TRUE`.
#' @param newOrder A character vector specifying a custom order for the clustering levels. This argument is used only if `alphabetic` is set to `FALSE`. If `alphabetic` is `TRUE`, `newOrder` is ignored.
#'
#' @return A `CYTdata` object with the reordered clustering levels.
#'
#' @details
#' This function allows you to reorder the clustering levels of a specified column in the `cellClustering` dataframe. You can reorder the levels alphabetically (using the `alphabetic` argument) or specify a custom order using the `newOrder` argument.
#'
#' If `alphabetic` is set to `TRUE`, the clustering levels will be reordered alphabetically. If `alphabetic` is set to `FALSE`, the `newOrder` argument must be provided, and the function will reorder the clustering levels according to the values in `newOrder`. The `newOrder` argument must contain values that are present in the existing clustering levels and must not contain duplicates.
#'
#' The function performs several checks to ensure the validity of the inputs, including checking if the specified clustering column exists, if the provided custom order is valid, and if the levels contain duplicates.
#'
#' @examples
#' # Example 1: Reorder clustering levels alphabetically
#' CYTdata_corrected <- reorderClustering(CYTdata, clustering = "clusterColumn", alphabetic = TRUE)
#'
#' # Example 2: Reorder clustering levels based on a custom order
#' CYTdata_corrected <- reorderClustering(CYTdata, clustering = "clusterColumn", alphabetic = FALSE, newOrder = c("clusterB", "clusterA", "clusterC"))
#'
#' @seealso
#' \code{\link{gtools::mixedsort}} for sorting the levels alphabetically while preserving the numeric order.
#'
#' @import checkmate
#' @import gtools
#' @export

reorderClustering <- function(CYTdata, clustering, alphabetic = TRUE, newOrder = NULL){

  CYTdata = checkValidity(CYTdata, mode = "error")
  checkmate::qassert(clustering, "S1")
  checkmate::qassert(alphabetic, "B1")
  checkmate::qassert(newOrder, c(0, "S*"))

  if (!clustering %in% colnames(CYTdata@clusteringData@cellClustering)) {
    stop("Error : 'clustering' argument (", clustering, ") is not a clustering column names")
  }

  oldOrder = levels(CYTdata@clusteringData@cellClustering[,clustering])
  if (alphabetic) { newOrder = gtools::mixedsort(oldOrder) }
  else {
    if (is.null(newOrder)) { stop("Error : 'alphabetic' argument is set to FALSE but 'newOrder' argument is NULL.") }

    clErr = setdiff(newOrder, oldOrder)
    if (length(clErr)>0) {
      stop("Error : 'newOrder' argument contains cluster IDs not present among", clustering, "'s cluster IDs (", paste0(clErr, collapse=", "), ".")
    }
    dupnewOrder = unique(newOrder[duplicated(newOrder)])
    if (length(dupnewOrder)>0) { stop("Error : Argument 'newOrder' contains duplicated", clustering, "'s cluster IDs :", paste0(dupnewOrder, collapse=", ")) }
    CYTdata@clusteringData@cellClustering[,clustering] = factor(CYTdata@clusteringData@cellClustering[,clustering], levels = gtools::mixedsort(oldOrder))
  }
  CYTdata@clusteringData@cellClustering[,clustering] = factor(CYTdata@clusteringData@cellClustering[,clustering], levels = newOrder)
  CYTdata = checkValidity(CYTdata, mode = "warning")
  return(CYTdata)
}

#' @title Check Metaclustering Membership for Clusters
#'
#' @description
#' This function checks the membership of clusters in a given metaclustering in the `CYTdata` object. It ensures that each cluster has cells that belong to only one metacluster and returns a data frame showing the metacluster membership for each cluster.
#'
#' @param CYTdata A `CYTdata` object containing the clustering and metaclustering data.
#' @param clustering A string specifying the name of the clustering column in the `cellClustering` data frame.
#' @param metaclustering A string specifying the name of the metaclustering column in the `cellClustering` data frame.
#'
#' @return A data frame with two columns:
#' - `metaclusters`: The metacluster labels corresponding to each cluster.
#' - `clusters`: The cluster labels.
#'
#' @details
#' This function checks if each cluster from the specified `clustering` column belongs exclusively to one metacluster from the `metaclustering` column. If any cluster contains cells from multiple metaclusters, an error is raised. The function then returns a data frame mapping each cluster to its corresponding metacluster.
#'
#' If the input `clustering` and `metaclustering` columns have multiple metacluster memberships for a cluster, the function stops and reports which clusters have cells assigned to distinct metaclusters.
#'
#' @examples
#' # Example: Get metaclustering membership for clusters in CYTdata
#' membership = checkgetMetaclusteringMembership(CYTdata, clustering = "clusterColumn", metaclustering = "metaClusterColumn")
#'
#' @seealso
#' \code{\link{table}} for generating contingency tables, which is used in this function to check memberships.
#'
#' @import checkmate
#' @export

checkgetMetaclusteringMembership <- function(CYTdata,
                                             clustering,
                                             metaclustering) {
  CYTdata = checkValidity(CYTdata, mode = "error")

  if (!clustering %in% colnames(CYTdata@clusteringData@cellClustering)) {
    stop("Error : 'clustering' argument (", clustering, ") is not a clustering column names")
  }
  if (!metaclustering %in% colnames(CYTdata@clusteringData@cellClustering)) {
    stop("Error : 'metaclustering' argument (", metaclustering, ") is not a clustering column names")
  }

  tableMembership = table(CYTdata@clusteringData@cellClustering[,clustering], CYTdata@clusteringData@cellClustering[,metaclustering])

  isHierarchised = apply(tableMembership, 1, function(x){return(sum(x!=0)==1)})
  if (sum(!isHierarchised)>0) {
    stop("Error : Following clusters (from ", clustering, " clustering) have cells belonging to distinct metaclusters (from ", metaclustering, " clustering) :", paste0(names(isHierarchised)[isHierarchised==0], collapse=", "))
  }
  membership = apply(tableMembership, 1, function(x){return(colnames(tableMembership)[x!=0])})
  membership = data.frame("metaclusters" = factor(as.vector(membership), levels = levels(CYTdata@clusteringData@cellClustering[,metaclustering])),
                          "clusters" = factor(names(membership), levels = levels(CYTdata@clusteringData@cellClustering[,clustering])))
  colnames(membership) = c(metaclustering, clustering)
  return(membership)
}




#' @title Compute MSI (Marker or DimRed) for Clusters
#'
#' @description
#' This function computes the MSI (Median/Mean of Scaled Expression) for each cluster in a specified `clustering` column using either marker expression data or dimensionality reduction (DimRed) data. The function performs scaling of the data and then computes either the median or mean for each cluster's MSI.
#'
#' @param CYTdata A `CYTdata` object containing the clustering and expression data.
#' @param clustering A string specifying the name of the clustering column in the `cellClustering` data frame.
#' @param clusters A vector of cluster labels (optional). If `NULL`, all clusters from `clustering` will be used.
#' @param samples A vector of sample identifiers (optional). If `NULL`, all samples will be used.
#' @param computeWith A string specifying whether to compute MSI based on "markers" or "DimRed". Defaults to "markers".
#' @param markers A vector of marker names (optional). Used only when `computeWith = "markers"`.
#' @param DimRed A string specifying the name of the dimensionality reduction method (optional). Used only when `computeWith = "DimRed"`.
#' @param typeScaling A string specifying the type of scaling to apply. Options include:
#'   - `"none"`: No scaling
#'   - `"center"`: Center the data
#'   - `"reduced"`: Apply reduced scaling
#'   - `"center_reduced"`: Apply center and reduced scaling
#'   - `"rescale_min_max"`: Rescale data to min-max scale
#'   Defaults to `"none"`.
#' @param rescaleQt A numeric vector of length 2 specifying the quantiles for rescaling. Defaults to `c(0.05, 0.95)`.
#' @param rescaleBorders A numeric vector of length 2 specifying the borders for rescaling. Defaults to `c(0, 1)`.
#' @param scaleByMetadata A string specifying a sample metadata column to scale by. Defaults to `NULL`, indicating no scaling by metadata.
#' @param typeMSI A string specifying whether to compute MSI using the `"mean"` or `"median"` method. Defaults to `"median"`.
#'
#' @return A data frame with MSI values for each cluster, with rows representing clusters and columns representing features.
#'
#' @details
#' The function computes the MSI (either mean or median) for each cluster based on the scaling of either marker expression data or dimensionality reduction (DimRed) data. If scaling by metadata is specified, the data is scaled separately for each metadata group before computing the MSI.
#'
#' The MSI values are computed for each cluster using the specified `clustering` column. Scaling options allow for different preprocessing steps to ensure that the data is appropriately transformed before computing MSI.
#'
#' @examples
#' # Example: Compute MSI for clusters using marker data
#' msi_values = computeMSI(CYTdata, clustering = "clusterColumn", computeWith = "markers", markers = c("marker1", "marker2"))
#'
#' # Example: Compute MSI for clusters using DimRed data
#' msi_values = computeMSI(CYTdata, clustering = "clusterColumn", computeWith = "DimRed", DimRed = "PCA")
#'
#' @seealso
#' \code{\link{exprsScaling}} for scaling the expression data, which is used in this function to preprocess data.
#'
#' @import checkmate
#' @export

computeMSI <- function(CYTdata,
                       clustering,
                       clusters = NULL,
                       samples = NULL,

                       computeWith = c("markers", "DimRed"),
                       markers = NULL,
                       DimRed = NULL,

                       typeScaling = c("none", "center", "reduced", "center_reduced", "rescale_min_max"),
                       rescaleQt = c(0.05, 0.95), rescaleBorders = c(0,1),
                       scaleByMetadata = NULL,
                       typeMSI = c("median", "mean")) {


  CYTdata = checkValidity(CYTdata, mode = "error")

  typeMSI = match.arg(typeMSI)
  checkmate::qassert(typeMSI, "S1")

  checkmate::qassert(samples, c(0, "S*"))
  samples = checkorderSamples(CYTdata, samples = samples, order = TRUE, checkDuplicates = TRUE)
  subIdx = CYTdata@sampleData@cellSample$SampleID %in% samples

  checkmate::qassert(clustering, "S1")
  if (!clustering %in% colnames(CYTdata@clusteringData@cellClustering)) {
    stop("Error : 'clustering' argument (", clustering, ") is not a clustering column name.")
  }
  clusters = tryCatch(
    { checkorderClustering(CYTdata, clustering = clustering, clusters = clusters, order = TRUE, checkDuplicates = TRUE) },
    error = function(e) {
      message("Error : 'clustering' argument is set to", clustering, "but, for this clustering, checkorderClustering function returns error with 'clusters' argument (",
              paste0(clusters, collapse=", "), ").\n checkorderClustering's error : ", e$message)
      return(NULL)
    })
  subIdx = subIdx & (CYTdata@clusteringData@cellClustering[,clustering] %in% clusters)

  computeWith = match.arg(computeWith)
  checkmate::qassert(computeWith, "S1")
  if (computeWith == "markers") {
    markers = tryCatch(
      { checkorderMarkers(CYTdata, markers, order=TRUE, checkDuplicates=TRUE) },
      error = function(e) {
        message("Error : 'computeWith' argument is set to 'markers' but checkorderMarkers function returns error with 'markers' argument (",
                paste0(markers, collapse=", "), ").\n checkorderMarkers's error : ", e$message)
        return(NULL)
      }
    )
    data = CYTdata@cellData@cellExprs[subIdx, markers, drop=FALSE]
  } else {
    if (is.null(DimRed) || (!DimRed %in% names(CYTdata@cellData@cellDimRed))) {
      stop("Error : 'computeWith' argument is set to 'DimRed' but 'DimRed' argument (", DimRed, ") is NULL or not a DimRed's list name")
    }
    data = CYTdata@cellData@cellDimRed[[DimRed]][subIdx,]
  }

  if (is.null(scaleByMetadata)) {
    data = exprsScaling(data, typeScaling = typeScaling, rescaleQt = rescaleQt, rescaleBorders = rescaleBorders)
  } else {
    if (!scaleByMetadata %in% colnames(CYTdata@sampleData@sampleMetadata)){
      stop("Error : 'scaleByMetadata' argument is not a sampleData slot's sampleMetadata column name.")
    }
    scaleByMetadataID = CYTdata@sampleData@sampleMetadata[CYTdata@sampleData@cellSample[rownames(data), "SampleID"], scaleByMetadata]
    SBMdata = data.frame()
    for (label in unique(scaleByMetadataID)) {
      SBMdataLab = data %>% subset(scaleByMetadataID == label) %>% exprsScaling(typeScaling = typeScaling, rescaleQt = rescaleQt, rescaleBorders = rescaleBorders)
      SBMdata = rbind.data.frame(SBMdata, SBMdataLab)
    }
    data = SBMdata
  }

  medoid = if (typeMSI == "median") median else mean
  data = data %>%
    cbind.data.frame("clustering" = CYTdata@clusteringData@cellClustering[rownames(data), clustering]) %>%
    group_by(clustering) %>%
    summarise(across(everything(), medoid)) %>% remove_rownames() %>% column_to_rownames("clustering")
  return(data)
}

#' @title Perform Clustering on CYTdata
#'
#' @description
#' This function performs clustering on the given CYTdata object using a variety of clustering algorithms. It allows clustering based on either markers or dimensionality reduction (DimRed) data and supports multiple clustering methods such as SOM, FlowSOM, KMeans, DBSCAN, and more. The function also provides options to preprocess the data using scaling and MSI computation before performing clustering.
#'
#' @param CYTdata A `CYTdata` object containing the data to be clustered, including cell expression data, clustering data, and sample data.
#' @param name A string specifying the name to assign to the resulting clustering.
#' @param addPrefix A string prefix to prepend to each cluster label. Defaults to `"C"`.
#' @param clusterWith A string specifying whether to perform clustering based on "markers" or "DimRed" (dimensionality reduction). Defaults to "markers".
#' @param markers A vector of marker names used for clustering when `clusterWith = "markers"`. This is optional.
#' @param DimRed A string specifying the name of the dimensionality reduction method used for clustering when `clusterWith = "DimRed"`. This is optional.
#' @param clusterBy A string specifying whether to cluster by "cell" or by existing "clustering". Defaults to "cell".
#' @param clustering A string specifying the clustering column to use when `clusterBy = "clustering"`. This is optional.
#' @param typeScaling A string specifying the type of scaling to apply to the data before clustering. Options include:
#'   - `"none"`: No scaling
#'   - `"center"`: Center the data
#'   - `"reduced"`: Apply reduced scaling
#'   - `"center_reduced"`: Apply center and reduced scaling
#'   - `"rescale_min_max"`: Rescale data to min-max scale
#'   Defaults to `"none"`.
#' @param rescaleQt A numeric vector of length 2 specifying the quantiles for rescaling. Defaults to `c(0.05, 0.95)`.
#' @param rescaleBorders A numeric vector of length 2 specifying the borders for rescaling. Defaults to `c(0, 1)`.
#' @param scaleByMetadata A string specifying a sample metadata column to scale by. Defaults to `NULL`, indicating no scaling by metadata.
#' @param typeMSI A string specifying whether to compute MSI using the `"mean"` or `"median"` method. Defaults to `"median"`.
#' @param type A string specifying the clustering algorithm to use. Options include:
#'   - `"SOM"`
#'   - `"FlowSOM"`
#'   - `"Phenograph"`
#'   - `"Hierarchical"`
#'   - `"ConsensusHierarchical"`
#'   - `"Spade"`
#'   - `"DBSCAN"`
#'   - `"Kmeans"`
#'   - `"Kmedoids"`
#'   - `"CLARA"`
#'   - `"PhenoGMM"`
#'   - `"Mean_shift"`
#'   - `"flowMeans"`
#'   Defaults to `"SOM"`.
#' @param checkOverwrite A boolean indicating whether to check if the clustering already exists. If `TRUE`, will ask the user whether to overwrite the clustering. Defaults to `TRUE`.
#' @param seed A numeric value to set the random seed for reproducibility. Defaults to 42.
#' @param ... Additional arguments passed to the clustering algorithms. These arguments are passed to the specific clustering functions (e.g., `kohonen::som()`, `FlowSOM::FlowSOM()`, `stats::kmeans()`, etc.).
#'   - **For SOM (`"SOM"`)**, additional arguments for the function `kohonen::som` include:
#'     - `grid`: A `som.grid` object specifying the grid structure (e.g., a 2D grid).
#'     - `rlen`: The number of learning iterations.
#'     - `toroidal`: A boolean to specify whether the map should be toroidal.
#'     - `alpha`: The learning rate.
#'     - `radius`: The radius of influence of the neighboring units in the map.
#'     - **Documentation**: [kohonen::som](https://cran.r-project.org/web/packages/kohonen/kohonen.pdf)
#'
#'   - **For FlowSOM (`"FlowSOM"`)**, additional arguments for the function `FlowSOM::FlowSOM` include:
#'     - `colsToUse`: A vector of column names to use for clustering.
#'     - `xdim`, `ydim`: The dimensions of the grid (e.g., number of rows and columns).
#'     - `maxNumberOfClusters`: The maximum number of clusters to generate.
#'     - `seed`: The seed for random number generation.
#'     - `normalize`: Logical value indicating if normalization should be applied before clustering.
#'     - **Documentation**: [FlowSOM::FlowSOM](https://cran.r-project.org/web/packages/FlowSOM/FlowSOM.pdf)
#'
#'   - **For K-means (`"Kmeans"`)**, additional arguments for the function `stats::kmeans` include:
#'     - `centers`: The number of clusters or initial centroids.
#'     - `nstart`: The number of random starts to use.
#'     - `iter.max`: The maximum number of iterations.
#'     - `algorithm`: The algorithm to use (e.g., `"Lloyd"`, `"Forgy"`).
#'     - **Documentation**: [stats::kmeans](https://www.rdocumentation.org/packages/stats)
#'
#'   - **For DBSCAN (`"DBSCAN"`)**, additional arguments for the function `dbscan::dbscan` include:
#'     - `eps`: The maximum distance between two points for them to be considered as in the same neighborhood.
#'     - `minPts`: The minimum number of points required to form a dense region.
#'     - **Documentation**: [dbscan::dbscan](https://cran.r-project.org/web/packages/dbscan/dbscan.pdf)
#'
#'   - **For Phenograph (`"Phenograph"`)**, additional arguments for the function `cytofkit2::Rphenograph` include:
#'     - `k`: The number of nearest neighbors to use.
#'     - `ncpus`: The number of CPUs to use for parallel processing.
#'     - **Documentation**: [cytofkit2::Rphenograph](https://cran.r-project.org/web/packages/cytofkit2/cytofkit2.pdf)
#'
#'   - **For CLARA (`"CLARA"`)**, additional arguments for the function `cluster::clara` include:
#'     - `k`: The number of clusters to form.
#'     - `samples`: The number of samples to randomly select from the data (used for faster clustering).
#'     - `metric`: The distance measure used to calculate distances between observations.
#'     - **Documentation**: [cluster::clara](https://cran.r-project.org/web/packages/cluster/cluster.pdf)
#'
#'   - **For Kmedoids (`"Kmedoids"`)**, additional arguments for the function `cluster::pam` (used for K-medoids) include:
#'     - `k`: The number of clusters to form.
#'     - `metric`: The distance metric to use.
#'     - `clustering`: A logical value indicating whether to compute the clustering.
#'     - **Documentation**: [cluster::pam](https://cran.r-project.org/web/packages/cluster/cluster.pdf)
#'
#' @return A list containing:
#'   - `"CYTdata"`: The updated `CYTdata` object with the new clustering added.
#'   - `"opts"`: A list of options and objects related to the clustering process.
#'
#' @details
#' This function performs clustering on the data contained in the `CYTdata` object. The user can choose to cluster based on either marker expression or dimensionality reduction data. It supports a variety of clustering algorithms and provides options for scaling and MSI computation. The function also checks if a clustering with the specified name already exists and prompts the user to overwrite it if necessary.
#'
#' Clustering algorithms available include:
#' - SOM (Self-Organizing Maps)
#' - FlowSOM
#' - Phenograph
#' - Hierarchical
#' - DBSCAN (Density-Based Spatial Clustering)
#' - K-means
#' - K-medoids (PAM)
#' - CLARA (Clustering Large Applications)
#' - and others.
#'
#' @examples
#' # Example: Perform clustering using K-means on markers data
#' result = runClustering(CYTdata, name = "Kmeans_clustering", type = "Kmeans", clusterWith = "markers", markers = c("marker1", "marker2"))
#'
#' # Example: Perform clustering using DBSCAN on DimRed data
#' result = runClustering(CYTdata, name = "DBSCAN_clustering", type = "DBSCAN", clusterWith = "DimRed", DimRed = "PCA")
#'
#' @seealso
#' \code{\link{computeMSI}} for computing MSI before clustering.
#'
#' @import checkmate
#' @export

runClustering <- function(CYTdata,
                          name,
                          addPrefix = "C",

                          clusterWith = c("markers", "DimRed"),
                          markers = NULL,
                          DimRed = NULL,
                          clusterBy = c("cell", "clustering"),
                          clustering = NULL,

                          typeScaling = c("none", "center", "reduced", "center_reduced", "rescale_min_max"),
                          rescaleQt = c(0.05, 0.95), rescaleBorders = c(0,1),
                          scaleByMetadata = NULL,
                          typeMSI = c("median", "mean"),

                          # Clustering params
                          type = c("SOM", "FlowSOM", "Phenograph", "Hierarchical", "ConsensusHierarchical", "Spade", "DBSCAN", "Kmeans", "Kmedoids", "CLARA", "PhenoGMM", "Mean_shift", "flowMeans"),
                          checkOverwrite = TRUE,
                          seed = 42,
                          ...) {

  CYTdata = checkValidity(CYTdata, mode = "error")

  checkmate::qassert(seed, "N1")
  type = match.arg(type)
  checkmate::qassert(type, "S1")

  if (checkOverwrite && (name %in% colnames(CYTdata@clusteringData@cellClustering))) {
    reply <- readline(prompt=paste("A clustering named", name, "is already performed, do you still want to continue and overwrite (yes or no): "))
    while (!tolower(reply) %in% c("yes", "y", "no", "n")) { reply <- readline(prompt="Reply not standard, please answer yes or no: ") }
    if (tolower(reply) %in% c("no", "n")){ stop("Function stopped, CYTdata object unchanged") }
  }
  checkmate::qassert(addPrefix, c(0, "S1"))

  if (clusterBy == "cell") {
    if (clusterWith == "markers") {
      markers = tryCatch(
        { checkorderMarkers(CYTdata, markers, order=TRUE, checkDuplicates=TRUE) },
        error = function(e) {
          message("Error : 'expression' argument is set to 'markers' but checkorderMarkers function returns error with markers' argument (",
                  paste0(markers, collapse=", "), ").\n checkorderMarkers's error : ", e$message)
          return(NULL)
        })
      data = CYTdata@cellData@cellExprs[,markers,drop=FALSE]
      cat("\nClustering is performed, using ", type, " on cell's markers expression (", paste0(markers, collapse=", "), ")")
    } else {
      if (is.null(DimRed) || (!DimRed %in% names(CYTdata@cellData@cellDimRed))) {
        stop("Error : 'clusterWith' argument is set to 'DimRed' but 'DimRed' argument (", DimRed, ") is NULL or not a DimRed list name")
      }
      data = CYTdata@cellData@cellDimRed[[DimRed]]
      cat("\nClustering is performed, using ", type, " on cell's reduced expression (DimRed name : ", DimRed, ")")
    }
  } else {
    data = computeMSI(CYTdata,
                      samples = NULL, clustering = clustering, clusters = NULL,
                      computeWith = clusterWith, markers = markers, DimRed = DimRed,
                      typeScaling = typeScaling, rescaleQt = rescaleQt, rescaleBorders = rescaleBorders, scaleByMetadata = scaleByMetadata, typeMSI = typeMSI)

    if (clusterWith == "markers") {
      cat("\nClustering is performed, using ", type, " on '", clustering, "' cluster's", typeMSI, "markers expression (", paste0(markers, collapse=", "), ")")
    } else {
      cat("\nClustering is performed, using ", type, " on '", clustering, "' cluster's", typeMSI, "reduced expression (DimRed name : ", DimRed, ")")
    }
  }

  print(data)

  switch (type,
          SOM = {
            t = system.time(somObject <- kohonen::som(as.matrix(data), ...))
            clusters <- somObject$unit.classif
            opts = list(somObject)
          },
          FlowSOM = {
            t = system.time(fSOMObject <- FlowSOM::FlowSOM(as.matrix(data), ...))
            clusters <- FlowSOM::GetClusters(fSOMObject)
            metaclusters <- FlowSOM::GetMetaclusters(fSOMObject)
            opts = list(fSOMObject)
          },
          Phenograph = {
            t = system.time(phenographObject <- cytofkit2::Rphenograph(data, ...))# phenograph_homemade(data, ...))
            clusters = phenographObject$membership
            opts = list(phenographObject)
          },
          Louvain = {
            stop()
          },
          Leiden = {
            stop()
          },
          Hierarchical = {
            t = system.time(hierarchicalObject <- doHierarchical(data, ...))
            clusters = hierarchicalObject$ClustersBelong
            opts = list(hierarchicalObject)
          },
          ConsensusHierarchical = {
            stop()
          },
          Spade = {
            t = system.time(spadeObject <- spade_homemade(data, ...))
            clusters = spadeObject$clustering
            opts = list(spadeObject)
          },
          DBSCAN = {
            t = system.time(clusters <- dbscan::dbscan(data, ...)$cluster)
          },
          Kmeans = {
            t = system.time(clusters <- stats::kmeans(data, ...)$cluster)
          },
          Kmedoids = {
            t = system.time(clusters <- cluster::pam(data, ...)$clustering)
          },
          clara = {
            t = system.time(clusters <- cluster::clara(data, ...)$clustering)
          },
          PhenoGMM = {
            stop()
          },
          Mean_shift = {
            stop()
          },
          flowMeans = {
            stop()
          })

  cat("\n", type, "completed in", round(t[3], 3) ,"seconds.\n")
  if( length(clusters) != nrow(data) ){
    stop("Error : Clustering is not complete. Some events are not assigned. Please try other cluster method(s).")
  }

  cat("Updating CYTdata's clusteringData slot.. \n")

  clusters = as.vector(clusters)
  if (!is.null(addPrefix)) { clusters = paste(addPrefix, clusters, sep="") }

  if (type == "FlowSOM") {
    if (clusterBy == "clustering") {
      clusters = structure(clusters, names = rownames(data))
      clusters = as.vector(clusters[CYTdata@clusteringData@cellClustering[,clustering]])
    }
    nameCL = paste(name, "FlowSOM_SOMclusters", sep="_")
    CYTdata@clusteringData@cellClustering[,nameCL] = factor(clusters, levels = gtools::mixedsort(unique(clusters)))
    CYTdata@clusteringData@clusteringPalette[[nameCL]] = structure(rainbow(nlevels(CYTdata@clusteringData@cellClustering[,nameCL])),
                                                                             names = levels(CYTdata@clusteringData@cellClustering[,nameCL]))


    metaclusters = as.vector(metaclusters)
    if (!is.null(addPrefix)) { metaclusters = paste(paste("Meta", addPrefix, sep=""), metaclusters, sep="") }

    if (clusterBy == "clustering") {
      metaclusters = structure(metaclusters, names = rownames(data))
      metaclusters = as.vector(metaclusters[CYTdata@clusteringData@cellClustering[,clustering]])
    }
    nameMC = paste(name, "FlowSOM_metaclusters", sep="_")
    CYTdata@clusteringData@cellClustering[,nameMC] = factor(metaclusters, levels = gtools::mixedsort(unique(metaclusters)))
    CYTdata@clusteringData@clusteringPalette[[nameMC]] = structure(rainbow(nlevels(CYTdata@clusteringData@cellClustering[,nameMC])),
                                                                             names = levels(CYTdata@clusteringData@cellClustering[,nameMC]))
  } else {
    if (clusterBy == "clustering") {
      clusters = structure(clusters, names = rownames(data))
      clusters2 = clusters
      clusters = as.vector(clusters[CYTdata@clusteringData@cellClustering[,clustering]])
    }
    CYTdata@clusteringData@cellClustering[,name] = factor(clusters, levels = gtools::mixedsort(unique(clusters)))
    CYTdata@clusteringData@clusteringPalette[[name]] = structure(rainbow(nlevels(CYTdata@clusteringData@cellClustering[,name])),
                                                                           names = levels(CYTdata@clusteringData@cellClustering[,name]))

    if (type=="Hierarchical" && clusterBy == "clustering") {
      dfClusters = data.frame("metaclusters" = as.vector(clusters2) %>% factor(levels = levels(CYTdata@clusteringData@cellClustering[,name])),
                              "clusters" = rownames(data)) %>% remove_rownames() %>% column_to_rownames("clusters")
      colnames(dfClusters) = name
      opts$HM = ComplexHeatmap::Heatmap(t(data),
                              cluster_columns = as.dendrogram(hierarchicalObject$hclustEvents),
                              column_split = nlevels(CYTdata@clusteringData@cellClustering[,name]),
                              cluster_rows =  as.dendrogram(hierarchicalObject$hclustFeatures),
                              top_annotation = ComplexHeatmap::HeatmapAnnotation(df = dfClusters, col = CYTdata@clusteringData@clusteringPalette))
    }
  }

  CYTdata = checkValidity(CYTdata, mode = "warning")
  return(list("CYTdata" = CYTdata, "opts" = opts))
}


#' @title Add Manual Clustering to CYTdata
#'
#' @description
#' This function allows the user to manually add a new clustering to the CYTdata object. The new clustering is stored in the \code{cellClustering} slot, and an optional color palette can be provided for visualization. The function also checks if the specified clustering name already exists and prompts the user for confirmation if overwrite is allowed.
#'
#' @param CYTdata An object of class \code{CYTdata}. It should contain clustering data in the \code{cellClustering} slot.
#' @param name A character string specifying the name of the new clustering to be added. This will be used as the column name in the \code{cellClustering} slot.
#' @param clusteringIDs A vector of clustering identifiers (IDs) to assign to cells. The length of this vector should match the number of cells in the \code{CYTdata} object.
#' @param additionalPalette (optional) A vector specifying custom colors for the new clustering. The number of colors should match the number of unique clustering IDs. If \code{NULL}, the default color palette will be applied using \code{rainbow}.
#' @param checkOverwrite A logical value indicating whether to prompt the user for confirmation if a clustering with the same name already exists in the \code{cellClustering} slot. Default is \code{TRUE}.
#'
#' @return An updated \code{CYTdata} object with the new clustering added to the \code{cellClustering} slot. The associated color palette for the clustering is also stored in the \code{clusteringPalette} list.
#'
#' @details
#' The function checks if the specified clustering name already exists in the \code{cellClustering} slot of the \code{CYTdata} object. If the name exists and \code{checkOverwrite} is \code{TRUE}, the user will be asked if they want to overwrite the existing clustering. If the user confirms, the existing clustering will be replaced with the new one.
#'
#' If no custom palette is provided, the function will generate a default color palette using the \code{rainbow} function.
#'
#' @examples
#' # Example usage
#' CYTdata <- addManualClustering(CYTdata = CYTdata,
#'                                 name = "newClustering",
#'                                 clusteringIDs = c(1, 2, 1, 3, 2, 1),
#'                                 checkOverwrite = TRUE)
#'
#' # Using custom color palette
#' customPalette <- c("red", "green", "blue")
#' CYTdata <- addManualClustering(CYTdata = CYTdata,
#'                                 name = "customClustering",
#'                                 clusteringIDs = c(1, 2, 1, 3, 2, 1),
#'                                 additionalPalette = customPalette,
#'                                 checkOverwrite = TRUE)
#'
#' @seealso \code{\link{checkValidity}}, \code{\link{qassert}}
#'
#' @import checkmate
#' @importFrom utils readline
#'
#' @export

addManualClustering <- function(CYTdata,
                                name,
                                clusteringIDs,
                                additionalPalette = NULL,
                                checkOverwrite = TRUE) {

  CYTdata = checkValidity(CYTdata, mode = "error")

  if (checkOverwrite && (name %in% colnames(CYTdata@clusteringData@cellClustering))) {
    reply <- readline(prompt=paste("A clustering named", name, "is already performed, do you still want to continue and overwrite (yes or no): "))
    while (!tolower(reply) %in% c("yes", "y", "no", "n")) { reply <- readline(prompt="Reply not standard, please answer yes or no: ") }
    if (tolower(reply) %in% c("no", "n")){ stop("Function stopped, CYTdata object unchanged") }
  }

  checkmate::qassert(clusteringIDs, "F*")
  CYTdata@clusteringData@cellClustering[,name] = clusteringIDs
  if (is.null(additionalPalette)) {
    CYTdata@clusteringData@clusteringPalette[[name]] = structure(rainbow(nlevels(CYTdata@clusteringData@cellClustering[,name])),
                                                                 names = levels(CYTdata@clusteringData@cellClustering[,name]))
  } else {
    CYTdata@clusteringData@clusteringPalette[[name]] = additionalPalette
  }

  CYTdata = checkValidity(CYTdata, mode = "warning")
  return(CYTdata)
}

#' @title Plot Population Counts by Clustering and Metadata
#'
#' @description
#' This function generates a stacked bar plot showing the number of cells per cluster, with colors representing the proportion of a metadata variable. It calculates the number of cells for each cluster, optionally groups them by metadata (e.g., sample groups), and visualizes the results.
#'
#' @param CYTdata An object of class \code{CYTdata}. This should contain cell clustering data in the \code{cellClustering} slot and sample metadata in the \code{sampleMetadata} slot.
#' @param clustering A character string specifying the name of the clustering to use. This should match a column name in the \code{cellClustering} slot.
#' @param metadata A character string specifying the metadata column name from the \code{sampleMetadata} slot that will be used for coloring the plot.
#' @param clusters A vector of cluster identifiers to be included in the plot. If \code{NULL}, all clusters are considered.
#' @param sort A logical value indicating whether the clusters should be sorted by their total cell count. Default is \code{TRUE}.
#' @param samples A vector of sample identifiers to filter the data. If \code{NULL}, all samples are considered.
#'
#' @return A \code{ggplot} object representing the stacked bar plot.
#'
#' @details
#' The function counts the number of cells per cluster for each sample. It then merges this count with the provided metadata and generates a stacked bar plot showing the count of cells for each cluster, colored by the metadata values. If the \code{sort} parameter is set to \code{TRUE}, the clusters are ordered by their total cell count in descending order.
#'
#' The plot provides a visual representation of the population distribution of cells within each cluster and how it relates to the metadata variable.
#'
#' @examples
#' # Example usage
#' plot <- plotPopulationCounts(CYTdata = CYTdata,
#'                              clustering = "SampleClustering",
#'                              metadata = "TreatmentGroup",
#'                              sort = TRUE)
#' print(plot)
#'
#' @import dplyr
#' @import ggplot2
#' @import viridis
#' @import reshape2
#'
#' @export

plotPopulationCounts <- function(CYTdata,
                                 clustering,
                                 metadata,
                                 clusters = NULL,
                                 sort = TRUE,
                                 samples = NULL) {

  CYTdata = checkValidity(CYTdata, mode = "error")

  cellcount = getCellCount(CYTdata = CYTdata,
                           clustering = clustering, clusters = clusters,
                           samples = samples,
                           type = "cellcount")


  checkmate::qassert(sort, "B1")
  checkmate::qassert(metadata, "S1")
  if (!metadata %in% colnames(CYTdata@sampleData@sampleMetadata)) {
    stop("Error : 'metadata' argument is not a sampleMetadata column name (available sampleMetadata column name/metadata : ",
         paste0(colnames(CYTdata@sampleData@sampleMetadata), collapse=","), ")")
  }

  medoid = if (typeMSI == "median") median else mean
  data = data %>%
    cbind.data.frame("clustering" = CYTdata@clusteringData@cellClustering[rownames(data), clustering]) %>%
    group_by(clustering) %>%
    summarise(across(everything(), medoid)) %>% remove_rownames() %>% column_to_rownames("clustering")


  sumCount = cellcount %>%
    t() %>%
    merge(CYTdata@sampleData@sampleMetadata, by = "row.names") %>%
    group_by(across(metadata)) %>%
    summarise(across(rownames(cellcount)), sum) %>%
    reshape2::melt(id = metadata)
  colnames(sumCount) = c("md", "clusters", "count")

  if (sort) { sumCount$clusters = factor(sumCount$clusters, levels = rownames(cellcount)[order(apply(cellcount,1,sum), decreasing=TRUE)]) }
  else { sumCount$clusters = factor(sumCount$clusters, levels = rownames(cellcount)) }

  plot <- ggplot2::ggplot(data = sumCount,
                          ggplot2::aes(x = clusters,
                                       y = count,
                                       fill = md)) +
    ggplot2::geom_bar(stat="identity", position = "stack", color="black") +
    viridis::scale_fill_viridis(option = "turbo", discrete = TRUE) +
    ggplot2::ggtitle(paste(clustering, "headcount and", metadata, "proportion")) +
    ggplot2::xlab(clustering) + ggplot2::ylab("Number of cells") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                   legend.key = ggplot2::element_blank())

  return(plot)
}

######################################################## Utils ########################################################

updateClusteringLevels <- function(CYTdata, andPalette = TRUE) {
  mdFrame = CYTdata@clusteringData@cellClustering %>% mutate(across(where(is.factor), droplevels))
  palList = CYTdata@clusteringData@clusteringPalette
  for (md in colnames(mdFrame)) { palList[[md]] = palList[[md]][levels(mdFrame[,md])] }
  CYTdata@clusteringData@cellClustering = mdFrame
  CYTdata@clusteringData@clusteringPalette = palList[colnames(mdFrame)]
  return(CYTdata)
}

checkorderClustering <- function(CYTdata, clustering, clusters = NULL, order = TRUE, checkDuplicates = TRUE){

  #CYTdata = checkValidity(CYTdata, mode = "error")
  checkmate::qassert(clustering, "S1")
  checkmate::qassert(clusters, c(0, "S*"))
  checkmate::qassert(order, "B1")
  checkmate::qassert(checkDuplicates, "B1")

  if (nrow(CYTdata@clusteringData@cellClustering)==0 && ncol(CYTdata@clusteringData@cellClustering)==0) { stop("Error : clusteringData slot is empty") }

  if (!clustering %in% colnames(CYTdata@clusteringData@cellClustering)) {
    stop("Error : 'clustering' argument is not a clustering column in cellClustering dataframe (clustering columns :", paste0(colnames(CYTdata@clusteringData@cellClustering), collapse = ", "), ").")
  }

  clusterID = levels(CYTdata@clusteringData@cellClustering[,clustering])
  if (is.null(clusters)) { clusters = clusterID }

  if (checkDuplicates) {
    clustersdup = clusters[duplicated(clusters)]
    if (length(clustersdup)>0) { stop("Error : 'clusters' argument contain duplicated values ( ", paste0(clustersdup, collapse = ", "), " ). It must be vector of unique clusters of" , clustering, "column.") }
  }
  splErr = setdiff(clusters, clusterID)
  if (length(splErr)>0) { stop("Error in CYTdata object, clusteringData slot: 'clusters' argument providing clusters", paste0(splErr, collapse=", "), "not present in", clustering, "column of cellClustering dataframe") }
  if (order) { clusters = clusterID[clusterID %in% clusters] }
  return(clusters)
}

exprsScaling <- function(data,
                         typeScaling = c("none", "center", "reduced", "center_reduced", "rescale_min_max"),
                         rescaleQt = c(0.05, 0.95), rescaleBorders = c(0,1)) {
  typeScaling = match.arg(typeScaling)
  checkmate::qassert(typeScaling, "S1")
  checkmate::qassert(rescaleBorders, "N2")
  checkmate::qassert(rescaleQt, "N2")
  if (any(rescaleQt<=0) || any(rescaleQt>=1)) {
    stop("Error : rescaleQt argument must be a vector a 2 integer between 0 and 1.")
  }

  sca <- function(exp) {
    switch(typeScaling,
           none = {exp = exp},
           center = {exp = scale(exp, center=T, scale=F)},
           reduced = {exp = scale(exp, center=F, scale=T)},
           center_reduced = {exp = scale(exp, center=T, scale=T)},
           rescale_min_max = {
             quantiles = stats::quantile(exp, probs = rescaleQt)
             low = quantiles[1]
             high = quantiles[2]
             exp[exp<low] = low
             exp[exp>high] = high
             exp = scales::rescale(exp, from = c(low, high), to = rescaleBorders)
           })
    return(exp)
  }
  data = data %>% mutate_all(sca)
  return(data)
}

doPhenograph <- function(data) {
  stop()
}

doHierarchical <- function(data,
                           N_Hierarchical,
                           criterion = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")) {
  checkmate::qassert(N_Hierarchical, "N1")
  if (N_Hierarchical>=nrow(data) | N_Hierarchical<2) {
    stop("Error : N_Hierarchical argument of doHierarchical function is set to", N_Hierarchical,
         " but must be an integer between 2 and the number of events to cluster (", nrow(data), ").")
  }
  criterion = match.arg(criterion)
  checkmate::qassert(criterion, "S1")

  ### Hclust for events (clusters, cells)
  dataDist = stats::dist(data)
  dataDist[dataDist==0] = 0.000001
  hclustEvents = stats::hclust(dataDist, method = criterion)
  isoMDSord = MASS::isoMDS(dataDist, k = 1)$points
  hclustEvents = dendextend::rotate(hclustEvents, rownames(isoMDSord)[order(isoMDSord)])
  labelEvents = hclustEvents$labels[hclustEvents$order]

  cutreeClusters = stats::cutree(hclustEvents, k = N_Hierarchical)
  cutreeClusters = structure(match(cutreeClusters, unique(cutreeClusters[hclustEvents$labels[hclustEvents$order]])), names = names(cutreeClusters))
  ClustersBelong = cutreeClusters[rownames(data)]

  ### Hclust for features
  dataDist = stats::dist(t(data))
  dataDist[dataDist==0] = 0.000001
  hclustFeatures = stats::hclust(dataDist, method = criterion)

  return(list("ClustersBelong" = ClustersBelong, "hclustEvents" = hclustEvents, "hclustFeatures" = hclustFeatures))
}

