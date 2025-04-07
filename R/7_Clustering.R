#' @title Renames markers within a CYTdata object
#'
#' @description This function aims to rename cell markers stored within a CYTdata object.
#'
#' This function is interesting to remove the names of the fluorochromes or metals recorded during the acquisition process.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#'
#' @return a S4 object of class 'CYTdata'
#'
#'@export
#'

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


#' @title Renames Clustering within a CYTdata object
#'
#' @description This function aims to rename Clustering stored within a CYTdata object.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param from a character vector providing the sample names to replace. By default, all sample names are replaced
#' @param to a character vector providing the new sample names to use
#'
#' @return a S4 object of class 'CYTdata'
#'
#'@export
#'

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



# pasteClustering <- function(CYTdata, clusterings, newClustering, clusteringSep = "_", clusteringPalette){
#
#   CYTdata = checkValidity(CYTdata, mode = "error")
#   if (nrow(CYTdata@clusteringData@cellClustering)==0 && ncol(CYTdata@clusteringData@cellClustering)==0) { stop("Error : clusteringData slot is empty") }
#
#   checkmate::qassert(clusterings, "S*")
#   checkmate::qassert(newClustering, "S1")
#   checkmate::qassert(clusteringSep, "S1")
#
#   if (length(clusterings)<2) { stop() }
#
#   clErr = setdiff(clusterings, colnames(CYTdata@clusteringData@cellClustering))
#   if (length(clErr)>0) {
#     stop("Error : 'clusterings' argument (", paste0(clusterings, collapse=", "),
#          ") contains names not being clustering column names(", paste0(clErr, collapse=", "), ")")
#   }
#   if (newClustering %in% colnames(CYTdata@clusteringData@cellClustering)) {
#     stop("Error : 'newClustering' argument provide name already being a clustering column name (", newClustering, ")")
#   }
#   if (!clusteringPalette %in% clusterings) {
#     stop("Error : 'clusteringPalette' argument provide name (", clusteringPalette,
#          ") not being a clustering column name contained in 'clusterings' argument (", paste0(clusterings, collapse=", "), ").")
#   }
#
#   cat("\n\nCurrent clustering column names of cellClustering dataframe :", paste0(clusterings, collapse=", "), "\n")
#   cat("\n\nWill be pasted, with the seperation operator '", clusteringSep,
#       "', into the following clustering column name : '", newClustering,
#       "'. And the clustering palette will be heritated from : '", clusteringPalette, "'.")
#
#   CYTdata@clusteringData@cellClustering[,newClustering] = CYTdata@clusteringData@cellClustering[,clusterings] %>%
#     apply(1, function(x) paste(x, collapse = clusteringSep))
#
#   #CYTdata@clusteringData@clusteringPalette[clusteringPalette] CYTdata@clusteringData@cellClustering[,clusteringPalette]
#
#   if (is.null(clustering)) { # Rename clustering columns
#     if (is.null(from)) { from = colnames(CYTdata@clusteringData@cellClustering) }
#     else {
#       clErr = setdiff(from, colnames(object@clusteringData@cellClustering))
#       if (length(clErr)>0) { stop("Error : 'from' argument contains names not present among cellClustering's clustering column names (", paste0(clErr, collapse=", "), ".") }
#     }
#     if (length(from)!=length(to)) { stop("Error : Length of argument 'from' (", length(from), ") and argument 'to' (", length(to), ") must be equal") }
#     dupfrom = unique(from[duplicated(from)])
#     if (length(dupfrom)>0) { stop("Error : Argument 'from' contains duplicated clustering column names :", paste0(dupfrom, collapse=", ")) }
#     dupto = unique(to[duplicated(to)])
#     if (length(dupto)>0) { stop("Error : Argument 'to' contains duplicated names :", paste0(dupto, collapse=", ")) }
#
#     cat("\n\nCurrent clustering column names of cellClustering dataframe are (in the order) :", paste0(colnames(CYTdata@clusteringData@cellClustering), collapse = ", "), "\n")
#     cat("\n\nThe following values :")
#     cat("\n - ", paste0(from, collapse=", "))
#     cat("\n\nwill be renamed, in the order, by :")
#     cat("\n - ", paste0(to, collapse=", "))
#     colnames(CYTdata@clusteringData@cellClustering)[match(colnames(CYTdata@clusteringData@cellClustering), from)] = to
#   }
#   else { # Rename clusters names
#     if (!clustering %in% colnames(CYTdata@clusteringData@cellClustering)) {
#       stop("Error : 'clustering' argument (", clustering, ") is not NULL or a clustering column names")
#     }
#     oldclulev = levels(CYTdata@clusteringData@cellClustering[,clustering])
#     if (is.null(from)) { from = oldclulev }
#     else {
#       clErr = setdiff(from, oldclulev)
#       if (length(clErr)>0) { stop("Error : 'from' argument contains clustering values not present among", clustering, "clustering values (", paste0(clErr, collapse=", "), ".") }
#       dupfrom = unique(from[duplicated(from)])
#       if (length(dupfrom)>0) { stop("Error : Argument 'from' contains duplicated", clustering, "clustering values :", paste0(dupfrom, collapse=", ")) }
#     }
#     if (length(from)!=length(to)) { stop("Error : Length of argument 'from' (", length(from), ") and argument 'to' (", length(to), ") must be equal") }
#     if (length(to) == 1) { to = rep(to, length(from)) }
#     newclulev = oldclulev
#     newclulev[match(from, newclulev)] = to
#
#     levdup = unique(newclulev[duplicated(newclulev)])
#     if (length(levdup)>0){
#       message("After renaming, several ", clustering, "'s cluster IDs have the same name (",  paste0(levdup, collapse=", "), "), either by renaming to an already existing and unchanged clustering value, or by duplicate in the 'to' argument, or both.")
#       if (merge) { message("\nWarning : 'merge' argument is set to TRUE. After renaming,", clustering, " column of cellClustering dataframe get its duplicated levels merged.") }
#       else { stop("\nError : 'merge' argument is set to FALSE. CYTdata unchanged.") }
#     }
#     cat("\n\nCurrent", clustering, "'s cluster IDs (levels) are (in the order) :", paste0(oldclulev, collapse = ", "), "\n")
#     cat("\n\nThe following values :")
#     cat("\n - ", paste0(from, collapse=", "))
#     cat("\n\nwill be renamed, in the order, by :")
#     cat("\n - ", paste0(to, collapse=", "))
#     CYTdata@clusteringData@cellClustering[,clustering] = plyr::mapvalues(CYTdata@clusteringData@cellClustering[,clustering], from = from, to = to)
#   }
#   CYTdata = checkValidity(CYTdata, mode = "warning")
#   return(CYTdata)
# }


#' @title Assigns meta-information about biological samples
#'
#' @description This function aims to attach meta-information to each biological sample.
#'
#' Especially, the following meta-information of each sample can be specified for subsequent analyses.
#' - The biological individual
#' - The biological condition (groups, vaccines, studies, etc.)
#' - The timepoint
#' Timepoint and Individual data must be specified.
#'
#' @param CYTdata a CYTdata object
#' @param clustering a dataframe containing meta-information about the biological samples.
#' The columns must contain, at least, a column named "Timepoint" and an other named "Individual".
#' The rownames have to be the biological samples, thus the number rows has to be equal to the number of samples.
#'
#' @return a S4 object of class 'CYTdata'
#'
#'
#' @export
#'

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

#' @title Assigns meta-information about biological samples
#'
#' @description This function aims to attach meta-information to each biological sample.
#'
#' Especially, the following meta-information of each sample can be specified for subsequent analyses.
#' - The biological individual
#' - The biological condition (groups, vaccines, studies, etc.)
#' - The timepoint
#' Timepoint and Individual data must be specified.
#'
#' @param CYTdata a CYTdata object
#' @param clustering a dataframe containing meta-information about the biological samples.
#' The columns must contain, at least, a column named "Timepoint" and an other named "Individual".
#' The rownames have to be the biological samples, thus the number rows has to be equal to the number of samples.
#'
#' @return a S4 object of class 'CYTdata'
#'
#'
#' @export
#'

checkgetMetaclutseringMembership <- function(CYTdata,
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


# @title Internal - Rescales marker expression by quantile
#
# @description This function is used internally to rescale the marker expression by the quantile method.
#
# @param exprs a vector containing one marker expression
# @param quant.low a numeric value providing the value of the first quantile
# @param quant.high a numeric value providing the value of the last quantile
#
# @return a data.frame containing quantile rescale marker expressions
#

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

# @title Internal - Rescales marker expression by quantile
#
# @description This function is used internally to rescale the marker expression by the quantile method.
#
# @param exprs a vector containing one marker expression
# @param quant.low a numeric value providing the value of the first quantile
# @param quant.high a numeric value providing the value of the last quantile
#
# @return a data.frame containing quantile rescale marker expressions
#

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

#' @title Identify cell cluster of having similar marker expressions using clustering techniques
#'
#'
#' @description This function aims to generate coordinates of data in an 2D-reduced space, stored in a CYTdata object.
#'
#' The algorithm used available are :
#' - Phenograph, a graph-based clustering technique
#' - Self-Organizing Maps (SOM), a neural network-based clustering technique
#' The whole set of cell markers or specific cell markers can be used during the clustering process.
#' The cell clustering can be performed on the DimReduction.coordinates representation or based on marker expression.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param markers a character vector providing the cell markers to use to generate the cluster data. By default, all markers are used
#' @param type a character value containing the type of clustering method to use. Possible values are: "Phenograph", "SOM" (default = Phenograph)
#' @param seed a numeric value providing the random seed to use during stochastic operations
#' @param ... additional arguments passed on to method from R package.
#' For SOM, please refer to som method from kohonen package : https://cran.r-project.org/web/packages/kohonen/kohonen.pdf
#' For DBSCAN, please refer to dbscan method from dbscan package :
#' For Phenograph, please refer to cytof_cluster method from cytofkit2 package : Rphenograph_k = 30 (ne fonctionne pas en dessous de 15 ?)
#'
#' @return a S4 object of class 'CYTdata'
#'
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
                              top_annotation = HeatmapAnnotation(df = dfClusters, col = CYTdata@clusteringData@clusteringPalette))
    }
  }

  CYTdata = checkValidity(CYTdata, mode = "warning")
  return(list("CYTdata" = CYTdata, "opts" = opts))
}


#' @title Identify meta cluster of clusters & Heatmap
#'
#'
#' @description This function aims to identify meta clusters, which are groups of clusters having similar expressions for selected markers. It allows to visualize the cluster marker expressions using heatmap.
#'
#' The mean of median marker expressions is computed for each cluster, and marker expressions displayed using a categorical heatmap (5 categories are defined by default).
#' The range expression of each cell marker is discretized into several categories between bounds of marker expressions.
#' To hierarchical clustering, shown using dendrogramm, can be computed on both marker and cluster levels.
#' Clustering method performed is Hierarchical clustering
#' Several clustering criterion are available such as Wald, Medoid, ...
#'
#' @param data a S4 object of class 'CYTdata'
#' @param N_Hierarchical a numeric value providing the number of metaclusters expected
#' @param criterion a string value providing the agglomeration method to be use.
#' Possible values are: 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'
#' Please refer to the function 'hclust' of the 'stats' package
#'
#'
#' @return a S4 object of class 'CYTdata'
#'
#' @export
#'
#'

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

#' @title
#'
#'
#' @description
#'
#' @param
#'
#'
#' @return
#'
#' @export
#'
#'

doPhenograph <- function(data) {
  stop()
}


#' @title Identify meta cluster of clusters & Heatmap
#'
#'
#' @description This function aims to identify meta clusters, which are groups of clusters having similar expressions for selected markers. It allows to visualize the cluster marker expressions using heatmap.
#'
#' The mean of median marker expressions is computed for each cluster, and marker expressions displayed using a categorical heatmap (5 categories are defined by default).
#' The range expression of each cell marker is discretized into several categories between bounds of marker expressions.
#' To hierarchical clustering, shown using dendrogramm, can be computed on both marker and cluster levels.
#' Clustering method performed is Hierarchical clustering
#' Several clustering criterion are available such as Wald, Medoid, ...
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param N_Metaclusters a numeric value providing the number of metaclusters expected
#' @param N_NewClusters a numeric value providing the number of phenotypic NewClusters expected
#' @param markers a character vector providing the marker names to use. By default, all markers are used
#' @param medoid.param a string value providing the type of cluster's markers expression used ("median" or "mean", "median" by default)
#' @param rescale.bounds a numeric vector providing the bounds ( in [0,1]) for rescale step
#' @param criterion a string value providing the agglomeration method to be use.
#' Possible values are: 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'
#' Please refer to the function 'hclust' of the 'stats' package
#' @param seed.order.hclust a numeric value providing the random seed to use during stochastic operations
#' @param nb.cat.heatmap a numeric value specifying the number of categories to use
#'
#'
#' @return a S4 object of class 'CYTdata'
#'
#' @export
#'

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



#' @title Renames markers within a CYTdata object
#'
#' @description This function aims to rename cell markers stored within a CYTdata object.
#'
#' This function is interesting to remove the names of the fluorochromes or metals recorded during the acquisition process.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param from a character vector providing the marker names to replace. By default, markers names from markers slot are replaced, in the order
#' @param to a character vector providing the new marker names to use
#'
#' @return a S4 object of class 'CYTdata'
#'
#'@export
#'

updateClusteringLevels <- function(CYTdata, andPalette = TRUE) {
  mdFrame = CYTdata@clusteringData@cellClustering %>% mutate(across(where(is.factor), droplevels))
  palList = CYTdata@clusteringData@clusteringPalette
  for (md in colnames(mdFrame)) { palList[[md]] = palList[[md]][levels(mdFrame[,md])] }
  CYTdata@clusteringData@cellClustering = mdFrame
  CYTdata@clusteringData@clusteringPalette = palList[colnames(mdFrame)]
  return(CYTdata)
}







#' @title Plots the numbers of cells of each population (cluster or metacluster)
#'
#' @description This function aims to visualize the number of cells associated to each population. The populations can be clusters or metaclusters
#'
#' This representation displays the population in the X-axis and the total number of associated cells in the Y-axis.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param population a character vector containing the identifiers of the population to use. By default, all the population are used
#' @param level a character value indicating the type of population plotted. Possible values are: "clusters", "metaclusters". By default, 'clusters' are used.
#' @param sort a boolean value indicating if population must be sorted by the number associated cluster
#' @param color.metadata a character value specifying the metadat used to color the barplot. By default, color.metadata is set to NULL (the barplot color is uniform)
#'
#' @return a ggplot2 object
#'
#' @export
#'
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





























# plotSOMmap <- function(CYTdata){
#   useSOMweights = FALSE
#   checkmate::qassert(useSOMweights, "B1")
#   if (useSOMweights) {
#     if (level == "metaclusters") {
#       stop("Error : SOM weights must be used but metaclusters cannot be determined using SOM")
#     }
#     if (is.null(CYTdata@Clustering@optionnal$SOMweights)) {
#       stop("Error : No SOM weights available in Clustering slot.
#            Impossible to use it for metaclustering")
#     }
#     SOMweights = CYTdata@Clustering@optionnal$SOMweights
#     if (!identical(rownames(SOMweights),population)){
#       stop("Error : SOM weights are inconsistent with clusters identifiers.")
#     }
#     markErr = setdiff(colnames(SOMweights), markers)
#     if (length(markErr)>0) {
#       stop("Error : Some markers selected for metaclustering are not available
#            in SOMweights matrix (", paste0(markErr, collapse=" ,"), ").")
#     }
#     data = SOMweights[population,]
#   }
# }
#
# plotPhenograph <- function(CYTdata){
#
# }
#
# plotSpadeMST <- function(CYTdata){
#
# }

#' BuildMST
#'
#' Build Minimal Spanning Tree
#'
#' Add minimal spanning tree description to the FlowSOM object
#'
#' @param fsom   FlowSOM object, as generated by \code{\link{BuildSOM}}
#' @param silent If \code{TRUE}, no progress updates will be printed
#' @param tSNE   If \code{TRUE}, an alternative t-SNE layout is computed as well
#'
#' @return FlowSOM object containing MST description
#'
#' @seealso \code{\link{BuildSOM}}, \code{\link{PlotStars}}
#'
#' @examples
#' # Read from file, build self-organizing map
#' fileName <- system.file("extdata", "68983.fcs", package="FlowSOM")
#' flowSOM.res <- ReadInput(fileName, compensate=TRUE, transform = TRUE,
#'                          scale = TRUE)
#' flowSOM.res <- BuildSOM(flowSOM.res, colsToUse = c(9, 12, 14:18))
#'
#' # Build the Minimal Spanning Tree
#' flowSOM.res <- BuildMST(flowSOM.res)
#'
#' @importFrom Rtsne Rtsne
#'
#' @export
BuildMST <- function(fsom, silent = FALSE, tSNE = FALSE){

  fsom$MST <- list()
  if(!silent) message("Building MST\n")

  adjacency <- as.matrix(stats::dist(fsom$map$codes, method = "euclidean"))
  fullGraph <- igraph::graph.adjacency(adjacency,
                                       mode = "undirected",
                                       weighted = TRUE)
  fsom$MST$graph <- igraph::minimum.spanning.tree(fullGraph)
  ws <- igraph::edge.attributes(fsom$MST$graph)$weight
  #normalize edge weights to match the grid size in coords (see below)
  ws <- ws / mean(ws)
  igraph::edge.attributes(fsom$MST$graph)$weight <- ws
  fsom$MST$l <- igraph::layout.kamada.kawai(
    coords = as.matrix(fsom$map$grid),
    fsom$MST$graph)


  if(tSNE){
    fsom$MST$l2 <- Rtsne::Rtsne(fsom$map$codes)$Y
    #library(RDRToolbox)
    #fsom$MST$l2 <- Isomap(fsom$map$codes, dims = 2, k = 3)[[1]]
  }

  library(vegan)
  data(BCI, BCI.env)
  BCI
  vegan::diversity(BCI, index = "simpson")

  return(fsom)
}









#' @title Plot compare cluster
#'
#' @description This function aims to visualise the expression of a given marker in each cluster.
#'
#' @param CYTdata a CYTdata object
#' @param clusters1 a character vector containing the identifier of the cluster to use
#' @param clusters2 a character vector containing the identifier of the cluster to use
#'
#' @return a ggplot2 object
#'
#' @export
#'
plotCompareClusters <- function(CYTdata,
                                clusters1,
                                clusters2) {

  checkmate::qassert(clusters1, "S1")
  checkmate::qassert(clusters2, "S1")

  matrix.exp <- CYTdata@matrix.expression
  cluster <- CYTdata@identify.clusters

  proj <- cbind(matrix.exp, cluster)

  parameters <- colnames(proj)
  markers <- parameters[!parameters %in% c("cluster")]

  exp.values.clusters1 <- proj[proj$cluster %in% clusters1, ]
  exp.values.clusters2 <- proj[proj$cluster %in% clusters2, ]

  exp.values.clusters <- rbind(exp.values.clusters1, exp.values.clusters2)

  melt.matrix <- reshape::melt(exp.values.clusters, id.vars = c("cluster"))

  plot <- ggplot2::ggplot() +
    ggplot2::geom_density(data = melt.matrix,
                          mapping = ggplot2::aes_string(x = "value",
                                                        fill = "cluster"),
                          size = 0.5,
                          alpha = 0.8) +
    ggplot2::facet_wrap(~variable)

  plot <- plot +
    ggplot2::scale_fill_manual(values = c("#ff0000", "#0000ff")) +
    ggplot2::scale_x_continuous(limits = c(-2, 6), breaks = c(-2:6)) +
    ggplot2::scale_y_continuous(labels = function(x) {format(round(x, 1),nsmall = 1)})

  plot <- plot +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 10, face = "bold"),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.key.width = grid::unit(0.35, "cm"))

  return(plot)
}

#' @title Plots of a distogram of marker co-expression
#'
#' @description This function aims to visualize the pairwise co-expression between all markers using a distogram representation.
#' Each tile corresponds to the co-expression between two markers and is gradient-colored based on the Pearson or Spearman correlation
#'
#' @param CYTdata a CYTdata object
#' @param clusters a character vector containing the identifier of the cluster to use
#' @param method a character value indicating the name of correlation method to use
#'
#' @return a ggplot2 object
#'
#' @export
#'
plotDistogram <- function(CYTdata,
                          markers = NULL,
                          population = NULL,
                          level = c("clusters", "metaclusters"),
                          samples = NULL,
                          method = c("pearson","spearman"),
                          palette = NULL) { #c("yellow", "orange", "red", "brown")

  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { validObject(CYTdata) }

  markers = checkorderMarkers(CYTdata, markers, order=FALSE, checkDuplicates=TRUE)

  level = match.arg(level)
  checkmate::qassert(level, "S1")
  population = checkorderPopulation(CYTdata, population=population, level=level,
                                    order=TRUE, checkDuplicates=TRUE)
  samples = checkorderSamples(CYTdata, samples,
                              order=TRUE, checkDuplicates=TRUE)

  method = match.arg(method)
  checkmate::qassert(method, "S1")

  if(level == "clusters") { popId = CYTdata@Clustering@clusters }
  else { popId = CYTdata@Metaclustering@metaclusters }

  data = subset(CYTdata@matrix.expression[,markers],
                popId %in% population & CYTdata@samples %in% samples)

  corMatrix = round(stats::cor(data, method = method), 3)
  dist = stats::as.dist(1 - corMatrix)
  hc = stats::hclust(dist)
  corMatrix = corMatrix[hc$order, hc$order]
  corMatrix[upper.tri(corMatrix, diag = TRUE)] = NA

  print(corMatrix)

  dimnames(corMatrix) = NULL
  corMatrixmelted = reshape2::melt(corMatrix)

  print(corMatrixmelted)

  plot <- ggplot2::ggplot(data = corMatrixmelted,
                          ggplot2::aes_string(x = "Var1",
                                              y = "Var2",
                                              fill = "value")) +
    ggplot2::ggtitle("Distogram") +
    ggplot2::geom_tile(color = "white")

  checkmate::qassert(palette, c("0", "S*"))
  if(!is.null(palette)){
    if(!all(areColors(palette))){
      stop("Error : 'palette' argument (", paste0(palette, collapse = ","),
           ") does not contain only hexadecimal color.)")
    }
    plot <- plot + ggplot2::scale_fill_gradientn(colours = palette, na.value = "white", name = paste(method, "correlation"))
  }
  else {
    plot <- plot + ggplot2::scale_fill_gradient2(low = "green", high = "red", mid = "black",
                                                 midpoint = 0, limit = c(-1, 1), na.value = "white",
                                                 name = paste(method, "correlation"),
                                                 oob = scales::oob_squish)
  }

  plot <- plot +
    ggplot2::annotate(geom = "text",
                      x = seq(1, length(markers)),
                      y = seq(1, length(markers)),
                      angle = -45,
                      size = 4,
                      label = markers,
                      hjust = 1) +
    ggplot2::coord_fixed() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.background = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.4, 0.7),
      legend.direction = "horizontal",
      legend.key = ggplot2::element_blank()) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = 7,
                                                   barheight = 1,
                                                   title.position = "top",
                                                   title.hjust = 0.5))

  rownames(corMatrix) = markers
  colnames(corMatrix) = markers
  return(list("distogram" = plot,
              "cor" = corMatrix))
}










#' @title  Computes the percentage of clusters with sufficient number of cells
#'
#' @description This function aims to :
#'
#'  - Compute and show cell clusters having a number of associated cells greater than a specific threshold
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' level = c("clusters", "metaclusters")
#' @param thSize a numeric value providing the minimum number of cells needed for a cluster to be considered a small cluster (default value = 50)
#'
#' @return a list containing QC size features
#'
#' @export
#'
#'
#'# @param sort a boolean value indicating if clusters must be sorted by the number associated cells

computePopulation_SizeQC <- function(CYTdata,
                                     population = NULL,
                                     level = c("clusters", "metaclusters"),
                                     thSize = 50,
                                     generateHeatmap = FALSE){

  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }

  level = match.arg(level)
  checkmate::qassert(level, "S1")
  population = checkorderPopulation(CYTdata=CYTdata, population=population, level=level, order=TRUE, checkDuplicates=TRUE)

  checkmate::qassert(thSize, "N1")
  if (thSize<=0) { stop("Error : 'thSize'  argument must be positive integer.") }
  checkmate::qassert(generateHeatmap, "B1")

  if (level == "clusters"){ cellcount = CYTdata@Clustering@cellcount[population,] }
  else { cellcount = CYTdata@Metaclustering@cellcount[population,] }

  effectif = apply(cellcount, 1, sum)
  vectorScores = as.logical(effectif > thSize)
  score = mean(vectorScores)*100

  res = list("cellcount" = cellcount,
             "vectorScores" = vectorScores,
             "thSize" = thSize,
             "score" = score)

  if (generateHeatmap){
    cat("\n\n Generating heatmap of QCsize ..")
    data = cbind.data.frame(cellcount, "total.cells" = vectorScores)
    data = reshape2::melt(data.matrix(data))
    colnames(data) = c("groupId", "samples", "count")
    data$passed = TRUE
    data$passed[as.logical(data$count < thSize)] = FALSE
    data$count = NULL
    title = paste("Quality control on", level, "size", sep = " ")
    subtitle = paste("Percentage of", level, "having a sufficient number of cells =",
                     format(round(score), nsmall = 2), "%", sep=" ")

    heatmap <- ggplot2::ggplot() +
      ggplot2::geom_tile(data = data,
                         ggplot2::aes_string(x = "samples", y = "groupId", fill = "passed"),
                         colour = "black") +
      ggplot2::scale_fill_manual(values = c("TRUE" = "green", "FALSE" = "red")) +
      ggplot2::geom_vline(xintercept = (length(unique(data$samples)) - 0.5), colour = "black", size = 2) +
      ggplot2::scale_x_discrete(expand = c(0, 0)) + ggplot2::scale_y_discrete(expand = c(0, 0)) +
      ggplot2::xlab("samples") + ggplot2::ylab("clusters") +
      ggplot2::ggtitle(bquote(atop(.(title), atop(italic(.(subtitle)), "")))) +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.background = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                     axis.line = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     panel.border = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_blank(),
                     legend.title = ggplot2::element_blank(),
                     legend.key = ggplot2::element_blank())
    res$heatmap = heatmap
  }

  cat("\n\n - Percentage of ", level, " having at least ", thSize ," cells = ",
      format(round(res$score, 2), nsmall = 2), "%\n\n")

  return(res)
}

#' @title  Computes the percentage of clusters with uniform phenotypes
#'
#' @description This function aims to :
#'
#'  - Identify and show cell clusters having a uniform phenotype (an unimodal expression and a low spread of expression for selected markers)
#'
#'  @details
#'
#' - Check unimodal expression : testing unimodal distribution of markers with a Hartigans test and check p-value
#'
#' - Check low spread of expression : check if a distribution of markers is not below the IQR threshold (interquantile range)
#'
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param thPvalue a numeric value providing the p-value threshold of the Hartigan's dip test (unimodal if pvalue > thPvalue)
#' @param th.IQR a numeric value providing the IQR (interquartile range) threshold to assume a distribution as uniform
#' @param th.good.cluster a numeric value providing the proportion threshold of markers expression distribution needed to pass the QC to assume a cluster as 'good'
#' @param th.good.clustering a numeric value providing the proportion threshold of 'good' clusters needed to assume the clustering passing the QC
#' @param ... arguments to pass onto "modetest" function from multimode package (https://www.rdocumentation.org/packages/multimode/versions/1.5/topics/modetest)
#' method silvermann recommend parallelize
#' @return a list containing QC uniform features
#'
#' @export
#'
#'

computeQCuniform <- function(CYTdata,
                             markers = NULL,
                             population = NULL,
                             level = c("clusters", "metaclusters"),
                             aggregate = FALSE,
                             thPvalue = 0.05,
                             thIQR = 2,
                             thGood = 0.75,
                             method = c("Hartigan", "Silvermann"),
                             NbSim = NULL,
                             generateHeatmap = TRUE,
                             sortQCbyMarkers = TRUE,
                             generateEcdf = TRUE, # empirical cumulative distribution
                             parallelize = FALSE,
                             nbCores = NULL) {

  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }

  checkmate::qassert(thPvalue, "N1")
  if (!(thPvalue>=0 && thPvalue<=1)) {
    stop("Error : 'thPvalue' argument is a p-value threshold and must be a positive numeric between 0 and 1.")
  }
  checkmate::qassert(thIQR, "N1")
  if (thIQR<=0) { stop("Error : 'thIQR' argument must be positive numerical value.") }
  checkmate::qassert(thGood, "N1")
  if (!(thGood>=0 && thGood<=1)) {
    stop("Error : 'thPvalue' argument is a proportion threshold and must be a positive numeric between 0 and 1.")
  }
  method = match.arg(method)
  checkmate::qassert(method, "S1")
  checkmate::qassert(NbSim, c("0","N1"))
  if (!is.null(NbSim) && NbSim<=0) { stop("Error : 'NbSim' argument must be positive integer.") }
  checkmate::qassert(generateHeatmap, "B1")
  checkmate::qassert(sortQCbyMarkers, "B1")
  checkmate::qassert(generateEcdf, "B1")
  checkmate::qassert(parallelize, "B1")
  checkmate::qassert(nbCores, c("0","N1"))
  if (!is.null(nbCores) && nbCores<=0) { stop("Error : 'nbCores' argument must be positive integer.") }

  level = match.arg(level)
  checkmate::qassert(level, "S1")
  population = checkorderPopulation(CYTdata, population=population, level=level,
                                    order=TRUE, checkDuplicates=TRUE)
  markers = checkorderMarkers(CYTdata, markers=markers, order=TRUE, checkDuplicates=TRUE)

  if (level == "clusters"){ popId = CYTdata@Clustering@clusters }
  else { popId = CYTdata@Metaclustering@metaclusters }

  switch(method,
         Hartigan = {
           if (is.null(NbSim)) NbSim = 2000
           getPvalue <- function(x) return(diptest::dip.test(x, B=NbSim)$p.value)
         },
         Silvermann = {
           if (is.null(NbSim)) NbSim = 20
           getPvalue <- function(x) return(multimode::modetest(x, method="SI", B=NbSim)$p.value)
         })

  if (length(population)==1) {
    cat("'population' argument has length equal to 1, representing a single population.
        \n Remark 1 : 'parallelize' argument is ignored here, no need to parallelize
        computation for a single population.
        \n Remark 2 : 'generateEcdf' argument is ignored also.")
    data = CYTdata@matrix.expression[popId %in% population, markers]
    IQR = apply(data, 2,
                function(z){
                  quantiles = stats::quantile(z)
                  return(quantiles[4] - quantiles[2])})
    names(IQR) = markers
    pvalues = apply(data, 2, getPvalue)
    names(pvalues) = markers

    uniform = pvalues <= thPvalue & IQR <= thIQR
    names(uniform) = markers

    popScore = mean(uniform)*100
    score = ifelse(popScore>=thGood*100, 100, 0)

    res <- list("Pvalue" = IQR,
                "IQR" =  pvalues,
                "Uniform" = uniform,
                "goodScores" = popScore,
                "Score" = score,
                "thPvalue" = thPvalue,
                "thIQR" = thIQR,
                "thGood" = thGood)

    if (generateHeatmap){
      data = data.frame("passed" = as.vector(uniform),
                        "markers" = names(uniform))
      res$heatmap <- ggplot2::ggplot() +
        ggplot2::geom_tile(ggplot2::aes_string(x = "markers",
                                               y = 1,
                                               fill = "passed"),
                           colour = "black", size = 0.25) +
        ggplot2::ggtitle(label = paste("Uniform", level, "quality control", sep=" "),
                         subtitle = paste("Percentage of markers passing the QC (for ",
                                          level, " ", population, " ) : ", format(round(score, 2), nsmall = 2), "%")) +
        ggplot2::scale_fill_manual(values = c("TRUE" = "green", "FALSE" = "red")) +
        ggplot2::scale_x_discrete(expand = c(0, 0)) +
        ggplot2::scale_y_discrete(expand = c(0, 0)) +
        ggplot2::ylab(level) + ggplot2::xlab("Markers") +
        ggplot2::theme_minimal() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                       plot.background  = ggplot2::element_blank(),
                       axis.text.x      = ggplot2::element_text(angle = 90, hjust = 1, vjust = 1),
                       axis.line        = ggplot2::element_blank(),
                       panel.background = ggplot2::element_blank(),
                       panel.border     = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_blank())
    }
  }
  else {

    data = cbind.data.frame("popId" = popId, CYTdata@matrix.expression[,markers])

    IQR = plyr::ddply(subset(data, popId %in% population),
                      "popId",
                      function(y){
                        apply(y[,-1], 2, function(z){
                          quantiles = stats::quantile(z)
                          return(quantiles[4] - quantiles[2])
                        })
                      })
    rownames(IQR) = population
    IQR$popId = NULL
    print(IQR)

    if(parallelize){
      cat("'parallelize' argument set to TRUE, creating core clusters..")
      maxCores = parallel::detectCores() - 1
      if (is.null(nbCores)) nbCores = maxCores
      if (nbCores > maxCores) stop("Error : nbCores argument require more cores than available.")

      message("\n- Parallelization activated : Register the parallel backend using ",
              nbCores, " cores...")
      cl = parallel:makeCluster(nbCores)
      doParallel::registerDoParallel(cl)

      Pvalue = foreach(i = 1:length(population), .combine = 'rbind.data.frame') %dopar% {
        pop = population[i]
        subdata = data[data$popId == pop, -1]
        pvalues = apply(subdata, 2, FUN = getPvalue)
      }
    }
    else {
      Pvalue = data.frame()
      for (pop in population) {
        subdata = data[data$popId == pop, -1]
        pvalues = apply(subdata, 2, FUN = getPvalue)
        Pvalue = rbind.data.frame(Pvalue, pvalues)
      }
    }

    rownames(Pvalue) = population
    colnames(Pvalue) = colnames(data)[-1]
    print(Pvalue)
    Uniform = data.frame(Pvalue <= thPvalue & IQR <= thIQR,
                         check.names = FALSE)
    print(Uniform)
    rownames(Uniform) = population
    print(Uniform)

    popScores = apply(Uniform, 1, mean)*100
    overallPassed = as.logical(as.vector(popScores) >= thGood*100)
    score = mean(overallPassed)*100
    res <- list("Pvalue" = Pvalue,
                "IQR" =  IQR,
                "Uniform" = Uniform,
                "goodScores" = popScores,
                "Score" = score,
                "thPvalue" = thPvalue,
                "thIQR" = thIQR,
                "thGood" = thGood)

    if (generateHeatmap) {

      cat("\n\n Generating heatmap of QCuniform ..")

      data = cbind.data.frame(Uniform, "overallPassed" = overallPassed)

      if (sortQCbyMarkers) {
        markersQC = apply(data[,-ncol(data)], 2, mean)*100
        lev = c(names(markersQC)[order(markersQC, decreasing=TRUE)], "overallPassed")
      }
      else {
        lev = colnames(Uniform)
      }

      data = reshape2::melt(data.matrix(data))
      colnames(data) = c(level, "markers", "passed") # row must be cluster
      data$passed = as.logical(data$passed)
      data$markers = factor(data$markers, levels = lev)

      res$heatmap = ggplot2::ggplot(data = data) +
        ggplot2::geom_tile(ggplot2::aes_string(x = "markers",
                                               y = level,
                                               fill = "passed"),
                           colour = "black", size = 0.25) +
        ggplot2::ggtitle(label = paste("Uniform", level, "quality control", sep=" "),
                         subtitle = paste("Percentage of", level, "being considered as 'good'",
                                          level, ": ", format(round(score, 2), nsmall = 2), "%")) +
        ggplot2::scale_fill_manual(values = c("TRUE" = "green", "FALSE" = "red")) +
        ggplot2::scale_x_discrete(expand = c(0, 0)) +
        ggplot2::scale_y_discrete(expand = c(0, 0)) +
        ggplot2::ylab(level) + ggplot2::xlab("Markers") +
        ggplot2::theme_minimal() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                       plot.background  = ggplot2::element_blank(),
                       axis.text.x      = ggplot2::element_text(angle = 90, hjust = 1, vjust = 1),
                       axis.line        = ggplot2::element_blank(),
                       panel.background = ggplot2::element_blank(),
                       panel.border     = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_blank())
    }

    if (generateEcdf) {

      cat("\n\n Generating empirical cumulative distribution..")

      res$ecdf = ggplot2::ggplot(data.frame("x" = popScores),
                                 ggplot2::aes(x)) +
        ggplot2::stat_ecdf(geom = "step", pad = FALSE) +
        ggplot2::geom_vline(xintercept = thGood, colour = "red", size = 1, linetype = "longdash") +
        ggplot2::geom_text(aes(x = thGood,
                               y = 1,
                               label = thGood),
                           hjust = -1, colour = "red", size = 5) +
        ggplot2::geom_hline(yintercept = score, colour = "red", size = 1) +
        ggplot2::geom_text(aes(y = score,
                               x = 1,
                               label = score),
                           vjust = -1, colour = "red", size = 5) +
        ggplot2::xlab(paste("Good", level, "threshold value", sep=" ")) +
        ggplot2::ylab(paste("Percentage of good", level, sep=" ")) +
        ggplot2::ggtitle(paste("Percentage of good", level, "over good",
                               level, "threshold value", sep=" ")) +
        ggplot2::theme_minimal() +
        ggplot2::theme(plot.title = ggplot2::element_text(size=15, hjust = 0.5),
                       axis.title = ggplot2::element_text(size = 18),
                       axis.text = ggplot2::element_text(size = 16))
    }
  }

  cat("\n\n - Percentage of ", level, " having at least",
      res$thGood ,"% of marker expression distribution passing the QC =",
      format(round(res$Score, 2), nsmall = 2), "% \n\n")

  return(res)
}






#' @title  Computes the percentage of clusters with sufficient number of cells
#'
#' @description This function aims to :
#'
#'  - Compute and show cell clusters having a number of associated cells greater than a specific threshold
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' level = c("clusters", "metaclusters")
#' @param THsize a numeric value providing the minimum number of cells needed for a cluster to be considered a small cluster (default value = 50)
#'
#' @return a list containing QC size features
#'
#' @export
#'
#'
#'# @param sort a boolean value indicating if clusters must be sorted by the number associated cells

computeQCsize <- function(CYTdata,
                          level = c("clusters", "metaclusters"),
                          THsize = 50,
                          generateHeatmap = FALSE){

  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.")}
  level = match.arg(level)
  checkmate::qassert(level, "S1")
  checkmate::qassert(THsize, "N1")
  checkmate::qassert(generateHeatmap, "B1")

  if (level == "clusters"){
    subObject = CYTdata@Clustering
    if (length(CYTdata@Clustering@clusters)==0) {
      stop("Error : 'level' argument require clustering step to be performed. Please perform it")
    }
  }
  else {
    subObject = CYTdata@Metaclustering
    if (length(CYTdata@Metaclustering@metaclusters)==0) {
      stop("Error : 'level' argument require metaclustering step to be performed. Please perform it")
    }
  }

  effectif = apply(subObject@cellcount, 1, sum)
  vectorScores = as.logical(effectif > THsize)
  score = mean(vectorScores)*100

  res = list("cellcount" = subObject@cellcount,
             "vectorScores" = vectorScores,
             "THsize" = THsize,
             "score" = score)

  if (generateHeatmap){
    data.heatmap = cbind.data.frame(subObject@cellcount, "total.cells" = clusters.effectif)
    res$heatmap = heatmapQCsize(data.heatmap, THsize, score, level)
  }

  cat("\n\n - Percentage of ", level, " having at least ", res$THsize ," cells = ",
      format(round(res$score, 2), nsmall = 2), "%\n\n")

  return(res)
}

# Heatmap of size QC for clusters/metaclusters :
# - Visualize the number of cells associated to each cluster, by samples
# - Coloration using QC results
# return a ggplot2 object

heatmapQCsize <- function(data, THsize, score, level) {

  cat("\n\n Generating heatmap of QCsize ..")

  data = reshape2::melt(data.matrix(data))
  colnames(data) = c("groupId", "samples", "count")
  data$passed = TRUE
  data$passed[as.logical(data$count < THsize)] = FALSE
  data$count = NULL
  title = paste("Quality control on", level, "size", sep = " ")
  subtitle = paste("Percentage of", level, "having a sufficient number of cells =",
                   format(round(score), nsmall = 2), "%", sep=" ")

  heatmap <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = data,
                       ggplot2::aes_string(x = "samples", y = "groupId", fill = "passed"),
                       colour = "black") +
    ggplot2::scale_fill_manual(values = c("TRUE" = "green", "FALSE" = "red")) +
    ggplot2::geom_vline(xintercept = (length(unique(data$samples)) - 0.5), colour = "black", size = 2) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) + ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::xlab("samples") + ggplot2::ylab("clusters") +
    ggplot2::ggtitle(bquote(atop(.(title), atop(italic(.(subtitle)), "")))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.background = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                   axis.line = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   legend.key = ggplot2::element_blank())

  return(heatmap)
}

lineplotQCuniform <- function(vectorScores, score, THgood, level) {

  cat("\n\n Generating lineplot of QCuniform ..")

  dataScores = data.frame("scores" = vectorScores,
                          "threshold" = as.numeric(names(vectorScores)))

  lineplot <- ggplot2::ggplot(dataScores, aes(x = threshold, y = scores)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::geom_vline(xintercept = THgood, colour = "red", size = 1, linetype = "longdash") +
    ggplot2::geom_text(aes(x = THgood,
                           y = min(dataScores$scores),
                           label = THgood),
                       hjust = -1, colour = "red", size = 5) +
    ggplot2::geom_hline(yintercept = score, colour = "red", size = 1) +
    ggplot2::geom_text(aes(y = score,
                           x = (mean(dataScores$threshold) + min(dataScores$threshold))/2,
                           label = score),
                       vjust = -1, colour = "red", size = 5) +
    ggplot2::xlab(paste("Good", level, "threshold value", sep=" ")) +
    ggplot2::ylab(paste("Percentage of good", level, sep=" ")) +
    ggplot2::ggtitle(paste("Percentage of good", level, "over good",
                           level, "threshold value", sep=" ")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(size=15, hjust = 0.5),
                   axis.title = ggplot2::element_text(size = 18),
                   axis.text = ggplot2::element_text(size = 16))

  return(lineplot)
}

# Heatmap of QC for uniform phenotype cluster :
# - Visualize the QC uniform results of clusters to each markers
# - Coloration using QC results
# return a ggplot2 object

heatmapQCuniform <- function(score, dataHeatmap, level) {

  cat("\n\n Generating heatmap of QCuniform ..")

  data = reshape2::melt(data.matrix(dataHeatmap))
  colnames(data) = c(level, "markers", "passed") # row must be cluster
  data$passed = as.logical(data$passed)

  markersQC = apply(dataHeatmap[,-ncol(dataHeatmap)], 2, mean)*100
  lev = c(names(markersQC)[order(markersQC, decreasing=TRUE)], "QC_passed")
  data$markers = factor(data$markers, levels = lev)

  heatmap <- ggplot2::ggplot(data = data) +
    ggplot2::geom_tile(ggplot2::aes_string(x = "markers",
                                           y = level,
                                           fill = "passed"),
                       colour = "black", size = 0.25) +
    ggplot2::ggtitle(label = paste("Uniform", level, "quality control", sep=" "),
                     subtitle = paste("Percentage of", level, "being considered as 'good'",
                                      level, ": ", format(round(score, 2), nsmall = 2), "%")) +
    ggplot2::scale_fill_manual(values = c("TRUE" = "green", "FALSE" = "red")) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::ylab(level) + ggplot2::xlab("Markers") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.background  = ggplot2::element_blank(),
                   axis.text.x      = ggplot2::element_text(angle = 90, hjust = 1, vjust = 1),
                   axis.line        = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.border     = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank())

  return(heatmap)
}


