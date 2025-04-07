#' @title
#'
#' @description
#'
#' @param
#'
#' @return
#'
#' @export
#'

checkorderMarkers <- function(CYTdata,
                              markers,
                              cellSlot = c("cellExprs", "cellAdditionalexprs"),
                              order = TRUE, checkDuplicates = TRUE){

  #CYTdata = checkValidity(CYTdata, mode = "error")
  checkmate::qassert(markers, c("0","S*"))
  cellSlot = match.arg(cellSlot)
  checkmate::qassert(cellSlot, "S1")
  checkmate::qassert(order, "B1")

  if (!is.null(markers) && length(markers)==0) { stop("Error : markers argument is an empty vector (length=0, but non NULL).") }

  checkmate::qassert(checkDuplicates, "B1")
  if (checkDuplicates) {
    markersdup = markers[duplicated(markers)]
    if (length(markersdup)>0) {
      stop("Error : markers argument contain duplicated values ( ", paste0(markersdup, collapse = ", "), " ). It must be vector of unique markers.")
    }
  }

  if (cellSlot == "cellExprs") { markersId = colnames(CYTdata@cellData@cellExprs) }
  else if (cellSlot == "cellAdditionalexprs") { markersId = colnames(CYTdata@cellData@cellAdditionalexprs) }

  if (is.null(markers)) { markers = markersId }
  else {
    markErr = setdiff(markers, markersId)
    if (length(markErr)>0) { stop("Error : 'markers' argument providing markers not present in", cellSlot, "'s colum names (", paste0(markErr, collapse=", "), ")") }
    if (order) { markers = markersId[markersId %in% markers] }
  }

  return(markers)
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

renameMarkers <- function(CYTdata,
                          cellSlot = c("cellExprs", "cellAdditionalexprs"),
                          from=NULL, to=NULL,
                          removeConjugate=FALSE, sepConjugate="_"){

  CYTdata = checkValidity(CYTdata, mode = "error")
  cellSlot = match.arg(cellSlot)
  checkmate::qassert(cellSlot, "S1")
  from = checkorderMarkers(CYTdata, markers=from, cellSlot=cellSlot, order=FALSE, checkDuplicates=TRUE)

  checkmate::qassert(removeConjugate, "B1")
  if (removeConjugate) {
    checkmate::qassert(sepConjugate, "S1")
    message("Warning : 'removeConjugate' argument is set to TRUE, 'to' argument is ignored. The first part (before underscore) of feature names given in 'from' argument are removed and features are renamed." )
    to = sapply(stringr::str_split(from, sepConjugate), function(x){
      if (length(x)>1) { return(paste(x[-1], collapse=sepConjugate)) }
      else { return(x[1]) }
    })
  }
  else {
    checkmate::qassert(to, "S*")
    if (length(from)!=length(to)) { stop("Error : Length of argument 'from' (", length(from), ") and argument 'to' (", length(to), ") must be equal") }
  }

  if (cellSlot == "cellExprs") {
    markersId = colnames(CYTdata@cellData@cellExprs)
    colnames(CYTdata@cellData@cellExprs)[match(from, markersId)] = to
  }
  else if (cellSlot == "cellAdditionalexprs") {
    markersId = colnames(CYTdata@cellData@cellAdditionalexprs)
    colnames(CYTdata@cellData@cellAdditionalexprs)[match(from, markersId)] = to
  }

  cat("\n\nCurrent marker names are (in the order) :", paste0(markersId, collapse = ", "))
  cat("\n\n\nThe following markers :")
  cat("\n\n - ", paste0(from, collapse=", "))
  cat("\n\n\nwill be renamed, in the order, by :")
  cat("\n - ", paste0(to, collapse=", "))

  CYTdata = checkValidity(CYTdata, mode = "warning")
  return(CYTdata)
}


#' @title
#'
#' @description
#'
#' @param
#' @return a S4 object of class 'CYTdata'
#'
#' @importFrom checkmate qassert
#' @importFrom flowCore read.FCS
#'
#' @export
#'

transformCYTdata <- function(CYTdata, markers = NULL, mode = c("arcsinh", "logicle", "biexponential"), cellSlot = c("cellExprs", "cellAdditionalexprs"), a=0, b=0.2, c=0, d=1, f=0, w=0){

  CYTdata = checkValidity(CYTdata, mode = "error")
  cellSlot = match.arg(cellSlot)
  checkmate::qassert(cellSlot, "S1")
  mode = match.arg(mode)
  checkmate::qassert(mode, "S1")

  markers = checkorderMarkers(CYTdata, markers=markers, cellSlot=cellSlot, order=FALSE, checkDuplicates=TRUE)

  switch(mode,
         arcsinh = { trans <- function(exp) {
           exp = sapply(exp, function(e){return(ifelse(is.na(e), NA, asinh(a+b*e)+c))})
           return(exp)
         }},
         logicle = { stop() },
         biexponential = { trans <- function(exp) {
           exp = sapply(exp, function(e){return(ifelse(is.na(e), NA, a*exp(b*(.-w))-c*exp(-d*(.-w))+f))})
           return(exp)
         }})

  if (cellSlot == "cellExprs") {
    CYTdata@cellData@cellExprs = CYTdata@cellData@cellExprs %>% mutate_at(markers, trans)
  }
  else if (cellSlot == "cellAdditionalexprs") {
    CYTdata@cellData@cellAdditionalexprs = CYTdata@cellData@cellAdditionalexprs %>% mutate_at(markers, trans)
  }

  CYTdata = checkValidity(CYTdata, mode = "warning")
  return(CYTdata)
}



# removeDuplicates <- function(CYTdata, markers = NULL) {
#
#   if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
#   else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
#
#   markers = checkorderMarkers(CYTdata, markers = markers, order = TRUE, checkDuplicates = TRUE)
#   if (length(markers)!=ncol(CYTdata@matrix.expression)) {
#     message("\nCells will be removed if their expression duplicated for combination of following markers : ", paste0(markers, collapse=", "))
#   }
#
#   idxDuplicates = duplicated(CYTdata@matrix.expression[,markers])
#   if (sum(idxDuplicates)==0) {
#     message("\nNo duplicated data found, identical CYTdata object returned")
#     return(CYTdata)
#   }
#   else {
#     message("\n", sum(idxDuplicates), " cell(s) is/are duplicated. Removing it from CYTdata object..")
#     gatedIdx = !idxDuplicates
#     newmatrix.expression = subset(CYTdata@matrix.expression, gatedIdx)
#
#     newsamples = CYTdata@samples[gatedIdx]
#     newsamples = droplevels(newsamples)
#     removedSpls = setdiff(levels(CYTdata@samples), levels(newsamples))
#     if(length(removedSpls)>0){
#       cat("\n\n - Gating operation remove following samples :", paste0(removedSpls, collapse = ", "))
#     }
#
#     cat("\n\n - Creating new CYTdata object..")
#     newCYTdata = methods::new("CYTdata",
#                               samples = newsamples,
#                               matrix.expression = newmatrix.expression)
#     if (ncol(CYTdata@raw.matrix.expression)>0) {
#       newCYTdata@raw.matrix.expression = subset(CYTdata@raw.matrix.expression, gatedIdx)
#     }
#
#     if (nrow(CYTdata@metadata)>0) {
#       newCYTdata@metadata = subset(CYTdata@metadata, rownames(CYTdata@metadata) %in% levels(newCYTdata@samples))
#     }
#
#     if(length(CYTdata@Clustering@clusters)>0){
#       cat("\n\n - Updating Clustering slot")
#       newClusters = CYTdata@Clustering@clusters[gatedIdx]
#       newClusters = droplevels(newClusters)
#       newpalette = CYTdata@Clustering@palette[levels(newClusters)]
#       cat("\ncomputing new cell cluster count, abundance matrix...")
#       newcellcount = compute.cellcount(newClusters, newsamples)
#       newabundance = compute.abundance(newcellcount)
#       cat("\nCreating new Clustering object")
#       newCYTdata@Clustering = methods::new("Clustering",
#                                            clusters = newClusters,
#                                            cellcount = newcellcount,
#                                            abundance = newabundance,
#                                            palette = newpalette,
#                                            optional_parameters = CYTdata@Clustering@optional_parameters)
#       message("\n\nRemark : Clustering results are preserved during gating operation but it is recommended to the user to run a
#         new clustering step with parameters adapted to the gated dataset")
#     }
#     if(length(CYTdata@Metaclustering@metaclusters)>0){
#       cat("\n\n - Updating Metaclustering slot")
#       newMetaclusters = CYTdata@Metaclustering@metaclusters[gatedIdx]
#       newMetaclusters = droplevels(newMetaclusters)
#       newpalette = CYTdata@Metaclustering@palette[levels(newMetaclusters)]
#       cat("\ncomputing new cell metacluster count, abundance matrix...")
#       newcellcount = compute.cellcount(newMetaclusters, newsamples)
#       newabundance = compute.abundance(newcellcount)
#       cat("\nCreating new Metaclustering object")
#       newCYTdata@Metaclustering = methods::new("Metaclustering",
#                                                metaclusters = newMetaclusters,
#                                                cellcount = newcellcount,
#                                                abundance = newabundance,
#                                                palette = newpalette,
#                                                optional_parameters = CYTdata@Metaclustering@optional_parameters)
#     }
#     if(nrow(CYTdata@DimReduction@coordinates)>0){
#       cat("\n\n - Updating DimReduction slot")
#       newcoordinates = subset(CYTdata@DimReduction@coordinates, gatedIdx)
#       newCYTdata@DimReduction = methods::new("DimReduction",
#                                              coordinates = newcoordinates,
#                                              optional_parameters = CYTdata@DimReduction@optional_parameters)
#     }
#     validObject(newCYTdata)
#     return(newCYTdata)
#   }
# }
