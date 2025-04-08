#' @title Rename Markers in CYTdata Object
#'
#' @description
#' This function renames the markers (channels) in a CYTdata object. It allows users to specify a set of markers to be renamed and the corresponding new names.
#' The function can also remove conjugate labels (e.g., parts of marker names separated by an underscore) before renaming.
#'
#' @param CYTdata A CYTdata object containing cytometry data. The object must be valid before and after renaming markers.
#' @param cellSlot A string specifying which data slot to rename markers in. It can be either:
#'  - "cellExprs" (default): The main expression data for the cells.
#'  - "cellAdditionalexprs": Additional expression data associated with the cells.
#' @param from A character vector of marker (channel) names to be renamed. The length of `from` should match the length of `to` (if not removing conjugates).
#' @param to A character vector of new marker (channel) names corresponding to `from` (element-wise). The length of `to` must match the length of `from`.
#' @param removeConjugate A boolean indicating whether to remove conjugate labels (e.g., part of the name before an underscore). Default is `FALSE`.
#' @param sepConjugate A string defining the separator between the conjugate part and the main part of the marker name. Default is "_".
#'
#' @return A modified CYTdata object with the renamed markers in the specified `cellSlot`.
#' If `removeConjugate` is `TRUE`, the conjugate part of the marker names will be removed before renaming.
#'
#' @details
#' The function renames the markers in the specified data slot (`cellExprs` or `cellAdditionalexprs`) by matching the marker names provided in the `from` argument
#' and replacing them with the corresponding names in the `to` argument. If the `removeConjugate` argument is set to `TRUE`, the function will ignore the `to` argument
#' and remove the conjugate part of the marker names (i.e., everything before the specified separator, typically an underscore "_").
#'
#' If the lengths of `from` and `to` are mismatched (when `removeConjugate` is `FALSE`), the function will raise an error. The function also validates the `CYTdata` object
#' before and after renaming to ensure its integrity.
#'
#' @examples
#' # Example: Rename markers "FSC-A" and "SSC-A" to "FSC-A-New" and "SSC-A-New" in the 'cellExprs' slot
#' renamed_data <- renameMarkers(CYTdata = cyt_data, markers = c("FSC-A", "SSC-A"),
#'                               cellSlot = "cellExprs", from = c("FSC-A", "SSC-A"), to = c("FSC-A-New", "SSC-A-New"))
#'
#' # Example: Remove conjugate labels from marker names in the 'cellExprs' slot
#' renamed_data <- renameMarkers(CYTdata = cyt_data, markers = c("CD3_A", "CD4_A"),
#'                               cellSlot = "cellExprs", from = c("CD3_A", "CD4_A"), removeConjugate = TRUE)
#'
#' @seealso
#' \code{\link{checkValidity}} for validating the CYTdata object.
#'
#' @import checkmate
#' @import dplyr
#' @import stringr
#' @export

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


#' @title Transform CYTdata Object
#'
#' @description
#' This function transforms the expression data in a CYTdata object using a specified transformation mode (e.g., arcsinh, logicle, or biexponential).
#' The transformation is applied to the specified markers (channels) in either the `cellExprs` or `cellAdditionalexprs` slot of the CYTdata object.
#' The function allows the user to define transformation parameters such as scaling and offsets.
#'
#' @param CYTdata A CYTdata object containing cytometry data. The object must be valid before and after transformation.
#' @param markers A character vector of marker (channel) names to apply the transformation to. If NULL, all markers will be transformed.
#' @param mode A string specifying the transformation mode. It can be one of the following:
#'  - "arcsinh" (default): Applies the arcsinh transformation to the data.
#'  - "logicle": Placeholder for future implementation of the logicle transformation.
#'  - "biexponential": Applies the biexponential transformation to the data.
#' @param cellSlot A string specifying which data slot to transform. It can be either:
#'  - "cellExprs" (default): The main expression data for the cells.
#'  - "cellAdditionalexprs": Additional expression data associated with the cells.
#' @param a A numeric parameter for the transformation. Default is 0.
#' @param b A numeric parameter for the transformation. Default is 0.2.
#' @param c A numeric parameter for the transformation. Default is 0.
#' @param d A numeric parameter for the transformation. Default is 1.
#' @param f A numeric parameter for the transformation. Default is 0.
#' @param w A numeric parameter for the transformation. Default is 0.
#'
#' @return A transformed CYTdata object with the expression data in the specified `cellSlot` transformed according to the chosen mode and parameters.
#' The resulting object contains the transformed expression data.
#'
#' @details
#' The function applies different transformations to the data depending on the specified `mode`:
#' - **"arcsinh"**: Applies the inverse hyperbolic sine transformation, defined as `arcsinh(a + b * x) + c`, where `x` is the expression value and `a`, `b`, `c` are parameters.
#' - **"logicle"**: A placeholder for the logicle transformation (currently not implemented).
#' - **"biexponential"**: Applies a biexponential transformation defined as `a * exp(b * (x - w)) - c * exp(-d * (x - w)) + f`, where `x` is the expression value and `a`, `b`, `c`, `d`, `f`, and `w` are parameters.
#'
#' The function uses `mutate_at` from `dplyr` to apply the transformation to the specified markers.
#'
#' Note that the function performs validation on the `CYTdata` object before and after transformation using the `checkValidity` function.
#' This ensures that the object is valid before processing and that any issues are flagged after the transformation.
#'
#' @examples
#' # Example: Apply arcsinh transformation to markers "FSC-A" and "SSC-A" in the 'cellExprs' slot
#' transformed_data <- transformCYTdata(CYTdata = cyt_data, markers = c("FSC-A", "SSC-A"),
#'                                      mode = "arcsinh", cellSlot = "cellExprs", a = 0, b = 0.2, c = 0)
#'
#' # Example: Apply biexponential transformation with custom parameters
#' transformed_data <- transformCYTdata(CYTdata = cyt_data, markers = c("CD3", "CD4"),
#'                                      mode = "biexponential", cellSlot = "cellExprs", a = 1, b = 0.2, c = 0, d = 1, f = 0, w = 0)
#'
#' @seealso
#' \code{\link{checkValidity}} for validating the CYTdata object.
#'
#' @import checkmate
#' @import dplyr
#' @importFrom methods new
#' @export

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

######################################################## Utils ########################################################

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
