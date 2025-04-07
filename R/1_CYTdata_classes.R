
########################################################################################################################################################################################
#################################################################################### cellData ####################################################################################
########################################################################################################################################################################################

#' @title ---------
#'
#' @description ---------
#'
#' @slot ---------
#'
#' @name cellData-class
#' @rdname cellData-class
#' @exportClass cellData

cellData <- methods::setClass("cellData",
                              slots = c(cellExprs = "data.frame",
                                        cellDimRed = "list",
                                        cellAdditionalexprs = "data.frame"))


#' @title ---------
#'
#' @description ---------
#'
#' @param ---------
#'
#' @return ---------
#'
#'
#' @export
#'

import_cellData <- function(CYTdata, imported_cellData_Object, checkOverwrite=TRUE){
  CYTdata = checkValidity(CYTdata, mode = "error")

  if (class(imported_cellData_Object)!="cellData") { stop("Error : argument 'imported_cellData_Object' must be a S4 object of class 'cellData'.") }
  if (checkOverwrite){
    reply <- readline(prompt="cellData slot is not empty, do you still want to continue and overwrite (yes or no): ")
    while (!tolower(reply) %in% c("yes", "y", "no", "n")) { reply <- readline(prompt="Reply not standard, please answer yes or no: ") }
    if (tolower(reply) %in% c("no", "n")){ stop("Function stopped, CYTdata object unchanged") }
  }

  message("\nUpdating CYTdata object, cellData slot...")
  CYTdata@cellData = imported_cellData_Object
  CYTdata = checkValidity(CYTdata, mode = "warning")
  return(CYTdata)
}

#' @title ---------
#'
#' @description ---------
#'
#' @param ---------
#'
#' @return ---------
#'
#'
#' @export
#'

remove_cellData <- function(CYTdata){
  CYTdata = checkValidity(CYTdata, mode = "error")
  message("\nRemoving cellData slot in CYTdata object, but keeping cellExprs dataframe (can not be empty)...")
  CYTdata@cellData =methods::new("cellData", slots = c(cellData = CYTdata@cellData@cellExprs))
  CYTdata = checkValidity(CYTdata, mode = "warning")
  return(CYTdata)
}

#' @title ---------
#'
#' @description ---------
#'
#' @param ---------
#'
#' @return ---------
#'
#'
#' @export
#'

export_cellData <- function(CYTdata){
  CYTdata = checkValidity(CYTdata, mode = "error")
  return(CYTdata@cellData)
}


########################################################################################################################################################################################
#################################################################################### sampleData ####################################################################################
########################################################################################################################################################################################


#' @title ---------
#'
#' @description ---------
#'
#' @slot ---------
#'
#' @name sampleData-class
#' @rdname sampleData-class
#' @exportClass sampleData

sampleData <- methods::setClass("sampleData",
                                slots = c(cellSample = "data.frame",
                                          sampleMetadata = "data.frame",
                                          metadataPalette = "list",
                                          sampleAdditionaldata = "data.frame"))

#' @title ---------
#'
#' @description ---------
#'
#' @param ---------
#'
#' @return ---------
#'
#'
#' @export
#'

import_sampleData <- function(CYTdata, imported_sampleData_Object, checkOverwrite=TRUE){
  CYTdata = checkValidity(CYTdata, mode = "error")

  if (class(imported_sampleData_Object)!="sampleData") { stop("Error : argument 'imported_sampleData_Object' must be a S4 object of class 'sampleData'.") }
  if (checkOverwrite){
    reply <- readline(prompt="sampleData slot is not empty, do you still want to continue and overwrite (yes or no): ")
    while (!tolower(reply) %in% c("yes", "y", "no", "n")) { reply <- readline(prompt="Reply not standard, please answer yes or no: ") }
    if (tolower(reply) %in% c("no", "n")){ stop("Function stopped, CYTdata object unchanged") }
  }

  message("\nUpdating CYTdata object, sampleData slot...")
  CYTdata@sampleData = imported_sampleData_Object
  CYTdata = checkValidity(CYTdata, mode = "warning")
  return(CYTdata)
}

#' @title ---------
#'
#' @description ---------
#'
#' @param ---------
#'
#' @return ---------
#'
#'
#' @export
#'

remove_sampleData <- function(CYTdata){
  CYTdata = checkValidity(CYTdata, mode = "error")
  message("\nRemoving sampleData slot in CYTdata object, but keeping cellSample dataframe (can not be empty)...")
  CYTdata@sampleData =methods::new("sampleData", slots = c(sampleData = CYTdata@sampleData@cellSample))
  CYTdata = checkValidity(CYTdata, mode = "warning")
  return(CYTdata)
}

#' @title ---------
#'
#' @description ---------
#'
#' @param ---------
#'
#' @return ---------
#'
#' @export
#'

export_sampleData <- function(CYTdata){
  CYTdata = checkValidity(CYTdata, mode = "error")
  return(CYTdata@sampleData)
}

########################################################################################################################################################################################
#################################################################################### clusteringData ####################################################################################
########################################################################################################################################################################################

#' @title ---------
#'
#' @description ---------
#'
#' @slot ---------
#'
#' @name clusteringData-class
#' @rdname clusteringData-class
#' @exportClass clusteringData

clusteringData <- methods::setClass("clusteringData",
                                    slots = c(cellClustering = "data.frame",
                                              clusteringPalette = "list"))

#' @title ---------
#'
#' @description ---------
#'
#' @param ---------
#'
#' @return ---------
#'
#' @export
#'

import_clusteringData <- function(CYTdata, imported_clusteringData_Object, checkOverwrite=TRUE){
  CYTdata = checkValidity(CYTdata, mode = "error")

  if (class(imported_clusteringData_Object)!="clusteringData") { stop("Error : argument 'imported_clusteringData_Object' must be a S4 object of class 'clusteringData'.") }
  if (checkOverwrite && (nrow(CYTdata@clusteringData@cellClustering)>0 || ncol(CYTdata@clusteringData@cellClustering)>0)) {
    reply <- readline(prompt="clusteringData slot is not empty, do you still want to continue and overwrite (yes or no): ")
    while (!tolower(reply) %in% c("yes", "y", "no", "n")) { reply <- readline(prompt="Reply not standard, please answer yes or no: ") }
    if (tolower(reply) %in% c("no", "n")){ stop("Function stopped, CYTdata object unchanged") }
  }

  message("\nUpdating CYTdata object, clusteringData slot...")
  CYTdata@clusteringData = imported_clusteringData_Object
  CYTdata = checkValidity(CYTdata, mode = "warning")
  return(CYTdata)
}

#' @title ---------
#'
#' @description ---------
#'
#' @param ---------
#'
#' @return ---------
#'
#' @export
#'

remove_clusteringData <- function(CYTdata){
  CYTdata = checkValidity(CYTdata, mode = "error")
  message("\nRemoving clusteringData slot in CYTdata object...")
  CYTdata@clusteringData =methods::new("clusteringData")
  CYTdata = checkValidity(CYTdata, mode = "warning")
  return(CYTdata)
}

#' @title ---------
#'
#' @description ---------
#'
#' @param ---------
#'
#' @return ---------
#'
#' @export
#'

export_clusteringData <- function(CYTdata){
  CYTdata = checkValidity(CYTdata, mode = "error")
  return(CYTdata@clusteringData)
}


########################################################################################################################################################################################
############################################################## differentialTesting, Kinetic, TrajectoryInference #####################################################################
########################################################################################################################################################################################

#' @title ---------
#'
#' @description ---------
#'
#' @slot ---------
#'
#' @name differentialTesting-class
#' @rdname differentialTesting-class
#' @exportClass differentialTesting

differentialTesting <- methods::setClass("differentialTesting", slots = c(data = "data.frame"))

#' @title ---------
#'
#' @description ---------
#'
#' @slot ---------
#'
#' @name Kinetic-class
#' @rdname Kinetic-class
#' @exportClass Kinetic

Kinetic <- methods::setClass("Kinetic", slots = c(data = "list"))

#' @title ---------
#'
#' @description ---------
#'
#' @slot ---------
#'
#' @name TrajectoryInference-class
#' @rdname TrajectoryInference-class
#' @exportClass TrajectoryInference

TrajectoryInference <- methods::setClass("TrajectoryInference", slots = c(root.cells = "vector", leaf.cells = "vector",  network = "list", walk = "list", diff.traj = "list"))


########################################################################################################################################################################################
####################################################################################### CYTdata ########################################################################################
########################################################################################################################################################################################

CYTdata <- methods::setClass("CYTdata",
                             slots = c(cellData = "cellData",
                                       sampleData = "sampleData",
                                       clusteringData = "clusteringData",
                                       differentialTesting = "differentialTesting",
                                       Kinetic = "Kinetic",
                                       TrajectoryInference = "TrajectoryInference"))


#' @title ---------
#'
#' @description ---------
#'
#' @param ---------
#'
#' @return ---------
#'
#'
#' @export
#'

checkValidity = function(object, mode = c("error", "warning")) {

  areColors <- function(x) {
    sapply(x, function(X) {
      tryCatch(is.matrix(col2rgb(X)),
               error = function(e) FALSE)
    })
  }

  switch(mode,
         error = { printFunction <- function(text) { stop("Error ", text) } },
         warning = { printFunction <- function(text) { message("Warning ", text) } }
  )

  if (class(object)!="CYTdata") { paste("in CYTdata object : must be of class 'CYTdata'.")
  }
  if (class(object@cellData)!="cellData") {
    paste("in CYTdata object : cellData slot must be of class 'cellData'.")
  }
  if (class(object@sampleData)!="sampleData") {
    paste("in CYTdata object : sampleData slot must be of class 'sampleData'.")
  }
  if (class(object@clusteringData)!="clusteringData") {
    paste("in CYTdata object : clusteringData slot must be of class 'clusteringData'.")
  }
  if (class(object@differentialTesting)!="differentialTesting") {
    paste("in CYTdata object : differentialTesting slot must be of class 'differentialTesting'.")
  }
  if (class(object@Kinetic)!="Kinetic") {
    paste("in CYTdata object : Kinetic slot must be of class 'Kinetic'.")
  }
  if (class(object@TrajectoryInference)!="TrajectoryInference") {
    paste("in CYTdata object : TrajectoryInference slot must be of class 'TrajectoryInference'.")
  }

  ########################### Check cellData slot ###########################

  ###### Check cellExprs, markers
  if (ncol(object@cellData@cellExprs)==0) {
    paste("in CYTdata object, cellData slot : cellExprs dataframe is empty.") %>% printFunction
  }
  if (is.null(rownames(object@cellData@cellExprs))) {
    paste("in CYTdata object, cellData slot : events ID (rownames) are missing in cellExprs dataframe.") %>% printFunction
  }
  if(any(is.na(colnames(object@cellData@cellExprs)))){
    paste("in CYTdata object, cellData slot : Column names (markers) of cellExprs dataframe (", paste0(colnames(object@cellData@cellExprs), collapse = ","), ") contain NA values.") %>% printFunction
  }
  if(any(nchar(colnames(object@cellData@cellExprs))==0)){
    paste("in CYTdata object, cellData slot : Column names (markers) of cellExprs dataframe (", paste0(colnames(object@cellData@cellExprs), collapse = ","), ") contain empty values. Please rename.") %>% printFunction
  }
  if(any(is.na(object@cellData@cellExprs))){
    paste("in CYTdata object, cellData slot : cellExprs dataframe contains NA values.") %>% printFunction
  }
  doubleMarkers = duplicated(colnames(object@cellData@cellExprs))
  if(any(doubleMarkers)){
    paste("in CYTdata object, cellData slot : The column names (markers) of cellExprs dataframe (", paste0(colnames(object@cellData@cellExprs), collapse = ","),
         ") contain duplicated values (", paste0(colnames(object@cellData@cellExprs)[doubleMarkers], collapse = ","),").") %>% printFunction
  }

  classes = sapply(object@cellData@cellExprs, class)
  names(classes) = colnames(object@cellData@cellExprs)
  if(any(classes!="numeric")) {
    paste("in CYTdata object, cellData slot : cellExprs columns must be numeric [Current classes :", paste0(paste(names(classes), " (", classes, ")", sep = ""), collapse = ", "), "].") %>% printFunction
  }

  ###### Check cellDimRed
  if (length(object@cellData@cellDimRed)>0) {
    if (is.null(names(object@cellData@cellDimRed)) || any(is.na(names(object@cellData@cellDimRed)))) {
      paste("in CYTdata object, cellData slot : cellDimRed list name missing (cellDimRed names : ", paste0(names(object@cellData@cellDimRed), collapse=", "), ").") %>% printFunction
    }
    for (i in seq_along(object@cellData@cellDimRed)) {
      if (!identical(rownames(object@cellData@cellExprs), rownames(object@cellData@cellDimRed[[i]]))) {
        paste("in CYTdata object, cellData slot : Cell IDs (rownames) of the cellExprs dataframe are not identical (order included)
              with the Cell IDs (rownames) of the cellDimRed element named", names(object@cellData@cellDimRed)[i], ".") %>% printFunction
      }
    }
  }

  ###### Check cellAdditionalexprs
  if (!identical(rownames(object@cellData@cellExprs), rownames(object@cellData@cellAdditionalexprs))) {
    paste("in CYTdata object, cellData slot : Cell IDs (rownames) of the cellExprs dataframe are not identical (order included) with the Cell IDs (rownames) of the cellAdditionalexprs dataframe.
          They must be identical even if cellAdditionalexprs dataframe is empty, with zero columns.") %>% printFunction
  }
  if (ncol(object@cellData@cellAdditionalexprs)>0) {
    if (length(unique(colnames(object@cellData@cellAdditionalexprs)))!=ncol(object@cellData@cellAdditionalexprs)) {
      paste("in CYTdata object, cellData slot : Colnames of the cellAdditionalexprs are not unique. They must be unique because they represent different biological conditions") %>% printFunction
    }
    classes = sapply(object@cellData@cellAdditionalexprs, class)
    names(classes) = colnames(object@cellData@cellAdditionalexprs)
    if(any(classes!="numeric")) {
      paste("in CYTdata object, cellData slot : cellAdditionalexprs columns must be numeric [Current classes  :",
            paste0(paste(names(classes), " (", classes, ")", sep = ""), collapse = ", "), "]") %>% printFunction
    }
  }

  ########################### Check sampleData slot ###########################

  ###### Check cellSample
  if (ncol(object@sampleData@cellSample)==0 || nrow(object@sampleData@cellSample)==0) {
    paste("in CYTdata object, sampleData slot : cellSample dataframe is empty.") %>% printFunction
  }
  if (nrow(object@sampleData@cellSample)!=length(rownames(object@cellData@cellExprs))) {
    paste("in CYTdata object, sampleData slot : cellSample dataframe doesn't have the same number of rows (",
         nrow(object@sampleData@cellSample), ") as the number of cells (", length(rownames(object@cellData@cellExprs)), ").") %>% printFunction
  }
  if (!identical(colnames(object@sampleData@cellSample), c("SampleID", "CellID"))) {
    paste("in CYTdata object, sampleData slot : cellSample column names must be 'SampleID' and 'CellID' (in this order) but are",
         paste0(colnames(object@sampleData@cellSample), collapse=", "), ".") %>% printFunction
  }
  if (!identical(rownames(object@cellData@cellExprs), rownames(object@sampleData@cellSample))) {
    paste("in CYTdata object, sampleData slot : events IDs (rownames) of cellExprs dataframe are not identical (order included) with the events IDs (rownames) of the cellSample dataframe.") %>% printFunction
  }


  if(class(object@sampleData@cellSample$SampleID)!="factor") {
    paste("in CYTdata object, sampleData slot : SampleID column of cellSample dataframe must be factors (currently", class(object@sampleData@cellSample$SampleID), ").") %>% printFunction
  }
  if(any(is.na(object@sampleData@cellSample$SampleID))){
    paste("in CYTdata object, sampleData slot : SampleID column of cellSample dataframe contains NA values.") %>% printFunction
  }
  errLev = setdiff(levels(object@sampleData@cellSample$SampleID), unique(object@sampleData@cellSample$SampleID))
  if (length(errLev)>0) {
    paste("in CYTdata object, sampleData slot : SampleIDcolumn of cellSample dataframe is factor but contains levels not present in the vector (",
          paste0(errLev, collapse = ", "), "). Please drop absent levels using droplevels function.") %>% printFunction
  }

  if(class(object@sampleData@cellSample$CellID)!="character") {
    paste("in CYTdata object, sampleData slot : CellID column of cellSample dataframe must be characters (currently", class(object@sampleData@cellSample$CellID), ").") %>% printFunction
  }
  if(any(is.na(object@sampleData@cellSample$CellID))){
    paste("in CYTdata object, sampleData slot : CellID column of cellSample dataframe contains NA values.") %>% printFunction
  }
  errLev = setdiff(levels(object@sampleData@cellSample$CellID), unique(object@sampleData@cellSample$CellID))
  if (length(errLev)>0) {
    paste("in CYTdata object, sampleData slot : CellIDcolumn of cellSample dataframe is character but contains levels not present in the vector (",
          paste0(errLev, collapse = ", "), "). Please drop absent levels using droplevels function.") %>% printFunction
  }


  ###### Check sampleMetadata, metadataPalette

  if (nrow(object@sampleData@sampleMetadata)>0){
    if (!identical(levels(object@sampleData@cellSample$SampleID), rownames(object@sampleData@sampleMetadata))){ # levels and rownames are unique
      paste("in CYTdata object, sampleData slot : Sample IDs (rownames) of the sampleMetadata dataframe (",
           paste0(rownames(object@sampleData@sampleMetadata), collapse = ", "),
           ") are not identical (order included) with the levels of cellSample's SampleID column (",
           paste0(levels(object@sampleData@cellSample$SampleID), collapse = ", "), ")") %>% printFunction
    }

    if (ncol(object@sampleData@sampleMetadata)!=length(object@sampleData@metadataPalette)) {
      paste("in CYTdata object, sampleData slot : sampleMetadata dataframe doesn't have the same number of metadata columns (",
           ncol(object@sampleData@sampleMetadata), ") as the length of metadataPalette list (", length(object@sampleData@metadataPalette), ").") %>% printFunction
    }
    if (!identical(colnames(object@sampleData@sampleMetadata), names(object@sampleData@metadataPalette))) {
      paste("in CYTdata object, sampleData slot : sampleMetadata dataframe doesn't have the same metadata columns (",
           paste0(colnames(object@sampleData@sampleMetadata), collapse = ","), ") (order included) as the names of metadataPalette list (",
           paste0(names(object@sampleData@metadataPalette), collapse = ","), ").") %>% printFunction
    }
    iscolor = areColors(object@sampleData@metadataPalette)
    if(!all(iscolor)){
      paste("in CYTdata object, sampleData slot : Following elements of metadataPalette list (", paste0(names(iscolor==FALSE), collapse = ","), "), does not contain hexadecimal color.)") %>% printFunction
    }

    classes = sapply(object@sampleData@sampleMetadata, class)
    names(classes) = colnames(object@sampleData@sampleMetadata)
    if(any(classes!="factor")) {
      paste("in CYTdata object, sampleData slot : sampleMetadata dataframe contains non factor metadata columns [Current classes  :",
           paste0(paste(names(classes), " (", classes, ")", sep = ""), collapse = ", "), "]") %>% printFunction
    }
    mdErr = colnames(object@sampleData@sampleMetadata)[duplicated(colnames(object@sampleData@sampleMetadata))]
    if (length(mdErr)>0) {
      paste("in CYTdata object, sampleData slot : sampleMetadata dataframe must have unique metadata column names, but", paste0(mdErr, collapse=", "), " are the names of several metadata columns") %>% printFunction
    }
    if (!all(c("Timepoint", "Individual") %in% colnames(object@sampleData@sampleMetadata))) {
      paste("in CYTdata object, sampleData slot : sampleMetadata dataframe must include 'Timepoint' and 'Individual' metadata columns (Current metadata columns :",
            paste0(colnames(object@sampleData@sampleMetadata), collapse = ","), ").") %>% printFunction
    }


    for (col in colnames(object@sampleData@sampleMetadata)) {
      column = object@sampleData@sampleMetadata[,col]
      if(any(is.na(column))){
        paste("in CYTdata object, sampleData slot : ", col, " metadata column of sampleMetadata dataframe contains NA values.") %>% printFunction
      }
      errLev = setdiff(levels(column), unique(column))
      if (length(errLev)>0) {
        paste("in CYTdata object, sampleData slot : ", col, " metadata column of sampleMetadata dataframe is factor but contains levels not present in the values (",
             paste0(errLev, collapse = ", "), "). Please drop absent levels using droplevels function.") %>% printFunction
      }
      if (!identical(levels(column), names(object@sampleData@metadataPalette[[col]]))){
        paste("in CYTdata object, sampleData slot : ", col, "metadata column of sampleMetadata dataframe has values (", paste0(levels(column), collapse = ","),
             ") not identical (order included) to the metadata values stored in correspond palette in metadataPalette list (",
             paste0(names(object@sampleData@metadataPalette[[col]]), collapse = ","), ")") %>% printFunction
      }
    }
  }

  ###### Check sampleAdditionaldata

  if (nrow(object@sampleData@sampleAdditionaldata)>0){
    if (!identical(levels(object@sampleData@cellSample$SampleID), rownames(object@sampleData@sampleAdditionaldata))){ # levels and rownames are unique
      paste("in CYTdata object, sampleData slot : Sample names (rownames) of the sampleAdditionaldata dataframe (",
           paste0(rownames(object@sampleData@sampleAdditionaldata), collapse = ", "),
           ") are not identical (order included) with the levels of cellSample's SampleID column (",
           paste0(levels(object@sampleData@cellSample$SampleID), collapse = ", "), ")") %>% printFunction
    }
    addErr = colnames(object@sampleData@sampleAdditionaldata)[duplicated(colnames(object@sampleData@sampleAdditionaldata))]
    if (length(addErr)>0) {
      paste("in CYTdata object, sampleData slot : sampleAdditionaldata dataframe must have unique column names, but", paste0(addErr, collapse=", "), " are the names of several columns") %>% printFunction
    }
  }


  ########################### Check clusteringData slot ###########################


  ###### Check cellClustering
  if (nrow(object@clusteringData@cellClustering)!=length(rownames(object@cellData@cellExprs))) {
    paste("in CYTdata object, clusteringData slot : cellClustering dataframe doesn't have the same number of rows (",
          nrow(object@clusteringData@cellClustering), ") as the number of cells (", length(rownames(object@cellData@cellExprs)), ").
          They must be the same even if cellClustering dataframe is empty, with zero columns (default setting)") %>% printFunction
  }
  if (!identical(rownames(object@cellData@cellExprs), rownames(object@clusteringData@cellClustering))) {
    paste("in CYTdata object, clusteringData slot : Event IDs (rownames) of the cellExprs dataframe are not identical (order included) with the Event IDs (rownames) of the cellClustering dataframe.
          They must be identical even if cellClustering dataframe is empty, with zero columns (default setting)") %>% printFunction
  }

  ###### Check cellClustering, clusteringPalette
  if (ncol(object@clusteringData@cellClustering)>0) {

    popErr = colnames(object@clusteringData@cellClustering)[duplicated(colnames(object@clusteringData@cellClustering))]
    if (length(mdErr)>0) {
      paste("in CYTdata object, clusteringData slot : cellClustering dataframe must have unique clustering column names, but", paste0(popErr, collapse=", "), " are the names of several clustering columns") %>% printFunction
    }

    if (ncol(object@clusteringData@cellClustering)!=length(object@clusteringData@clusteringPalette)) {
      paste("in CYTdata object, clusteringData slot : cellClustering dataframe doesn't have the same number of clustering columns (",
           ncol(object@clusteringData@cellClustering), ") as the length of clusteringPalette list (", length(object@clusteringData@clusteringPalette), ").") %>% printFunction
    }
    if (!identical(colnames(object@clusteringData@cellClustering), names(object@clusteringData@clusteringPalette))) {
      paste("in CYTdata object, clusteringData slot : cellClustering dataframe doesn't have the same clustering columns (",
           paste0(colnames(object@clusteringData@cellClustering), collapse = ","), ") (order included) as the names of clusteringPalette list (",
           paste0(names(object@clusteringData@clusteringPalette), collapse = ","), ").") %>% printFunction
    }
    iscolor = areColors(object@clusteringData@clusteringPalette)
    if(!all(iscolor)){
      paste("in CYTdata object, clusteringData slot :  Following elements of clusteringPalette list (", paste0(names(iscolor==FALSE), collapse = ","), "), does not contain hexadecimal color.)") %>% printFunction
    }


    classes = sapply(object@clusteringData@cellClustering, class)
    names(classes) = colnames(object@clusteringData@cellClustering)
    if(any(classes!="factor")) {
      paste("in CYTdata object, clusteringData slot : cellClustering dataframe contains non factor metadata column [Current classes  :",
           paste0(paste(names(classes), "(", classes, ")", sep = ""),collapse = ", "), "]") %>% printFunction
    }
    for (col in colnames(object@clusteringData@cellClustering)) {
      column = object@clusteringData@cellClustering[,col]
      if(any(is.na(column))){
        paste("in CYTdata object, clusteringData slot :", col, "clustering column of cellClustering dataframe contain NA values.") %>% printFunction
      }
      errLev = setdiff(levels(column), unique(column))
      if (length(errLev)>0) {
        paste("in CYTdata object, clusteringData slot :", col, "clustering column of cellClustering dataframe is factor but contains levels not present in the value (",
             paste0(errLev, collapse = ", "), "). Please drop absent levels using droplevels function.") %>% printFunction
      }
      if (!identical(levels(column), names(object@clusteringData@clusteringPalette[[col]]))){
        paste("in CYTdata object, clusteringData slot :", col, "clustering column of cellClustering dataframe has values (", paste0(levels(column), collapse = ","),
             ") not identical (order included) to the clustering values stored in correspond palette stored in clusteringPalette list (",
             paste0(names(object@clusteringData@clusteringPalette[[col]]), collapse = ","), ")") %>% printFunction
      }
    }
  } else {
    if (length(object@clusteringData@clusteringPalette)>0) { paste("in CYTdata object, clusteringData slot : clusteringPalette list is given but cellClustering dataframe is empty.") %>% printFunction }
  }
  return(object)
}
