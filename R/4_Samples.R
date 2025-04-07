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

checkorderSamples <- function(CYTdata, samples, order = TRUE, checkDuplicates = TRUE) {
  #CYTdata = checkValidity(CYTdata, mode = "error")
  checkmate::qassert(samples, c("0","S*"))
  checkmate::qassert(order, "B1")
  checkmate::qassert(checkDuplicates, "B1")
  samplesId = levels(CYTdata@sampleData@cellSample$SampleID)
  if (is.null(samples)) {
    samples = samplesId
  } else {
    if (length(samples)==0) { stop("Error : 'samples' argument is an empty vector (length=0, but non NULL).") }
    if (checkDuplicates) {
      spddup = unique(samples[duplicated(samples)])
      if (length(spddup)>0) { stop("Error : 'samples' argument contain duplicated values ( ", paste0(spddup, collapse = ", "), " ). It must be vector of unique samples levels.") }
    }
    splErr = setdiff(samples, samplesId)
    if (length(splErr)>0) { stop("Error in CYTdata object, sampleData slot: 'samples' argument providing samples", paste0(splErr, collapse=", "), "not present in cellSample's SampleID column") }
    if (order) { samples = samplesId[samplesId %in% samples] }
  }
  return(samples)
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

updateMetadataLevels <- function(CYTdata, andPalette = TRUE) {
  mdFrame = CYTdata@sampleData@sampleMetadata %>% mutate(across(where(is.factor), droplevels))
  palList = CYTdata@sampleData@metadataPalette
  for (md in colnames(mdFrame)) { palList[[md]] = palList[[md]][levels(mdFrame[,md])] }
  CYTdata@sampleData@sampleMetadata = mdFrame
  CYTdata@sampleData@metadataPalette = palList[colnames(mdFrame)]
  return(CYTdata)
}


#' @title Renames samples within a CYTdata object
#'
#' @description This function aims to rename samples stored within a CYTdata object.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param from a character vector providing the sample names to replace. By default, all sample names are replaced
#' @param to a character vector providing the new sample names to use
#'
#' @return a S4 object of class 'CYTdata'
#'
#'@export
#'

renameSamples <- function(CYTdata, merge = FALSE, from = NULL, to){
  CYTdata = checkValidity(CYTdata, mode = "error")
  from = checkorderSamples(CYTdata, samples = from, order=FALSE, checkDuplicates=TRUE)
  checkmate::qassert(to, "S*")
  checkmate::qassert(merge, "B1")

  if (length(to) == 1) { to = rep(to, length(from)) }
  if (length(from)!=length(to)) { stop("Error : Length of argument 'from' (", length(from), ") and argument 'to' (", length(to), ") must be equal") }

  oldspllev = levels(CYTdata@sampleData@cellSample$SampleID)
  newspllev = oldspllev
  newspllev[match(from, newspllev)] = to

  levdup = unique(newspllev[duplicated(newspllev)])
  if (length(levdup)>0){
    message("After renaming, several samples have the same name (",  paste0(levdup, collapse=", "), "), either by renaming to an already existing and unchanged sample name, or by duplicate in the 'to' argument, or both.")
    if (merge) { message("\nWarning : 'merge' argument is set to TRUE. After renaming, sampleID column of cellSample dataframe get their duplicated levels merged.") }
    else { stop("\nError : Common names between samples, after renaming, but 'merge' argument is set to FALSE. CYTdata unchanged.") }
  }

  cat("\n\nCurrent sample names are (in the order) :", paste0(oldspllev, collapse = ", "), "\n")
  cat("\n\nThe following samples :")
  cat("\n - ", paste0(from, collapse=", "))
  cat("\n\nwill be renamed, in the order, by :")
  cat("\n - ", paste0(to, collapse=", "))
  CYTdata@sampleData@cellSample$SampleID = plyr::mapvalues(CYTdata@sampleData@cellSample$SampleID, from = from, to = to)

  oldmd = CYTdata@sampleData@sampleMetadata
  if (nrow(oldmd)>0){
    message("After renaming, sampleID column of cellSample dataframe contained duplicated levels, which were merged. As it concerns metadata : \n
              - If the name of merged samples was already associated to an existing sample name before renaming, the metadata of the existing sample are the new metadata of newly renamed sample(s).\n
              - If the name of merged samples was a new one (argument 'to' contain this name several times). The metadata of the resulting sample is the one of the sample contained in 'from' argument
                and which was ordered first (in 'from' argument) among all the samples renamed to this name.")
    newmd = data.frame()
    for (spl in levels(CYTdata@sampleData@cellSample$SampleID)) {
      fromspl = oldspllev[newspllev == spl]
      if (length(fromspl)>1){ # if duplicated
        message("\n\n\nWarning : ", paste0(fromspl, collapse=", "), " samples have been merged into ", spl, " sample.")
        if (spl %in% fromspl) {
          message("Warning :", spl, "sample original metadata has been conserved in sample Metadata frame.")
          metadataRow = oldmd[spl,]
        } # sample name already present
        else {
          spl2 = from[to==spl][1]
          message("Warning :", spl, "sample wasn't already present so no original metadata has been conserved. The metadata of the first sample, in the 'from' argument renamed to", spl, " (", spl2, ") has been conferred to", spl, ".")
          metadataRow = oldmd[spl2,]
        } # sample name not present, first sample
      }
      else {
        if (spl %in% to) {
          message("\n\n\nWarning : ", from[to==spl], " sample has been renamed to ", spl, ", but no merging happened. So ", from[to==spl], " metadata values have been given to new sample ", spl, ".")
          metadataRow = oldmd[from[to==spl],]
        }
        else { metadataRow = oldmd[spl,] }
      }
      newmd = rbind.data.frame(newmd, metadataRow)
    }
    message("Updating sampleMetadata dataframe levels (using droplevels) and updating metadataPalette list (if droplevels used)..")
    rownames(newmd) = levels(CYTdata@sampleData@cellSample$SampleID)
    CYTdata@sampleData@sampleMetadata = newmd
    CYTdata = updateMetadataLevels(CYTdata)
  }
  CYTdata = checkValidity(CYTdata, mode = "warning")
  return(CYTdata)
}

#' @title Renames samples within a CYTdata object
#'
#' @description This function aims to rename samples stored within a CYTdata object.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param from a character vector providing the sample names to replace. By default, all sample names are replaced
#' @param to a character vector providing the new sample names to use
#'
#' @return a S4 object of class 'CYTdata'
#'
#'@export
#'
#'

reorderSamples <- function(CYTdata, alphabetic = TRUE, newOrder = NULL){
  CYTdata = checkValidity(CYTdata, mode = "error")
  checkmate::qassert(alphabetic, "B1")
  if (alphabetic) {
    CYTdata@sampleData@cellSample$SampleID = factor(CYTdata@sampleData@cellSample$SampleID, levels = gtools::mixedsort(levels(CYTdata@sampleData@cellSample$SampleID)))
  } else {
    if (is.null(newOrder)) { stop("Error : 'alphabetic' argument is set to FALSE but 'newOrder' argument is NULL.") }
    newOrder = checkorderSamples(CYTdata, samples = newOrder, order = FALSE, checkDuplicates = TRUE)
    CYTdata@sampleData@cellSample$SampleID = factor(CYTdata@sampleData@cellSample$SampleID, levels = newOrder)
  }
  CYTdata@sampleData@sampleMetadata = CYTdata@sampleData@sampleMetadata[levels(CYTdata@sampleData@cellSample$SampleID),]
  CYTdata = checkValidity(CYTdata, mode = "warning")
}

#' @title Numbers of cells for each sample
#'
#' @description This function aims to visualize the number of cells associated to each sample.This representation displays the samples in the X-axis and the number of associated cells in the Y-axis.
#' Several statistics can be computed and shown.
#'
#' @details
#' The following statistic can be computed:
#' -'min' corresponds to the lowest number of cells within a data set
#' -'q25' corresponds to the number of cells separates the quantiles 25% within data set
#' -'median' corresponds to the number of cells separates the lower half from the upper half within data set
#' -'mean' corresponds to the number of cells quantity shared within data set
#' -'q75' corresponds to the number of cells separates the quantiles 75% within data set
#' -'max' corresponds to the largest number of cells within a data set
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param stats a character vector providing the statistics to display. Possible values are: 'min', 'median', 'mean', 'q75', 'max'
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param sort a boolean value indicating if samples must be sorted by the number of cells
#'
#' @return a ggplot2 object
#'
#' @export
#'

samplesCounts <- function(CYTdata,
                          samples = NULL,
                          sort = TRUE) {
  CYTdata = checkValidity(CYTdata, mode = "error")
  checkmate::qassert(sort, "B1")
  samples = checkorderSamples(CYTdata, samples = samples, order = FALSE, checkDuplicates = TRUE)
  data = CYTdata@sampleData@cellSample
  if (!is.null(samples)) { data = subset(data, SampleID %in% samples)}
  nbCells = plyr::ddply(data, "SampleID", nrow)
  colnames(nbCells) = c("SampleID", "counts")
  if (sort == TRUE) { nbCells = nbCells[order(nbCells$counts, decreasing = TRUE),] }
  nbCells  = nbCells %>% mutate(SampleID = factor(SampleID))
  mean = mean(nbCells$counts)
  median = stats::median(nbCells$counts)
  min = min(nbCells$counts)
  max = max(nbCells$counts)
  q25 = stats::quantile(nbCells$counts, probs = 0.25)
  q75 = stats::quantile(nbCells$counts, probs = 0.75)
  xscale.middle <- nrow(nbCells)/2
  if (xscale.middle%%2==0) { xscale.middle <- xscale.middle + 0.5}
  plot <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = nbCells, ggplot2::aes_string(x="samples", xend="samples", yend="counts"), y=0, color = "gray50") +
    ggplot2::geom_point(data = nbCells, ggplot2::aes_string(x="samples", y="counts"), size=2) +
    ggplot2::geom_hline(yintercept = mean, color="green") +
    ggplot2::geom_text(data = nbCells, ggplot2::aes_string(x="xscale.middle", y="mean"), label=paste0("mean: ", format(mean, scientific = FALSE, big.mark= "," )),color="green", vjust=0) +
    ggplot2::geom_hline(yintercept = min, color="blue") +
    ggplot2::geom_text(data = nbCells, ggplot2::aes_string(x="xscale.middle", y="min"), label=paste0("min: ", format(min, scientific = FALSE, big.mark = ",")), color="blue", vjust=0) +
    ggplot2::geom_hline(yintercept = q25, color="purple") +
    ggplot2::geom_text(data = nbCells, ggplot2::aes_string(x="xscale.middle", y="q25"), label=paste0("q25: ", format(q25, scientific = FALSE, big.mark = ",")), color="purple", vjust=0) +
    ggplot2::geom_hline(yintercept = median, color="red") +
    ggplot2::geom_text(data = nbCells, ggplot2::aes_string(x="xscale.middle", y="median"),label=paste0("median: ", format(median, scientific = FALSE, big.mark = ",")), color="red", vjust=0) +
    ggplot2::geom_hline(yintercept = q75, color="purple") +
    ggplot2::geom_text(data = nbCells, ggplot2::aes_string(x="xscale.middle", y="q75"), label=paste0("q75: ", format(q75, scientific = FALSE, big.mark = ",")), color="purple", vjust=0) +
    ggplot2::geom_hline(yintercept = max, color="blue") +
    ggplot2::geom_text(data = nbCells, ggplot2::aes_string(x="xscale.middle", y="max"), label=paste0("max: ", format(max, scientific = FALSE, big.mark = ",")), color="blue", vjust=0) +
    ggplot2::scale_y_continuous(labels = function(x){format(x, scientific = FALSE, big.mark=",")}) +
    ggplot2::ylab("Number of cells") +
    ggplot2::xlab("Samples") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic"),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15.5),
                   axis.title.x = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill=NA),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   legend.position = "none")
  return(plot)
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
#' @param metadata a dataframe containing meta-information about the biological samples.
#' The columns must contain, at least, a column named "Timepoint" and an other named "Individual".
#' The rownames have to be the biological samples, thus the number rows has to be equal to the number of samples.
#'
#' @return a S4 object of class 'CYTdata'
#'
#'
#' @export
#'

importsampleMetadata <- function(CYTdata, sampleMetadata_df, optionnalPalette = NULL, checkOverwrite = TRUE){
  CYTdata = checkValidity(CYTdata, mode = "error")
  checkmate::qassert(checkOverwrite, "B1")
  checkmate::qassert(sampleMetadata_df, "D*")
  if (checkOverwrite && ncol(CYTdata@sampleData@sampleMetadata)!=0) {
    reply = readline(prompt="sampleMetadata already present, do you still want to continue and overwrite (yes or no): ")
    while (!tolower(reply) %in% c("yes", "y", "no", "n")) { reply <- readline(prompt="Reply not standard, please answer yes or no: ") }
    if (tolower(reply) %in% c("no", "n")){ stop("Function stopped, CYTdata object unchanged") }
  }
  cat("\n\nConverting sampleMetadata's metadata columns to factor..")
  sampleMetadata_df = sampleMetadata_df %>% dplyr::mutate_if(function(x) { return(!is.factor(x)) }, factor)
  cat("\n\nUpdating sampleMetadata dataframe...")
  CYTdata@sampleData@sampleMetadata = sampleMetadata_df
  cat("\n\nUpdating metadataPalette list...")
  if(is.null(optionnalPalette)) {
    metadataPalette = sampleMetadata_df %>% lapply(function(col) { structure(rainbow(nlevels(col)), names = levels(col)) })
    names(metadataPalette) = colnames(sampleMetadata_df)
    CYTdata@sampleData@metadataPalette = metadataPalette
  } else {
    CYTdata@sampleData@metadataPalette = optionnalPalette
  }
  CYTdata = checkValidity(CYTdata, mode = "error")
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
#' @param metadata a dataframe containing meta-information about the biological samples.
#' The columns must contain, at least, a column named "Timepoint" and an other named "Individual".
#' The rownames have to be the biological samples, thus the number rows has to be equal to the number of samples.
#'
#' @return a S4 object of class 'CYTdata'
#'
#'
#' @export
#'

importsampleAdditionaldata <- function(CYTdata, sampleAdditionaldata_df, checkOverwrite = TRUE){
  CYTdata = checkValidity(CYTdata, mode = "error")
  checkmate::qassert(checkOverwrite, "B1")
  checkmate::qassert(sampleAdditionaldata_df, "D*")
  if (checkOverwrite && ncol(CYTdata@sampleData@sampleAdditionaldata)!=0) {
    reply = readline(prompt="sampleAdditionaldata already present, do you still want to continue and overwrite (yes or no): ")
    while (!tolower(reply) %in% c("yes", "y", "no", "n")) { reply <- readline(prompt="Reply not standard, please answer yes or no: ") }
    if (tolower(reply) %in% c("no", "n")){ stop("Function stopped, CYTdata object unchanged") }
  }
  cat("\n\nUpdating sampleAdditionaldata dataframe...")
  CYTdata@sampleData@sampleAdditionaldata = sampleAdditionaldata_df
  CYTdata = checkValidity(CYTdata, mode = "error")
  return(CYTdata)
}


#' #' @title Assigns meta-information about biological samples
#' #'
#' #' @description This function aims to attach meta-information to each biological sample.
#' #'
#' #' Especially, the following meta-information of each sample can be specified for subsequent analyses.
#' #' - The biological individual
#' #' - The biological condition (groups, vaccines, studies, etc.)
#' #' - The timepoint
#' #' Timepoint and Individual data must be specified.
#' #'
#' #' @param CYTdata a CYTdata object
#' #' @param metadata a dataframe containing meta-information about the biological samples.
#' #' The columns must contain, at least, a column named "Timepoint" and an other named "Individual".
#' #' The rownames have to be the biological samples, thus the number rows has to be equal to the number of samples.
#' #'
#' #' @return a S4 object of class 'CYTdata'
#' #'
#' #'
#' #' @export
#' #'
#'
#' renameMetadata <- function(CYTdata, metadata = NULL, from = NULL, to){
#'   CYTdata = checkValidity(CYTdata, mode = "error")
#'   if (nrow(CYTdata@sampleData@sampleMetadata)==0 && ncol(CYTdata@sampleData@sampleMetadata)==0) {
#'     stop("Error : sampleData's sampleMetadata slot is empty")
#'   }
#'
#'   checkmate::qassert(metadata, c(0, "S*"))
#'   checkmate::qassert(to, "S*")
#'   checkmate::qassert(merge, "B1")
#'
#'   if (is.null(metadata)) { # Rename columns
#'     if (is.null(from)) { from = colnames(CYTdata@sampleData@sampleMetadata) }
#'     else {
#'       mdErr = setdiff(from, colnames(object@sampleData@sampleMetadata))
#'       if (length(mdErr)>0) { stop("Error : 'from' argument contains names not present among sampleMetadata's metadata column names (", paste0(mdErr, collapse=", "), ".") }
#'     }
#'
#'     if (length(from)!=length(to)) { stop("Error : Length of argument 'from' (", length(from), ") and argument 'to' (", length(to), ") must be equal") }
#'     dupfrom = unique(from[duplicated(from)])
#'     if (length(dupfrom)>0) { stop("Error : Argument 'from' contains duplicated metadata column names :", paste0(dupfrom, collapse=", ")) }
#'     dupto = unique(to[duplicated(to)])
#'     if (length(dupto)>0) { stop("Error : Argument 'to' containes duplicated names :", paste0(dupto, collapse=", ")) }
#'
#'     cat("\n\nCurrent sampleMetadata's metadata column names are (in the order) :", paste0(colnames(CYTdata@sampleData@sampleMetadata), collapse = ", "), "\n\n\nThe following columns :\n - ",
#'         paste0(from, collapse=", "), "\n\nwill be renamed, in the order, by :\n - ", paste0(to, collapse=", "))
#'
#'     colnames(CYTdata@sampleData@sampleMetadata)[match(colnames(CYTdata@sampleData@sampleMetadata), from)] = to
#'     names(CYTdata@sampleData@metadataPalette)[match(names(CYTdata@sampleData@metadataPalette), from)] = to
#'     massage("Check unicity of meadata column names..")
#'     CYTdata = checkValidity(CYTdata, mode = "warning")
#'   }
#'   else { # Rename metadata
#'     if (!metadata %in% colnames(CYTdata@sampleData@sampleMetadata)) { stop("Error : 'metadata' argument (", metadata, ") is not NULL or a metadata column names") }
#'     oldmdlev = levels(CYTdata@sampleData@sampleMetadata[,metadata])
#'
#'     if (is.null(from)) { from = oldmdlev }
#'     else {
#'       mdErr = setdiff(from, oldmdlev)
#'       if (length(mdErr)>0) { stop("Error : 'from' argument contains metadata values not present among", metadata, "metadata values (", paste0(mdErr, collapse=", "), ".") }
#'       dupfrom = unique(from[duplicated(from)])
#'       if (length(dupfrom)>0) { stop("Error : Argument 'from' contains duplicated", metadata, "metadata values :", paste0(dupfrom, collapse=", ")) }
#'     }
#'
#'     if (length(from)!=length(to)) { stop("Error : Length of argument 'from' (", length(from), ") and argument 'to' (", length(to), ") must be equal") }
#'     if (length(to) == 1) { to = rep(to, length(from)) }
#'     newmdlev = oldmdlev
#'     newmdlev[match(from, newmdlev)] = to
#'
#'     levdup = unique(newmdlev[duplicated(newmdlev)])
#'     if (length(levdup)>0){
#'       message("After renaming, several ", metadata, "metadata values have the same name (",  paste0(levdup, collapse=", "), "), either by renaming to an already existing and unchanged metadata value, or by duplicate in the 'to' argument, or both.")
#'       if (merge) { message("\nWarning : 'merge' argument is set to TRUE. After renaming,", metadata, " column of sampleMetadata dataframe get its duplicated levels merged.") }
#'       else { stop("\nError : 'merge' argument is set to FALSE. CYTdata unchanged.") }
#'     }
#'     cat("\n\nCurrent", metadata, "values (levels) are (in the order) :", paste0(oldmdlev, collapse = ", "), "\n\n\nThe following values :\n - ", paste0(from, collapse=", "),
#'         "\n\nwill be renamed, in the order, by :\n - ", paste0(to, collapse=", "))
#'     CYTdata@sampleData@sampleMetadata[,metadata] = plyr::mapvalues(CYTdata@sampleData@sampleMetadata[,metadata], from = from, to = to)
#'
#'     message("After renaming,", metadata, "column of sampleMetadata dataframe contained duplicated value levels, which were merged. As it concerns metadataPalette : \n
#'               - If the name of merged samples was already associated to an existing sample name before renaming, the metadata of the existing sample are the new metadata of newly renamed sample(s).\n
#'               - If the name of merged samples was a new one (argument 'to' contain this name several times). The metadata of the resulting sample is the one of the sample contained in 'from' argument and which was ordered first (in 'from' argument) among all the samples renamed to this name.")
#'     pal = CYTdata@sampleData@metadataPalette[metadata]
#'     palnames = names(pal)
#'     palnames[match(palnames, from)] = to
#'     CYTdata@sampleData@metadataPalette[metadata] = structure(pal[!duplicated(palnames)], names = palnames[!duplicated(palnames)])
#'     !duplicated(palnames)
#'
#'
#'
#'
#'       newmd = data.frame()
#'       for (spl in levels(CYTdata@sampleData@cellSample$SampleID)) {
#'         fromspl = oldspllev[newspllev == spl]
#'
#'         if (length(fromspl)>1){ # if duplicated
#'           message("Warning :", paste0(fromspl, collapse=", "), "samples have been merged into", spl, "sample.")
#'           if (spl %in% fromspl) {
#'             message("Warning :", spl, "sample original metadata has been conserved.")
#'             metadataRow = oldmd[spl,]
#'           } # sample name already present
#'           else {
#'             spl2 = from[to==spl][1]
#'             message("Warning :", spl, "sample wasn't already present so no original metadata has been conserved. The metadata of the first sample, in the 'from' argument (", spl2, "), renamed to", spl, "has been conferred to", spl, ".")
#'             metadataRow = oldmd[spl2,]
#'           } # sample name not present, first sample
#'         }
#'         else {
#'           if (spl %in% to) { metadataRow = oldmd[from[to==spl],] }
#'           else { metadataRow = oldmd[spl,] }
#'         }
#'         newmd = rbind.data.frame(newmd, metadataRow)
#'       }
#'       rownames(newmd) = levels(CYTdata@sampleData@cellSample$SampleID)
#'       CYTdata@sampleData@sampleMetadata = newmd
#'     }
#'
#'     CYTdata = checkValidity(CYTdata, mode = "warning")
#'     return(CYTdata)
#'   }

#' #' @title Assigns meta-information about biological samples
#' #'
#' #' @description This function aims to attach meta-information to each biological sample.
#' #'
#' #' Especially, the following meta-information of each sample can be specified for subsequent analyses.
#' #' - The biological individual
#' #' - The biological condition (groups, vaccines, studies, etc.)
#' #' - The timepoint
#' #' Timepoint and Individual data must be specified.
#' #'
#' #' @param CYTdata a CYTdata object
#' #' @param metadata a dataframe containing meta-information about the biological samples.
#' #' The columns must contain, at least, a column named "Timepoint" and an other named "Individual".
#' #' The rownames have to be the biological samples, thus the number rows has to be equal to the number of samples.
#' #'
#' #' @return a S4 object of class 'CYTdata'
#' #'
#' #'
#' #' @export
#' #'
#'
#' reorderMetadata <- function(CYTdata, metadata, alphabetic = TRUE, newOrder = NULL){
#'
#'   CYTdata = checkValidity(CYTdata, mode = "error")
#'   checkmate::qassert(alphabetic, "B1")
#'
#'   oldOrder = levels(CYTdata@sampleData@sampleMetadata[,metadata])
#'   if (alphabetic) { newOrder = gtools::mixedsort(oldOrder) }
#'   else {
#'
#'     if (is.null(newOrder)) { stop("Error : 'alphabetic' argument is set to FALSE but 'newOrder' argument is NULL.") }
#'     mdErr = setdiff(newOrder, oldOrder)
#'     if (length(mdErr)>0) { stop("Error : 'newOrder' argument contains metadata values not present among", metadata, "metadata values (Not present : ", paste0(mdErr, collapse=", "), ".") }
#'     mdErr = setdiff(oldOrder, newOrder)
#'     if (length(mdErr)>0) { stop("Error : 'newOrder' argument does not contain all metadata values present among", metadata, "metadata values (Missing : ", paste0(mdErr, collapse=", "), ".") }
#'     dupnewOrder = unique(newOrder[duplicated(newOrder)])
#'     if (length(dupnewOrder)>0) { stop("Error : Argument 'newOrder' contains duplicated", metadata, "metadata values :", paste0(dupnewOrder, collapse=", ")) }
#'
#'
#'     CYTdata@sampleData@sampleMetadata[,metadata] = factor(CYTdata@sampleData@sampleMetadata[,metadata], levels = gtools::mixedsort(oldOrder))
#'   }
#'   CYTdata@sampleData@sampleMetadata[,metadata] = factor(CYTdata@sampleData@sampleMetadata[,metadata], levels = newOrder)
#'   CYTdata = checkValidity(CYTdata, mode = "warning")
#'   return(CYTdata)
#' }

#'
#'
#'
#'
#' getSamplesMetadata <- function(CYTdata, metadataCondition){
#'
#'   if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
#'   else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }
#'
#'   checkmate::qassert(metadataCondition, "S1")
#'
#'   if (ncol(CYTdata@metadata)==0) {
#'     stop("Error : No metadata slot present in CYTdata object")
#'   }
#'   if (!metadataCondition %in% colnames(CYTdata@metadata)){
#'     stop("Error : 'metadataCondition' argument metadataCondition is not the
#'          name of a column present in metadata (", paste0(colnames(CYTdata@metadata), collapse=","), ")")
#'   }
#'
#'   metacondVect = CYTdata@metadata[[metadataCondition]]
#'   res = lapply(levels(metacondVect),
#'                FUN = function(co){
#'                  spls = subset(CYTdata@metadata, metacondVect == co)
#'                  return (rownames(spls))
#'                })
#'   names(res) = levels(metacondVect)
#'   return(res)
#' }


