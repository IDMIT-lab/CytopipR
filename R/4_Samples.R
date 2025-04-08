#' @title Rename Samples in CYTdata Object
#'
#' @description
#' This function allows users to rename samples in a CYTdata object. The function can handle duplicate sample names, merge samples if necessary, and update the associated metadata.
#'
#' @param CYTdata A CYTdata object that will be updated with the renamed samples.
#' @param merge A logical value indicating whether to merge samples with duplicate names after renaming. Default is `FALSE`.
#' @param from A character vector of current sample names to be renamed.
#' @param to A character vector of new sample names corresponding to the names in the `from` argument.
#' @return A modified CYTdata object with the samples renamed and metadata updated.
#'
#' @details
#' This function renames samples in the `cellSample` dataframe of the CYTdata object. It also updates the associated sample metadata in the `sampleMetadata` slot.
#' If there are duplicate sample names after renaming, the function can either stop execution or merge the samples based on the `merge` argument.
#' If `merge` is `TRUE`, the samples with duplicate names are merged, and their metadata is handled accordingly:
#' - If the name of the merged samples existed previously, the metadata of the existing sample is retained.
#' - If the name of the merged sample was new, the metadata from the first renamed sample is retained.
#'
#' The function also checks that the `from` and `to` arguments are the same length. If they are not, an error is raised.
#'
#' @examples
#' # Example 1: Rename samples without merging
#' updated_CYTdata <- renameSamples(CYTdata = cyt_data, from = c("Sample1", "Sample2"), to = c("NewSample1", "NewSample2"))
#'
#' # Example 2: Rename samples and merge duplicate samples
#' updated_CYTdata <- renameSamples(CYTdata = cyt_data, from = c("Sample1", "Sample1"), to = c("MergedSample1", "MergedSample1"), merge = TRUE)
#'
#' @seealso
#' \code{\link{checkValidity}} for validating the CYTdata object before and after renaming samples.
#'
#' @import checkmate
#' @import dplyr
#' @import plyr
#' @export

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

#' @title Reorder Samples levels in CYTdata Object
#'
#' @description
#' This function allows users to reorder the samples in a CYTdata object either alphabetically or based on a user-provided order.
#' The reordering will apply to both the sample identifiers (`SampleID`) and the associated metadata in the `sampleMetadata` slot.
#'
#' @param CYTdata A CYTdata object containing the cytometry data. The object must be valid before and after reordering the samples.
#' @param alphabetic A logical value indicating whether the samples should be reordered alphabetically. If `TRUE`, the samples will be sorted alphabetically by `SampleID`.
#' If `FALSE`, the user must provide a custom ordering in `newOrder`.
#' @param newOrder A character vector specifying the custom order of the samples. This argument is required if `alphabetic` is `FALSE`.
#'
#' @return A modified CYTdata object with reordered samples. The sample order will be reflected in both the `SampleID` column of `cellSample` and the rows of `sampleMetadata`.
#'
#' @details
#' This function allows users to control the order of samples in the `CYTdata` object. If `alphabetic` is set to `TRUE`, the samples are sorted in lexicographical order
#' based on their `SampleID` values. If `alphabetic` is set to `FALSE`, the function expects the `newOrder` argument to provide a custom order for the samples.
#'
#' The function also updates the `sampleMetadata` to reflect the new sample order.
#'
#' If `newOrder` is provided but does not match the available samples, the function will raise an error.
#'
#' @examples
#' # Example: Reorder samples alphabetically
#' reordered_data <- reorderSamples(CYTdata = cyt_data, alphabetic = TRUE)
#'
#' # Example: Reorder samples based on a custom order
#' custom_order <- c("SampleB", "SampleA", "SampleC")
#' reordered_data <- reorderSamples(CYTdata = cyt_data, alphabetic = FALSE, newOrder = custom_order)
#'
#' @seealso
#' \code{\link{checkValidity}} for validating the CYTdata object before and after reordering.
#'
#' @import checkmate
#' @import gtools
#' @export

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

#' @title Sample Counts and Plot for CYTdata Object
#'
#' @description
#' This function computes the number of cells per sample in a CYTdata object and generates a plot with key statistical indicators such as mean, median, quantiles, and min/max values. It can optionally filter and sort the samples.
#'
#' @param CYTdata A CYTdata object containing sample data.
#' @param samples A character vector of sample names to include in the analysis. If `NULL`, all samples are used.
#' @param sort A logical value indicating whether the resulting data should be sorted in descending order by the number of cells (default is `TRUE`).
#' @return A ggplot2 plot visualizing the number of cells per sample along with key statistical metrics.
#'
#' @details
#' This function calculates the number of cells in each sample and then creates a plot to visualize the cell counts. The plot includes the following statistical lines:
#' - **Mean**: Green line
#' - **Min**: Blue line
#' - **25th Percentile (q25)**: Purple line
#' - **Median**: Red line
#' - **75th Percentile (q75)**: Purple line
#' - **Max**: Blue line
#'
#' The function also includes data labels with the actual values of these metrics.
#'
#' If `samples` is provided, only the specified samples are considered. If `sort` is set to `TRUE`, the resulting data will be sorted by the number of cells in descending order.
#'
#' @examples
#' # Example 1: Plot the cell counts for all samples
#' cell_count_plot <- samplesCounts(CYTdata = cyt_data)
#'
#' # Example 2: Plot the cell counts for specific samples
#' cell_count_plot <- samplesCounts(CYTdata = cyt_data, samples = c("Sample1", "Sample2"))
#'
#' # Example 3: Plot the cell counts without sorting
#' cell_count_plot <- samplesCounts(CYTdata = cyt_data, sort = FALSE)
#'
#' @seealso
#' \code{\link{checkValidity}} for validating the CYTdata object before generating the plot.
#'
#' @import checkmate
#' @import plyr
#' @import dplyr
#' @import ggplot2
#' @export

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


#' @title Import Sample Metadata into CYTdata Object
#'
#' @description
#' This function allows users to import a sample metadata dataframe (`sampleMetadata_df`) into a CYTdata object.
#' It also provides the option to overwrite existing sample metadata, convert columns to factors, and optionally specify a palette for metadata visualization.
#'
#' @param CYTdata A CYTdata object that will be updated with the new sample metadata.
#' @param sampleMetadata_df A dataframe containing sample metadata. The columns of this dataframe should represent metadata variables associated with the samples.
#' @param optionnalPalette An optional list providing custom palettes for metadata visualization. If `NULL`, a default palette using `rainbow` is applied to each metadata variable.
#' @param checkOverwrite A logical value indicating whether to prompt the user to confirm overwriting existing sample metadata in the CYTdata object if it already exists. Default is `TRUE`.
#'
#' @return A modified CYTdata object with the updated sample metadata and associated metadata palette.
#'
#' @details
#' This function updates the `sampleMetadata` slot of the `sampleData` component in the CYTdata object with the provided `sampleMetadata_df`.
#' All columns in the `sampleMetadata_df` are converted to factors (if they are not already).
#' If the `checkOverwrite` parameter is `TRUE` and the `sampleMetadata` slot is not empty, the user is prompted to confirm whether to overwrite the existing metadata.
#' The function also updates the `metadataPalette` slot with a default palette or an optional custom palette provided by the user.
#'
#' @examples
#' # Example: Import sample metadata with default palette and overwrite confirmation
#' updated_CYTdata <- importsampleMetadata(CYTdata = cyt_data, sampleMetadata_df = metadata_df)
#'
#' # Example: Import sample metadata with a custom palette
#' custom_palette <- list(SampleType = c("Type1" = "red", "Type2" = "blue"))
#' updated_CYTdata <- importsampleMetadata(CYTdata = cyt_data, sampleMetadata_df = metadata_df, optionnalPalette = custom_palette)
#'
#' @seealso
#' \code{\link{checkValidity}} for validating the CYTdata object before and after importing metadata.
#'
#' @import checkmate
#' @import dplyr
#' @export

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

######################################################## Utils ########################################################

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

updateMetadataLevels <- function(CYTdata) {
  mdFrame = CYTdata@sampleData@sampleMetadata %>% mutate(across(where(is.factor), droplevels))
  palList = CYTdata@sampleData@metadataPalette
  for (md in colnames(mdFrame)) { palList[[md]] = palList[[md]][levels(mdFrame[,md])] }
  CYTdata@sampleData@sampleMetadata = mdFrame
  CYTdata@sampleData@metadataPalette = palList[colnames(mdFrame)]
  return(CYTdata)
}
