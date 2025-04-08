#' @title Change the Color Palette for Metadata or Clustering
#'
#' @description
#' This function allows users to change the color palette used for metadata or clustering visualizations in the given `CYTdata` object. The palette can either be automatically generated using a rainbow color scale or defined by the user with a custom color palette.
#'
#' @param CYTdata An object of class \code{CYTdata} containing cytometry data and associated metadata or clustering information.
#' @param type A character string specifying whether to change the palette for "metadata" or "clustering". Must be one of "metadata" or "clustering".
#' @param metadata A character string specifying the metadata column name (if \code{type = "metadata"}). Ignored if \code{type = "clustering"}.
#' @param clustering A character string specifying the clustering column name (if \code{type = "clustering"}). Ignored if \code{type = "metadata"}.
#' @param autoColorRainbow A logical value indicating whether to automatically generate a rainbow palette (default is \code{TRUE}). If \code{FALSE}, the user must provide a custom palette via the \code{homemadePalette} parameter.
#' @param homemadePalette A character vector of colors to be used as a custom palette (only used if \code{autoColorRainbow = FALSE}). Ignored if \code{autoColorRainbow = TRUE}.
#'
#' @return The modified \code{CYTdata} object with the updated color palette for the specified metadata or clustering.
#'
#' @details
#' This function allows users to either:
#' - Automatically generate a rainbow color palette for metadata or clustering using the `rainbow()` function, or
#' - Provide a custom color palette by supplying a vector of color names or hexadecimal color codes in the `homemadePalette` argument.
#'
#' The function checks that the given metadata or clustering exists in the respective data slots of the `CYTdata` object and modifies the corresponding palette accordingly.
#'
#' @examples
#' # Example usage:
#' # Automatically generate a rainbow palette for metadata
#' changePalette(CYTdata, type = "metadata", metadata = "cellType", autoColorRainbow = TRUE)
#'
#' # Provide a custom color palette for clustering
#' changePalette(CYTdata, type = "clustering", clustering = "cluster1", autoColorRainbow = FALSE, homemadePalette = c("red", "blue", "green"))
#'
#' @import checkmate
#' @importFrom methods stop
#'
#' @export

changePalette <- function(CYTdata,
                          type = c("metadata", "clustering"),
                          metadata = NULL,
                          clustering = NULL,
                          autoColorRainbow = TRUE,
                          homemadePalette = NULL){

  CYTdata = checkValidity(CYTdata, mode = "error")
  type = match.arg(type)
  checkmate::qassert(type, "S1")

  if (type=="metadata") {
    if (is.null(metadata)) { stop() }
    checkmate::qassert(metadata, "S1")
    if (!metadata %in% colnames(CYTdata@sampleData@sampleMetadata)) { stop() }

    checkmate::qassert(autoColorRainbow, "B1")
    if (autoColorRainbow) {
      message("'autoColorRainbow' set to TRUE, automated rainbow palette built for ", metadata)
      CYTdata@sampleData@metadataPalette[[metadata]] = structure(rainbow(nlevels(CYTdata@sampleData@sampleMetadata[,metadata])),
                                                                 names = levels(CYTdata@sampleData@sampleMetadata[,metadata]))
    } else {
      if (is.null(homemadePalette)) { stop() }
      checkmate::qassert(homemadePalette, "S*")
      CYTdata@sampleData@metadataPalette[[metadata]] = homemadePalette
    }
  }
  else {
    if (is.null(clustering)) { stop() }
    checkmate::qassert(clustering, "S1")
    if (!clustering %in% colnames(CYTdata@clusteringData@cellClustering)){ stop() }

    checkmate::qassert(autoColorRainbow, "B1")
    if (autoColorRainbow) {
      message("'autoColorRainbow' set to TRUE, automated rainbow palette built for ", clustering)
      CYTdata@clusteringData@clusteringPalette[[clustering]] = structure(rainbow(nlevels(CYTdata@clusteringData@cellClustering[,clustering])),
                                                                         names = levels(CYTdata@clusteringData@cellClustering[,clustering]))
    } else {
      if (is.null(homemadePalette)) { stop() }
      checkmate::qassert(homemadePalette, "S*")
      CYTdata@clusteringData@clusteringPalette[[clustering]] = homemadePalette
    }
  }

  CYTdata = checkValidity(CYTdata, mode = "warning")
}


#' @title Show Color Palettes for Metadata or Clustering
#'
#' @description
#' This function visualizes the color palette associated with either metadata or clustering information in the given `CYTdata` object. It allows users to view the color assignments for categories in metadata or clustering results and provides a plot to show the color distribution.
#'
#' @param CYTdata An object of class \code{CYTdata} containing cytometry data and associated metadata or clustering information.
#' @param type A character string specifying whether to show the palette for "metadata" or "clustering". Must be one of "metadata" or "clustering".
#' @param metadata A character string specifying the metadata column name (if \code{type = "metadata"}). Ignored if \code{type = "clustering"}.
#' @param clustering A character string specifying the clustering column name (if \code{type = "clustering"}). Ignored if \code{type = "metadata"}.
#' @param ncol An integer indicating the number of columns to display the palette in. Default is 10.
#' @param labelSize A numeric value specifying the size of the labels displayed in the palette plot. Default is 5.
#'
#' @return A \code{ggplot} object displaying the color palette for the chosen metadata or clustering.
#'
#' @details
#' This function extracts the color palette either from the metadata or the clustering information stored in the `CYTdata` object. It then arranges the colors in a grid format and labels each color with the corresponding category name. The plot shows the colors with the labels, allowing the user to see how each category is represented by color.
#'
#' The `ncol` parameter controls the number of columns in the palette plot, and `labelSize` adjusts the font size of the category labels on the plot.
#'
#' @examples
#' # Example usage:
#' # Display the color palette for a metadata column
#' showPalettes(CYTdata, type = "metadata", metadata = "cellType", ncol = 5, labelSize = 6)
#'
#' # Display the color palette for a clustering result
#' showPalettes(CYTdata, type = "clustering", clustering = "cluster1", ncol = 8, labelSize = 5)
#'
#' @import checkmate
#' @import ggplot2
#' @importFrom dplyr %>%
#' @importFrom methods stop
#'
#' @export

showPalettes <- function(CYTdata,
                         type = c("metadata", "clustering"),
                         metadata = NULL,
                         clustering = NULL,
                         ncol = 10,
                         labelSize = 5) {

  CYTdata = checkValidity(CYTdata, mode = "error")
  type = match.arg(type)
  checkmate::qassert(type, "S1")

  checkmate::qassert(labelSize, "N1")
  checkmate::qassert(ncol, "N1")
  if (ncol<=0) { stop("Error : 'ncol' argument must be a positive integer") }
  if (labelSize<1) { stop("Error : argument 'labelSize' must be positive integer") }

  if (type=="metadata") {
    if (is.null(metadata)) { stop() }
    checkmate::qassert(metadata, "S1")
    if (!metadata %in% colnames(CYTdata@sampleData@sampleMetadata)) { stop() }
    pal = CYTdata@sampleData@metadataPalette[[metadata]]

  }
  else {
    if (is.null(clustering)) { stop() }
    checkmate::qassert(clustering, "S1")
    if (!clustering %in% colnames(CYTdata@clusteringData@cellClustering)){ stop() }
    pal = CYTdata@clusteringData@clusteringPalette[[clustering]]
  }

  X = (1:length(pal)) %% ncol
  X[X==0] = ncol
  Y = ((1:length(pal)-1) %/% ncol)+1


  plot = data.frame("pop" = names(pal), "X" = X, "Y" = Y) %>%
    ggplot2::ggplot(ggplot2::aes(X, Y, fill= pop)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::geom_label(ggplot2::aes(X, Y, label= pop), fill = "white", size = labelSize) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.title = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   legend.position = "none")
  plot(plot)
}

############################################# Utils

areColors <- function(x) {
  sapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)),
             error = function(e) FALSE)
  })
}
