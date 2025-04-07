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


#' @title Get the samples associated to each value of a specific condition displayed in metadata data.frame
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param metadataCondition a string specifying the metadata condition chosen. The string has to be the name of the column (in metadata data.frame) containing the condition chosen
#'
#' @return a list containing the samples for each value of condition chosen
#'
#' @export
#'
#'

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

#' @title Get the samples associated to each value of a specific condition displayed in metadata data.frame
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param metadataCondition a string specifying the metadata condition chosen. The string has to be the name of the column (in metadata data.frame) containing the condition chosen
#'
#' @return a list containing the samples for each value of condition chosen
#'
#' @export
#'
#'

areColors <- function(x) {
  sapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)),
             error = function(e) FALSE)
  })
}
