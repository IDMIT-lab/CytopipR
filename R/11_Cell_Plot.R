#' @title Plot Cells in 2D with Various Colorings
#'
#' @description
#' This function creates a 2D plot of cells from cytometry data. It supports visualization of the data in multiple ways, including dimensionality reduction (e.g., PCA or t-SNE), marker expression comparisons, clustering results, density plots, and metadata coloring. The function also allows for downsampling and various customization options such as linear regression plots.
#'
#' @param CYTdata An object of class \code{CYTdata} containing cytometry data.
#' @param expression A character string specifying the type of expression data to plot. Options are "DimRed" for dimensionality reduction or "markers" for specific marker expressions.
#' @param DimRed A string specifying the name of the dimensionality reduction method (e.g., PCA, t-SNE). Used if \code{expression = "DimRed"}.
#' @param coupleMarkers A vector of two marker names to compare. Used if \code{expression = "markers"}.
#' @param colorBy A character string specifying how to color the plot. Options are "clustering", "marker", "density", "metadata", or "none".
#' @param clustering A string specifying the clustering method (used if \code{colorBy = "clustering"}).
#' @param clusters A vector of cluster names or indices to include in the plot (used if \code{colorBy = "clustering"}).
#' @param printClustering A logical value indicating whether to print the centroids of each cluster on the plot (used if \code{colorBy = "clustering"}).
#' @param marker A string specifying the marker to use for coloring the plot (used if \code{colorBy = "marker"}).
#' @param markerPalette A vector of colors to be used for marker expression coloring (used if \code{colorBy = "marker"}).
#' @param densityPalette A vector of colors for density plot coloring (used if \code{colorBy = "density"}).
#' @param metadata A string specifying the metadata column name to use for coloring the plot (used if \code{colorBy = "metadata"}).
#' @param samples A vector of sample IDs to include in the plot. Optional.
#' @param filteringMethod A string specifying the filtering method to apply. Options are "remove" or "color".
#' @param filteringColor The color to use for points that are filtered out (used if \code{filteringMethod = "color"}).
#' @param noneColor The color to use for the plot when no coloring is applied (used if \code{colorBy = "none"}).
#' @param percentDownsampling A numeric value between 0 and 1 to downsample the data by a given percentage.
#' @param bounds A vector of two numeric values representing the quantile bounds for marker expression (used if \code{colorBy = "marker"}).
#' @param plotLinearRegression A logical value indicating whether to plot a linear regression line (used if \code{expression = "markers"}).
#' @param ylim The y-axis limits for the plot.
#' @param xlim The x-axis limits for the plot.
#'
#' @return A \code{ggplot} object representing the 2D plot of cells.
#'
#' @details
#' This function provides several ways to visualize the relationship between cells and different features:
#' - **DimRed**: Displays a dimensionality reduction (e.g., PCA, t-SNE) of the cell data.
#' - **markers**: Plots the expression of two markers and can color by a selected feature (clustering, metadata, etc.).
#' - **colorBy options**:
#'   - **"clustering"**: Colors the points by cluster label, and optionally adds centroids for each cluster.
#'   - **"marker"**: Colors the points based on marker expression, with optional quantile clipping.
#'   - **"density"**: Colors the points based on cell density.
#'   - **"metadata"**: Colors the points based on metadata values.
#'   - **"none"**: Displays the data with no color applied.
#'
#' Users can customize the appearance of the plot, downsample the data for better performance, and apply different filtering methods to exclude points.
#'
#' @examples
#' # Example usage for DimRed expression with clustering coloring:
#' plotCells(CYTdata, expression = "DimRed", DimRed = "PCA1", colorBy = "clustering", clustering = "kmeans")
#'
#' # Example usage for marker expression comparison with density coloring:
#' plotCells(CYTdata, expression = "markers", coupleMarkers = c("CD4", "CD8"), colorBy = "density")
#'
#' # Example usage for marker expression with custom color palette:
#' plotCells(CYTdata, expression = "markers", coupleMarkers = c("CD4", "CD8"), colorBy = "marker", markerPalette = c("blue", "red"))
#'
#' @import ggplot2
#' @import checkmate
#' @import ggnewscale
#' @import ggpointdensity
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom viridis scale_color_viridis
#' @import dplyr
#' @import ggrepel
#'
#' @export

plotCells <- function(CYTdata,

                      expression = c("DimRed", "markers"),
                      DimRed = NULL,
                      coupleMarkers = NULL,

                      colorBy = c("clustering", "marker", "density", "metadata", "none"),

                      clustering = NULL,
                      clusters = NULL,
                      printClustering = TRUE,

                      marker = NULL,
                      markerPalette = NULL,

                      densityPalette = NULL,

                      metadata = "Individual",
                      samples = NULL,

                      filteringMethod = c("remove", "color"),
                      filteringColor = "grey",

                      noneColor = "red",

                      percentDownsampling = NULL,
                      bounds = c(0.05, 0.95),
                      plotLinearRegression = FALSE,
                      ylim = NULL, xlim = NULL) {

  CYTdata = checkValidity(CYTdata, mode = "error")

  if (expression == "DimRed") {
    if (is.null(DimRed) || (!DimRed %in% names(CYTdata@cellData@cellDimRed))) {
      stop("Error : 'DimRed' argument (", DimRed, ") is NULL or not a DimRed's list name")
    }
    data = CYTdata@cellData@cellDimRed[[DimRed]]
    title = "Dimensionality reduction view"
  } else {
    checkmate::qassert(coupleMarkers, "S2")
    coupleMarkers = checkorderMarkers(CYTdata, coupleMarkers, order=FALSE, checkDuplicates=FALSE)
    data = CYTdata@cellData@cellExprs[,coupleMarkers,drop=FALSE]
    title = paste(coupleMarkers[1], "vs", coupleMarkers[2], "expression comparison", sep=" ")
  }

  if (!is.null(percentDownsampling)) {
    if (percentDownsampling>=1 || percentDownsampling<=0) {
      stop("Error : 'percentDownsampling' argument must be a positive integer between 0 and 1.")
    }
    data = data[sample(1:nrow(data), round(nrow(data)*percentDownsampling), replace=F),]
  }
  labs = colnames(data)

  checkmate::qassert(samples, c(0, "S*"))
  samples = checkorderSamples(CYTdata, samples = samples, order = TRUE, checkDuplicates = TRUE)
  subIdx = CYTdata@sampleData@cellSample[rownames(data), "SampleID"] %in% samples

  if (!is.null(clusters)) {
    checkmate::qassert(clustering, c(0, "S1"))
    if (is.null(clustering) || (!clustering %in% colnames(CYTdata@clusteringData@cellClustering))) {
      stop("Error : 'clusters' argument is not null but 'clustering' argument (", clustering, ") is NULL or not a clustering column name.")
    }
    clusters = tryCatch(
      { checkorderClustering(CYTdata, clustering = clustering, clusters = clusters, order = TRUE, checkDuplicates = TRUE) },
      error = function(e) {
        message("Error : 'clustering' argument is set to", clustering, "but, for this clustering, checkorderClustering function returns error with 'clusters' argument (",
                paste0(clusters, collapse=", "), ").\n checkorderClustering's error : ", e$message)
        return(NULL)
      })
    subIdx = subIdx & CYTdata@clusteringData@cellClustering[rownames(data),clustering] %in% clusters
  }


  plot <- ggplot2::ggplot()
  if (sum(!subIdx)>0) {
    filteringMethod = match.arg(filteringMethod)
    checkmate::qassert(filteringMethod, "S1")
    if (filteringMethod!="remove") {
      checkmate::qassert(filteringColor, "S1")
      if (!areColors(filteringColor)) { stop("Error : 'filteringColor' argument is not a hexadecimal color.") }
      dataGrey = data[!subIdx,]
      if (nrow(dataGrey)>0) {
        plot <- ggplot2::ggplot() +
          ggplot2::geom_point(data = dataGrey,
                              aes_string(x = labs[1], y = labs[2]),
                              color = filteringColor, size = 0.001) +
          ggnewscale::new_scale_color()
      }
    }
  }

  data = data[subIdx,]

  if (colorBy=="clustering") {

    data$clustering = CYTdata@clusteringData@cellClustering[rownames(data), clustering]

    cat("\nPlotting 2D-reduced data and using palette slot to color the clustering", clustering,
        ": \n - ", paste0(levels(data$clustering)[levels(data$clustering) %in% unique(data$clustering)], collapse = ", "))

    plot <- plot +
      ggplot2::geom_point(data = data,
                          ggplot2::aes_string(x = labs[1], y = labs[2], col = "clustering"),
                          size = 0.01) +
      ggplot2::scale_color_manual(values = CYTdata@clusteringData@clusteringPalette[[clustering]]) +
      ggplot2::labs(title = paste0(title, "with coloration by clustering (", clustering, ")"), col = clustering) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 1.5), ncol = 1))

    checkmate::qassert(printClustering, "B1")
    if (printClustering) {
      centroids = data %>%
        group_by(clustering) %>%
        summarise(across(everything(), median))

      plot <- plot +
        ggrepel::geom_text_repel(data = centroids,
                                 ggplot2::aes_string(x = labs[1], y = labs[2], label ="clustering"),
                                 color = "black", size = 8, force = 5)
    }

  } else if (colorBy=="marker") {

    checkmate::qassert(marker, "S1")
    if (is.null(marker)) { stop() }
    marker = checkorderMarkers(CYTdata, marker, order=TRUE, checkDuplicates=TRUE)

    checkmate::qassert(bounds, "N2")
    if (any(bounds<=0) || any(bounds>=1)) {
      stop("Error : 'bounds' argument is a vector of quantile bounds and must be two positive integer between 0 and 1")
    }

    checkmate::qassert(markerPalette, c("0", "S*"))
    if (!is.null(markerPalette)){
      areColors <- function(x) { sapply(x, function(X) {
        tryCatch(is.matrix(col2rgb(X)), error = function(e) FALSE) }) }
      if(!all(areColors(markerPalette))){
        stop("Error : 'markerPalette' argument (", paste0(markerPalette, collapse = ","), "), does not contain only hexadecimal color.)")
      }
    }
    else { markerPalette = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9,"RdBu")))(85) }

    data$value = CYTdata@cellData@cellExprs[rownames(data), marker]
    limits = stats::quantile(data$value, probs = bounds)
    data$value[data$value < limits[1]] = limits[1]
    data$value[data$value > limits[2]] = limits[2]

    plot <- plot +
      ggplot2::geom_point(data = data,
                          ggplot2::aes_string(x = labs[1], y = labs[2], col = "value"),
                          size = 0.001) +
      ggplot2::scale_color_gradientn(colours = markerPalette) +
      ggplot2::labs(title = paste0(title, "with", marker, "'s gradient representation"), col = marker)

  } else if (colorBy=="density") {
    if (sum(subIdx)<nrow(CYTdata@sampleData@cellSample)) {
      message("Warning : 'colorBy' argument is set to 'density' and 'samples' and 'clusters' arguments are not set to NULL
              which means that density is not computed on the whole dataset but the subsetted dataset data are available for density computation.")
    }
    checkmate::qassert(densityPalette, c("0", "S*"))
    if(!is.null(densityPalette)){
      if(!all(areColors(densityPalette))){ stop("Error : 'densityPalette' argument is a color palette (",
                                                paste0(densityPalette, collapse = ","),
                                                "), does not contain only hexadecimal color.") }
    }
    else { densityPalette = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11,'Spectral')))(11) }
    plot <- plot +
      ggpointdensity::geom_pointdensity(data = data,
                                        aes_string(x = labs[1], y = labs[2]),
                                        size=0.001) +
      ggplot2::scale_color_manual(values = densityPalette) +
      viridis::scale_color_viridis() +
      ggplot2::labs(color = "Cell density", title = paste0(title, "with coloration by density"))

  } else if (colorBy=="metadata") {
    checkmate::qassert(metadata, "S1")
    if (!metadata %in% colnames(CYTdata@sampleData@sampleMetadata)) { stop() }
    data = cbind.data.frame(data,
                            "metadata" = CYTdata@sampleData@sampleMetadata[CYTdata@sampleData@cellSample[rownames(data),"SampleID"],metadata])
    plot <- plot +
      ggplot2::geom_point(data = data,
                          ggplot2::aes_string(x = labs[1],
                                              y = labs[2],
                                              color = "metadata"), size = 0.001) +
      ggplot2::scale_color_manual(values = CYTdata@sampleData@metadataPalette[[metadata]]) +
      ggplot2::labs(color = metadata, title = paste0(title, " with coloration by metadata (", metadata, ").")) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 2), ncol = 1))

  } else {
    checkmate::qassert(noneColor, "S1")
    if (!areColors(noneColor)) { stop("Error : 'noneColor' argument is not a hexadecimal color.") }
    plot <- plot +
      ggplot2::geom_point(data = data,
                          aes_string(x = labs[1],
                                     y = labs[2]),
                          size = 0.001, color = noneColor) +
      ggplot2::ggtitle(title)
  }

  if (expression == "markers") {
    if (!is.null(ylim)){ plot <- plot + ggplot2::ylim(ylim[1], ylim[2]) }
    if (!is.null(xlim)){ plot <- plot + ggplot2::xlim(xlim[1], xlim[2]) }
    if (plotLinearRegression){
      plot <- plot +
        ggplot2::geom_smooth(data, ggplot2::aes_string(x = labs[1], y = labs[2], colour="linear", fill="linear"),
                             method="lm", formula=labs[1]~labs[2])
    }
  }

  plot <- plot +
    ggplot2::xlab(labs[1]) + ggplot2::ylab(labs[2]) +
    ggplot2::scale_x_continuous(expand = c(0.01,0.01)) +
    ggplot2::scale_y_continuous(expand = c(0.01,0.01)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.title = ggplot2::element_text(size=10,face="bold"),
                   #axis.text = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black"),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   aspect.ratio = 1,
                   legend.position = "right",
                   legend.title = ggplot2::element_text(hjust = 0.5, vjust = 1, face = 'bold'))

  return(plot)
}

