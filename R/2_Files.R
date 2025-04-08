#' Check if files have the correct FCS extension, and extract and visualize the channels of FCS files
#'
#' This function :
#' - checks whether the provided files have the `.fcs` extension. If any files do not have the correct extension, the error message displayed them.
#' - extracts the marker names (channels) from multiple FCS files, compares them across all files, and generates a heatmap to visualize the presence/absence across the files
#'
#' @param files A character vector containing the paths to the FCS files. Each file should be in FCS format.
#' @param metadata A data frame containing metadata to annotate the heatmap.
#'
#' @return A list containing:
#' \item{tab}{A matrix of 0s and 1s indicating the presence (1) or absence (0) of each channel in each file.}
#' \item{Heatmap}{A heatmap visualizing the presence/absence of channels across files. The heatmap uses grey shades to indicate presence and absence of channels.}
#' \item{completeChannels}{A character vector containing the channels that are present in all files.}
#' \item{partialChannels}{A character vector containing the channels that are not present in all files.}
#'
#' @examples
#' # Example with FCS files and metadata
#' files <- c("data/sample1.fcs", "data/sample2.fcs")
#' metadata <- data.frame(Condition = c("A", "B"))
#' rownames(metadata) <- gsub(pattern = ".fcs", replacement = "", x = basename(files))
#' result <- checkFCS_getChannels(files, metadata)
#'
#' # Accessing the results
#' result$tab  # Presence/absence matrix
#' result$completeChannels  # Channels present in all files
#' result$partialChannels  # Channels present in some but not all files
#' # Example with invalid files
#' files <- c("data/sample1.csv", "data/sample2.fcs")
#' result <- checkFCS_getChannels(files, metadata)  # This will throw an error for "sample1.csv"

checkFCS_getChannels <- function(files, metadata) {

  checkmate::qassert(files, "S*")
  checkmate::qassert(metadata, "D*")

  files_extensions = files %>%
    basename() %>%
    strsplit(split="\\.") %>%
    lapply(FUN = function(x) return(x[-1])) %>%
    unlist()
  if(any(files_extensions != "fcs")){
    stop("Error : Files", paste0(basename(files)[files_extensions != "fcs"], collapse = ","), "are not FCS files. Please convert the file(s) into FCS format before import.")
  }

  channels_list = list()
  markers_list = list()
  for (f in files) {
    fcs = flowCore::read.FCS(f, which.lines = 1:2, truncate_max_range = FALSE)
    channels_list = append(channels_list, list(flowCore::markernames(fcs)))
  }
  names(channels_list) = sapply(files, basename)

  channels = unique(unlist(channels_list))
  completeChannels = Reduce(intersect, channels_list)
  partialChannels = setdiff(channels, completeChannels)

  tab = matrix(0, nrow = length(channels_list), ncol = length(channels), dimnames = list(names(channels_list), channels))
  for (i in seq_along(channels_list)) { tab[i, channels_list[[i]]] = 1 }
  if (length(unique(as.vector(tab)))==2) { col=c("white", "grey") } else { col="grey" } # Couleurs : blanc pour 0, noir pour 1
  heatmap = ComplexHeatmap::Heatmap(t(tab), name = "Presence", cluster_columns = F, col = col,border = TRUE,
                                    row_names_side = "left", column_names_side = "top", heatmap_legend_param = list(title = "Presence"),
                                    width = ncol(t(tab))*unit(1, "cm"), height = nrow(t(tab))*unit(0.8, "cm"),
                                    cell_fun = function(j, i, x, y, width, height, fill) { grid::grid.rect(x, y, width, height, gp = grid::gpar(fill = fill, col = "black")) },
                                    top_annotation = ComplexHeatmap::HeatmapAnnotation(df = metadata))
  return(list("tab" = tab, "Heatmap" = ComplexHeatmap::draw(heatmap),
              "completeChannels" = completeChannels,
              "partialChannels" = partialChannels))
}


#' @title Create CYTdata Object
#'
#' @description
#' This function processes cytometry data from FCS files or a data frame to create a CYTdata object,
#' which is a standardized format used for downstream analysis, such as clustering or statistical tests.
#' The function allows users to specify channels to retain, perform optional downsampling, and validate the input data.
#'
#' @param files A character vector containing paths to FCS files. Required if the format is set to 'fcs'.
#' @param data A data.frame containing the cytometry data. Required if the format is set to 'data.frame'.
#' @param format A string specifying the format of the input data. Can either be "fcs" (for FCS files) or "data.frame" (for an existing data frame). Default is 'fcs'.
#' @param channels A character vector specifying the channel names to retain from the input data.
#' @param Ndownsampling An integer specifying the number of events (cells) to downsample to for each FCS file. If NULL, no downsampling is performed. Default is NULL.
#' @param ... Additional arguments passed to the `flowCore::read.FCS` function when importing FCS files.
#'
#' @return A CYTdata object containing the processed cytometry data and metadata.
#' The object contains:
#' - `cellData`: Expression data for each cell.
#' - `sampleData`: Metadata for each sample.
#' - `clusteringData`: Empty data frame (reserved for future clustering results).
#'
#' @details
#' If the `format` argument is set to `"fcs"`, the function reads FCS files and processes them based on the specified channels.
#' It also performs optional downsampling to a specified number of events per file. Missing channels are filled with NA values,
#' and warnings are issued when channels are missing.
#'
#' If the `format` argument is set to `"data.frame"`, the function processes the provided data frame, ensuring that it contains
#' the required columns (`SampleID` and `CellID`). It samples up to 1000 cells per sample and splits the data into expression
#' values (`cellExprs`) and sample metadata (`cellSample`).
#'
#' The resulting CYTdata object can be used for further analysis, such as clustering or statistical tests.
#'
#' @examples
#' # Example using FCS files
#' fcs_files <- c("sample1.fcs", "sample2.fcs")
#' channels_to_keep <- c("FSC-A", "SSC-A", "CD3")
#' cyt_data <- createCYTdata(files = fcs_files, format = "fcs", channels = channels_to_keep, Ndownsampling = 500)
#'
#' # Example using a data frame
#' df <- data.frame(SampleID = rep(c("Sample1", "Sample2"), each = 1000),
#'                  CellID = rep(1:1000, 2),
#'                  FSC_A = rnorm(2000), SSC_A = rnorm(2000), CD3 = rnorm(2000))
#' cyt_data <- createCYTdata(data = df, format = "data.frame", channels = c("FSC_A", "SSC_A", "CD3"))
#'
#' @seealso
#' \code{\link[flowCore]{read.FCS}} for reading FCS files.
#' \code{\link{checkValidity}} for validating the CYTdata object.
#'
#' @import checkmate
#' @import flowCore
#' @import dplyr
#' @importFrom methods new
#' @importFrom dplyr select_if mutate_at group_by sample_n ungroup column_to_rownames
#' @export

createCYTdata <- function(files = NULL, data = NULL, format = c("fcs", "data.frame"), channels, Ndownsampling = NULL, ...){

  checkmate::qassert(Ndownsampling, c(0, "N1"))
    if (!is.null(Ndownsampling)) {
      if (Ndownsampling<=0) { stop() }
    }

  format = match.arg(format)
  checkmate::qassert(format, "S1")

  if (format=="fcs") {
    if (is.null(files)) { stop("Error : 'format' argument set to 'fcs' but files argument is equal to NULL")}
    checkmate::qassert(files, "S*")
    checkmate::qassert(channels, "S*")

    files_extensions = files %>%
      basename() %>%
      strsplit(split="\\.") %>%
      lapply(FUN = function(x) return(x[-1])) %>%
      unlist()
    if(any(files_extensions != "fcs")){
      stop("Error : Files", paste0(basename(files)[files_extensions != "fcs"], collapse = ","), "are not FCS files. Please convert the file(s) into FCS format before import.")
    }

    cat("Starting import of FCS files..")
    exprs = list() # Final matrix.expression data.frame
    sample = list()
    for (i in seq_along(files)) {
      f = files[[i]]
      cat(paste0("\n\nImporting ", basename(f)," file.. (", i, "/", length(files), ")"))
      fcs = flowCore::read.FCS(f)#, ...)

      channelsSub = intersect(channels, flowCore::markernames(fcs))
      if (length(channelsSub)==0) { stop("\nError : All the Channels that should be kept are not present in ", basename(f), " file. Please reset 'channels' argument.") }
      channelsAdd = setdiff(channels, channelsSub)
      if (length(channelsAdd)>0) { message("\nWarning : Some channels (", paste0(channelsAdd, collapse=", "), ") are not present in ", basename(f),
                                           " file. They will be filled with NA values for ", basename(f), " file and place in cellAdditionalexprs dataframe.") }

      df = as.data.frame(flowCore::exprs(fcs))[,names(flowCore::markernames(fcs))[match(channelsSub, flowCore::markernames(fcs))], drop = FALSE]
      cat("\n - Number of events for sample", basename(f), ":", nrow(df))
      colnames(df) = channelsSub
      df[,channelsAdd] = NA

      sampleID = gsub(pattern = ".fcs", replacement = "", x = basename(f))
      cellID = paste("Cell", 1:nrow(df), sep="")
      df2 = data.frame("SampleID" = sampleID, "CellID" = cellID)
      rownames(df) = paste(sampleID, cellID, sep=" - ")
      rownames(df2) = rownames(df)

      if (!is.null(Ndownsampling)) {
        if (Ndownsampling>nrow(df)) {
          cat("\n\nNdownsampling argument is bigger (", Ndownsampling, ") than the number of events for this sample. No dowsmapling perfromed for this sample")
        }
        else {
          idx = sample(1:nrow(df), Ndownsampling, replace = F)
          cat("\n\nEvents from this sample are dowsampled to ", Ndownsampling, " cells.")
          df = df[idx,,drop = FALSE]
          df2 = df2[idx,,drop = FALSE]
        }
      }

      exprs = append(exprs, list(df[,channels, drop = FALSE]))
      sample = append(sample, list(df2))
    }

    cellExprs = do.call(rbind.data.frame, exprs)
    cellAdditionalexprs = cellExprs %>% select_if(~ any(is.na(.)))
    cellExprs = cellExprs %>% select_if(~ !any(is.na(.)))
    if(ncol(cellExprs) == 0) {
      stop("Error : No channels is present with the same name for all the sample. Resulting cellExprs dataframe is empty.")
    }
    cellSample = do.call(rbind.data.frame, sample)
    cellSample = cellSample %>% mutate_at("SampleID", as.factor)

    cat("\nTotal number of events :", nrow(cellExprs))

  } else {
    if (is.null(data)) { stop("Error : 'format' argument set to 'data.frame' but data argument is equal to NULL")}
    checkmate::qassert(data, "D*")

    if (!"SampleID" %in% colnames(data)) { stop("Error : 'data' argument is a data.frame but does not contain a column 'SampleID'.")}
    if (!"CellID" %in% colnames(data)) { stop("Error : 'data' argument is a data.frame but does not contain a column 'CellID'.")}
    if (ncol(data)<3) { stop("Error : 'data' argument is a data.frame but contains only columns 'CellID' et 'SampleID'.")}

    data = data %>%
      group_by(SampleID) %>%
      sample_n(min(1000, n())) %>%
      ungroup() %>%
      mutate(ID = paste(SampleID, CellID, sep=" - ")) %>% column_to_rownames("ID")

    cellExprs = data %>% dplyr::select(-c("SampleID", "CellID"))
    cellSample = data %>% dplyr::select(c("SampleID", "CellID"))
    cellAdditionalexprs = cellExprs %>% select_if(~ any(is.na(.)))
    cellExprs = cellExprs %>% select_if(~ !any(is.na(.)))
  }

  cellClustering = data.frame(matrix(ncol = 0, nrow = nrow(cellSample)))
  rownames(cellClustering) = rownames(cellSample)

  cellData <- methods::new("cellData", cellExprs = cellExprs, cellAdditionalexprs = cellAdditionalexprs)
  sampleData <- methods::new("sampleData", cellSample = cellSample)
  clusteringData <- methods::new("clusteringData", cellClustering = cellClustering)
  CYTdata <- methods::new("CYTdata", cellData = cellData, sampleData = sampleData, clusteringData = clusteringData)

  CYTdata = checkValidity(CYTdata, mode = "warning")
  return(CYTdata)
}
