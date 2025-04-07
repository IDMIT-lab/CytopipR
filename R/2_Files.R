checkFiles <- function(files) {
  files_extensions = files %>%
    basename() %>%
    strsplit(split="\\.") %>%
    lapply(FUN = function(x) return(x[-1])) %>%
    unlist()
  if(any(files_extensions != "fcs")){
    stop("Error : Files", paste0(basename(files)[files_extensions != "fcs"], collapse = ","), "are not FCS files. Please convert the file(s) into FCS format before import.")
  }
}

#' @title Get the markers (and channels associated) common to a set of FCS or TXT files
#'
#' @description This function aims to import acquired cell events into a CYTdata object.
#'
#' Input files must be FCS or TXT files.
#' Only common channels and markers are returned
#'
#' @param files a character vector specifying the path of the tab-separated or FCS files to load
#'
#' @return a named character vector
#'
#' @export
#'

getChannels <- function(files, metadata) {

  checkmate::qassert(files, "S*")
  checkmate::qassert(metadata, "D*")
  checkFiles(files)

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
  heatmap = Heatmap(t(tab), name = "Presence", cluster_columns = F, col = col,border = TRUE,
                    row_names_side = "left", column_names_side = "top", heatmap_legend_param = list(title = "Presence"),
                    cell_fun = function(j, i, x, y, width, height, fill) { grid.rect(x, y, width, height, gp = gpar(fill = fill, col = "black")) },
                    top_annotation = HeatmapAnnotation(df = metadata))
  return(list("tab" = tab, "Heatmap" = draw(heatmap),
              "completeChannels" = completeChannels,
              "partialChannels" = partialChannels))
}


#' @title Imports of cell expression profiles from FCS files
#'
#' @description This function aims to import acquired cell events into a CYTdata object.
#'
#' Input files must be FCS files.
#' Different transformations can be applied such as logicle, arcsinh. If yes, parameters are added "..."
#' Cell marker having technical or biological biaises can be excluded during the import.
#'
#' @param files a character vector specifying the path of the tab-separated or FCS files to load
#' @param channels a character vector providing the channels to import from FCS/text files. By default, only de common channels to all the FCS files are kept
#' @param transformList an S4 object of class transformList with features integrated (transform function, channels to transform, parameters). See Flowcore package : https://www.bioconductor.org/packages/release/bioc/manuals/flowCore/man/flowCore.pdf
#' @param rawData boolean specifying if raW data should be also imported and stored into CYTdata's slot "rawExpression".
#' @param ... argument to pass onto method read.FCS from Flowcore package
#' @return a S4 object of class 'CYTdata'
#'
#' @importFrom checkmate qassert
#' @importFrom flowCore read.FCS
#'
#' @export
#'

createCYTdata <- function(files, channels, Ndownsampling = NULL, ...){

  checkmate::qassert(files, "S*")
  checkmate::qassert(channels, "S*")
  checkmate::qassert(Ndownsampling, c(0, "N1"))
  if (!is.null(Ndownsampling)) {
    if (Ndownsampling<=0) { stop() }
  }
  checkFiles(files)

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

  cellClustering = data.frame(matrix(ncol = 0, nrow = nrow(cellSample)))
  rownames(cellClustering) = rownames(cellSample)

  cat("\nTotal number of events :", nrow(cellExprs))

  cellData <- methods::new("cellData", cellExprs = cellExprs, cellAdditionalexprs = cellAdditionalexprs)
  sampleData <- methods::new("sampleData", cellSample = cellSample)
  clusteringData <- methods::new("clusteringData", cellClustering = cellClustering)
  CYTdata <- methods::new("CYTdata", cellData = cellData, sampleData = sampleData, clusteringData = clusteringData)


  CYTdata = checkValidity(CYTdata, mode = "warning")
  return(CYTdata)
}
