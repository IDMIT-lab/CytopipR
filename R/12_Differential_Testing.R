#' @title Get Cell Count or Abundance
#' @description
#' This function calculates either the raw cell counts or the cell abundance (as percentages) for a specified clustering,
#' clusters, and samples from the CYTdata object. The result can be either a contingency table (cell count) or a relative abundance table.
#'
#' @param CYTdata A Cytometry data object containing the cell clustering and sample data.
#' @param clustering A string specifying the name of the clustering column to use for analysis.
#' @param clusters A vector of cluster labels to include. If NULL, all clusters will be used.
#' @param samples A vector of sample IDs to include in the analysis. If NULL, all samples will be used.
#' @param type A string specifying the type of output. Either "cellcount" for raw counts or "abundance" for relative percentages.
#'
#' @return A table (data.frame) containing either the cell counts or the relative abundance for each cluster-sample combination.
#' If `type = "cellcount"`, the result will be a table of raw counts. If `type = "abundance"`, the result will be a table of relative abundances (percentages).
#'
#' @examples
#' # Get raw cell counts for specific clusters and samples
#' cellcount_result <- getCellCount(CYTdata, clustering = "ClusterType", clusters = c("Cluster1", "Cluster2"), samples = c("SampleA", "SampleB"), type = "cellcount")
#'
#' # Get abundance (percentage) for the same clusters and samples
#' abundance_result <- getCellCount(CYTdata, clustering = "ClusterType", clusters = c("Cluster1", "Cluster2"), samples = c("SampleA", "SampleB"), type = "abundance")
#'
#' @seealso \code{\link{checkorderClustering}} for clustering order checking function.
#' @export

getCellCount <- function(CYTdata,

                         clustering,
                         clusters = NULL,
                         samples = NULL,

                         type = c("cellcount", "abundance")) {

  #CYTdata = checkValidity(CYTdata, mode = "error")

  type = match.arg(type)
  checkmate::qassert(type, "S1")

  checkmate::qassert(clustering, "S1")
  if (!clustering %in% colnames(CYTdata@clusteringData@cellClustering)) {
    stop("Error : 'clustering' argument (", clustering, ") is not a clustering column name.")
  }
  clusters = tryCatch(
    { checkorderClustering(CYTdata, clustering = clustering, clusters = clusters, order = TRUE, checkDuplicates = TRUE) },
    error = function(e) {
      message("Error : 'clustering' argument is set to", clustering, "but, for this clustering, checkorderClustering function returns error with 'clusters' argument (",
              paste0(clusters, collapse=", "), ").\n checkorderClustering's error : ", e$message)
      return(NULL)
    })
  gatedIdxClusters = CYTdata@clusteringData@cellClustering[,clustering] %in% clusters
  levClusters = CYTdata@clusteringData@cellClustering[gatedIdxClusters,clustering] %>% droplevels() %>% levels()

  checkmate::qassert(samples, c(0, "S*"))
  samples = checkorderSamples(CYTdata, samples = samples, order = TRUE, checkDuplicates = TRUE)
  gatedIdxSpl = CYTdata@sampleData@cellSample$SampleID %in% samples
  levSpl = CYTdata@sampleData@cellSample$SampleID[gatedIdxSpl] %>% droplevels() %>% levels()

  cellcount = table(CYTdata@clusteringData@cellClustering[gatedIdxSpl&gatedIdxClusters,clustering] %>% factor(levels = levClusters),
                    CYTdata@sampleData@cellSample$SampleID[gatedIdxSpl&gatedIdxClusters] %>% factor(levels = levSpl)) %>% as.data.frame.matrix()  # clusters full of 0 -> NA donc changer

  if (type=="cellcount") { return(cellcount) }
  else {
    abundance = apply(cellcount, 2, function(df){df / sum(df)}) * 100
    if(nrow(cellcount)==1){
      abundance = t(abundance)
      rownames(abundance) = rownames(cellcount)
    }
    abundance = data.frame(abundance, check.names = FALSE)
    return(abundance[clusters,,drop=FALSE])
  }
}

#' @title Run Differential Testing
#' @description
#' This function performs differential testing for cell clustering data based on sample group comparisons.
#' It supports both two-group and multi-group tests, with options for various statistical tests and post-hoc analyses.
#' The function computes fold change (FC) and adjusts p-values for multiple testing.
#' The result is stored in the differential testing data within the CYTdata object.
#'
#' @param CYTdata A Cytometry data object containing clustering and sample data.
#' @param clustering A string specifying the name of the clustering column to be analyzed.
#' @param clusters A vector of cluster labels for the analysis. If NULL, all clusters are used.
#' @param comparisonSpl A list of sample groups to compare against each other. Must contain at least two elements.
#' @param referenceSplName A string specifying the name of the reference sample group. If NULL, the first group is used.
#' @param differentialTestingTitle A string for the title of the differential testing results.
#' @param variable A string indicating the type of variable to analyze. Can be one of "cellcount", "abundance", or "marker".
#' @param NFSValues A named vector of normalization factor values for visualizations. Used only if variable is "abundance".
#' @param marker A string specifying the marker to analyze. Used only if variable is "marker".
#' @param comparisonSplFCref A list of sample groups for fold change reference.
#' @param FCrefMethod A string specifying the method for calculating fold change. Options are "ratio" or "percIncrease".
#' @param eps A small value to avoid division by zero when calculating fold change.
#' @param twogroupTest A string specifying the statistical test for two-group comparisons. Options are "wilcox.test", "t.test", or "permutation.test".
#' @param moregroupTest A string specifying the statistical test for multi-group comparisons. Options are "kruskal", "ANOVA", "ANOVA.MR", or "friedman".
#' @param posthocTest A string specifying the post-hoc test to apply. Options are "none", "dunn", or "tukey".
#' @param p.adjust A string specifying the method for adjusting p-values. Options are "none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", or "fdr".
#' @param paired A logical indicating whether to perform a paired test. Defaults to FALSE.
#' @param verbose A logical indicating whether to print verbose messages.
#' @param checkOverwrite A logical indicating whether to prompt for overwriting existing differential testing results.
#' @param ... Additional arguments passed to other functions.
#'
#' @return The updated CYTdata object with differential testing results stored in the `differentialTesting` slot.
#'         The results include p-values, fold changes, and significance values for each cluster comparison.
#'
#' @examples
#' # Example for performing differential testing with cell counts
#' CYTdata <- runDifferentialTesting(CYTdata,
#'                                   clustering = "ClusterType",
#'                                   comparisonSpl = list("GroupA" = c("Sample1", "Sample2"),
#'                                                        "GroupB" = c("Sample3", "Sample4")),
#'                                   differentialTestingTitle = "Test1",
#'                                   variable = "cellcount")
#'
#' # Example for performing differential testing with abundance values
#' CYTdata <- runDifferentialTesting(CYTdata,
#'                                   clustering = "ClusterType",
#'                                   comparisonSpl = list("GroupA" = c("Sample1", "Sample2"),
#'                                                        "GroupB" = c("Sample3", "Sample4")),
#'                                   differentialTestingTitle = "Test2",
#'                                   variable = "abundance",
#'                                   NFSValues = c(Sample1 = 1.2, Sample2 = 1.1, Sample3 = 0.9, Sample4 = 1.0))
#'
#' # Example for performing differential testing with a marker
#' CYTdata <- runDifferentialTesting(CYTdata,
#'                                   clustering = "ClusterType",
#'                                   comparisonSpl = list("GroupA" = c("Sample1", "Sample2"),
#'                                                        "GroupB" = c("Sample3", "Sample4")),
#'                                   differentialTestingTitle = "Test3",
#'                                   variable = "marker",
#'                                   marker = "Marker1")
#'
#' @seealso \code{\link{checkorderClustering}} for ordering and checking the clustering data.
#' @export

runDifferentialTesting <- function(CYTdata,

                                   clustering,
                                   clusters = NULL,

                                   comparisonSpl,
                                   referenceSplName = NULL,
                                   differentialTestingTitle = "cmp",

                                   variable = c("cellcount", "abundance", "marker"),
                                   NFSValues = NULL,
                                   marker = NULL,
                                   comparisonSplFCref = NULL,
                                   FCrefMethod = c("ratio", "percIncrease"),
                                   eps = 0.001,

                                   twogroupTest = c("wilcox.test", "t.test", "permutation.test"),
                                   moregroupTest = c("kruskal", "ANOVA", "ANOVA.MR", "friedman"),
                                   posthocTest = c("none", "dunn", "tukey"),
                                   p.adjust = c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),
                                   paired = FALSE,

                                   verbose = FALSE,
                                   checkOverwrite = TRUE, ...) {

  CYTdata = checkValidity(CYTdata, mode = "error")

  checkmate::qassert(checkOverwrite, "B1")
  if (checkOverwrite && differentialTestingTitle %in% unique(CYTdata@differentialTesting@data$title)) {
    reply <- readline(prompt="\nComparison already performed, do you still want to continue and overwrite (yes or no): ")
    while (!tolower(reply) %in% c("yes", "y", "no", "n")) { reply <- readline(prompt="Reply not standard, please answer yes or no: ") }
    if (tolower(reply) %in% c("no", "n")){ stop("Function stopped, CYTdata object unchanged") }
    CYTdata@differentialTesting@data = CYTdata@differentialTesting@data %>% dplyr::filter(title != differentialTestingTitle)
  }

  checkmate::qassert(clustering, "S1")
  if (!clustering %in% colnames(CYTdata@clusteringData@cellClustering)) {
    stop("Error : 'clustering' argument (", clustering, ") is not a clustering column name.")
  }
  clusters = tryCatch(
    { checkorderClustering(CYTdata, clustering = clustering, clusters = clusters, order = TRUE, checkDuplicates = TRUE) },
    error = function(e) {
      message("Error : 'clustering' argument is set to", clustering, "but, for this clustering, checkorderClustering function returns error with 'clusters' argument (",
              paste0(clusters, collapse=", "), ").\n checkorderClustering's error : ", e$message)
      return(NULL)
    })

  checkmate::qassert(comparisonSpl, "L*")
  checkmate::qassert(referenceSplName, c(0, "S1"))

  if (length(comparisonSpl)>1) {
    namesDup = names(comparisonSpl)[duplicated(names(comparisonSpl))]
    if (length(namesDup)>0) {
      stop("Error : 'comparisonSpl' argument list contain duplicated names (", paste0(namesDup, collapse=", "), ").")
    }
    for (i in seq_along(comparisonSpl)) {
      if (is.null(names(comparisonSpl)[i]) || is.na(names(comparisonSpl)[i]) || (names(comparisonSpl)[i] == "")) {
        stop("Error : Name of element", i, "of 'comparisonSpl' argument list contain NULL, NA or empty names.")
      }
      comparisonSpl[[i]] = tryCatch(
        { checkorderSamples(CYTdata, samples = comparisonSpl[[i]], order = TRUE, checkDuplicates = TRUE) },
        error = function(e) {
          message("Error : Element", i, "of 'comparisonSpl' argument list (", names(comparisonSpl)[i], ") is getting an error with checkorderSamples function.\n checkorderSamples's error : ", e$message)
          return(NULL)
        })
    }
    if (!is.null(referenceSplName)) {
      if (!referenceSplName %in% names(comparisonSpl)) {
        stop("Error : 'referenceSplName' argument (", referenceSplName, ") is not a comparisonSpl element name.")
      }
    }
  } else { stop("Error : 'comparisonSpl' argument must be a list of length >2") }

  checkmate::qassert(comparisonSplFCref, c(0, "L*"))
  if (!is.null(comparisonSplFCref)) {
    if (length(comparisonSplFCref)>1) {
      namesDup = names(comparisonSplFCref)[duplicated(names(comparisonSplFCref))]
      if (length(namesDup)>0) {
        stop("Error : 'comparisonSplFCref' argument list contain duplicated names (", paste0(namesDup, collapse=", "), ").")
      }
      if (!setequal(names(comparisonSplFCref), names(comparisonSpl))) {
        stop("Error : 'comparisonSplFCref' argument list contain differents names than 'comparisonSpl' argument list.")
      }
      comparisonSplFCref = comparisonSplFCref[names(comparisonSpl)]
      for (i in seq_along(comparisonSplFCref)) {
        comparisonSplFCref[[i]] = tryCatch(
          { checkorderSamples(CYTdata, samples = comparisonSplFCref[[i]], order = TRUE, checkDuplicates = TRUE) },
          error = function(e) {
            message("Error : Element", i, "of 'comparisonSplFCref' argument list (", names(comparisonSplFCref)[i], ") is getting an error with checkorderSamples function.\n checkorderSamples's error : ", e$message)
            return(NULL)
          })
        if (length(comparisonSplFCref[[i]])!=length(comparisonSpl[[i]])) {
          stop("Error : Element", i, "of 'comparisonSplFCref' argument list is a samples vector with not the same number of samples than the element", i, "of 'comparisonSpl' argument.")
        }
      }
      splList = append(comparisonSpl, comparisonSplFCref)
    } else { stop("Error : 'comparisonSplFCref' argument must be a list of length >2") }

    FCrefMethod = match.arg(FCrefMethod)
    checkmate::qassert(eps, "N1")
    switch(FCrefMethod,
           ratio = {
             getValue <- function(num, denom) {
               if (denom==0) { return(num/(denom+eps)) }
               else { return(num/denom) }
             }
           },
           percIncrease = {
             getValue <- function(num, denom) {
               if (denom==0) { return((num-denom)/(denom+eps)) }
               else { return((num-denom)/denom) }
             }
           })
    prefix = paste(FCrefMethod, "_", sep="")
  } else {
    splList = comparisonSpl
    prefix = ""
  }

  splVec = as.vector(unlist(splList))
  splDup = splVec[duplicated(splVec)]
  if (length(splDup)>0) {
    message("Warning : 'comparisonSpl' and 'comparisonSplFCref' arguments contain duplicated samples among all the samples vectors stored in it (",
            paste0(splDup, collapse=", "), ").")
  }
  splVec = unique(splVec)

  variable = match.arg(variable)
  if (variable == "marker") {
    checkmate::qassert(marker, "S1")
    if (is.null(marker)) { stop("Error : 'variable' argument is set to marker but 'marker' argument is NULL") }
    marker = tryCatch(
      { checkorderMarkers(CYTdata, markers = marker, cellSlot = "cellExprs", order = TRUE, checkDuplicates = TRUE) },
      error = function(e) {
        message("Error : 'variable' argument is set to marker but checkorderMarkers function returns error with 'marker' argument (", marker, ").\n checkorderMarkers's error : ", e$message)
        return(NULL)
      })

    df = computeMSI(CYTdata, samples = splVec,
                    clustering = clustering, clusters = clusters,
                    computeWith = "markers", marker = marker, ...)
    valueVar = paste(prefix, marker, sep="")
  } else {
    df = getCellCount(CYTdata, samples = splVec,
                      clustering = clustering, clusters = clusters,
                      type = variable)
    valueVar = paste(prefix, variable, sep="")
  }
  Ylab = paste(valueVar, "level")

  if (!is.null(NFSValues)) {
    if (variable != "abundance") {
      message("Warning : NFSValues argument is given but variable is set to", variable, "and not to abundancce. NFSValues argument is not taken into account.")
    }
    else {
      if (is.null(names(NFSValues))) { stop("Error : NFSValues argument is not a named vector, please see documentation.") }
      if (sum(is.na(NFSValues))>0) { stop("Error : NFSValues does not contain some samples (", paste0(names(NFSValues)[is.na(NFSValues)], collapse=", "), ") as named elements.") }

      spls_NFS = colnames(df)
      splErr = setdiff(spls_NFS, names(NFSValues))
      if (length(splErr)>0) { stop("Error : NFSValues does not contain needed samples (", paste0(splErr, collapse=", "), ") for visualisation.") }

      df = as.matrix(df) %*% diag(as.vector(NFSValues[spls_NFS]))
      df = as.data.frame(df)
      colnames(df) = spls_NFS
    }
  }

  df = df %>%
    data.matrix() %>%
    reshape2::melt() %>%
    dplyr::rename("clusters" = "Var1", "samples" = "Var2")

  df2 = do.call(rbind, lapply(names(comparisonSpl), function(splName) { data.frame(samples = comparisonSpl[[splName]], idx = 1:length(comparisonSpl[[splName]]), group = splName, ratio = "num") }))
  if (!is.null(comparisonSplFCref)) {
    df2 = rbind(df2, do.call(rbind, lapply(names(comparisonSplFCref), function(splName) { data.frame(samples = comparisonSplFCref[[splName]], idx = 1:length(comparisonSplFCref[[splName]]), group = splName, ratio = "denom") })))
    df = df %>% merge(df2, by = "samples") %>% dplyr::group_by(clusters, idx, group) %>% summarise(value = getValue(value[ratio=="num"], value[ratio=="denom"])) %>% arrange(idx)
  } else {
    df = df %>% merge(df2, by = "samples") %>% arrange(idx)
  }

  # df = merge(df, CYTdata@sampleData@sampleMetadata %>%
  #              rownames_to_column("samples"), by = "samples")
  #
  # plotlist = list()
  # for (i in 1:length(unique(df$clusters))) {
  #   plotlist[[i]] = df %>%
  #     dplyr::filter(clusters==unique(df$clusters)[i]) %>%
  #     ggplot(aes(x=group, y=value, fill=group)) +
  #     ggplot2::geom_point(ggplot2::aes(group = group, shape = Study),
  #                           position = ggplot2::position_dodge(0.5), color = "black", size = 5)
  #     geom_boxplot()
  #
  #     ggboxplot(, x = "group", y = "value",
  #                             title = unique(df$clusters)[i],
  #                             color = "group", palette = "jco",
  #                             add = "jitter") + stat_compare_means() + stat_compare_means(method = twogroupTest)
  # }
  # library(gridExtra)
  # n <- length(plotlist)
  # nCol <- floor(sqrt(n))
  # do.call("grid.arrange", c(plotlist, ncol=nCol))
  #
  # ggarrange(plotlist=plotlist, ncol = 3) %>%
  #   ggexport(filename = paste(differentialTestingTitle, ".pdf", sep=""))

  p.adjust = match.arg(p.adjust)
  checkmate::qassert(p.adjust, "S1")
  checkmate::qassert(paired, "B1")
  if (length(comparisonSpl)==2) {
    twogroupTest = match.arg(twogroupTest)
    checkmate::qassert(twogroupTest, "S1")

    if (is.null(referenceSplName)) { referenceSplName = names(comparisonSpl)[1] }
    comparisonSplName = setdiff(names(comparisonSpl), referenceSplName)

    switch(twogroupTest,
           wilcox.test = { test_fct <- rstatix::wilcox_test },
           t.test = { test_fct <- rstatix::t_test },
           permut.test = { test_fct <- coin::oneway_test(distribution="exact") })

    print(df)

    statsTest = df %>%
      dplyr::group_by(clusters) %>%
      arrange(idx) %>%
      test_fct(value ~ group, paired = paired) %>%
      mutate(.y. = valueVar, group1 = comparisonSplName, group2 = referenceSplName) %>%
      rstatix::add_significance("p") %>%
      rstatix::adjust_pvalue(p.col = "p", output.col = "p.adj", method = p.adjust) %>%
      rstatix::add_significance("p.adj")

    #print(statsTest)

    foldchange = df %>%
      mutate(group = ifelse(group == referenceSplName, "referenceSpl", "comparisonSpl")) %>%
      dplyr::group_by(clusters) %>%
      dplyr::summarise(`log2FC(g1/g2)` = log(mean(value[group == "comparisonSpl"])/
                                               mean(value[group == "referenceSpl"]), 2))

    stats = merge(foldchange, statsTest, by="clusters") %>% arrange(clusters)
    stats$title = rep(differentialTestingTitle, nrow(stats))
    stats$clustering = rep(clustering, nrow(stats))
    if (is.null(stats$p.adj)) { stats$p.adj = rep(NA, nrow(stats)) }

  } else if (length(comparisonSpl)>2) {

    moregroupTest = match.arg(moregroupTest)
    checkmate::qassert(moregroupTest, "S1")
    posthocTest = match.arg(posthocTest)
    checkmate::qassert(posthocTest, "S1")

  } else { stop() }

  if (nrow(CYTdata@differentialTesting@data) == 0) {
    if(verbose) cat("\n\nCreating differentialTesting data.frame with ", differentialTestingTitle, " statistical results..")
    CYTdata@differentialTesting@data = stats
  } else {
    if(verbose) cat("\n\nUpdating differentialTesting data.frame with ", differentialTestingTitle, " statistical results..")
    CYTdata@differentialTesting@data = rbind.data.frame(CYTdata@differentialTesting@data, stats)
  }
  CYTdata@differentialTesting@data = CYTdata@differentialTesting@data[,c("title", "clustering", "clusters",
                                                                         ".y.", "group1", "n1", "group2", "n2", "log2FC(g1/g2)",
                                                                         "statistic", "p", "p.signif", "p.adj", "p.adj.signif")]

  CYTdata = checkValidity(CYTdata, mode = "warning")
  return(CYTdata)
}


#' @title Volcano Plot for Differential Testing Results
#'
#' @description
#' Generates a volcano plot based on the differential testing results from `CYTdata`. The plot visualizes the log2 fold change (log2FC)
#' between two groups (`g1` and `g2`) and the significance based on p-values or adjusted p-values. The points are color-coded based on
#' significance or phenotype clusters.
#'
#' @param CYTdata An object containing the results of differential testing (class: CYTdata).
#' @param differentialTestingTitle The title of the differential testing result to visualize (string).
#' @param thPval The p-value threshold for significance (numeric, default is 0.05).
#' @param thFC The fold change threshold for significance (numeric, default is 2).
#' @param displayAdjust A logical value indicating whether adjusted p-values should be displayed (default is FALSE).
#' @param errorNA A logical value specifying whether to stop execution when there are missing adjusted p-values (default is FALSE).
#' @param phenoColor A logical value indicating whether the points should be color-coded by phenotype clusters (default is TRUE).
#' @param signifOnly A logical value indicating whether to only plot significant points (default is TRUE).
#'
#' @return A ggplot2 object representing the volcano plot.
#' @details
#' The function will filter the data based on the given title (`differentialTestingTitle`). It computes the log-transformed fold
#' change (log2FC) and p-values (`log10Pval`). Points in the volcano plot are color-coded based on significance thresholds
#' for p-value and fold change.
#'
#' The plot can show all points, but when `signifOnly` is set to TRUE, only the significant points are shown (colored
#' based on whether they are up-regulated or down-regulated).
#'
#' The plot also includes dashed lines for the thresholds on the x-axis (log2FC) and y-axis (-log10(p-value)).
#'
#' @examples
#' # Example of creating a volcano plot
#' volcano_plot <- plotVolcano(
#'   CYTdata = CYTdata,
#'   differentialTestingTitle = "MyDifferentialTestingResult",
#'   thPval = 0.05,
#'   thFC = 2,
#'   displayAdjust = TRUE,
#'   phenoColor = TRUE,
#'   signifOnly = TRUE
#' )
#' print(volcano_plot)
#'
#' @seealso
#' \code{\link{CYTdata}}, \code{\link{ggplot2}}, \code{\link{ggrepel}}
#'

plotVolcano <- function(CYTdata,
                        differentialTestingTitle,
                        thPval = 0.05, thFC = 2,
                        displayAdjust = FALSE, errorNA = FALSE,
                        phenoColor = TRUE, signifOnly = TRUE) {

  CYTdata = checkValidity(CYTdata, mode = "error")

  stats = CYTdata@differentialTesting@data

  checkmate::qassert(differentialTestingTitle, "S1")
  stats = stats %>% dplyr::filter(title == differentialTestingTitle)
  if (nrow(stats)==0){
    stop("Error : No differential abundant analysis has been performed with ", differentialTestingTitle, " title. Nothing can be vizualized")
  }
  g1 = unique(stats$group1)
  if (length(g1)>1){ stop() }
  g2 = unique(stats$group2)
  if (length(g2)>1){ stop() }
  clustering = unique(stats$clustering)
  if (length(clustering)>1){ stop() }

  checkmate::qassert(thPval, "N1")
  if (!(thPval>=0 && thPval<=1)) {
    stop("Error : 'thPval' argument is a p-value threshold and must be a positive numeric between 0 and 1.")
  }
  checkmate::qassert(thFC, "N1")
  if (thFC<=1) { stop("Error : 'thFC' argument must be a numerical value greater than 1.") }

  checkmate::qassert(errorNA, "B1")
  checkmate::qassert(displayAdjust, "B1")
  if (displayAdjust) {
    if (any(is.na(stats$p.adj))) {
      if (errorNA) {
        stop("Error : 'displayAdjust' argument set to TRUE, 'errorNA' argument set to TRUE but missing adjusted p-value computed for ", differentialTestingTitle, ".")
      }
      else {
        message("Warning : 'displayAdjust' argument set to TRUE, 'errorNA' argument set to FALSE and missing adjusted p-value computed for ", differentialTestingTitle,
                ". So Following clusters are removed from volcano plot : ", paste0(stats$clusters[is.na(stats$p.adj)], collapse=", "))
        stats = stats %>% drop_na(p.adj)
      }
    }
    stats$log10Pval = -log10(stats$p.adj)
    ylab = "-log10(pvalue-adj)"
  } else {
    if (any(is.na(stats$p.adj))) {
      if (errorNA) {
        stop("Error : 'errorNA' argument set to TRUE but missing p-value computed for ", differentialTestingTitle, ".")
      }
      else {
        message("Warning : 'errorNA' argument set to FALSE and missing p-value computed for ", differentialTestingTitle,
                ". So Following clusters are removed from volcano plot : ", paste0(stats$clusters[is.na(stats$p.adj)], collapse=", "))
        stats = stats %>% drop_na(p)
      }
    }
    stats$log10Pval = -log10(stats$p)
    ylab = "-log10(pvalue)"
  }

  stats$log2FC = stats$`log2FC(g1/g2)`
  stats$dir = "ns"
  stats$dir[stats$log10Pval > -log10(thPval) & stats$log2FC > log2(thFC)] = "up" #green
  stats$dir[stats$log10Pval > -log10(thPval) & stats$log2FC < -log2(thFC)] = "down" #red

  checkmate::qassert(phenoColor, "B1")
  checkmate::qassert(signifOnly, "B1")
  if (phenoColor) {
    palette = CYTdata@clusteringData@clusteringPalette[[clustering]]
    lev = levels(CYTdata@clusteringData@cellClustering[,clustering])
    stats$clusters = factor(stats$clusters, levels = lev[lev %in% stats$clusters])

    if (signifOnly) {
      signifStats = stats %>% dplyr::filter(dir!="ns")
      nosignifStats = stats %>% dplyr::filter(dir=="ns")
      plot <- ggplot2::ggplot() +
        ggplot2::geom_point(data = signifStats,
                            ggplot2::aes_string(x="log2FC", y="log10Pval", fill="clusters"),
                            shape = 21, stroke = 0, size=8) +
        ggplot2::scale_fill_manual(values = palette) +
        ggrepel::geom_text_repel(data = signifStats,
                                 ggplot2::aes_string(x = "log2FC", y = "log10Pval", label = "clusters"),
                                 color = "black", size = 6, force = 4) +
        ggplot2::geom_point(data = nosignifStats,
                            ggplot2::aes_string(x="log2FC", y="log10Pval"),
                            shape = 21, stroke = 0, fill="grey", size=5)
    }
    else {
      stats = dplyr::mutate_at(stats, "clusters", as.character)
      plot <- ggplot2::ggplot() +
        ggplot2::geom_point(data = stats,
                            ggplot2::aes_string(x="log2FC", y="log10Pval", fill="clusters"),
                            shape = 21, stroke = 0, size=8) +
        ggplot2::scale_fill_manual(values = palette) +
        ggrepel::geom_text_repel(data = stats,
                                 ggplot2::aes_string(x = "log2FC", y = "log10Pval", label = "clusters"),
                                 color = "black", size = 6, force = 4)
    }
  } else {
    signifStats = stats %>%
      dplyr::filter(dir != "ns") %>%
      dplyr::mutate_at("popId", as.character)
    plot <- ggplot2::ggplot() +
      ggplot2::geom_point(data = stats,
                          ggplot2::aes_string(x="log2FC", y="log10Pval", fill="dir"),
                          shape = 21, stroke = 0, size=8) +
      ggplot2::scale_fill_manual(values = c("up" = "green", "down" = "red", "ns" = "grey")) +
      ggrepel::geom_text_repel(data = signifStats,
                               ggplot2::aes_string(x = "log2FC", y = "log10Pval", label = "clustering"),
                               color = "black", size = 6, force = 4)
  }

  max.fc = max(abs(stats$log2FC[is.finite(stats$log2FC)]))
  plot <- plot +
    ggplot2::geom_hline(yintercept = -log10(thPval), linetype = "dashed") +
    ggplot2::geom_vline(xintercept = c(log2(thFC), -log2(thFC)), linetype = "dashed") +
    ggplot2::labs(title = paste(differentialTestingTitle, "(Volcano plot representation)", sep=" "),
                  subtitle = paste0(c(g1, "= group1 (right) and", g2, "= group2 (left)"), collapse=" ")) +
    ggplot2::xlab("log2FC (g1/g2)") + ggplot2::ylab(ylab) +
    ggplot2::scale_x_continuous(limits = c(-max.fc, max.fc), breaks = seq(-10, 10, 1)) +
    ggplot2::scale_y_continuous(breaks = c(seq(0, 10, 1), -log10(thPval))) +
    ggplot2::guides(color=guide_legend(title=clustering)) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic"),
                   aspect.ratio = 1,
                   legend.position = "bottom")
  return(plot)
}


#' @title Plot Piechart for Cell Abundance by Clustering
#'
#' @description
#' Generates a pie chart to visualize the cell abundance across different clusters, samples, and metadata categories.
#' The pie chart is stratified based on metadata values for X-axis and Y-axis, and the data can be optionally weighted
#' by a set of custom NFS values.
#'
#' @param CYTdata An object containing cytometry data (class: CYTdata).
#' @param clustering The name of the clustering column used for grouping cells in the plot (string).
#' @param clusters A vector of clusters to include in the plot. If `NULL`, all clusters will be included (default is `NULL`).
#' @param samples A vector of sample names to include in the plot. If `NULL`, all samples are used (default is `NULL`).
#' @param XaxisMetadata The metadata variable to be used for the X-axis (string, default is "Timepoint").
#' @param YaxisMetadata The metadata variable to be used for the Y-axis (string, default is "Individual").
#' @param NFSValues A named vector of numeric values for scaling cell counts (optional, default is `NULL`).
#' @param pieSize A logical value indicating whether the size of the pie chart slices should be proportional to the total value (default is `TRUE`).
#' @param colorBar A string representing the color for the border of the pie chart slices (hexadecimal color, default is "#000000").
#'
#' @return A ggplot2 object representing the pie chart.
#' @details
#' The pie chart is constructed by first calculating the sum of cell counts for each combination of X-axis and Y-axis metadata values.
#' The chart is faceted by the levels of the X-axis and Y-axis metadata, with each slice representing the abundance of cells in a given group.
#' If `NFSValues` is provided, the cell counts are weighted by these values. The color of the slices is determined by the clustering of cells.
#'
#' If `pieSize` is set to `TRUE`, the size of each pie slice will be proportional to the sum of the cell counts for the group.
#'
#' @examples
#' # Example of creating a pie chart
#' piechart_plot <- plotPiechart(
#'   CYTdata = CYTdata,
#'   clustering = "Cluster1",
#'   samples = c("Sample1", "Sample2"),
#'   XaxisMetadata = "Timepoint",
#'   YaxisMetadata = "Individual",
#'   pieSize = TRUE,
#'   colorBar = "#FF5733"
#' )
#' print(piechart_plot)
#'
#' @seealso
#' \code{\link{CYTdata}}, \code{\link{ggplot2}}
#'

plotPiechart <- function(CYTdata,

                         clustering,
                         clusters = NULL,

                         samples = NULL,
                         XaxisMetadata = "Timepoint",
                         YaxisMetadata = "Individual",

                         NFSValues = NULL,

                         pieSize = TRUE,
                         colorBar = "#000000") {

  CYTdata = checkValidity(CYTdata, mode = "error")

  checkmate::qassert(samples, c(0, "S*"))
  samples = checkorderSamples(CYTdata, samples = samples, order = TRUE, checkDuplicates = TRUE)

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
  }


  if (nrow(CYTdata@sampleData@sampleMetadata)==0){
    stop("Error in CYTdata object, sampleMetadata dataframe from sampleData slot is empty but necessary.")
  }
  checkmate::qassert(XaxisMetadata, "S1")
  if (!XaxisMetadata %in% colnames(CYTdata@sampleData@sampleMetadata)) { stop() }
  checkmate::qassert(YaxisMetadata, "S1")
  if (!YaxisMetadata %in% colnames(CYTdata@sampleData@sampleMetadata)) { stop() }
  if (XaxisMetadata == YaxisMetadata){ stop("Error : 'XaxisMetadata' and 'YaxisMetadata' arguments are the same metadata. Must be different.") }

  metadata = CYTdata@sampleData@sampleMetadata[rownames(CYTdata@sampleData@sampleMetadata) %in% samples,, drop=FALSE]
  XMetadatalevels = levels(CYTdata@sampleData@sampleMetadata[[XaxisMetadata]])
  XMetadatalevels = XMetadatalevels[XMetadatalevels %in% metadata[[XaxisMetadata]]]
  YMetadatalevels = levels(CYTdata@sampleData@sampleMetadata[[YaxisMetadata]])
  YMetadatalevels = YMetadatalevels[YMetadatalevels %in% metadata[[YaxisMetadata]]]
  comblevels = expand.grid(XMetadatalevels, YMetadatalevels)

  checkmate::qassert(NFSValues, c("0", "N*"))
  if (!is.null(NFSValues)) {
    df = getCellCount(CYTdata, samples = samples, clustering = clustering, clusters = clusters, type = "abundance")

    if (is.null(names(NFSValues))) { stop("Error : NFSValues argument is not a named vector, please see documentation.") }
    if (sum(is.na(NFSValues))>0) { stop("Error : NFSValues does not contain some samples (", paste0(names(NFSValues)[is.na(NFSValues)], collapse=", "), ") as named elements.") }

    spls_NFS = colnames(df)
    splErr = setdiff(spls_NFS, names(NFSValues))
    if (length(splErr)>0) { stop("Error : NFSValues argument's names not providing all the samples needed by samples argument (", paste0(splErr, collapse=", "), ").") }

    df = as.matrix(df) %*% diag(as.vector(NFSValues[spls_NFS]))
    df = as.data.frame(df)
    colnames(df) = spls_NFS
  }
  else {
    df = getCellCount(CYTdata, samples = samples, clustering = clustering, clusters = clusters, type = "cellcount")
  }

  checkmate::qassert(pieSize, "B1")
  checkmate::qassert(colorBar, "S1")
  if (!areColors(colorBar)) { stop("Error : 'colorBar' argument must be a hexadecimal color (ex : #FFFFFF).") }

  data = data.frame()
  for (i in 1:nrow(comblevels)) {
    spls = metadata %>% subset((metadata[[XaxisMetadata]] == comblevels[i,1]) & (metadata[[YaxisMetadata]] == comblevels[i,2])) %>% rownames()
    valueSpl = apply(df[,spls, drop=FALSE], 1, sum)
    dataPie = data.frame("value" = valueSpl, "group" = rownames(df),
                         "XaxisMetadata" = rep(comblevels[i,1], length(valueSpl)),
                         "YaxisMetadata" = rep(comblevels[i,2], length(valueSpl)))
    if (pieSize) { dataPie$size = rep(sum(dataPie$value), nrow(dataPie)) }
    else { dataPie$size = rep(1, nrow(dataPie)) }
    data = rbind.data.frame(data, dataPie)
  }
  data$group = factor(data$group, level = rownames(df))

  plot = ggplot2::ggplot(data, ggplot2::aes(x= size/2, y = value, fill = group, width = size)) +
    ggplot2::geom_bar(position="fill", stat="identity", color = colorBar, width = 1) +
    ggplot2::scale_fill_manual(values = CYTdata@clusteringData@clusteringPalette[[clustering]]) +
    ggplot2::facet_grid(factor(YaxisMetadata, levels=YMetadatalevels) ~ factor(XaxisMetadata, levels=XMetadatalevels)) +
    ggplot2::coord_polar("y", start = 0) +
    ggplot2::ylab("") + ggplot2::xlab("") +
    ggplot2::guides(fill=guide_legend(title=clustering)) +
    ggplot2::theme_classic() +
    ggplot2::theme(strip.text.x = ggplot2::element_text(size=12, face="bold.italic"),
                   strip.text.y = ggplot2::element_text(size=12, face="bold.italic"),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic"),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.grid  = ggplot2::element_blank(),
                   aspect.ratio = 1,
                   legend.position = "bottom")

  return(plot)
}
