
#Represent pval interaction


# @title Internal - Computes the number of cells for each population
#
# @description This function is used internally to computes the number of cells for each population.
#
# @param population a character vector providing the population to count the cells belonging to
# @param samples a character vector providing for each cell the associated biological sample
#
# @return a dataframe containing the numbers of cells associated for each population for each sample (rownames = population / colnames = samples)
#

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

#' @title Computes differential analysis statistics for cell clusters
#'
#' @description This function aims to identify of differentially abundant clusters.
#'
#' Such clusters correspond to cell clusters having abundances statistically different between two biological comparisonSpls.
#' The statistical test used for the comparisons can be defined by users.
#' For each cluster, the p-value, log2 fold-change and effect size relative to the reference comparisonSpl are computed.
#' Statistical comparison can be performed in a paired and unpaired manner.
#' Computed p-values can be corrected for multiple testing.
#'
#' @param CYTdata a CYTdata object
#' @param comparisonSpl a character value providing the name of the comparisonSpl to be compared
#' @param referenceSpl a character value providing the name of reference comparisonSpl
#' @param comparison a character value providing the name of comparison
#' @param test.statistics a character value providing the type of statistical test to use. Possible values are: 'wilcoxon' or 't-test'
#' @param p.adjust a character value providing the type of p-value adjustment  to use.
#' Possible values are: 'none', 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr'
#' By default, no p-value adjustement made
#' @param paired a boolean value indicating if individuals are paired
#'
#' @return a S4 object of class 'CYTdata'
#'
#' @export
#' @import rstatix
#'

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
    CYTdata@differentialTesting@data = CYTdata@differentialTesting@data %>% filter(title != differentialTestingTitle)
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
    rename("clusters" = "Var1", "samples" = "Var2")

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
  #     filter(clusters==unique(df$clusters)[i]) %>%
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


#' @title Plots of a volcano plot of statistical analysis
#'
#' @description This function aims to visualize the results of a differential abundant analysis using a Volcano plot
#' In such representation, each in dot corresponds to a pheotpic family and dots are positioned in two dimensional space where the X-axis represents the log2(fold-change) and the Y-axis represents the -log10 of the p-value.
#' Un horizontal line is displayed accordingly to the p-value threshold and to vertical lines are displayed accordingly to the fold-change threshold.
#' The metaclusters are plotted and colored according to their color in the metaclustering heatmap
#'
#' @param CYTdata a CYTdata object
#' @param comparison a character value containing the comparison to study
#' @param thPval a numeric value containing the p-value threshold to use
#' @param thFC a numeric value containing the fold-change threshold to use
#'
#'
#' @return a ggplot2 object
#'
#' @export
#'

plotVolcano <- function(CYTdata,
                        differentialTestingTitle,
                        thPval = 0.05, thFC = 2,
                        displayAdjust = FALSE, errorNA = FALSE,
                        phenoColor = TRUE, signifOnly = TRUE) {

  CYTdata = checkValidity(CYTdata, mode = "error")

  stats = CYTdata@differentialTesting@data

  checkmate::qassert(differentialTestingTitle, "S1")
  stats = stats %>% filter(title == differentialTestingTitle)
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
      signifStats = stats %>% filter(dir!="ns")
      nosignifStats = stats %>% filter(dir=="ns")
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


#' @title Plot pie chart of abundace
#
#' @description This function aims to visualize the composition of differents cell subset according to metadata information by plotting pie chart with each part being a metaluster
#'
#' It is possible to plot the abundance profile of each metacluster within each kinetic family or
#' plot the mean abundance profile of metaclusters within each kinetic family
#'
#' @param CYTdata a CYTdata object
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param supermetaclusters a character vector containing the supermetaclusters to plot in the pie chart. By default, "all the supermetaclusters are used
#' @param color.palette color palette if want to change (one color by part, same length than number of metaclusters)
#' @param Xaxis a string being the metadata used to plot pie charts on the horizontal line. By default, "Timepoint" metadata is selected
#' @param Yaxis a string being the metadata used to plot pie charts on the vertical line. By default, "Individual" metadata is selected
#'
#' @return a ggplot2 object
#'
#' @export
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


























#' @title Plots cell cluster abundances using a boxplot representation
#'
#' @description This function aims to visualize and compare the cell cluster abundances for each biological condition using boxplot and jitter representations.
#'
#' The abundance of a specific cell cluster or a set of cell clusters can be displayed.
#' The representation can be restricted to a specific set of samples.
#' Moreover, boxplot can be constructed based on sample meta information.
#' Statistic can be computed for all comparisons.
#'
#' @param CYTdata a CYTdata object
#' @param population a character vector containing the identifiers of the population to use
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param observationMetadata a character value containing the biological condition to use.
#' The name of the condition must be the one of a metadat condition (column names of CYTdata's metadat slot). By default, metadata "Timepoint" is used
#' @param Yvalue a character value containing the parameters to use. Possible value are percentage (default) or absolute.
#' @param group.by a character value containing the biological condition to group boxplot by. If set to NULL (default), simple boxplot is generated
#' @param test.statistics a character value providing the type of statistical test to use. Possible values are: 'wilcox.test' or 't.test'
#' @param paired a boolean value indicating if a paired or unpaired comparison should be applied
#' @param hideNS a boolean value indicating if non-significant p-value must be hidden
#'
#' @return a ggplot2 object
#'
#' @export
#'

plotBoxplot <- function(CYTdata,
                        population,
                        level = c("clusters", "metaclusters"),
                        samples = NULL,
                        observationMetadata = "Timepoint",
                        groupMetadata = NULL,
                        pointMetadata = NULL,
                        Yvalue = c("percentage", "absolute"),
                        computePval = FALSE,
                        test.statistics = c("wilcox.test", "t.test"),
                        paired = FALSE,
                        hideNS = TRUE) {

  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }

  if (ncol(CYTdata@metadata)==0) {
    stop("Error : Missing metadata slot for 'CYTdata' argument.")
  }
  checkmate::qassert(observationMetadata, "S1")
  checkmate::qassert(groupMetadata, c("0", "S1"))
  checkmate::qassert(pointMetadata, c("0", "S1"))
  if (!observationMetadata %in% colnames(CYTdata@metadata)){
    stop("Error : 'observationMetadata' argument is not a metadata (ex : 'Timepoint', 'Individual')")
  }
  if (!is.null(groupMetadata)){
    if (!groupMetadata %in% colnames(CYTdata@metadata)){
      stop("Error : 'groupMetadata' argument is not a metadata (ex : 'Timepoint', 'Individual')")
    }
    if (groupMetadata==observationMetadata){
      stop("Error : Arguments 'observationMetadata' and 'groupMetadata' are the same metadata")
    }
  }
  if (!is.null(pointMetadata)){
    if (!pointMetadata %in% colnames(CYTdata@metadata)){
      stop("Error : 'pointMetadata' argument is not a metadata (ex : 'Timepoint', 'Individual')")
    }
    if (pointMetadata==observationMetadata){
      stop("Error : Arguments 'observationMetadata' and 'pointMetadata' are the same metadata")
    }
  }

  level = match.arg(level)
  checkmate::qassert(level, "S1")
  Yvalue = match.arg(Yvalue)
  checkmate::qassert(Yvalue, "S1")
  data = getRelativeAbundance(CYTdata, population, level, samples, Yvalue)

  if (Yvalue == "percentage") { Ylab = "Abundance of population" }
  else { Ylab = "Absolute count of population" }

  checkmate::qassert(computePval, "B1")
  test.statistics = match.arg(test.statistics)
  checkmate::qassert(test.statistics, "S1")
  checkmate::qassert(paired, "B1")
  checkmate::qassert(hideNS, "B1")
  switch(test.statistics,
         wilcox.test = { test_fct <- rstatix::wilcox_test },
         t.test = { test_fct <- rstatix::t_test })

  data[[observationMetadata]] = droplevels(data[[observationMetadata]])
  colnames(data)[colnames(data)==observationMetadata] = "observationMetadata"
  if (!is.null(groupMetadata)){
    data[[groupMetadata]] = droplevels(data[[groupMetadata]])
    colnames(data)[colnames(data)==groupMetadata] = "groupMetadata"
  }
  if (!is.null(pointMetadata)){
    data[[pointMetadata]] = droplevels(data[[pointMetadata]])
    colnames(data)[colnames(data)==pointMetadata] = "pointMetadata"
    if (length(unique(data[[pointMetadata]]))>6) {
      message("Warning : 'pointMetadata' argument's corresponding metadata has more
            than 6 different values which is greater than the number of point shapes
            available in ggplot2/geom_point. pointMetadata not taken into account and
            argument set to NULL.")
      pointMetadata = NULL
    }
  }

  if (!is.null(groupMetadata)){

    plot <- ggplot2::ggplot(data, ggplot2::aes(x = observationMetadata, y = value)) +
      ggplot2::geom_boxplot(ggplot2::aes(color = groupMetadata),
                            size=0.1, fatten=1, linewidth=1, outlier.shape=NA, width=0.3,
                            position=position_dodge(0.5)) +
      ggplot2::labs(fill = groupMetadata) +
      ggplot2::guides(color = ggplot2::guide_legend(title = groupMetadata, override.aes = list(size=3, shape=NA)))

    if (computePval) {
      statGroup = data %>%
        group_by(observationMetadata) %>%
        test_fct(value ~ groupMetadata, paired = paired) %>%
        rstatix::add_xy_position(x = "observationMetadata", dodge = 0.8) %>%
        rstatix::add_significance("p")
      #print(statGroup)
      statGroup = statGroup %>%
        filter(p.signif != "ns") %>%
        filter(!is.na(p.signif))
      plot <- plot + ggpubr::stat_pvalue_manual(data = statGroup,
                                                label = "p = {p.signif}",
                                                hide.ns = hideNS, color = "black", size = 5, tip.length = 0)
    }

    if (!is.null(pointMetadata)){
      plot <- plot + ggplot2::geom_point(ggplot2::aes(group = groupMetadata, shape = pointMetadata),
                                         position = ggplot2::position_dodge(0.5), color = "black", size = 5) +
        ggplot2::scale_shape_manual(values=1:nlevels(data$pointMetadata)) +
        ggplot2::labs(shape = pointMetadata)
    } else {
      plot <- plot + ggplot2::geom_point(ggplot2::aes(group = groupMetadata),
                                         position = ggplot2::position_dodge(0.5), color = "black", shape = 16, size = 5)
    }
  } else {
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = observationMetadata, y = value, color = observationMetadata)) +
      ggplot2::geom_boxplot(ggplot2::aes(color = observationMetadata), outlier.shape = NA, linewidth=1, size = 0.2, fatten = 1, width = 0.3, position=position_dodge(0.5)) +
      ggplot2::guides(color = ggplot2::guide_legend(title = observationMetadata, override.aes = list(size=3, shape=NA))) +

      if (!is.null(pointMetadata)){
        plot <- plot + ggplot2::geom_point(ggplot2::aes(group = observationMetadata),
                                           position = ggplot2::position_dodge(0.5), color = "black", size = 5) +
          ggplot2::labs(shape = pointMetadata)
      } else {
        plot <- plot + ggplot2::geom_point(ggplot2::aes(group = observationMetadata),
                                           position = ggplot2::position_dodge(0.5), color = "black", shape = 16, size = 5)
      }
  }

  if (computePval) {
    statObservation = data %>%
      test_fct(value ~ observationMetadata, paired = paired) %>%
      rstatix::add_xy_position(x = "observationMetadata") %>%
      rstatix::add_significance("p") %>%
      filter(p.signif != "ns") %>%
      filter(!is.na(p.signif))

    plot <- plot + ggpubr::stat_pvalue_manual(data = statObservation,
                                              label = "p = {p.signif}",
                                              hide.ns=hideNS, color="black", size=5, tip.length=0.02, step.increase=0.05)
  }

  plot <- plot +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
    ggplot2::labs(title = "Boxplot representation",
                  subtitle = paste0(level, ": ", paste0(population, collapse = ", "))) +
    ggplot2::ylab(Ylab) +
    ggplot2::xlab(observationMetadata) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic"),
                   aspect.ratio = 1,
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   legend.position = "right")

  return(plot)
}

#' @title Plots cell cluster abundances using a boxplot representation
#'
#' @description This function aims to visualize and compare the cell cluster abundances for each biological condition using boxplot and jitter representations.
#'
#' The abundance of a specific cell cluster or a set of cell clusters can be displayed.
#' The representation can be restricted to a specific set of samples.
#' Moreover, boxplot can be constructed based on sample meta information.
#' Statistic can be computed for all comparisons.
#'
#' @param CYTdata a CYTdata object
#' @param population a character vector containing the identifiers of the population to use
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param observationMetadata a character value containing the biological condition to use.
#' The name of the condition must be the one of a metadat condition (column names of CYTdata's metadat slot). By default, metadata "Timepoint" is used
#' @param Yvalue a character value containing the parameters to use. Possible value are percentage (default) or absolute.
#' @param group.by a character value containing the biological condition to group boxplot by. If set to NULL (default), simple boxplot is generated
#' @param test.statistics a character value providing the type of statistical test to use. Possible values are: 'wilcox.test' or 't.test'
#' @param paired a boolean value indicating if a paired or unpaired comparison should be applied
#' @param hideNS a boolean value indicating if non-significant p-value must be hidden
#'
#' @return a ggplot2 object
#'
#' @export
#'

plotBoxLinedPlot <- function(CYTdata,
                             population,
                             level = c("clusters", "metaclusters"),
                             samples = NULL,
                             observationMetadata = "Timepoint",
                             Yvalue = c("percentage", "absolute"),
                             computeStats = TRUE,
                             test.statistics = c("wilcox.test", "t.test"),
                             paired = FALSE,
                             hideNS = TRUE) {

  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }

  level = match.arg(level)
  checkmate::qassert(level, "S1")
  Yvalue = match.arg(Yvalue)
  checkmate::qassert(Yvalue, "S1")
  data = getRelativeAbundance(CYTdata, population, level, samples, Yvalue)
  if (Yvalue == "percentage") { Ylab = "% abundance of population" }
  else { Ylab = "absolute abundance of population" }


  if (ncol(CYTdata@metadata)==0) {
    stop("Error : Missing metadata slot for 'CYTdata' argument.")
  }
  checkmate::qassert(observationMetadata, "S1")
  if (!observationMetadata %in% colnames(CYTdata@metadata)){
    stop("Error : 'observationMetadata' argument is not a metadata (ex : 'Timepoint', 'Individual')")
  }
  if (observationMetadata=="Individual"){
    stop("Error : 'observationMetadata' argument can not be 'Individual' metadata.")
  }
  data[[observationMetadata]] = droplevels(data[[observationMetadata]])
  colnames(data)[colnames(data)==observationMetadata] = "observationMetadata"

  plot <- ggplot2::ggplot(data,
                          ggplot2::aes(x = "observationMetadata",
                                       y = "value",
                                       group = "Individual")) +
    ggplot2::geom_line(ggplot2::aes(color="Individual")) +
    ggplot2::geom_point(ggplot2::aes(shape = "Individual", color="Individual")) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = "observationMetadata"),
                          outlier.shape = NA, size = 0.2, fatten = 1) +
    ggplot2::labs(fill = observationMetadata) +
    ggplot2::labs(group = "Individual") +
    ggplot2::labs(color = "Individual") +
    ggplot2::labs(shape = "Individual")

  if (computeStats) {
    test.statistics = match.arg(test.statistics)
    checkmate::qassert(test.statistics, "S1")
    checkmate::qassert(paired, "B1")
    checkmate::qassert(hideNS, "B1")

    switch(test.statistics,
           wilcox.test = { test_fct <- rstatix::wilcox_test },
           t.test = { test_fct <- rstatix::t_test })
    statObservation = data %>%
      test_fct(value ~ observationMetadata, paired = paired) %>%
      rstatix::add_xy_position(x = "observationMetadata") %>%
      rstatix::add_significance("p") %>%
      filter(p.signif != "ns")

    plot <- plot +
      ggpubr::stat_pvalue_manual(data = statObservation,
                                 label = "p = {p.signif}",
                                 hide.ns = hideNS,
                                 color = "black", size = 5, tip.length = 0.02, step.increase = 0.05)
  }

  plot <- plot +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
    ggplot2::labs(title = "Boxplot representation",
                  subtitle = paste0(level, ": ", paste0(population, collapse = ", "))) +
    ggplot2::ylab(Ylab) +
    ggplot2::xlab(observationMetadata) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic"),
                   aspect.ratio = 1,
                   panel.grid.minor =  ggplot2::element_blank(),
                   panel.grid.major =  ggplot2::element_blank(),
                   legend.position = "right")
  return(plot)
}

#' @title Plots cell cluster abundances using a boxplot representation
#'
#' @description This function aims to visualize and compare the cell cluster abundances for each biological condition using boxplot and jitter representations.
#'
#' The abundance of a specific cell cluster or a set of cell clusters can be displayed.
#' The representation can be restricted to a specific set of samples.
#' Moreover, boxplot can be constructed based on sample meta information.
#' Statistic can be computed for all comparisons.
#'
#' @param CYTdata a CYTdata object
#' @param population a character vector containing the identifiers of the population to use
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param observationMetadata a character value containing the biological condition to use.
#' The name of the condition must be the one of a metadata condition (column names of CYTdata's metadat slot). By default, metadata "Timepoint" is used
#' @param Yvalue a character value containing the parameters to use. Possible value are percentage (default) or absolute.
#' @param group.by a character value containing the biological condition to group boxplot by. If set to NULL (default), simple boxplot is generated
#' @param test.statistics a character value providing the type of statistical test to use. Possible values are: 'wilcox.test' or 't.test'
#' @param paired a boolean value indicating if a paired or unpaired comparison should be applied
#' @param hideNS a boolean value indicating if non-significant p-value must be hidden
#'
#' @return a ggplot2 object
#'
#' @export
#'
#'
plotLineplot2 <- function(CYTdata,
                          population = NULL,
                          level = c("clusters", "metaclusters"),
                          samples = NULL,
                          observationMetadata = "Timepoint",
                          colorMetadata = NULL,
                          Yvalue = c("percentage", "absolute"),
                          stats = TRUE,
                          test.statistics = c("wilcox.test", "t.test"),
                          paired = FALSE,
                          hideNS = TRUE) {

  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }

  if (ncol(CYTdata@metadata)==0) {
    stop("Error : Missing metadata slot for 'CYTdata' argument.")
  }
  checkmate::qassert(observationMetadata, "S1")
  checkmate::qassert(colorMetadata, c("0", "S1"))
  if (!observationMetadata %in% colnames(CYTdata@metadata)){
    stop("Error : 'observationMetadata' argument is not a metadata (ex : 'Timepoint', 'Individual')")
  }
  if (!is.null(colorMetadata)){
    if (!colorMetadata %in% colnames(CYTdata@metadata)){
      stop("Error : 'colorMetadata' argument is not a metadata (ex : 'Timepoint', 'Individual')")
    }
    if (colorMetadata==observationMetadata){
      stop("Error : Arguments 'observationMetadata' and 'colorMetadata' are the same metadata")
    }
  }

  level = match.arg(level)
  checkmate::qassert(level, "S1")
  Yvalue = match.arg(Yvalue)
  checkmate::qassert(Yvalue, "S1")
  data = getRelativeAbundance(CYTdata, population, level, samples, Yvalue)

  if (Yvalue == "percentage") { Ylab = "Abundance of population" }
  else { Ylab = "Absolute count of population" }

  data[[observationMetadata]] = droplevels(data[[observationMetadata]])
  colnames(data)[colnames(data)==observationMetadata] = "observationMetadata"
  if (!is.null(colorMetadata)){
    data[[colorMetadata]] = droplevels(data[[colorMetadata]])
    colnames(data)[colnames(data)==colorMetadata] = "colorMetadata"
  }

  if (!is.null(colorMetadata)){


    plot <- ggplot2::ggplot(data,
                            ggplot2::aes(x = observationMetadata, y = value, group=colorMetadata)) +
      ggplot2::geom_line(size=2) +
      ggplot2::geom_point(size = 8) +
      ggplot2::labs(color = colorMetadata,
                    group = "Individual",
                    shape = colorMetadata)
  }
  else { }

  plot <- plot +
    ggplot2::labs(title = "Lineplot representation",
                  subtitle = paste0(level, ": ", paste0(population, collapse = ", "))) +
    ggplot2::ylab(Ylab) +
    ggplot2::xlab(observationMetadata) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic"),
                   #aspect.ratio = 1,
                   text = ggplot2::element_text(size = 35),
                   axis.title.y = element_text(size = 40, face = "bold"),
                   axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 35, face = "bold"),
                   axis.text.y = ggplot2::element_text(size = 35, face = "bold"),
                   panel.grid.minor =  ggplot2::element_blank(),
                   panel.grid.major =  ggplot2::element_blank(),
                   legend.position = "right")
  return(plot)
}


plotLineplot <- function(CYTdata,
                         population = NULL,
                         level = c("clusters", "metaclusters"),
                         samples = NULL,
                         observationMetadata = "Timepoint",
                         groupMetadata = NULL,
                         colorMetadata = NULL,
                         Yvalue = c("percentage", "absolute"),
                         addJitter = FALSE,
                         stats = TRUE,
                         test.statistics = c("wilcox.test", "t.test"),
                         paired = FALSE,
                         hideNS = TRUE) {

  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }

  if (ncol(CYTdata@metadata)==0) {
    stop("Error : Missing metadata slot for 'CYTdata' argument.")
  }
  checkmate::qassert(observationMetadata, "S1")
  checkmate::qassert(groupMetadata, c("0", "S1"))
  if (!observationMetadata %in% colnames(CYTdata@metadata)){
    stop("Error : 'observationMetadata' argument is not a metadata (ex : 'Timepoint', 'Individual')")
  }
  if (!is.null(groupMetadata)){
    if (!groupMetadata %in% colnames(CYTdata@metadata)){
      stop("Error : 'groupMetadata' argument is not a metadata (ex : 'Timepoint', 'Individual')")
    }
    if (groupMetadata==observationMetadata){
      stop("Error : Arguments 'observationMetadata' and 'groupMetadata' are the same metadata")
    }
  }

  level = match.arg(level)
  checkmate::qassert(level, "S1")
  Yvalue = match.arg(Yvalue)
  checkmate::qassert(Yvalue, "S1")

  checkmate::qassert(addJitter, "B1")

  data = getRelativeAbundance(CYTdata, population, level, samples, Yvalue)

  if (Yvalue == "percentage") { Ylab = "Abundance of population" }
  else { Ylab = "Absolute count of population" }

  data[[observationMetadata]] = droplevels(data[[observationMetadata]])
  colnames(data)[colnames(data)==observationMetadata] = "observationMetadata"
  if (!is.null(groupMetadata)){
    data[[groupMetadata]] = droplevels(data[[groupMetadata]])
    colnames(data)[colnames(data)==groupMetadata] = "groupMetadata"
  }

  if (!is.null(groupMetadata)){

    data = data %>%
      dplyr::group_by(observationMetadata, groupMetadata) %>%
      dplyr::summarize(mean = mean(value, na.rm=TRUE),
                       sd = sd(value, na.rm=TRUE))

    plot <- ggplot2::ggplot(data,
                            ggplot2::aes(x = observationMetadata, y = mean,
                                         group = groupMetadata)) +
      ggplot2::geom_line(ggplot2::aes(color = groupMetadata), size=2) +
      ggplot2::geom_point(ggplot2::aes(shape = groupMetadata, color = groupMetadata), size = 8) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=mean-sd, ymax=mean+sd, color = groupMetadata),
                             size = .8, width=.5, position=position_dodge(0)) +
      ggplot2::labs(color = groupMetadata,
                    group = groupMetadata,
                    shape = groupMetadata)
  }
  else {
    data = data %>%
      dplyr::group_by(observationMetadata) %>%
      dplyr::summarize(mean = mean(value, na.rm=TRUE),
                       sd = sd(value, na.rm=TRUE))
    plot <- ggplot2::ggplot(data,
                            ggplot2::aes(x = observationMetadata, y = mean, group=1)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_path() +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=mean-sd, ymax=mean+sd),
                             size = .1, width=.1, position=position_dodge(0))
  }

  plot <- plot +
    ggplot2::labs(title = "Lineplot representation",
                  subtitle = paste0(level, ": ", paste0(population, collapse = ", "))) +
    ggplot2::ylab(Ylab) +
    ggplot2::xlab(observationMetadata) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic"),
                   #aspect.ratio = 1,
                   text = ggplot2::element_text(size = 35),
                   axis.title.y = element_text(size = 40, face = "bold"),
                   axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 35, face = "bold"),
                   axis.text.y = ggplot2::element_text(size = 35, face = "bold"),
                   panel.grid.minor =  ggplot2::element_blank(),
                   panel.grid.major =  ggplot2::element_blank(),
                   legend.position = "right")
  return(plot)
}





#' @title Plots a proportional stacked area chart of clusters evolution along the samples
#'
#' @description This function aims to visualize a proportional stacked area chart of clusters evolution along the samples
#'
#'
#' @param CYTdata a CYTOF.analysis object
#' @param samples a character vector containing the names of biological samples to use
#' @param clusters a character vector containing the names of clusters to use
#'
#' @return a ggplot2 object
#'
#' @export
#'
#'


plotStackedBarplot <- function(CYTdata,
                               population = NULL,
                               level = c("clusters", "metaclusters"),
                               samples = NULL,
                               observationMetadata = "Timepoint",
                               groupMetadata = NULL,
                               Yvalue = c("ABrelative", "CCabsolute")) {

  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }

  level = match.arg(level)
  checkmate::qassert(level, "S1")
  population = checkorderPopulation(CYTdata, population=population, level=level,
                                    order=TRUE, checkDuplicates=TRUE)
  samples = checkorderSamples(CYTdata, samples, order=TRUE, checkDuplicates=TRUE)

  if (ncol(CYTdata@metadata)==0) {
    stop("Error : Missing metadata slot for 'CYTdata' argument.")
  }
  checkmate::qassert(observationMetadata, "S1")
  checkmate::qassert(groupMetadata, c("0", "S1"))
  if (!observationMetadata %in% colnames(CYTdata@metadata)){
    stop("Error : 'observationMetadata' argument is not a metadata (ex : 'Timepoint', 'Individual')")
  }
  if (!is.null(groupMetadata)){
    if (!groupMetadata %in% colnames(CYTdata@metadata)){
      stop("Error : 'groupMetadata' argument is not a metadata (ex : 'Timepoint', 'Individual')")
    }
    if (groupMetadata==observationMetadata){
      stop("Error : Arguments 'observationMetadata' and 'groupMetadata' are the same metadata")
    }
  }

  Yvalue = match.arg(Yvalue)
  checkmate::qassert(Yvalue, "S1")

  if (level == "clusters"){
    data = CYTdata@Clustering@cellcount[population, samples, drop=FALSE]
    palette = CYTdata@Clustering@palette
    lev = levels(CYTdata@Clustering@clusters)
  }
  else {
    data = CYTdata@Metaclustering@cellcount[population, samples, drop=FALSE]
    palette = CYTdata@Metaclustering@palette
    lev = levels(CYTdata@Metaclustering@metaclusters)
  }

  data = suppressWarnings(reshape::melt(data.matrix(data)))
  colnames(data) = c("population", "samples", "value")
  md = CYTdata@metadata
  md$samples = rownames(md)
  data = merge(data, md, by = "samples")
  data[[observationMetadata]] = droplevels(data[[observationMetadata]])
  colnames(data)[colnames(data)==observationMetadata] = "observationMetadata"
  if (!is.null(groupMetadata)){
    data[[groupMetadata]] = droplevels(data[[groupMetadata]])
    colnames(data)[colnames(data)==groupMetadata] = "groupMetadata"
  }
  if (length(unique(data$observationMetadata))==1) {
    stop("Error : 'samples' argument procides only one level of 'observationMetadata' argument (",
         unique(data$observationMetadata), "). Impossible to build stacked area plot.")
  }
  data$population = factor(data$population, levels = lev)

  if (!is.null(groupMetadata)){
    data = data %>%
      dplyr::group_by(population, observationMetadata, groupMetadata) %>%
      dplyr::summarise(y = sum(value))

    if(Yvalue == "ABrelative") {
      data = data %>%
        dplyr::group_by(observationMetadata, groupMetadata) %>%
        dplyr::mutate(y = y / sum(y))
    }
    plot <- ggplot2::ggplot(data %>% arrange(observationMetadata),
                            ggplot2::aes(fill = population,
                                         y = y,
                                         x = groupMetadata))

    if(Yvalue == "ABrelative") {
      plot <- plot + ggplot2::geom_bar(position="fill", stat="identity", width = 0.75)
      Ylab = "Abundance (Relative)"
    }
    else {
      plot <- plot + ggplot2::geom_bar(position="stack", stat="identity", width = 0.75)
      Ylab =  "Absolute counts"
    }

    plot <- plot + facet_grid(~observationMetadata, switch = "x")
  }
  else {
    data = data %>%
      dplyr::group_by(population, observationMetadata) %>%
      dplyr::summarise(y = sum(value))
    if(Yvalue == "ABrelative") {
      data = data %>%
        dplyr::group_by(observationMetadata) %>%
        dplyr::mutate(y = y / sum(y))
    }
    plot <- ggplot2::ggplot(data,
                            ggplot2::aes(fill = population,
                                         y = y,
                                         x = observationMetadata))

    if(Yvalue == "ABrelative") {
      plot <- plot + ggplot2::geom_bar(position="fill", stat="identity")
      Ylab = "Abundance (Relative)"
    }
    else {
      plot <- plot + ggplot2::geom_bar(position="stack", stat="identity")
      Ylab =  "Absolute counts"
    }

  }

  plot <- plot +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(title = "Stacked Area representation",
                  fill = "Population") +
    ggplot2::ylab(Ylab) +
    #ggplot2::xlab(observationMetadata) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size=55, face="bold"),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic"),
                   panel.grid.minor =  ggplot2::element_blank(),
                   panel.grid.major =  ggplot2::element_blank(),
                   text = ggplot2::element_text(size = 35),
                   axis.title.y = element_text(size = 40, face = "bold"),
                   axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 35, face = "bold"),
                   axis.text.y = ggplot2::element_text(size = 35, face = "bold"),
                   legend.position = "right")
  if (!is.null(groupMetadata)){ plot <- plot + ggplot2::theme(strip.placement = "outside",
                                                              strip.text.x = element_text(size = 40, face = "bold"),
                                                              strip.background = element_rect(fill = NA, color = "white"),
                                                              panel.spacing = unit(1.4, "lines")) }

  return(plot)
}



#' @title Plots a proportional stacked area chart of clusters evolution along the samples
#'
#' @description This function aims to visualize a proportional stacked area chart of clusters evolution along the samples
#'
#'
#' @param CYTdata a CYTOF.analysis object
#' @param samples a character vector containing the names of biological samples to use
#' @param clusters a character vector containing the names of clusters to use
#'
#' @return a ggplot2 object
#'
#' @export
#'

plotStackedArea <- function(CYTdata,
                            population = NULL,
                            level = c("clusters", "metaclusters"),
                            samples = NULL,
                            observationMetadata = "Timepoint",
                            Yvalue = c("ABrelative", "CCabsolute")) {

  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.") }
  else { CYTdata = MakeValid(CYTdata, verbose = TRUE) }

  level = match.arg(level)
  checkmate::qassert(level, "S1")
  population = checkorderPopulation(CYTdata, population=population, level=level,
                                    order=TRUE, checkDuplicates=TRUE)
  samples = checkorderSamples(CYTdata, samples, order=TRUE, checkDuplicates=TRUE)

  if (ncol(CYTdata@metadata)==0) {
    stop("Error : Missing metadata slot for 'CYTdata' argument.")
  }
  checkmate::qassert(observationMetadata, "S1")
  if (!observationMetadata %in% colnames(CYTdata@metadata)){
    stop("Error : 'observationMetadata' argument is not a metadata (ex : 'Timepoint', 'Individual')")
  }

  Yvalue = match.arg(Yvalue)
  checkmate::qassert(Yvalue, "S1")

  if (level == "clusters"){
    data = CYTdata@Clustering@cellcount[population, samples, drop=FALSE]
    palette = CYTdata@Clustering@palette
    lev = levels(CYTdata@Clustering@clusters)
  }
  else {
    data = CYTdata@Metaclustering@cellcount[population, samples, drop=FALSE]
    palette = CYTdata@Metaclustering@palette
    lev = levels(CYTdata@Metaclustering@metaclusters)
  }

  data = suppressWarnings(reshape::melt(data.matrix(data)))
  colnames(data) = c("population", "samples", "value")
  md = CYTdata@metadata
  md$samples = rownames(md)
  data = merge(data, md, by = "samples")
  data[[observationMetadata]] = droplevels(data[[observationMetadata]])
  colnames(data)[colnames(data)==observationMetadata] = "observationMetadata"
  if (length(unique(data$observationMetadata))==1) {
    stop("Error : 'samples' argument procides only one level of 'observationMetadata' argument (",
         unique(data$observationMetadata), "). Impossible to build stacked area plot.")
  }
  data$population = factor(data$population, levels = lev)

  data = data %>%
    dplyr::group_by(population, observationMetadata) %>%
    dplyr::summarise(y = sum(value))
  if(Yvalue == "ABrelative") {
    data = data %>%
      dplyr::group_by(observationMetadata) %>%
      dplyr::mutate(y = y / sum(y))
    Ylab = "Abundance (Relative)"
  }
  else { Ylab =  "Absolute counts" }

  plot <- ggplot2::ggplot(data,
                          ggplot2::aes(x = observationMetadata,
                                       y = y)) +
    ggplot2::geom_area(ggplot2::aes(colour = population,
                                    group = population,
                                    fill = population),
                       alpha=0.6, size=1, colour="white") +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::labs(title = "Stacked Area representation",
                  fill = level) +
    ggplot2::ylab(Ylab) +
    ggplot2::xlab(observationMetadata) +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic"),
                   aspect.ratio = 1,
                   panel.grid.minor =  ggplot2::element_blank(),
                   panel.grid.major =  ggplot2::element_blank(),
                   legend.position = "right")
  return(plot)
}








#' @title Plots a PCA representation based cell cluster abundances
#'
#' @description This function aims to represent a Principal Component Analysis representation based on cell cluster abundances.
#' In such representation, clusters or samples are positioned based on computed principal components.
#' The representation can be displayed based on specific principal components.
#' The representation can be restricted to specific cell clusters and samples. In addition, it is possible to choose the levels displayed, clusters or samples.
#'
#' @param CYTdata a CYTdata object
#' @param levels a character value containing the variable to be displayed. Possible values are: 'both', 'clusters' or 'samples'
#' @param clusters a character vector containing the identifier of the cluster to use. By default, all clusters are used
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param components a numeric vector providing the components to display
#' @param condition.samples a character vector containing the variable to be studied for the samples. Possible values are: 'condition' or 'timepoint"
#' @param cor.radius.th a numeric value specifying the radius of the correlation plot radius
#' @param plot.text a boolean value specifying if adds text directly at the plot
#'
#' @return a ggplot2 object
#'
#' @export
#'
plotPCA <- function(CYTdata,
                    levels = c("both", "clusters", "samples"),
                    clusters = NULL,
                    samples = NULL,
                    components = c(1, 2),
                    condition.samples = c("condition", "timepoint"),
                    cor.radius.th = 0.6,
                    plot.text = TRUE) {

  levels <- match.arg(levels)
  condition.samples <- match.arg(condition.samples)

  checkmate::qassert(levels, "S1")
  checkmate::qassert(clusters, c("0", "S*"))
  checkmate::qassert(samples, c("0", "S*"))
  checkmate::qassert(components, "N2")
  checkmate::qassert(condition.samples, "S1")
  checkmate::qassert(cor.radius.th, "N1")
  checkmate::qassert(plot.text, "B1")

  if (levels != "both" && levels != "clusters" && levels != "samples") {
    stop("The levels name is invalid")
  }

  circleFun <- function(center = c(0, 0), diameter = 1, npoints = 100) {
    r <- diameter / 2
    t <- seq(0, 2 * pi, length.out = npoints)
    x <- center[1] + r * cos(t)
    y <- center[2] + r * sin(t)
    return(data.frame(x = x, y = y))
  }

  matrix.abundance <- CYTdata@matrix.abundance

  if (!is.null(clusters)) {
    matrix.abundance <- matrix.abundance[rownames(matrix.abundance) %in% clusters, ]
  }
  if (!is.null(samples)) {
    matrix.abundance <- matrix.abundance[, colnames(matrix.abundance) %in% samples]
  }

  res.PCA <- FactoMineR::PCA(t(matrix.abundance), graph = FALSE)
  var.explained <- res.PCA$eig[, 2]

  data.clusters <- data.frame(res.PCA$var$coord[, components])
  colnames(data.clusters) <- c("x", "y")
  data.clusters$id <- rownames(matrix.abundance)

  data.variables <- data.frame(res.PCA$ind$coord[, components])
  colnames(data.variables) <- c("x", "y")
  data.variables$id <- colnames(matrix.abundance)

  data.variables$x <- scales::rescale(data.variables$x, to = c(-1, 1))
  data.variables$y <- scales::rescale(data.variables$y, to = c(-1, 1))

  data.variables <- merge(data.variables, CYTdata@metadata, by = "row.names")
  data.variables$Row.names <- NULL

  circle1 <- circleFun(c(0, 0), 2, npoints = 100)
  circle2 <- circleFun(c(0, 0), 2 * cor.radius.th, npoints = 100)

  data.clusters <- data.clusters[sqrt(data.clusters$x^2 + data.clusters$y^2) > cor.radius.th, ]

  xlab <- paste0("PCA", components[1], " (", round(var.explained[components[1]], 2), "%)")
  ylab <- paste0("PCA", components[2], " (", round(var.explained[components[2]], 2), "%)")

  plot <- ggplot2::ggplot()

  if (levels == "clusters" || levels == "both") {
    plot <- plot +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", size = 0.2) +
      ggplot2::geom_path(data = circle1,
                         ggplot2::aes_string(x = "x", y = "y"),
                         size = 0.2, color = "gray") +
      ggplot2::geom_path(data = circle2,
                         ggplot2::aes_string(x = "x", y = "y"),
                         size = 0.2, linetype = "dashed", color = "gray") +
      ggiraph::geom_point_interactive(data = data.clusters,
                                      ggplot2::aes_string(x = "x", y = "y",
                                                          tooltip = "id",
                                                          data_id = "id"),
                                      shape = 21, size = 2,
                                      stroke = 0.1, fill = "black", color = "black")

    if (plot.text == TRUE) {
      plot <- plot +
        ggrepel::geom_text_repel(data = data.clusters,
                                 ggplot2::aes_string(x = "x",
                                                     y = "y",
                                                     label = "id"),
                                 size = 3,
                                 color = "black",
                                 max.overlaps = Inf,
                                 min.segment.length = 0,
                                 segment.color = NA,
                                 segment.size = 0.1)
    }
  }

  if (levels == "samples" || levels == "both") {
    plot <- plot +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", size = 0.2) +
      ggplot2::geom_path(data = circle1,
                         ggplot2::aes_string(x = "x", y = "y"),
                         size = 0.2, color = "gray") +
      ggplot2::geom_path(data = circle2,
                         ggplot2::aes_string(x = "x", y = "y"),
                         size = 0.2, linetype = "dashed", color = "gray") +
      ggiraph::geom_point_interactive(data = data.variables,
                                      ggplot2::aes_string(x = "x", y = "y",
                                                          fill = condition.samples,
                                                          tooltip = "id",
                                                          data_id = "individual"),
                                      shape = 21, size = 2,
                                      stroke = 0.1, color = "black")

    if (plot.text == TRUE) {
      plot <- plot +
        ggrepel::geom_text_repel(data = data.variables,
                                 ggplot2::aes_string(x = "x", y = "y",
                                                     label = "id"),
                                 size = 3,
                                 color = "black",
                                 max.overlaps = Inf,
                                 min.segment.length = 0,
                                 segment.color = NA,
                                 segment.size = 0.1)
    }
  }

  plot <- plot + ggplot2::labs(title = "PCA representation")

  plot <- plot +
    ggplot2::coord_fixed() +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab)

  plot <- plot +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      aspect.ratio = 1,
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(color = NA),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank())

  return(plot)
}


#' @title Plots a MDS representation based on cell cluster abundances/median expression or sample median expression
#'
#' @description This function aims to visualize the similarities between samples or clusters based on their abundances or median expression, using a Multidimensional Scaling representation.
#' Each dot represents a sample or a cluster and the distances between the dots are proportional to the Euclidean distance between these objects.
#' The representation can be restricted to specific cell clusters and samples. In addition, it is possible to choose the levels displayed, clusters or samples.
#'
#' @param CYTdata a S4 object of class 'CYTdata'
#' @param matrix a character vector containing the matrix to be studied. Possible values are: 'abundance' or 'expression' (default = expression)
#' @param levels a character value containing the variable to be displayed. Possible values are: 'metaclusters' or 'samples' (default = metaclusters)
#' @param markers a character vector containing the names of biological markers to use (if matrix set to 'expression'). By default, all markers are used
#' @param metaclusters a character vector containing the identifiers of the metaclusters to use. By default, all metaclusters are used. Only if levels set to 'metaclusters'
#' @param samples a character vector containing the names of biological samples to use. By default, all samples are used
#' @param color.metadata a character vector containing the metadata variable to color the samples points by (if levels set to 'samples')
#' @param printPolygon a boolean indicating if hull is plotted around each of the groups made by color.metadata argument (default = TRUE)
#'
#' @return a ggplot2 object
#'
#' @export
#'

plotMDS <- function(CYTdata,
                    matrix = c("expression", "abundance"),
                    levels = c("samples", "metaclusters"),
                    markers = NULL,
                    metaclusters = NULL,
                    samples = NULL,
                    color.metadata = NULL,
                    printPolygon = TRUE) {

  if (class(CYTdata)!="CYTdata") { stop("Error : argument 'CYTdata' a S4 object of class 'CYTdata'.")}
  matrix = match.arg(matrix)
  checkmate::qassert(matrix, "S1")
  levels = match.arg(levels)
  checkmate::qassert(levels, "S1")
  checkmate::qassert(metaclusters, c("0", "S+"))
  checkmate::qassert(samples, c("0", "S+"))
  checkmate::qassert(color.metadata, c("0", "S1"))
  checkmate::qassert(printPolygon, "B1")

  if (is.null(metaclusters) && length(CYTdata@Metaclustering@metaclusters)>0) { metaclusters = unique(CYTdata@Metaclustering@metaclusters) }
  samples = checkorderSamples(CYTdata, samples, order = TRUE, checkDuplicates = TRUE)
  markers = checkorderMarkers(CYTdata, markers, order = TRUE, checkDuplicates = TRUE)
  if (is.null(samples)) { samples = unique(CYTdata@samples) }
  if (is.null(markers)) { markers = cocCYTdata@markers }

  if (matrix == "abundance") {
    if (length(CYTdata@Metaclustering@metaclusters)==0) { stop("Error in plot.MDS function : 'matrix' arguments set to 'abundance' require
    non empty 'Metaclustering' slots (from CYTdata object)") }
    data.matrix = CYTdata@Metaclustering@abundance
    if (is.null(metaclusters)) { metaclusters = rownames(data.matrix) }
    if (is.null(samples)) { samples = colnames(data.matrix) }
    data.matrix = data.matrix[rownames(data.matrix) %in% metaclusters, colnames(data.matrix) %in% samples]  }
  else {
    if (levels == "metaclusters") {
      if (length(CYTdata@Metaclustering@metaclusters)==0) { stop("Error in plot.MDS function : 'levels' arguments set to 'metaclusters' require
      non empty 'Metaclustering' slots (from CYTdata object)") }
      if (is.null(metaclusters)) { metaclusters = unique(CYTdata@Metaclustering@metaclusters) }
      if (is.null(samples)) { samples = unique(CYTdata@samples) }
      data.matrix = cbind(CYTdata@matrix.expression[, markers],
                          "SPL" = CYTdata@samples,
                          "MC" = CYTdata@Metaclustering@metaclusters)
      data.matrix = subset(data.matrix, (MC %in% metaclusters) & (SPL %in% samples))
      data.matrix$SPL = NULL
      data.matrix = plyr::ddply(data.matrix, "MC",
                                function(x) {
                                  x$MC = NULL
                                  apply(x, 2, stats::median)
                                })
      rownames(data.matrix) = data.matrix$MC
      data.matrix$MC = NULL
    } else {
      data.matrix = cbind(CYTdata@matrix.expression[, markers], "SPL" = CYTdata@samples)
      if (length(CYTdata@Metaclustering@metaclusters)==0){
        if (!is.null(metaclusters)){ message("Warning : Metaclustering slot empty but arguments 'matrix' set to 'expression'
        and 'levels' set to 'samples'. So given 'metaclusters' argument no taken into account") }
      } else{
        if (is.null(metaclusters)) { metaclusters = unique(CYTdata@Metaclustering@metaclusters) }
        data.matrix = cbind(data.matrix, "MC" = CYTdata@Metaclustering@metaclusters)
        data.matrix = subset(data.matrix, MC %in% metaclusters)
        data.matrix$MC = NULL
      }
      data.matrix = subset(data.matrix, SPL %in% samples)
      data.matrix = plyr::ddply(data.matrix, "SPL",
                                function(x) {
                                  x$SPL = NULL
                                  apply(x, 2, stats::median)
                                })
      rownames(data.matrix) = data.matrix$SPL
      data.matrix$SPL = NULL
      data.matrix = t(data.matrix)
    }
  }

  plot <- ggplot2::ggplot()

  if (levels == "metaclusters") {
    colors.metaclusters = CYTdata@Metaclustering@heatmap.object$metaclusters.colors

    fit1 = data.matrix %>% stats::dist() %>% MASS::isoMDS(k = 2, trace = FALSE)
    x1 = fit1$points[, 1]
    y1 = fit1$points[, 2]
    min.lim = min(min(x1), min(y1)) * 1.1
    max.lim = max(max(x1), max(y1)) * 1.1

    proj.metaclusters = data.frame(x = x1, y = y1)
    proj.metaclusters$MC = rownames(data.matrix)

    plot <- plot + ggplot2::geom_point(data = proj.metaclusters,
                                       ggplot2::aes_string(x = "x", y = "y", col = "MC"), size = 1) +
      ggplot2::labs(title = "Multidimensional Scaling",
                    subtitle = paste0("Kruskal's stress = ", round(fit1$stress, 2))) +
      ggplot2::geom_hline(yintercept = (min.lim + max.lim) / 2, linetype = "dashed") +
      ggplot2::geom_vline(xintercept = (min.lim + max.lim) / 2, linetype = "dashed") +
      ggrepel::geom_text_repel(data = proj.metaclusters,
                               ggplot2::aes_string(x = "x", y = "y",
                                                   label = "MC", color = "MC"),
                               size = 3, force = 3) +
      ggplot2::scale_color_manual(values = colors.metaclusters)
  }
  else {

    if (!is.null(color.metadata)){
      if (!color.metadata %in% colnames(CYTdata@metadata)){ stop("Error : 'color.metadata' argument
                                                                  is not a biological condition present in metadata (",
                                                                 paste0(colnames(CYTdata@metadata), collapse=","), ")") }
    } else { color.metadata = "samples" }

    fit2 = t(data.matrix) %>% stats::dist() %>% MASS::isoMDS(k = 2, trace = FALSE)
    x2 = fit2$points[, 1]
    y2 = fit2$points[, 2]
    min.lim = min(min(x2), min(y2)) * 1.1
    max.lim = max(max(x2), max(y2)) * 1.1

    proj.samples = data.frame(x = x2, y = y2)
    proj.samples$samples = colnames(data.matrix)

    md = CYTdata@metadata
    md$samples = rownames(md)
    proj.samples = merge(proj.samples, md, by = "samples")

    hull.data = c()
    for (grp in unique(proj.samples[[color.metadata]])){
      proj.grp = proj.samples[proj.samples[[color.metadata]] == grp,]
      hull.coord = grDevices::chull(proj.grp[,c("x","y")])
      hull.grp = proj.grp[hull.coord,]
      hull.data = rbind(hull.data, hull.grp)
    }
    plot <- plot +
      ggplot2::geom_point(data = proj.samples,
                          ggplot2::aes_string(x = "x", y = "y",
                                              col = color.metadata), size = 1) +
      ggplot2::geom_polygon(data = hull.data,
                            ggplot2::aes_string(x = "x", y = "y",
                                                fill = color.metadata,
                                                group = color.metadata), alpha=0.30) +
      ggplot2::labs(title = "Multidimensional Scaling",
                    subtitle = paste0("Kruskal's stress = ", round(fit2$stress, 2))) +
      ggplot2::geom_hline(yintercept = (min.lim + max.lim) / 2, linetype = "dashed") +
      ggplot2::geom_vline(xintercept = (min.lim + max.lim) / 2, linetype = "dashed") +
      ggrepel::geom_text_repel(data = proj.samples,
                               ggplot2::aes_string(x = "x", y = "y",
                                                   label = "samples", col = color.metadata),
                               size = 3, force = 3)
  }

  plot <- plot +
    ggplot2::coord_fixed() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic"),
                   aspect.ratio = 1,
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill = NA),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank())

  return(plot)
}
