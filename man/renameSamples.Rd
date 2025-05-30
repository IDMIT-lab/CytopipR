% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4_Samples.R
\name{renameSamples}
\alias{renameSamples}
\title{Rename Samples in CYTdata Object}
\usage{
renameSamples(CYTdata, merge = FALSE, from = NULL, to)
}
\arguments{
\item{CYTdata}{A CYTdata object that will be updated with the renamed samples.}

\item{merge}{A logical value indicating whether to merge samples with duplicate names after renaming. Default is `FALSE`.}

\item{from}{A character vector of current sample names to be renamed.}

\item{to}{A character vector of new sample names corresponding to the names in the `from` argument.}
}
\value{
A modified CYTdata object with the samples renamed and metadata updated.
}
\description{
This function allows users to rename samples in a CYTdata object. The function can handle duplicate sample names, merge samples if necessary, and update the associated metadata.
}
\details{
This function renames samples in the `cellSample` dataframe of the CYTdata object. It also updates the associated sample metadata in the `sampleMetadata` slot.
If there are duplicate sample names after renaming, the function can either stop execution or merge the samples based on the `merge` argument.
If `merge` is `TRUE`, the samples with duplicate names are merged, and their metadata is handled accordingly:
- If the name of the merged samples existed previously, the metadata of the existing sample is retained.
- If the name of the merged sample was new, the metadata from the first renamed sample is retained.

The function also checks that the `from` and `to` arguments are the same length. If they are not, an error is raised.
}
\examples{
# Example 1: Rename samples without merging
updated_CYTdata <- renameSamples(CYTdata = cyt_data, from = c("Sample1", "Sample2"), to = c("NewSample1", "NewSample2"))

# Example 2: Rename samples and merge duplicate samples
updated_CYTdata <- renameSamples(CYTdata = cyt_data, from = c("Sample1", "Sample1"), to = c("MergedSample1", "MergedSample1"), merge = TRUE)

}
\seealso{
\code{\link{checkValidity}} for validating the CYTdata object before and after renaming samples.
}
