\name{GetScores}
\alias{GetScores}
\title{Compute Directional Pathway Activity Scores}
\usage{
GetScores(
  expression_profile,
  activation_geneset,
  inhibition_geneset,
  method = c("gsva", "ssgsea", "zscore", "plage"),
  min.sz = 1,
  max.sz = 10000,
  equity = F,
  weights = NULL
)
}
\arguments{
\item{expression_profile}{A matrix or data frame of expression values, with
genes as rows and samples as columns.}
\item{activation_geneset}{A named list of character vectors containing the
activation (positive-direction) gene sets.}
\item{inhibition_geneset}{A named list of character vectors containing the
inhibition (negative-direction) gene sets.}
\item{method}{Single-sample enrichment method: \code{"gsva"}, \code{"ssgsea"},
\code{"zscore"}, or \code{"plage"}.}
\item{min.sz}{Minimum gene-set size. Default \code{1}.}
\item{max.sz}{Maximum gene-set size. Default \code{10000}.}
\item{equity}{Logical. If \code{TRUE}, aggregate by the difference of mean
activation and inhibition scores; if \code{FALSE}, normalize by the total
number of gene sets. Default \code{FALSE}.}
\item{weights}{Optional named numeric vector of gene-set weights applied before
aggregation.}
}
\value{
A list with \code{final_activity_score} (a data frame with columns
\code{activation_score}, \code{inhibition_score},
\code{positive_aggregation_score}, and \code{negative_aggregation_score}) and
\code{single_score_matrix} (the standardized per-gene-set enrichment matrix).
}
\description{
Projects an expression matrix onto activation and inhibition gene sets using
single-sample enrichment methods, then combines the directional scores into
activation, inhibition, and signed aggregation scores.
}
\examples{
# See the package README for a complete simulated-data example.
}
