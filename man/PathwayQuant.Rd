\name{PathwayQuant}
\alias{PathwayQuant}
\title{Discover and Refine Directional Pathway Gene Signatures}
\usage{
PathwayQuant(
  expression_accession_vector,
  alternative = c("less", "greater"),
  denovo_genes = 1000,
  top_genes = NULL,
  top_threshold = NULL,
  knowledge_genesets = NULL,
  concordance = TRUE,
  knn = F,
  p_combine_method = c("geometric_mean", "rra", "invchisq", "logitp", "cct",
    "meanp", "meanz", "sumlog", "sumz", "sump", "votep", "wilkinsonp",
    "rankproduct"),
  na_ratio = 0.3,
  magnitude = F,
  similarity_metric = c("Jaccard", "Sorensen-Dice", "Hub-Promoted",
    "Hub-Depressed"),
  deduplication = c("network", "clustering"),
  linkage = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty",
    "median", "centroid"),
  similarity_threshold = 0.5,
  min_genes = 5
)
}
\arguments{
\item{expression_accession_vector}{Character vector of dataset names. Each name
must exist in the calling environment as an expression data frame (genes as
rows, samples as columns), together with a \code{<name>_G} group data frame
containing \code{Tag} and \code{group} columns.}
\item{alternative}{One-sided test direction: \code{"less"} (inhibition,
logFC < 0) or \code{"greater"} (activation, logFC > 0).}
\item{denovo_genes}{Number of top de novo genes retained per dataset. Default
\code{1000}.}
\item{top_genes}{Number of top genes retained in the final signature
(alternative to \code{top_threshold}).}
\item{top_threshold}{Critical-score threshold for the final signature
(alternative to \code{top_genes}).}
\item{knowledge_genesets}{A named list of prior gene sets used for
concordance-based refinement.}
\item{concordance}{Logical. If \code{TRUE}, intersect de novo gene sets with
knowledge gene sets pairwise; if \code{FALSE}, concatenate them. Default
\code{TRUE}.}
\item{knn}{Logical. If \code{TRUE}, impute missing values with KNN. Default
\code{FALSE}.}
\item{p_combine_method}{Cross-dataset p-value integration method.}
\item{na_ratio}{Maximum allowed proportion of missing values per gene. Default
\code{0.3}.}
\item{magnitude}{Logical. If \code{TRUE}, weight the p-value score by logFC
magnitude. Default \code{FALSE}.}
\item{similarity_metric}{Gene-set similarity coefficient: \code{"Jaccard"},
\code{"Sorensen-Dice"}, \code{"Hub-Promoted"}, or \code{"Hub-Depressed"}.}
\item{deduplication}{Redundancy-removal strategy: \code{"network"} or
\code{"clustering"}.}
\item{linkage}{Linkage method for hierarchical clustering (used when
\code{deduplication = "clustering"}).}
\item{similarity_threshold}{Similarity threshold above which gene sets are
merged. Default \code{0.5}.}
\item{min_genes}{Minimum number of genes required in a final gene set. Default
\code{5}.}
}
\value{
A list with \code{genesets} (final directional gene sets named \code{PAGS_*} or
\code{PIGS_*}), \code{pmatrix}, \code{critical_score}, and, when
\code{magnitude = TRUE}, \code{fcmatrix}.
}
\description{
Derives de novo pathway-associated gene signatures from multiple expression
datasets through one-sided differential analysis, cross-dataset p-value
integration, knowledge-based refinement, and redundancy removal, returning a
compact, directionally resolved signature panel.
}
\examples{
# See the package README for a complete simulated-data example.
}
