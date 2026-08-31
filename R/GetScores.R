#' Compute directional pathway activity scores
#'
#' @description
#' Projects an expression matrix onto activation and inhibition gene sets using
#' single-sample enrichment methods, then combines the directional scores into
#' activation, inhibition, and signed aggregation scores.
#'
#' @param expression_profile A matrix or data frame of expression values, with
#'   genes as rows and samples as columns.
#' @param activation_geneset A named list of character vectors containing the
#'   activation (positive-direction) gene sets.
#' @param inhibition_geneset A named list of character vectors containing the
#'   inhibition (negative-direction) gene sets.
#' @param method Single-sample enrichment method: `"gsva"`, `"ssgsea"`,
#'   `"zscore"`, or `"plage"`.
#' @param min.sz Minimum gene-set size. Default `1`.
#' @param max.sz Maximum gene-set size. Default `10000`.
#' @param equity Logical. If `TRUE`, aggregate by the difference of mean
#'   activation and inhibition scores; if `FALSE`, normalize by the total number
#'   of gene sets. Default `FALSE`.
#' @param weights Optional named numeric vector of gene-set weights applied
#'   before aggregation.
#'
#' @return A list with `final_activity_score` (columns `activation_score`,
#'   `inhibition_score`, `positive_aggregation_score`,
#'   `negative_aggregation_score`) and `single_score_matrix` (the standardized
#'   per-gene-set enrichment matrix).
#' @export



GetScores=function(expression_profile,
                   activation_geneset,
                   inhibition_geneset,
                   method=c("gsva","ssgsea","zscore","plage"),
                   min.sz=1,
                   max.sz=10000,
                   equity=F,
                   weights=NULL
){
  options(matrixStats.useNames = FALSE)
  library(GSVA)
  
  exp = expression_profile
  dimnames = list(rownames(exp), colnames(exp))
  exp = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
  exp = as.data.frame(exp)
  
  print("Start calculating geneset score profiles.")
  geneSets = c(activation_geneset, inhibition_geneset)
  
  if(method == "gsva"){
    param = gsvaParam(exprData = as.matrix(exp), geneSets = geneSets, 
                      minSize = min.sz, maxSize = max.sz, kcdf = "Gaussian")
    result = gsva(param)
  } else if(method == "ssgsea"){
    param = ssgseaParam(exprData = as.matrix(exp), geneSets = geneSets,
                        minSize = min.sz, maxSize = max.sz)
    result = gsva(param)
  } else if(method == "zscore"){
    param = zscoreParam(exprData = as.matrix(exp), geneSets = geneSets,
                        minSize = min.sz, maxSize = max.sz)
    result = gsva(param)
  } else if(method == "plage"){
    param = plageParam(exprData = as.matrix(exp), geneSets = geneSets,
                       minSize = min.sz, maxSize = max.sz)
    result = gsva(param)
  }
  
  result = as.data.frame(t(result))
  
  # Handle missing gene sets and all-NaN columns
  diff = setdiff(c(names(activation_geneset), names(inhibition_geneset)), colnames(result))
  nan_columns = sapply(result, function(x) all(is.nan(x) | is.na(x)))
  result = result[, !nan_columns, drop = FALSE]
  
  # Standardize the enrichment scores
  result = as.data.frame(scale(result))
  
  # Apply gene-set weights (if provided)
  if(!is.null(weights)){
    weights = weights[!nan_columns]
    weights = weights[setdiff(names(weights), diff)]
    result = result[, colnames(result) %in% names(weights), drop = FALSE]
    stopifnot(ncol(result) == length(weights))
    result = sweep(result, 2, weights, "*")
  }
  
  # Define the directional activity-score helper
  calculate_activity_scores = function(df,equity) {
    pos_cols = grep("pos", names(df), ignore.case = TRUE, value = TRUE)
    neg_cols = grep("neg", names(df), ignore.case = TRUE, value = TRUE)
    if(equity){
      if (length(pos_cols) > 0) {
        df$activation_score = rowMeans(df[, pos_cols, drop = FALSE], na.rm = TRUE)
      } else {
        df$activation_score = NA
      }
      
      if (length(neg_cols) > 0) {
        df$inhibition_score = rowMeans(df[, neg_cols, drop = FALSE], na.rm = TRUE)
      } else {
        df$inhibition_score = NA
      }
      
      df$positive_aggregation_score = ifelse(!is.na(df$activation_score) & !is.na(df$inhibition_score),
                                 df$activation_score - df$inhibition_score,
                                 NA)
      df$negative_aggregation_score = ifelse(!is.na(df$inhibition_score) & !is.na(df$activation_score),
                                 df$inhibition_score - df$activation_score,
                                 NA)
      
    }else{
      if (length(pos_cols) > 0) {
        df$activation_score = rowMeans(df[, pos_cols, drop = FALSE], na.rm = TRUE)
        df$activation_score1 = rowSums(df[, pos_cols, drop = FALSE], na.rm = TRUE)
      } else {
        df$activation_score = NA
      }
      
      if (length(neg_cols) > 0) {
        df$inhibition_score = rowMeans(df[, neg_cols, drop = FALSE], na.rm = TRUE)
        df$inhibition_score1 = rowSums(df[, neg_cols, drop = FALSE], na.rm = TRUE)
      } else {
        df$inhibition_score = NA
      }
      
      df$positive_aggregation_score = ifelse(!is.na(df$activation_score1) & !is.na(df$inhibition_score1),
                                 (df$activation_score1 - df$inhibition_score1)/(length(pos_cols)+length(neg_cols)),
                                 NA)
      df$negative_aggregation_score = ifelse(!is.na(df$inhibition_score1) & !is.na(df$activation_score1),
                                 (df$inhibition_score1 - df$activation_score1)/(length(pos_cols)+length(neg_cols)),
                                 NA)
      df=df[,c("activation_score","inhibition_score","positive_aggregation_score","negative_aggregation_score")]
    }
    return(df)
  }
  
  # Compute the final directional scores
  final_result = calculate_activity_scores(df=result,equity=equity)[, c("activation_score","inhibition_score","positive_aggregation_score","negative_aggregation_score")]
  final_result = list(final_result, result)
  names(final_result) = c("final_activity_score", "single_score_matrix")
  return(final_result)
}

























