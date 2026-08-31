#' Discover and refine directional pathway gene signatures
#'
#' @description
#' Derives de novo pathway-associated gene signatures from multiple expression
#' datasets through one-sided differential analysis, cross-dataset p-value
#' integration, knowledge-based refinement, and redundancy removal, returning a
#' compact, directionally resolved signature panel.
#'
#' @param expression_accession_vector Character vector of dataset names. Each
#'   name must exist in the calling environment as an expression data frame
#'   (genes as rows, samples as columns), together with a `<name>_G` group data
#'   frame containing `Tag` and `group` columns.
#' @param alternative One-sided test direction: `"less"` (inhibition, logFC < 0)
#'   or `"greater"` (activation, logFC > 0).
#' @param denovo_genes Number of top de novo genes retained per dataset.
#'   Default `1000`.
#' @param top_genes Number of top genes retained in the final signature
#'   (alternative to `top_threshold`).
#' @param top_threshold Critical-score threshold for the final signature
#'   (alternative to `top_genes`).
#' @param knowledge_genesets A named list of prior gene sets used for
#'   concordance-based refinement.
#' @param concordance Logical. If `TRUE`, intersect de novo gene sets with
#'   knowledge gene sets pairwise; if `FALSE`, concatenate them. Default `TRUE`.
#' @param knn Logical. If `TRUE`, impute missing values with KNN. Default
#'   `FALSE`.
#' @param p_combine_method Cross-dataset p-value integration method: one of
#'   `"geometric_mean"`, `"rra"`, `"invchisq"`, `"logitp"`, `"cct"`, `"meanp"`,
#'   `"meanz"`, `"sumlog"`, `"sumz"`, `"sump"`, `"votep"`, `"wilkinsonp"`, or
#'   `"rankproduct"`.
#' @param na_ratio Maximum allowed proportion of missing values per gene.
#'   Default `0.3`.
#' @param magnitude Logical. If `TRUE`, weight the p-value score by logFC
#'   magnitude. Default `FALSE`.
#' @param similarity_metric Gene-set similarity coefficient: `"Jaccard"`,
#'   `"Sorensen-Dice"`, `"Hub-Promoted"`, or `"Hub-Depressed"`.
#' @param deduplication Redundancy-removal strategy: `"network"` (graph
#'   components) or `"clustering"` (hierarchical clustering).
#' @param linkage Linkage method for hierarchical clustering (used when
#'   `deduplication = "clustering"`).
#' @param similarity_threshold Similarity threshold above which gene sets are
#'   merged. Default `0.5`.
#' @param min_genes Minimum number of genes required in a final gene set.
#'   Default `5`.
#'
#' @return A list with `genesets` (final directional gene sets named `PAGS_*` or
#'   `PIGS_*`), `pmatrix`, `critical_score`, and, when `magnitude = TRUE`,
#'   `fcmatrix`.
#' @export




PathwayQuant=function(expression_accession_vector,
                      alternative=c("less","greater"),
                      denovo_genes=1000,
                      top_genes=NULL,
                      top_threshold=NULL,
                      knowledge_genesets=NULL,
                      concordance=TRUE,
                      knn=F,
                      p_combine_method=c("geometric_mean","rra",
                                         "invchisq","logitp","cct",
                                         "meanp","meanz","sumlog","sumz",
                                         "sump","votep","wilkinsonp",
                                         "rankproduct"),
                      na_ratio=0.3,
                      magnitude=F,
                      similarity_metric=c("Jaccard","Sorensen-Dice","Hub-Promoted","Hub-Depressed"),
                      deduplication=c("network","clustering"),
                      linkage=c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid"),
                      similarity_threshold=0.5,
                      min_genes=5
){
  library(limma)
  library(metap)
  library(reshape2)
  library(igraph)
  library(tidyverse)
  library(RobustRankAggreg)
  library(RankProd)
  library(impute)
  ##########################################CTT: a p-value integration method
  #Liu Y, Xie J. Cauchy combination test: a powerful test with analytic p-value calculation under arbitrary dependency structures. J Am Stat Assoc. 2020;115(529):393-402. doi: 10.1080/01621459.2018.1554485. Epub 2019 Apr 25. PMID: 33012899; PMCID: PMC7531765.
  #Liu Y, Chen S, Li Z, Morrison AC, Boerwinkle E, Lin X. ACAT: A Fast and Powerful p Value Combination Method for Rare-Variant Analysis in Sequencing Studies. Am J Hum Genet. 2019 Mar 7;104(3):410-421. doi: 10.1016/j.ajhg.2019.01.002. PMID: 30849328; PMCID: PMC6407498.
  CCT = function(pvals, weights=NULL){
    if(sum(is.na(pvals)) > 0){
      stop("Cannot have NAs in the p-values!")
    }
    if((sum(pvals<0) + sum(pvals>1)) > 0){
      stop("All p-values must be between 0 and 1!")
    }
    is.zero = (sum(pvals==0)>=1)
    is.one = (sum(pvals==1)>=1)
    if(is.zero && is.one){
      zero_indices = which(pvals == 0)
      one_indices = which(pvals == 1)
      remove_count = min(length(zero_indices), length(one_indices))
      remove_zero = zero_indices[1:remove_count]
      remove_one = one_indices[1:remove_count]
      remove_indices = c(remove_zero, remove_one)
      pvals = pvals[-remove_indices]
      if (length(pvals) == 0){
        return(0)
      }
      if (!is.null(weights)) {
        weights = weights[-remove_indices]
      }
      warning("Removed ", remove_count, " p-value(s) of 0 and ", remove_count, " p-value(s) of 1 to avoid conflict.")
    }
    is.zero = (sum(pvals==0)>=1)
    is.one = (sum(pvals==1)>=1)
    if(is.zero){
      return(0)
    }
    if(is.one){
      warning("There are p-values that are exactly 1!")
      return(1)
    }
    if(is.null(weights)){
      weights = rep(1/length(pvals),length(pvals))
    }else if(length(weights)!=length(pvals)){
      stop("The length of weights should be the same as that of the p-values!")
    }else if(sum(weights < 0) > 0){
      stop("All the weights must be positive!")
    }else{
      weights = weights/sum(weights)
    }
    is.small = (pvals < 1e-16)
    if (sum(is.small) == 0){
      cct.stat = sum(weights*tan((0.5-pvals)*pi))
    }else{
      cct.stat = sum((weights[is.small]/pvals[is.small])/pi)
      cct.stat = cct.stat + sum(weights[!is.small]*tan((0.5-pvals[!is.small])*pi))
    }
    if(cct.stat > 1e+15){
      pval = (1/cct.stat)/pi
    }else{
      pval = 1-pcauchy(cct.stat)
    }
    return(pval)
  }
  ##########################################pvalue_integration_strategy: Parallel integration of multiple P-value integration methods.
  pvalue_integration_strategy=function(matrix,p_combine_method){
    if(length(colnames(matrix))>1){
      if(p_combine_method %in% c("invchisq")){
        p_combine_list=apply(matrix,1,function(row) invchisq(p=na.omit(row), k=length(na.omit(row))))
        p_score=c()
        for(pnum in 1:length(p_combine_list)){
          p_combine_list_single=p_combine_list[[pnum]]$p
          names(p_combine_list_single)=names(p_combine_list)[pnum]
          p_score=c(p_score,p_combine_list_single)
        }
        p_score=as.data.frame(p_score)
      }else if(p_combine_method=="geometric_mean"){
        p_score=apply(matrix,1,function(row) exp(mean(log(row),na.rm=T)))
        p_score=as.data.frame(p_score)
      }else if(p_combine_method=="cct"){
        p_score=apply(matrix,1,function(row) CCT(p=na.omit(row)))
        p_score=as.data.frame(p_score)
      }else if(p_combine_method=="logitp"){
        p_score=apply(matrix,1,function(row) as.numeric(logitp(na.omit(row))$p))
        p_score=as.data.frame(p_score)
      }else if(p_combine_method=="meanp"){
        p_score=apply(matrix,1,function(row) as.numeric(meanp(na.omit(row))$p))
        p_score=as.data.frame(p_score)
      }else if(p_combine_method=="meanz"){
        p_score=apply(matrix,1,function(row) as.numeric(meanz(na.omit(row))$p))
        p_score=as.data.frame(p_score)
      }else if(p_combine_method=="sumlog"){
        p_score=apply(matrix,1,function(row) as.numeric(sumlog(na.omit(row))$p))
        p_score=as.data.frame(p_score)
      }else if(p_combine_method=="sumz"){
        p_score=apply(matrix,1,function(row) as.numeric(sumz(na.omit(row))$p))
        p_score=as.data.frame(p_score)
      }else if(p_combine_method=="sump"){
        p_score=apply(matrix,1,function(row) as.numeric(sump(na.omit(row))$p))
        p_score=as.data.frame(p_score)
      }else if(p_combine_method=="votep"){
        p_score=apply(matrix,1,function(row) as.numeric(votep(na.omit(row))$p))
        p_score=as.data.frame(p_score)
      }else if(p_combine_method=="wilkinsonp"){
        p_score=apply(matrix,1,function(row) as.numeric(wilkinsonp(na.omit(row))$p))
        p_score=as.data.frame(p_score)
      }
    }else{
      p_score=matrix
    }
    p_score[is.na(p_score)]=1
    return(p_score)
  }
  ##########################################Limma one-sided test
  #Artmann S, Jung K, Bleckmann A, Beissbarth T. Detection of simultaneous group effects in microRNA expression and related target gene sets. PLoS One. 2012;7(6):e38365. doi: 10.1371/journal.pone.0038365. Epub 2012 Jun 19. PMID: 22723856; PMCID: PMC3378551.
  limma.one.sided = function(fit, lower = FALSE){
    se.coef=sqrt(fit$s2.post) * fit$stdev.unscaled
    df.total=fit$df.prior + fit$df.residual
    pt(fit$t, df = df.total, lower.tail = lower)
  }
  ##########################################
  vector=expression_accession_vector
  list_data=list()
  list_G_data=list()
  print("Start loading expression data and grouping information.")
  for (q in 1:length(vector)) {
    list_data[[q]]=get(vector[q])
    list_G_data[[q]]=get(paste0(vector[q],"_G"))
    names(list_data)[q]=paste0(vector[q])
    names(list_G_data)[q]=paste0(vector[q],"_G")
  }
  modified_vector=lapply(vector, function(x) {
    ifelse(grepl("_",x),substr(x,1,regexpr("_", x)-1),x)
  })
  modified_vector=unique(as.character(modified_vector))
  large_list=c()
  for (item in 1:length(modified_vector)) {
    assign(paste0("list_", modified_vector[item]),list_data[grep(paste0("^",modified_vector[item]), names(list_data))])
    large_list=c(large_list,paste0("list_", modified_vector[item]))
  }
  print("Start the differential analysis of genes.")
  for (j in 1:length(large_list)) {
    expList=get(large_list[j])
    for (i in 1:length(expList)){
      exp=expList[[i]]
      dimnames=list(rownames(exp),colnames(exp))
      exp=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
      exp=as.data.frame(exp)
      group=list_G_data[[paste0(names(expList)[i],"_G")]]
      group=as.data.frame(group)
      Get_limma=function(gene=gene,group=group,alternative_limma=alternative){
        design=model.matrix(~0+factor(group$group))
        colnames(design)=levels(factor(group$group))
        rownames(design)=colnames(gene)
        compare=makeContrasts(H-L, levels=design)
        fit=lmFit(gene, design)
        fit=contrasts.fit(fit, compare)
        fit=eBayes(fit)
        if(alternative_limma=="greater"){
          result=limma.one.sided(fit, lower=FALSE)
          result=as.data.frame(result)
          colnames(result)="P.Value"
          result$gene_id=rownames(result)
          topTable=topTable(fit, coef=1, number=Inf)
          topTable$gene_id=rownames(topTable)
          topTable=topTable[,c("gene_id","logFC")]
          result=merge(result,topTable,by="gene_id")
        }else if(alternative_limma=="less"){
          result=limma.one.sided(fit, lower=TRUE)
          result=as.data.frame(result)
          colnames(result)="P.Value"
          result$gene_id=rownames(result)
          topTable=topTable(fit, coef=1, number=Inf)
          topTable$gene_id=rownames(topTable)
          topTable=topTable[,c("gene_id","logFC")]
          result=merge(result,topTable,by="gene_id")
        }else{
          print("Please enter a correct method.(greater/less)")
        }
        return(result)
      }
      result=Get_limma(gene=exp,group=group,alternative_limma=alternative)
      print("Start screening genes.")
      result=result[,c("gene_id","logFC","P.Value")]
      assign(paste0(names(expList)[i],"_result"),result)
    }
    if (length(expList)>1){
      for (mp in 1:length(expList)){
        if (mp == 1) { 
          result_single=get(paste0(names(expList)[mp],"_result"))
          result_p=result_single[,c("gene_id","P.Value")]
          colnames(result_p)[2]=c("P.Value1")
          result_logFC=result_single[,c("gene_id","logFC")]
          colnames(result_logFC)[2]=c("logFC1")
          result_p_all=result_p
          result_logFC_all=result_logFC
        } else {
          result_single=get(paste0(names(expList)[mp],"_result"))
          result_p=result_single[,c("gene_id","P.Value")]
          colnames(result_p)[2]=paste0("P.Value",mp)
          result_logFC=result_single[,c("gene_id","logFC")]
          colnames(result_logFC)[2]=paste0("logFC",mp)
          result_p_all=merge(result_p_all,result_p,by="gene_id")
          result_logFC_all=merge(result_logFC_all,result_logFC,by="gene_id")
        }
      }
      rownames(result_p_all)=result_p_all[,1]
      result_p_all=result_p_all[,-1]
      rownames(result_logFC_all)=result_logFC_all[,1]
      result_logFC_all=result_logFC_all[,-1]
      if(p_combine_method%in%c("geometric_mean",
                               "invchisq","logitp","cct",
                               "meanp","meanz","sumlog","sumz",
                               "sump","votep","wilkinsonp")){
        result_p_all_combine=pvalue_integration_strategy(matrix=result_p_all,p_combine_method=p_combine_method)
      }
      if(p_combine_method=="rra"){
        rank_list = lapply(result_p_all, function(x) {
          if(is.numeric(x)) {
            ranks = rank(x, na.last = "keep", ties.method = "min")
            rownames(result_p_all)[order(ranks, na.last = NA)]
          }
        })
        aggregated_ranks = aggregateRanks(rank_list)
        result_p_all_combine = aggregated_ranks[,c("Score"),drop=F]
      }
      if(p_combine_method=="rankproduct"){
        data_matrix = as.matrix(result_p_all)
        cl = rep(1, ncol(data_matrix))
        rp_result = RP(data_matrix, cl = cl, num.perm = 100)
        result_p_all_combine = rp_result$pval[, 1, drop = FALSE]
        rownames(result_p_all_combine) = rownames(result_p_all)
      }
      result_p_all_combine=as.data.frame(result_p_all_combine)
      colnames(result_p_all_combine)="p_score"
      result_p_all_combine$p_score[result_p_all_combine$p_score == 0]=min(result_p_all_combine$p_score[result_p_all_combine$p_score!=0])
      result_p_all=result_p_all_combine
      result_p_all$rowname=rownames(result_p_all)
      result_p_all=result_p_all[,c("rowname","p_score")]
      colnames(result_p_all)=c("gene_id","P.Value")
      result_logFC_all$logFC_mean=rowMeans(result_logFC_all)
      result_logFC_all$gene_id=rownames(result_logFC_all)
      result_logFC_all=result_logFC_all[,c("gene_id","logFC_mean")]
      colnames(result_logFC_all)=c("gene_id","logFC")
      result=merge(result_p_all,result_logFC_all,by="gene_id")
      assign(paste0(gsub("_.*$","",names(expList)[1]),"_result"),result)
    }else{
      print("Genesets are not processed for merging.")
    }
  }
  print("Start screening genes.")
  for (mq in 1:length(modified_vector)) {
    result=get(paste0(modified_vector[mq],"_result"))
    if(alternative=="less"){
      downregulated_gene=subset(result,logFC<0)
      downregulated_gene=downregulated_gene[order(downregulated_gene$P.Value), ]
      downregulated_gene=downregulated_gene[1:denovo_genes, ]
      downregulated_gene=downregulated_gene$gene_id
      assign(paste0(modified_vector[mq],"_DN"),downregulated_gene)
    }else if(alternative=="greater"){
      upregulated_gene=subset(result,logFC>0)
      upregulated_gene=upregulated_gene[order(upregulated_gene$P.Value), ]
      upregulated_gene=upregulated_gene[1:denovo_genes, ]
      upregulated_gene=upregulated_gene$gene_id
      assign(paste0(modified_vector[mq],"_UP"),upregulated_gene)
    }else{
      print("Please enter a correct method.(greater/less)")
    }
  }
  if(alternative=="less"){
    combined_DN=list()
    for (mn in 1:length(modified_vector)) {
      combined_DN[[mn]]=get(paste0(modified_vector[mn],"_DN"))
      names(combined_DN)[mn]=paste0(modified_vector[mn], "_DN")
    }
    padded_list = list()
    print("Start pairwise intersection operations.")
    if(concordance){
      for (name1 in names(combined_DN)) {
        for (name2 in names(knowledge_genesets)) {
          intersect_set = intersect(combined_DN[[name1]], knowledge_genesets[[name2]])
          padded_list[[paste(name1, name2, "intersect", sep="_")]] = intersect_set
        }
      }
      padded_list = padded_list[sapply(padded_list, length) > 0]
      padded_list = padded_list[!duplicated(lapply(padded_list, sort))]
    }else{
      if(exists("knowledge_genesets")){
        padded_list = c(combined_DN, knowledge_genesets)
      }else{
        padded_list = combined_DN
      }
    }
  }else if(alternative=="greater"){
    combined_UP=list()
    for (mn in 1:length(modified_vector)) {
      combined_UP[[mn]]=get(paste0(modified_vector[mn],"_UP"))
      names(combined_UP)[mn]=paste0(modified_vector[mn], "_UP")
    }
    padded_list = list()
    print("Start pairwise intersection operations.")
    if(concordance){
      for (name1 in names(combined_UP)) {
        for (name2 in names(knowledge_genesets)) {
          intersect_set = intersect(combined_UP[[name1]], knowledge_genesets[[name2]])
          padded_list[[paste(name1, name2, "intersect", sep="_")]] = intersect_set
        }
      }
      padded_list = padded_list[sapply(padded_list, length) > 0]
      padded_list = padded_list[!duplicated(lapply(padded_list, sort))]
    }else{
      if(exists("knowledge_genesets")){
        padded_list = c(combined_UP,knowledge_genesets)
      }else{
        padded_list = combined_UP
      }
    }
  }else{
    print("Please enter a correct method.(greater/less)")
  }
  names(padded_list)=paste0("PGS",1:length(padded_list))
  for (i in 1:length(padded_list)) {
    padded_list[[i]]=unique(padded_list[[i]])
  }
  print("Check for duplicates within each gene set.")
  print("Start building the p-value matrix.")
  for (mn in 1:length(modified_vector)) {
    if (mn == 1){
      pmatrix=get(paste0(modified_vector[mn],"_result"))[,c("gene_id","P.Value")]
      colnames(pmatrix)[mn+1]=paste0(modified_vector[mn])
    }else{
      pmatrix=merge(pmatrix,get(paste0(modified_vector[mn],"_result"))[,c("gene_id","P.Value")],by="gene_id",all=T)
      colnames(pmatrix)[mn+1]=paste0(modified_vector[mn])
    }
  }
  rownames(pmatrix)=pmatrix[,1]
  pmatrix=pmatrix[,-1]
  pmatrix$na_count=rowSums(is.na(pmatrix))/ncol(pmatrix)
  pmatrix=pmatrix[pmatrix$na_count <= na_ratio,]
  pmatrix=pmatrix[,!colnames(pmatrix) %in% "na_count"]
  if(knn == TRUE){
    data_matrix = as.matrix(pmatrix)
    imputed_result = impute.knn(data_matrix, k = 10, rowmax = 1, colmax = 1)
    pmatrix = as.data.frame(imputed_result$data)
  }
  if(p_combine_method%in%c("geometric_mean",
                           "invchisq","logitp","cct",
                           "meanp","meanz","sumlog","sumz",
                           "sump","votep","wilkinsonp")){
    pscores=pvalue_integration_strategy(matrix=pmatrix,p_combine_method=p_combine_method)
    pmatrix=merge(pmatrix,pscores,by="row.names")
    rownames(pmatrix)=pmatrix[,1]
    pmatrix=pmatrix[,-1]
    pmatrix$p_score=-log10(pmatrix$p_score)
    pmatrix=pmatrix[order(pmatrix$p_score,decreasing = T),]
  }
  if(p_combine_method == "rra"){
    rank_list = lapply(pmatrix, function(x) {
      if(is.numeric(x)) {
        ranks = rank(x, na.last = "keep", ties.method = "min")
        rownames(pmatrix)[order(ranks, na.last = NA)]
      }
    })
    aggregated_ranks = aggregateRanks(rank_list)
    pmatrix$p_score = aggregated_ranks$Score[match(rownames(pmatrix), aggregated_ranks$Name)]
    pmatrix$p_score=-log10(pmatrix$p_score)
    pmatrix=pmatrix[order(pmatrix$p_score,decreasing = T),]
  }
  if(p_combine_method == "rankproduct"){
    data_matrix = as.matrix(pmatrix)
    cl = rep(1, ncol(data_matrix))
    rp_result = RP(data_matrix, cl = cl, num.perm = 1000)
    pmatrix$p_score = rp_result$pval[, 1]
    pmatrix$p_score=-log10(pmatrix$p_score)
    pmatrix=pmatrix[order(pmatrix$p_score,decreasing = T),]
  }
  ##############################################################################
  print("Start building the logFC matrix.")
  for (mn in 1:length(modified_vector)) {
    if (mn == 1){
      logFCmatrix=get(paste0(modified_vector[mn],"_result"))[,c("gene_id","logFC")]
      colnames(logFCmatrix)[mn+1]=paste0(modified_vector[mn])
    }else{
      logFCmatrix=merge(logFCmatrix,get(paste0(modified_vector[mn],"_result"))[,c("gene_id","logFC")],by="gene_id",all=T)
      colnames(logFCmatrix)[mn+1]=paste0(modified_vector[mn])
    }
  }
  rownames(logFCmatrix)=logFCmatrix[,1]
  logFCmatrix=logFCmatrix[,-1]
  logFCmatrix$na_count=rowSums(is.na(logFCmatrix))/ncol(logFCmatrix)
  logFCmatrix=logFCmatrix[logFCmatrix$na_count <= na_ratio,]
  logFCmatrix=logFCmatrix[,!colnames(logFCmatrix) %in% "na_count"]
  if(knn == TRUE){
    data_matrix = as.matrix(logFCmatrix)
    imputed_result = impute.knn(data_matrix, k = 10, rowmax = 1, colmax = 1)
    logFCmatrix = as.data.frame(imputed_result$data)
  }
  logFCmatrix$logFC_score=rowMeans(logFCmatrix,na.rm = T)
  logFCmatrix=logFCmatrix[order(logFCmatrix$logFC_score,decreasing = T),]
  ##############################################################################
  critical_score=pmatrix[,ncol(pmatrix),drop=F]
  if(magnitude){
    logFCmatrix0=logFCmatrix[,ncol(logFCmatrix),drop=F]
    critical_score=merge(critical_score,logFCmatrix0,by="row.names")
    rownames(critical_score)=critical_score[,1]
    critical_score=critical_score[,-1]
    critical_score$p_score = as.numeric(as.character(critical_score$p_score))
    critical_score$logFC_score = as.numeric(as.character(critical_score$logFC_score))
    if(alternative=="greater"){
      critical_score$critical_score=critical_score$logFC_score*critical_score$p_score
    }
    if(alternative=="less"){
      critical_score$critical_score=(-critical_score$logFC_score)*critical_score$p_score
    }
    critical_score=critical_score[order(critical_score$critical_score,decreasing = T),]
  }else{
    colnames(critical_score)="critical_score"
  }
  ##############################################################################
  if(alternative=="less"){
    if(!is.null(top_genes)){
      pmatrix_DN=critical_score[1:top_genes,,drop=F]
    }
    if(!is.null(top_threshold)){
      pmatrix_DN=critical_score[critical_score$critical_score>top_threshold,,drop=F]
    }
    Element_DN=as.character(rownames(pmatrix_DN))
    update_inhibition_geneset=list()
    inhibition_geneset=padded_list
    for (mi in 1:length(names(inhibition_geneset))){
      coldn=inhibition_geneset[[mi]]
      coldn=coldn[coldn != ""]
      col_intersect_dn=intersect(coldn,Element_DN)
      update_inhibition_geneset[[mi]]=col_intersect_dn
      names(update_inhibition_geneset)[mi]=paste0(names(inhibition_geneset)[mi],"_purification")
    }
    aim=update_inhibition_geneset
  }else if(alternative=="greater"){
    if(!is.null(top_genes)){
      pmatrix_UP=critical_score[1:top_genes,,drop=F]
    }
    if(!is.null(top_threshold)){
      pmatrix_UP=critical_score[critical_score$critical_score>top_threshold,,drop=F]
    }
    Element_UP=as.character(rownames(pmatrix_UP))
    update_activation_geneset=list()
    activation_geneset=padded_list
    for (ma in 1:length(names(activation_geneset))){
      colup=activation_geneset[[ma]]
      colup=colup[colup != ""]
      col_intersect_up=intersect(colup,Element_UP)
      update_activation_geneset[[ma]]=col_intersect_up
      names(update_activation_geneset)[ma]=paste0(names(activation_geneset)[ma],"_purification")
    }
    aim=update_activation_geneset
  }else{
    print("Please enter a correct method.(greater/less)")
  }
  GENE_SETS=aim[
    sapply(aim, function(x) {
      sum(x == "" | is.na(x) | x == " ") < length(x)
    })
  ]
  GENE_SETS = GENE_SETS[!duplicated(lapply(GENE_SETS, sort))]
  update_GENE_SETS=list()
  for (i in 1:length(names(GENE_SETS))){
    colup=GENE_SETS[[i]]
    colup=colup[colup != ""]
    update_GENE_SETS[[i]]=colup
    names(update_GENE_SETS)[i]=paste0(names(GENE_SETS)[i])
  }
  print("Start calculating correlations between genesets using integration coefficients.")
  if(similarity_metric=="Jaccard"){
    integration_ratio=function(set1, set2) {
      intersection=length(intersect(set1, set2))
      union=length(union(set1, set2))
      ratio=intersection/union
      return(ratio)
    }
  }else if(similarity_metric=="Sorensen-Dice"){
    integration_ratio=function(set1, set2) {
      intersection=length(intersect(set1, set2))
      Nset1=length(set1)
      Nset2=length(set2)
      ratio=(intersection*2)/(Nset1+Nset2)
      return(ratio)
    }
  }else if(similarity_metric=="Hub-Promoted"){
    integration_ratio=function(set1, set2) {
      intersection=length(intersect(set1, set2))
      Nset1=length(set1)
      Nset2=length(set2)
      ratio=intersection/min(Nset1,Nset2)
      return(ratio)
    }
  }else if(similarity_metric=="Hub-Depressed"){
    integration_ratio=function(set1, set2) {
      intersection=length(intersect(set1, set2))
      Nset1=length(set1)
      Nset2=length(set2)
      ratio=intersection/max(Nset1,Nset2)
      return(ratio)
    }
  }
  if(deduplication=="network"){
    print("Remove redundancy in the network.")
    matrix=outer(
      names(update_GENE_SETS), names(update_GENE_SETS), 
      Vectorize(function(x, y) integration_ratio(update_GENE_SETS[[x]], update_GENE_SETS[[y]]))
    )
    dimnames(matrix)=list(names(update_GENE_SETS), names(update_GENE_SETS))
    matrix[upper.tri(matrix)]=NA
    matrix=reshape2::melt(matrix)
    colnames(matrix)=c("rownames", "colnames", "value")
    matrix=matrix[matrix$rownames != matrix$colnames, ]
    matrix=na.omit(matrix)
    matrix0=matrix[matrix$value>=similarity_threshold,]
    if (is.null(matrix0)) {
      print("Low correlation between genesets, no need for joint genesets.")
    } else {
      g=graph_from_data_frame(matrix0, directed = FALSE)
      clusters=components(g)
      membership=as.data.frame(clusters$membership)
      membership$genesets=rownames(membership)
      colnames(membership)=c("cluster","genesets")
      num_cluster=unique(membership$cluster)
      update_GENE_SETS_disjunction=update_GENE_SETS[!names(update_GENE_SETS)%in%membership$genesets]
      for(z in 1:length(num_cluster)){
        membershipud=membership[membership$cluster==num_cluster[z],,drop=F]
        union_vector=c()
        for (name in membershipud$genesets) {
          union_vector=union(union_vector, update_GENE_SETS[[name]])
        }
        update_GENE_SETS_disjunction[[length(update_GENE_SETS_disjunction)+1]]=unique(union_vector)
        names(update_GENE_SETS_disjunction)[length(update_GENE_SETS_disjunction)]=paste0("Integration_",num_cluster[z],"_",alternative)
      }
    } 
  }
  if(deduplication=="clustering"){
    print("Remove redundancy by clustering.")
    aux.merge = function(gsets, minjac, linkage_method){
      integration_matrix = outer(1:length(gsets), 1:length(gsets), 
                                 Vectorize(function(x, y) integration_ratio(gsets[[x]], 
                                                                            gsets[[y]])))
      dimnames(integration_matrix) = list(names(gsets), names(gsets))
      integration_matrix[is.na(integration_matrix)] = 0
      integration_matrix[is.infinite(integration_matrix)] = 0
      integration_matrix = pmax(pmin(integration_matrix, 1), 0)
      dist_matrix = as.dist(1 - integration_matrix)
      hc = hclust(dist_matrix, method = linkage_method)
      if (is.unsorted(hc$height)) {
        ord = order(hc$height)
        hc$height = hc$height[ord]
        hc$merge = hc$merge[ord, ]
      }
      clust = cutree(hc, h = 1 - minjac)
      tapply(1:length(gsets), clust, function(i){sort(unique(unlist(gsets[i])))})	
    }
    update_GENE_SETS_disjunction = aux.merge(gsets = update_GENE_SETS, minjac = similarity_threshold, linkage_method=linkage)
    update_GENE_SETS_disjunction=update_GENE_SETS_disjunction[!duplicated(lapply(update_GENE_SETS_disjunction, sort))]
  }
  update_GENE_SETS_disjunction=update_GENE_SETS_disjunction[sapply(update_GENE_SETS_disjunction,length) >= min_genes]
  if(alternative=="greater"){
    names(update_GENE_SETS_disjunction)=paste0("PAGS", 1:length(update_GENE_SETS_disjunction), "_pos")
  }
  if(alternative=="less"){
    names(update_GENE_SETS_disjunction)=paste0("PIGS", 1:length(update_GENE_SETS_disjunction), "_neg")
  }
  if(magnitude){
    output=list(update_GENE_SETS_disjunction,pmatrix,logFCmatrix,critical_score)
    names(output)=c("genesets","pmatrix","fcmatrix","critical_score")
  }else{
    output=list(update_GENE_SETS_disjunction,pmatrix,critical_score)
    names(output)=c("genesets","pmatrix","critical_score")
  }
  return(output)
}