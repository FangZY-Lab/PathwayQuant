# PathwayQuant

A unified, direction-resolved framework for pathway activity quantification and
de novo gene-signature discovery from multi-cohort expression data.

## Highlights

- **Direction-aware pathway quantification.** Pathway activation and inhibition
  are modeled separately and fused into a single signed activity score,
  preserving the regulatory direction that undirected enrichment methods
  discard.
- **Cohort-agnostic single-sample scoring.** `GetScores` supports GSVA, ssGSEA,
  z-score, and PLAGE for sample-level, batch-robust pathway activity
  estimation.
- **Multi-cohort de novo discovery.** `PathwayQuant` couples one-sided
  differential analysis with thirteen cross-dataset p-value integration
  strategies—including robust rank aggregation and the Cauchy combination
  test—to synthesize evidence across heterogeneous cohorts.
- **Knowledge-guided refinement.** Prior gene sets are fused with de novo
  signatures through pairwise concordance filtering, bridging data-driven and
  hypothesis-driven discovery.
- **Redundancy-aware consolidation.** Jaccard, Sørensen-Dice, and hub-based
  similarity, combined with graph-component detection or hierarchical
  clustering, resolve biologically redundant gene sets.
- **Translation-ready signatures.** The resulting `PAGS` (activation) and
  `PIGS` (inhibition) panels are compact, interpretable, and directly usable
  for downstream scoring, validation, and biomarker translation.

## Pipeline overview

```mermaid
flowchart TB

    %% ================= INPUT =================
    subgraph IN["Input"]
        IN_EXP["Expression matrices<br/>(genes x samples)"]
        IN_GROUP["Group files<br/>(Tag, group: H / L)"]
        IN_POS["Activation gene sets<br/>(positive direction)"]
        IN_NEG["Inhibition gene sets<br/>(negative direction)"]
        IN_PRIOR["Prior knowledge gene sets"]
    end

    %% ================= GetScores =================
    subgraph GS["GetScores: directional activity scoring"]
        direction TB
        GS_METHOD["Single-sample enrichment"]
        GS_gsva["gsva"]
        GS_ssgsea["ssgsea"]
        GS_zscore["zscore"]
        GS_plage["plage"]
        GS_SCORES["Gene-set enrichment scores"]
        GS_NAN["Drop all-NaN columns"]
        GS_SCALE["Standardize (scale)"]
        GS_WT{"weights?"}
        GS_WT_YES["Apply gene-set weights"]
        GS_ACT["activation_score<br/>(positive gene sets)"]
        GS_INH["inhibition_score<br/>(negative gene sets)"]
        GS_AGG["positive / negative<br/>aggregation scores"]
        GS_OUT["final_activity_score +<br/>single_score_matrix"]
    end

    %% ================= PathwayQuant =================
    subgraph PQ["PathwayQuant: de novo signature discovery"]
        direction TB
        PQ_LOAD["Load data via get()"]
        PQ_GROUP["Group datasets by prefix (_)"]
        PQ_LIMMA["One-sided limma test<br/>(H vs L)"]
        PQ_FC["Per-gene logFC + P-value"]
        PQ_META{"Multiple<br/>sub-datasets?"}
        PQ_COMB["Cross-dataset<br/>p-value integration"]
        PQ_c1["geometric_mean"]
        PQ_c2["rra"]
        PQ_c3["invchisq"]
        PQ_c4["logitp"]
        PQ_c5["cct"]
        PQ_c6["meanp / meanz"]
        PQ_c7["sumlog / sumz"]
        PQ_c8["sump / votep"]
        PQ_c9["wilkinsonp"]
        PQ_c10["rankproduct"]
        PQ_SCREEN["Screen genes<br/>(denovo_genes)"]
        PQ_DIR{"alternative?"}
        PQ_UP["Activation<br/>(logFC > 0)"]
        PQ_DN["Inhibition<br/>(logFC < 0)"]
        PQ_KNOW["Knowledge integration"]
        PQ_CONC{"concordance?"}
        PQ_CONC_YES["Pairwise intersection"]
        PQ_CONC_NO["Concatenation"]
        PQ_PMAT["P-value matrix"]
        PQ_FCMAT["logFC matrix"]
        PQ_NAF["NA filter (na_ratio)"]
        PQ_KNN{"knn?"}
        PQ_KNN_YES["KNN imputation"]
        PQ_PSCORE["p_score = -log10(P)"]
        PQ_MAG{"magnitude?"}
        PQ_MAG_YES["critical_score =<br/>logFC x p_score"]
        PQ_MAG_NO["critical_score = p_score"]
        PQ_TOP["top_genes / top_threshold"]
        PQ_PURIFY["Purification (intersection)"]
        PQ_DEDUP{"deduplication?"}
        PQ_NET["network<br/>(graph components)"]
        PQ_CLUST["clustering<br/>(hclust + cutree)"]
        PQ_SIM["Similarity metric"]
        PQ_s1["Jaccard"]
        PQ_s2["Sorensen-Dice"]
        PQ_s3["Hub-Promoted"]
        PQ_s4["Hub-Depressed"]
        PQ_MIN["min_genes filter"]
        PQ_OUT["genesets (PAGS / PIGS)<br/>pmatrix + critical_score"]
    end

    %% ================= OUTPUT =================
    subgraph OUT["Output"]
        OUT_FINAL["Directional gene signatures"]
    end

    %% ================= EDGES: GetScores =================
    IN_EXP --> GS_METHOD
    IN_POS --> GS_METHOD
    IN_NEG --> GS_METHOD
    GS_METHOD --> GS_gsva
    GS_METHOD --> GS_ssgsea
    GS_METHOD --> GS_zscore
    GS_METHOD --> GS_plage
    GS_gsva --> GS_SCORES
    GS_ssgsea --> GS_SCORES
    GS_zscore --> GS_SCORES
    GS_plage --> GS_SCORES
    GS_SCORES --> GS_NAN --> GS_SCALE --> GS_WT
    GS_WT -->|yes| GS_WT_YES
    GS_WT -->|no| GS_ACT
    GS_WT_YES --> GS_ACT
    GS_ACT --> GS_AGG
    GS_INH --> GS_AGG
    GS_AGG --> GS_OUT
    GS_OUT --> OUT_FINAL

    %% ================= EDGES: PathwayQuant =================
    IN_EXP --> PQ_LOAD
    IN_GROUP --> PQ_LOAD
    PQ_LOAD --> PQ_GROUP --> PQ_LIMMA --> PQ_FC --> PQ_META
    PQ_META -->|yes| PQ_COMB
    PQ_META -->|no| PQ_SCREEN
    PQ_COMB --> PQ_c1
    PQ_COMB --> PQ_c2
    PQ_COMB --> PQ_c3
    PQ_COMB --> PQ_c4
    PQ_COMB --> PQ_c5
    PQ_COMB --> PQ_c6
    PQ_COMB --> PQ_c7
    PQ_COMB --> PQ_c8
    PQ_COMB --> PQ_c9
    PQ_COMB --> PQ_c10
    PQ_c1 --> PQ_SCREEN
    PQ_c2 --> PQ_SCREEN
    PQ_c3 --> PQ_SCREEN
    PQ_c4 --> PQ_SCREEN
    PQ_c5 --> PQ_SCREEN
    PQ_c6 --> PQ_SCREEN
    PQ_c7 --> PQ_SCREEN
    PQ_c8 --> PQ_SCREEN
    PQ_c9 --> PQ_SCREEN
    PQ_c10 --> PQ_SCREEN
    PQ_SCREEN --> PQ_DIR
    PQ_DIR -->|greater| PQ_UP
    PQ_DIR -->|less| PQ_DN
    PQ_UP --> PQ_KNOW
    PQ_DN --> PQ_KNOW
    IN_PRIOR --> PQ_KNOW
    PQ_KNOW --> PQ_CONC
    PQ_CONC -->|TRUE| PQ_CONC_YES
    PQ_CONC -->|FALSE| PQ_CONC_NO
    PQ_CONC_YES --> PQ_PMAT
    PQ_CONC_NO --> PQ_PMAT
    PQ_PMAT --> PQ_NAF
    PQ_FCMAT --> PQ_NAF
    PQ_NAF --> PQ_KNN
    PQ_KNN -->|TRUE| PQ_KNN_YES
    PQ_KNN -->|FALSE| PQ_PSCORE
    PQ_KNN_YES --> PQ_PSCORE
    PQ_PSCORE --> PQ_MAG
    PQ_MAG -->|TRUE| PQ_MAG_YES
    PQ_MAG -->|FALSE| PQ_MAG_NO
    PQ_MAG_YES --> PQ_TOP
    PQ_MAG_NO --> PQ_TOP
    PQ_TOP --> PQ_PURIFY --> PQ_DEDUP
    PQ_DEDUP -->|network| PQ_NET
    PQ_DEDUP -->|clustering| PQ_CLUST
    PQ_NET --> PQ_SIM
    PQ_CLUST --> PQ_SIM
    PQ_SIM --> PQ_s1
    PQ_SIM --> PQ_s2
    PQ_SIM --> PQ_s3
    PQ_SIM --> PQ_s4
    PQ_s1 --> PQ_MIN
    PQ_s2 --> PQ_MIN
    PQ_s3 --> PQ_MIN
    PQ_s4 --> PQ_MIN
    PQ_MIN --> PQ_OUT
    PQ_OUT --> OUT_FINAL
```

## Installation

Run the following in a fresh R session:

```r
install.packages("pak")
pak::pkg_install("FangZY-Lab/PathwayQuant")
```

`pak` installs all dependencies automatically. Alternatively, install them
explicitly:

```r
install.packages("BiocManager")
BiocManager::install(c("GSVA", "limma", "RankProd", "impute"))
install.packages(c("metap", "reshape2", "igraph", "RobustRankAggreg", "tidyverse"))

install.packages("remotes")
remotes::install_github("FangZY-Lab/PathwayQuant", upgrade = FALSE)
```

Then load the package:

```r
library(PathwayQuant)
```

## Quick start

The example below simulates ten expression datasets grouped into five studies
(`SimDataA`, `SimDataB`, `SimDataC`, `SimDataD`, `SimDataE`), where A, C, and D
each comprise several sub-datasets. `GetScores` computes sample-level
directional activity scores, and `PathwayQuant` discovers activation and
inhibition gene signatures.

```r
library(PathwayQuant)
set.seed(123)

make_dataset <- function(n_genes = 100, n_h = 30, n_l = 30, seed = 1) {
  set.seed(seed)
  n_samples <- n_h + n_l
  samples <- c(paste0("S", seq_len(n_h)), paste0("S", n_h + seq_len(n_l)))

  mat <- matrix(rnorm(n_samples * n_genes), nrow = n_genes, ncol = n_samples)
  # 10 activation genes: higher expression in the high-pathway (H) group.
  mat[1:10, seq_len(n_h)] <- mat[1:10, seq_len(n_h)] + 2
  # 10 inhibition genes: higher expression in the low-pathway (L) group.
  mat[11:20, n_h + seq_len(n_l)] <- mat[11:20, n_h + seq_len(n_l)] + 2

  expr <- as.data.frame(mat)
  rownames(expr) <- paste0("GENE", seq_len(n_genes))
  colnames(expr) <- samples

  group <- data.frame(
    Tag = samples,
    group = c(rep("H", n_h), rep("L", n_l)),
    stringsAsFactors = FALSE
  )
  list(expr = expr, group = group)
}

# Ten datasets across five studies. A, C, and D are split into sub-datasets.
dataset_seed <- c(
  SimDataA_a = 101, SimDataA_b = 102, SimDataA_c = 103,
  SimDataB   = 104,
  SimDataC_a = 105, SimDataC_b = 106, SimDataC_c = 107,
  SimDataD_a = 108, SimDataD_b = 109,
  SimDataE   = 110
)

for (nm in names(dataset_seed)) {
  d <- make_dataset(seed = dataset_seed[[nm]])
  assign(nm, d$expr)
  assign(paste0(nm, "_G"), d$group)
}

datasets <- names(dataset_seed)

# Twenty knowledge gene sets, each containing ten genes (drawn from GENE1-GENE100).
knowledge_genesets <- setNames(
  lapply(1:20, function(i) paste0("GENE", ((0:9) + (i - 1) * 5) %% 100 + 1)),
  paste0("pathway_", sprintf("%02d", 1:20))
)

# 1. Sample-level directional pathway activity scores.
scores <- GetScores(
  expression_profile = SimDataA_a,
  activation_geneset = list(
    pos_activation_A = paste0("GENE", 1:5),
    pos_activation_B = paste0("GENE", 6:10)
  ),
  inhibition_geneset = list(
    neg_inhibition_A = paste0("GENE", 11:15),
    neg_inhibition_B = paste0("GENE", 16:20)
  ),
  method = "ssgsea"
)
head(scores$final_activity_score)

# 2. Discover activation signatures (up-regulated genes).
activation <- PathwayQuant(
  expression_accession_vector = datasets,
  alternative = "greater",
  denovo_genes = 50,
  top_genes = 20,
  knowledge_genesets = knowledge_genesets,
  concordance = TRUE,
  knn = FALSE,
  p_combine_method = "geometric_mean",
  magnitude = FALSE,
  similarity_metric = "Jaccard",
  deduplication = "network",
  linkage = "ward.D",
  similarity_threshold = 0.5,
  min_genes = 5
)
str(activation$genesets)

# 3. Discover inhibition signatures (down-regulated genes).
inhibition <- PathwayQuant(
  expression_accession_vector = datasets,
  alternative = "less",
  denovo_genes = 50,
  top_genes = 20,
  knowledge_genesets = knowledge_genesets,
  concordance = TRUE,
  knn = FALSE,
  p_combine_method = "geometric_mean",
  magnitude = FALSE,
  similarity_metric = "Jaccard",
  deduplication = "network",
  linkage = "ward.D",
  similarity_threshold = 0.5,
  min_genes = 5
)
str(inhibition$genesets)
```

If the demo works, the activation signatures contain the planted activation
genes (`GENE1`-`GENE10`) and the inhibition signatures contain the planted
inhibition genes (`GENE11`-`GENE20`).

## Function reference

### `GetScores()`

Sample-level directional pathway activity scoring via single-sample enrichment.

```r
GetScores(
  expression_profile,
  activation_geneset,
  inhibition_geneset,
  method = "ssgsea",
  min.sz = 1,
  max.sz = 10000,
  equity = FALSE,
  weights = NULL
)
```

| Argument | Description | Default |
| --- | --- | --- |
| `expression_profile` | Matrix/data frame of expression values, genes as rows and samples as columns. | required |
| `activation_geneset` | Named list of activation (positive-direction) gene sets; names should contain `"pos"`. | required |
| `inhibition_geneset` | Named list of inhibition (negative-direction) gene sets; names should contain `"neg"`. | required |
| `method` | Enrichment method: `"gsva"`, `"ssgsea"`, `"zscore"`, or `"plage"`. | `"ssgsea"` |
| `min.sz` | Minimum gene-set size. | `1` |
| `max.sz` | Maximum gene-set size. | `10000` |
| `equity` | If `TRUE`, aggregate by the mean-score difference; otherwise sum-normalized. | `FALSE` |
| `weights` | Optional named numeric vector of gene-set weights. | `NULL` |

Returns `list(final_activity_score, single_score_matrix)`.

### `PathwayQuant()`

Multi-cohort de novo gene-signature discovery with knowledge-guided refinement
and redundancy resolution.

```r
PathwayQuant(
  expression_accession_vector,
  alternative = "greater",
  denovo_genes = 1000,
  top_genes = NULL,
  top_threshold = NULL,
  knowledge_genesets = NULL,
  concordance = TRUE,
  knn = FALSE,
  p_combine_method = "geometric_mean",
  na_ratio = 0.3,
  magnitude = FALSE,
  similarity_metric = "Jaccard",
  deduplication = "network",
  linkage = "ward.D",
  similarity_threshold = 0.5,
  min_genes = 5
)
```

| Argument | Description | Default |
| --- | --- | --- |
| `expression_accession_vector` | Dataset names; each must exist in the calling environment together with `<name>_G` (columns `Tag`, `group` with `"H"`/`"L"`). | required |
| `alternative` | Test direction: `"greater"` (activation) or `"less"` (inhibition). | `"greater"` |
| `denovo_genes` | Number of top de novo genes retained per dataset. | `1000` |
| `top_genes` | Number of top genes retained in the final signature. | `NULL` |
| `top_threshold` | Critical-score threshold for the final signature. | `NULL` |
| `knowledge_genesets` | Named list of prior gene sets. | `NULL` |
| `concordance` | If `TRUE`, intersect de novo and knowledge sets pairwise. | `TRUE` |
| `knn` | If `TRUE`, KNN-impute missing values. | `FALSE` |
| `p_combine_method` | Cross-dataset p-value integration: `"geometric_mean"`, `"rra"`, `"invchisq"`, `"logitp"`, `"cct"`, `"meanp"`, `"meanz"`, `"sumlog"`, `"sumz"`, `"sump"`, `"votep"`, `"wilkinsonp"`, `"rankproduct"`. | `"geometric_mean"` |
| `na_ratio` | Maximum missing-value proportion per gene. | `0.3` |
| `magnitude` | If `TRUE`, weight p-score by logFC magnitude. | `FALSE` |
| `similarity_metric` | `"Jaccard"`, `"Sorensen-Dice"`, `"Hub-Promoted"`, `"Hub-Depressed"`. | `"Jaccard"` |
| `deduplication` | `"network"` or `"clustering"`. | `"network"` |
| `linkage` | Hierarchical clustering linkage (when `deduplication = "clustering"`). | `"ward.D"` |
| `similarity_threshold` | Similarity threshold for merging gene sets. | `0.5` |
| `min_genes` | Minimum genes per final gene set. | `5` |

Returns `list(genesets, pmatrix, critical_score)` and, when
`magnitude = TRUE`, `fcmatrix`.

## Notes

1. `PathwayQuant()` reads expression and group data from the calling
   environment via `get()`. Create the objects in the same environment where
   you call the function.
2. Specify the method-selection arguments explicitly (for example `method`,
   `alternative`, `p_combine_method`, `similarity_metric`, `deduplication`,
   `linkage`), as shown in the examples.

## Author

Dingkang Zhao (赵定康)

Email: <dingkang.25@intl.zju.edu.cn>
