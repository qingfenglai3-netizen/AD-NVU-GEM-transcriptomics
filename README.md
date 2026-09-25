# AD-NVU-GEM-transcriptomics

Code and derived data for *GEM as a glia-linked vascular response candidate in Alzheimer's
disease: Multi-layer transcriptomic integration of the human temporal cortex* (PLOS ONE).

## Layout

| Path | Contents |
|---|---|
| `analysis/` | Analysis scripts, numbered by workflow step |
| `figures/` | Scripts that draw the main and supporting figures |
| `reproduce/` | Scripts that recompute the statistical results from `data/` |
| `data/` | Derived tables; `data/figure_source/` holds the plotted values and `MANIFEST.csv` |
| `config/paths.example.yml` | Path template for the analysis scripts |

## Workflow

| Step | Scripts | Output |
|---|---|---|
| 1 Atlas | `01_integrate_snRNA_atlas.R` | QC, integration and clustering of GSE188545 and GSE237718 |
| 2 Annotation | `02_annotate_major_cell_types.R` | Major cell types |
| 3 Bulk | `03_bulk_differential_expression.R`, `03b`, `03c` | GSE132903 differential expression; gene-level module features |
| 4 Subtypes | `04_subtype_annotation.R`, `04a_subtype_worker.R` | Vascular, astrocyte and microglial subtypes |
| 5 Abundance | `05a`–`05f` | MiloR differential abundance; donor APOE genotypes; covariate and APOE models of vascular GEM |
| 6 Communication | `06_liana_pathway_analysis.R`, `06a`, `06b` | LIANA ligand-receptor changes and filtering counts |
| 7 Trajectory | `07_sctour_state_continuity.py`, `07b` | scTour pseudotime and donor-level statistics |
| 8 Spatial | `08_spatial_rctd.R`, `08b` | RCTD deconvolution of GSE220442; spatial scores of prioritized axes |
| 9 Spatial communication | `09_holonet_spatial_communication.py`, `09a` | HoloNet scores |
| 10 Target nomination and perturbation | `10a`–`10c` | scTenifoldKnk GEM knockdown; NicheNet ligand activity and target priority; selection sensitivity |
| 11 Classification | `11a`–`11c` | 72 selector-classifier pipelines with 10-fold cross-validation |
| 12 External cohorts | `12a`–`12j`, `GSE125583/` | SEA-AD, GSE5281, GSE36980, GSE118553 and GSE125583 |

Steps 1–10 require the raw public data and large single-cell objects. These scripts are the
versions that were run; machine-specific paths appear as placeholders (`<ANALYSIS_ROOT>`,
`<WORKDIR>`, `<RAW_DATA>`, `/path/to/...`, R library paths) and must be set before use.
Their outputs are provided in `data/`.

## Reproducing the statistical results

Python 3.10 or later with the packages in `requirements.txt`; R 4.4 with edgeR and limma
(`R_sessionInfo.txt`). From the repository root:

```
python reproduce/discovery_72_pipelines.py
python reproduce/GSE36980_models.py
python reproduce/GSE118553_models.py
python reproduce/seaad_apoe.py
Rscript reproduce/milor_and_vascular_de.R
Rscript reproduce/discovery_apoe.R
python figures/figureS14_S15.py
```

Results are written to `results/`. GradientBoosting pipelines can differ slightly between operating systems (up to 0.011 AUC between Windows and Linux); all other pipelines reproduce the published values exactly. The R scripts read `GEM_REPO` (repository root) and
`GEM_OUTPUT` (output directory) if set.

The GSE125583 evaluation has its own instructions in
`analysis/12_external_validation/GSE125583/README.md`.

## Classifier settings

Input: the 18 module genes measurable on GSE132903, each the mean of all probes assigned to
the gene. Within each training fold, genes are standardized (`StandardScaler`) and, except
for the all-gene setting, the top k = 3, 5, 8, 12 or 16 genes are kept by `SelectKBest(f_classif)`.
Cross-validation is stratified 10-fold with seed 20260608; AUC is computed on the pooled
out-of-fold probabilities. Classifiers (scikit-learn 1.5.2, `random_state=20260608`):

| Model | Setting |
|---|---|
| Logistic_L2 / Logistic_L1 | `LogisticRegression`, liblinear, C = 1, max_iter = 2000 |
| ElasticNet_logistic | `LogisticRegression`, saga, l1_ratio = 0.5, C = 1, max_iter = 5000 |
| Linear_SVM / RBF_SVM | `SVC(probability=True)`, C = 1; RBF gamma = "scale" |
| RandomForest / ExtraTrees | 180 trees, class_weight = "balanced" |
| GradientBoosting | 100 stages, learning rate 0.1, depth 3 |
| DecisionTree | Gini, class_weight = "balanced" |
| KNN_3 / KNN_5 | k = 3 or 5 |
| GaussianNB | default |

The external model is Logistic_L2 on all 18 genes, trained once on GSE132903 and applied
without refitting.

## Data sources

| Accession | Data |
|---|---|
| GSE188545, GSE237718 | Temporal cortex single-nucleus RNA-seq |
| GSE132903 | Middle temporal gyrus microarray |
| GSE220442 | Visium spatial transcriptomics |
| GSE118553, GSE36980 | Temporal cortex microarray |
| GSE125583 | Fusiform gyrus RNA-seq |
| GSE5281 | Laser-captured neurons, microarray |
| SEA-AD | Vascular subsets via CELLxGENE; donor metadata in `data/SEAAD_official_donor_metadata.xlsx` |

Download locations and checksums are in `data/download_sources.json`.

## License

MIT (`LICENSE`). Public data remain under their original terms.
