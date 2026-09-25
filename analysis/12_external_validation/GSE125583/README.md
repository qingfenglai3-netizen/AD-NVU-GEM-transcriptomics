# GSE125583

Evaluation of the fixed 18-gene logistic L2 model in bulk RNA-seq of the fusiform gyrus
(219 AD and 70 control donors), with a comparison of GEM expression. Results correspond to
S16 Fig and the GSE125583 sheets of S10 Table.

## Run

From the repository root:

```sh
python analysis/12_external_validation/GSE125583/download_input.py
python analysis/12_external_validation/GSE125583/run_validation.py
```

`download_input.py` retrieves the expression archive (about 92 MB) from GEO and checks its
SHA-256. `run_validation.py` reads the discovery features and GSE118553 predictions from
`data/`, reproduces the GSE118553 predictions, and writes to `recomputed/` in this directory.
`--input-dir` points to an existing copy of the archive; `--output-dir` sets another output
directory.

## Analysis

The analysis plan is in `analysis_spec_frozen.json`. The model was trained on GSE132903 and
applied without refitting. The primary analysis uses the discovery scaling; a secondary
analysis standardizes each gene within GSE125583 without diagnostic labels; a subset analysis
excludes the 89 donors also present in GSE95587 (158 AD, 42 controls). AUC is computed on the
continuous decision score with 2,000 diagnosis-stratified bootstrap resamples. GEM is
compared on the log2(nRPKM + 1) scale with a Welch test and a linear model adjusted for age
and sex (HC3 standard errors).

## Results

| Analysis | AD / control | AUC (95% CI) |
|---|---|---|
| Discovery scaling | 219 / 70 | 0.799 (0.742–0.850) |
| Within-cohort standardization | 219 / 70 | 0.827 (0.778–0.874) |
| Discovery scaling, GSE95587 donors excluded | 158 / 42 | 0.802 (0.730–0.867) |
| Within-cohort standardization, GSE95587 donors excluded | 158 / 42 | 0.837 (0.780–0.892) |

GEM was higher in AD by 0.480 log2(nRPKM + 1) after adjustment for age and sex (95% CI
0.287–0.673). Estimates, per-donor predictions, donor metadata and gene coverage are in
`results/`; `results/validation_summary.pdf` is the source of S16 Fig.

Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125583
