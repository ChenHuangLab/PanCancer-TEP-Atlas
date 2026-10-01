# HSC, megakaryocyte and platelet single-cell integration

`HSC_Megakaryocyte_Platelet_altas_integration.R` is the author's original `2026_07_05_single_cell_atlas_construction.R`, uploaded without changes to its contents. The repository filename is retained to preserve existing links. The earlier adapted integration script has been replaced.

Source-file SHA-256: `e353eba6c91d8dfc6d48141f4e42ad0b89b24f980845bab71e9e9a47461e9c85`.

## Inputs and paths

The original script uses the author's analysis-environment paths. Set `base_dir` and the output paths for your own environment before running it. The previous command-line options (`--input-dir`, `--output-dir`, `--gene-map`, `--cluster-map`) are not part of this original script.

The input directory contains six prepared Seurat objects:

```text
MK_dataset_supplementary/
├── BloodCancer/first_qc_cleaned.rds
├── Malignant Lymphoma/first_qc_cleaned.rds
├── Multiple Myeloma/first_qc_cleaned.rds
├── Osteosarcoma/first_qc_cleaned.rds
├── other Cancer/first_qc_cleaned.rds
└── other disease/first_qc_cleaned.rds
```

The working directory also needs `mapped_genes.csv` (three headerless columns: original gene ID, Entrez ID and symbol) and `unmapped_genes.csv` (one headerless column). The input objects contain raw RNA counts and the sample/dataset metadata used in the script, including `sample_id`, `dataset_id`, `sample_source`, `disease`, `disease_state`, `tissue` and `marker_condition`.

`in-slico_FACS.R` is the separately released upstream cell-selection script. The integration script starts from the prepared `first_qc_cleaned.rds` objects and does not call FACS itself.

## Original analysis

The source merges the prepared objects, retains mapped genes, converts Ensembl IDs to symbols, aggregates duplicate symbols with Matrix.utils and excludes dotted gene names. Its effective cell filtering retains samples with at least 10 cells. It splits by `sample_id`, finds 2,000 HVGs per sample, selects integration features, normalizes, scales and runs PCA using the source defaults. Harmony uses `sample_id`, `theta = 6` and `lambda = 1`. Neighbors and the final UMAP use Harmony components 1–30; clustering uses `resolution = 0.1` and final UMAP uses `min.dist = 0.2`.

The original manual cluster-label block is preserved as supplied. It includes overlapping assignments for clusters 18 and 19, cluster IDs outside its initial 0–24 table, and trailing empty arguments in some `c()` calls; it needs author review before use. The previously described optional CSV annotation mechanism is not used by this original file. The source also overwrites its dataset-filtered object with a sample-filtered object and uses `save()` for one intermediate file with an `.rds` suffix; that intermediate is an RData file, not a `readRDS()` input.

## Requirements and verification

Install the packages loaded or called by the original script, including Seurat 5, Matrix, Matrix.utils, data.table, R.matlab, tidyverse, biomaRt, stringr, future, future.apply and harmony. Use the original analysis package versions and inputs when reproducing the atlas.

The uploaded file was checked against the supplied local source using SHA-256 and parsed with R 4.6.0. It has not been validated by a complete integration run; the upstream inputs and required packages were unavailable in the preparation environment. Syntax parsing does not check the runtime behavior of the original manual annotation block.

The deposited processed platelet object with pathway scores is described in the [repository data availability statement](../README.md#data-availability). The larger Harmony-integrated object is not included in that deposit.
