# HSC, megakaryocyte and platelet single-cell integration

`HSC_Megakaryocyte_Platelet_altas_integration.R` is a cleaned publication copy of `2026_07_05_single_cell_atlas_construction.R`. The previously published entry-point filename is retained so existing links keep working. The original local script was not changed.

## Input layout

The integration script starts from six **prepared** Seurat 5 objects. Each needs an RNA assay with raw counts and a `sample_id` metadata column. Place the inputs under `input/` next to the scripts:

```text
single-cell HSC MK platelet integration/
├── HSC_Megakaryocyte_Platelet_altas_integration.R
├── in-slico_FACS.R
└── input/
    ├── BloodCancer/first_qc_cleaned.rds
    ├── Malignant Lymphoma/first_qc_cleaned.rds
    ├── Multiple Myeloma/first_qc_cleaned.rds
    ├── Osteosarcoma/first_qc_cleaned.rds
    ├── other Cancer/first_qc_cleaned.rds
    ├── other disease/first_qc_cleaned.rds
    ├── mapped_genes.csv
    ├── unmapped_genes.csv              # optional diagnostic input
    └── cluster_celltype_map.csv        # optional reviewed annotations
```

`mapped_genes.csv` is the upstream MyGene.info table, with columns `original_id,entrezgene,symbol` (the original headerless export also works). The script retains IDs present in this table, maps Ensembl IDs to gene symbols, and sums duplicate symbols in a sparse matrix. `unmapped_genes.csv` is an optional one-column list used to report the fraction of raw UMIs from unmapped genes. The upstream CELLxGENE/GEO download, Table S1 sample metadata, and creation of these mapping files are outside this script.

`in-slico_FACS.R` is the separately released upstream cell-selection script, implementing the manuscript's MK/platelet and HSC marker rules. The new integration entry point reads `first_qc_cleaned.rds` from the six groups and does **not** call FACS itself. It checks the 250-UMI threshold and removes samples with fewer than 10 retained cells.

## Integration and optional annotations

For each sample, the script finds 2,000 variable genes and selects 2,000 shared integration features. It then log-normalizes at `scale.factor = 1e4`, scales selected genes, computes 30 PCs, runs Harmony across source group and sample ID with `theta = 6` and `lambda = 1`, builds neighbors from Harmony components 1–30, clusters at `resolution = 0.1`, and computes UMAP. Source group is included in the batch key so matching sample IDs in different inputs remain distinct.

The source script's manual cell-type assignments contain conflicting labels for clusters 18 and 19 and refer to clusters outside the initialized range. No cell-type labels are assigned by default. To add reviewed labels, supply `input/cluster_celltype_map.csv` with the columns below and one row for **every cluster observed in that run**:

```csv
ClusterID,celltype
0,Megakaryocyte
1,Hematopoietic stem cell
```

These two rows are **format examples only**, not validated labels. Without this file, `seurat_clusters` remains available and `celltype` is omitted.

## Run

Install Seurat 5, SeuratObject, Matrix, Harmony and their R dependencies:

```sh
Rscript HSC_Megakaryocyte_Platelet_altas_integration.R
```

Or supply another input/output location:

```sh
Rscript HSC_Megakaryocyte_Platelet_altas_integration.R --input-dir /path/to/input --output-dir /path/to/results
```

`--gene-map` and `--cluster-map` override the default mapping paths. The script writes `first_qc_final_harmony.rds`, `integration_features.csv`, `retained_cells_per_sample.csv`, and UMAP PDFs under `results/`. The Zenodo object named `first_qc_final_harmony_with_nuclear_score_and_MK_platelet_score.rds` includes additional downstream scores, so it is not the direct output of this script. Input RDS objects and results are ignored by Git.

This publication copy was syntax-checked with R 4.6.0. A full integration run requires the six upstream RDS inputs and Seurat/Harmony packages, which were unavailable in the preparation environment.
