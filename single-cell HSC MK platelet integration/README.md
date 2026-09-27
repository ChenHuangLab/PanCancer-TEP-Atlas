# HSC, megakaryocyte and platelet single-cell integration

These are cleaned publication copies of `in-slico_FACS.R` and `HSC_Megakaryocyte_Platelet_altas_integration.R`. The original filenames are retained for traceability. The source files outside this repository were not modified.

## Scope and inputs

The scripts start from five **prepared Seurat RDS objects**. Collection from CELLxGENE/GEO, annotation to gene symbols with MyGene.info against the study's version 70 annotation, and the normal/tumor/disease-control sample definitions from Table S1 are upstream preparation steps. Those source data and mapping tables are not included here, so these two scripts do not recreate that acquisition or mapping stage. Each input must have an `RNA` assay with raw `counts`, a nonempty `sample_id` metadata column, and gene-symbol rows including all 12 FACS markers. Other sample metadata such as `sample_source`, `tissue`, `disease`, and `disease_state` are preserved if present.

Place the inputs next to the scripts:

```text
single-cell HSC MK platelet integration/
├── HSC_Megakaryocyte_Platelet_altas_integration.R
├── in-slico_FACS.R
└── input/
    ├── BloodCancer/first_qc_Megakaryocyte.rds
    ├── Malignant Lymphoma/first_qc_Megakaryocyte.rds
    ├── Multiple Myeloma/first_qc_Megakaryocyte.rds
    ├── Osteosarcoma/first_qc_Megakaryocyte.rds
    └── other Cancer/first_qc_Megakaryocyte.rds
```

The analysis requires Seurat 5, SeuratObject, Harmony, and their dependencies. `in-slico_FACS.R` also supports a `SingleCellExperiment` input when `SummarizedExperiment` is installed, though the integration entry point uses Seurat objects.

## Methods implemented

1. Retain cells with **at least 250 raw UMIs**.
2. Label a cell as MK/platelet when **any one** of `PF4`, `PPBP`, `MAST1`, `ITGA2B`, `GP9`, or `GP1BA` has a raw count **greater than 1**. Label a cell as HSC when **more than three** of `SPINK2`, `CYTL1`, `EGFL7`, `GATA1`, `GATA2`, and `CD34` have positive raw counts. A cell can meet both rules. Missing marker genes cause an explicit error rather than silently weakening a rule.
3. Drop samples with fewer than **10** retained target cells.
4. Find 2,000 variable features per sample and select the top 2,000 shared integration features. Log-normalize at a scale factor of `1e4`, scale selected features, and run 30 PCs.
5. Run Harmony on a batch ID made from the source group and `sample_id` with `theta = 6` and `lambda = 1`. This keeps identically named samples from different input groups separate. Then build neighbors from Harmony dimensions 1–30, cluster at **resolution 0.1**, and compute UMAP from the same Harmony dimensions. Resolution belongs to `FindClusters`, not `FindNeighbors`.

The script first normalizes each sample for feature ranking, and then normalizes the merged object before scaling and PCA. This is the Seurat preprocessing sequence used in the release script.

## Run

```sh
Rscript HSC_Megakaryocyte_Platelet_altas_integration.R
```

`input/` and `results/` are resolved relative to the script file. To use another location:

```sh
Rscript HSC_Megakaryocyte_Platelet_altas_integration.R --input-dir /path/to/inputs --output-dir /path/to/results
```

The script writes `HSC_MK_Platelet_integrated.rds`, a retained-cell count CSV, and UMAP PDFs for metadata fields present in the inputs. The input RDS files and generated results are ignored by Git; this directory releases the analysis scripts and their input contract.

## References for function behavior

- [Seurat: SelectIntegrationFeatures](https://satijalab.org/seurat/reference/selectintegrationfeatures)
- [Seurat: NormalizeData](https://satijalab.org/seurat/reference/normalizedata)
- [Harmony: RunHarmony for Seurat](https://github.com/immunogenomics/harmony/blob/master/vignettes/Seurat.Rmd)
