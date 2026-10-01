# Biological characterization of PAS-associated regions

These scripts use a modified [AttriMIL](https://github.com/MedCAI/AttriMIL) model to locate WSI patches associated with transcriptomically defined PAS groups. The analysis is exploratory and uses the same cohort from which the model checkpoint was obtained. The saved class probabilities are **in-sample scores**, not evidence of PAS prediction in an independent cohort or a clinically deployable model. HoVer-Net and HistoTME outputs are used to describe the selected regions, not to validate a PAS classifier.

## Files and workflow

| File | Role |
| --- | --- |
| `trainer_attrimil_nsclc_all_sample.py` | Load `full_model.pt`, rank class-specific AttriMIL attribute scores, export the top 100 level-0 coordinates per slide, and optionally make WSI heatmaps. Despite its historical filename, this script does not train a model. |
| `models/AttriMIL.py`, `dataloader.py` | Modified upstream model and H5 feature reader required by the mining script. |
| `extract_all_patches_stratified.py` | Select slides within each PAS group using in-sample scores and extract their top five patches for downstream characterization. |
| `stat_hovernet_result.py`, `morphology_calculation.py` | Summarize externally generated HoVer-Net PanNuke nuclei JSON files. |
| `stat_histotme_result.py` | Match selected coordinates to externally generated HistoTME signatures, one cohort at a time. |
| `stat_histotme_aligned.py` | Repeat HistoTME summary on the slides that also have HoVer-Net output. |
| `full_model.pt` | Saved model weights supplied with this release. SHA-256: `48fe950922aff0305fb6c76e8d5043d924f9bde32937970892dd88d94f45c034`. |

The exact training script that produced `full_model.pt` was not available in the source material used to prepare this folder. The code here reproduces the **region-mining and biological characterization steps from that checkpoint**, not the initial training. No independent validation cohort is included.

## Inputs

The metadata CSV needs `slide_id`, `case_id`, and `group` (`Low`, `Medium`, `High`). A `Patho`, `source`, or `project` column identifying LUAD or LUSC is recommended. If absent, `dataloader.py` locates each slide in the two feature roots and rejects ambiguous matches. Each feature root must contain `features_uni_v1/<slide_id>.h5` with aligned `features` (1024 dimensions) and level-0 `coords` arrays. Heatmaps and patch extraction additionally need the corresponding `.svs` files. Patch size is 256 pixels in the original analysis.

Install the Python packages in `requirements.txt` and the OpenSlide native library appropriate for your platform. Supply the external feature files, slides, HoVer-Net JSON, and HistoTME CSV separately; they are not part of this code release.

## Example commands

Run commands from this folder, replacing every placeholder path:

```bash
python trainer_attrimil_nsclc_all_sample.py \
  --metadata-csv /data/NSCLC_wsi_with_group.csv \
  --luad-features /data/TCGA-LUAD/20x_256px_0px_overlap \
  --lusc-features /data/TCGA-LUSC/20x_256px_0px_overlap \
  --output-dir ./results \
  --luad-slides /data/TCGA-LUAD/slides \
  --lusc-slides /data/TCGA-LUSC/slides

python extract_all_patches_stratified.py \
  --mining-csv ./results/mining_results_top100_final.csv \
  --luad-slides /data/TCGA-LUAD/slides \
  --lusc-slides /data/TCGA-LUSC/slides \
  --output-dir ./selected_patches

python stat_hovernet_result.py \
  --hovernet-dir /data/hovernet_pannuke_output \
  --output-csv ./hovernet_cell_stats.csv

python morphology_calculation.py \
  --hovernet-dir /data/hovernet_pannuke_output \
  --output-csv ./morphology_per_slide.csv \
  --summary-csv ./morphology_group_summary.csv --mpp 0.5

python stat_histotme_result.py \
  --mining-csv ./results/mining_results_top100_final.csv \
  --histotme-dir /data/TCGA-LUAD/HistoTME_result \
  --matched-csv ./histotme_luad_matched.csv \
  --summary-csv ./histotme_luad_summary.csv

python stat_histotme_aligned.py \
  --hovernet-dir /data/hovernet_pannuke_output \
  --histotme-dir /data/TCGA-LUSC/HistoTME_result \
  --mining-csv ./results/mining_results_top100_final.csv \
  --matched-csv ./histotme_aligned_matched.csv \
  --summary-csv ./histotme_aligned_summary.csv
```

For `stat_histotme_result.py`, run separately for LUAD and LUSC with their corresponding HistoTME roots. For `stat_histotme_aligned.py`, likewise use the appropriate cohort root; its HoVer-Net directory must contain `High`, `Medium`, and `Low` subdirectories. These scripts expect per-slide HistoTME tables at `<root>/<slide_id>/uni/<slide_id>_5fold.csv` with `x` and `y` matching the mining coordinates.

## Provenance and license

`models/AttriMIL.py` and `dataloader.py` are adaptations of [MedCAI/AttriMIL](https://github.com/MedCAI/AttriMIL) (Cai et al., *Medical Image Analysis*, 2025), released under Apache-2.0. The adaptations add the three PAS groups, dropout settings, UNI feature loading, and cohort-path resolution. The mining and downstream characterization scripts are study-specific additions. The original Apache-2.0 license is included as `ATTRIMIL_LICENSE`; it applies to the adapted upstream files alongside the repository's other licensing information. HoVer-Net and HistoTME are external tools and are not redistributed here.
