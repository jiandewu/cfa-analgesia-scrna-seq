# CFA analgesia single-cell RNA-seq analysis

R workflows for quality control, integration, clustering, marker exploration, and treatment comparisons in a six-sample 10x Genomics single-cell RNA-seq experiment.

The experimental design combines two conditions (Saline and CFA) with three treatments (Vehicle, APAP, and 3DDA).

## Repository layout

- `analysis/01_quality_control.R`: initial six-sample QC and integration workflow.
- `analysis/02_compare_cfa_vs_saline.R`: CFA versus Saline under Vehicle.
- `analysis/03_compare_saline_treatments.R`: Vehicle, APAP, and 3DDA under Saline.
- `analysis/04_compare_cfa_treatments.R`: Vehicle, APAP, and 3DDA under CFA.
- `analysis/05_all_samples.R`: combined factorial analysis.
- `analysis/06_new_dataset_all_samples.R`: later six-sample dataset using generic sample labels.
- `environment.yml`: reproducible R environment.

These are exploratory analysis scripts preserved from the original project and cleaned for public release. They use an older Seurat-style workflow and should be reviewed before applying thresholds or biological labels to another dataset.

## Environment

Create and activate the Conda environment:

```bash
conda env create -f environment.yml
conda activate cfa-analgesia-scrna-seq
```

## Input layout

By default, scripts read Cell Ranger output from `data/`. Each sample should contain:

```text
data/
└── SAMPLE_NAME/
    └── outs/
        └── filtered_feature_bc_matrix/
            ├── barcodes.tsv.gz
            ├── features.tsv.gz
            └── matrix.mtx.gz
```

The named-sample workflows expect these directories:

```text
Saline_Veh_A1
CFA_Veh_B1
Saline_APAP_C1
CFA_APAP_D1
Saline_3DDA_E1
CFA_3DDA_F1
```

The new-dataset workflow expects `Sample1` through `Sample6`.

## Configuration

Input and output roots can be changed without editing the scripts:

```bash
export SCRNA_DATA_DIR=/path/to/cellranger-runs
export SCRNA_RESULTS_DIR=/path/to/results
Rscript analysis/05_all_samples.R
```

When unset, the scripts use `data/` and `results/` relative to the current working directory.

## Public-data policy

No count matrices, sample metadata, generated results, credentials, manuscripts, or drafts are distributed in this repository. Those paths and common generated formats are ignored by Git.

Before publishing changes, inspect exactly what will be committed:

```bash
git status
git diff --cached
```
