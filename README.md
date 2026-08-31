# WGCNA on Breast Cancer Transcriptomic Data (METABRIC)

Weighted Gene Co-expression Network Analysis (WGCNA) pipeline in R, applied to the METABRIC breast cancer dataset.

## Dataset

- **Source:** [Breast Cancer Gene Expression Profiles (METABRIC)](https://www.kaggle.com/datasets/raghadalharbi/breast-cancer-gene-expression-profiles-metabric) — Kaggle
- **Samples:** 1,904 breast cancer patients
- **Expression features:** 489 genes (mRNA z-scores)
- **Clinical traits used:** PAM50 molecular subtype, tumor size, tumor stage (30 clinical variables available in total)

> The raw CSV is not included in this repository due to size/licensing — download it directly from the Kaggle link above and place it as `data/METABRIC_RNA_Mutation.csv`.

## Requirements

```r
install.packages("BiocManager")
BiocManager::install("WGCNA")
install.packages(c("dplyr", "fastDummies"))
```

Tested in R on Google Colab (R runtime).

## Pipeline

| Step | Script section | Function(s) |
|---|---|---|
| 1. Load & split columns (clinical / mutation / gene expression) | `01_prepare_data.R` | base R |
| 2. Build sample × gene and sample × trait matrices | `01_prepare_data.R` | base R |
| 3. Data quality check | `02_run_wgcna.R` | `goodSamplesGenes()` |
| 4. Soft-thresholding power selection | `02_run_wgcna.R` | `pickSoftThreshold()` |
| 5. Network construction & module detection | `02_run_wgcna.R` | `blockwiseModules()` |
| 6. Module–trait relationships | `03_module_trait.R` | `moduleEigengenes()`, `cor()`, `corPvalueStudent()`, `labeledHeatmap()` |
| 7. Export results | `04_export_results.R` | `write.csv()`, `save()` |

Full annotated script: [`wgcna_metabric_analysis.R`](./wgcna_metabric_analysis.R)

## How to run

1. Download `METABRIC_RNA_Mutation.csv` from Kaggle and place it in `data/`.
2. Open `wgcna_metabric_analysis.R` in R / RStudio / Google Colab (R runtime).
3. Run top to bottom. Figures are written to `output/figures/`, tables to `output/tables/`.

## Results

- **6 gene modules** detected (plus one unassigned/grey module), at the soft-thresholding power selected via scale-free topology fit (R² ≈ 0.9 criterion).
- Module eigengenes were correlated against tumor size, tumor stage, and PAM50 subtype (one-hot encoded, as it is a nominal categorical variable).
- Notable pattern: the *blue* and *turquoise* modules show opposing correlation trends across PAM50 subtypes (e.g., Basal vs. Luminal A), consistent with known biology of these breast cancer subtypes.

See `output/figures/module_trait_relationship.png` and `output/figures/dendrogram_modul.png` for the full visualizations, and `output/tables/` for per-module gene lists and correlation/p-value tables.

## Repository structure

```
.
├── data/                          # place METABRIC_RNA_Mutation.csv here (not tracked)
├── wgcna_metabric_analysis.R       # full annotated pipeline
├── output/
│   ├── figures/
│   │   ├── dendrogram_modul.png
│   │   └── module_trait_relationship.png
│   └── tables/
│       ├── metabric_gene_modules.csv
│       ├── module_<color>_genes.csv   # one file per module
│       ├── module_trait_correlation.csv
│       └── module_trait_pvalue.csv
└── README.md
```

## Notes / lessons learned

- `patient_id` and `cohort` were excluded from trait correlation despite being numeric — they are identifiers/batch labels, not biological variables, and including them produces spurious correlations.
- `cancer_type` was excluded as a trait: 1,903 of 1,904 samples fall into a single category, leaving no meaningful variance to correlate against.
- PAM50 subtype (nominal, unordered) was one-hot encoded rather than integer/label-encoded, to avoid implying a false ordinal relationship between subtypes.
