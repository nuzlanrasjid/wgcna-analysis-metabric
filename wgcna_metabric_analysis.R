# ============================================================
# WGCNA on Breast Cancer Transcriptomic Data (METABRIC)
# ============================================================
# Weighted Gene Co-expression Network Analysis, applied to the
# METABRIC breast cancer dataset (Kaggle) as a training exercise
# ahead of applying the same pipeline to plant RNA-seq data.
#
# Dataset: https://www.kaggle.com/datasets/raghadalharbi/breast-cancer-gene-expression-profiles-metabric
# ============================================================

# ------------------------------------------------------------
# 1. Install & load packages
# ------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("WGCNA", update = FALSE, ask = FALSE)
install.packages("fastDummies")

library(WGCNA)
library(dplyr)
library(fastDummies)

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# ------------------------------------------------------------
# 2. Load dataset
# ------------------------------------------------------------
df <- read.csv("data/METABRIC_RNA_Mutation.csv")

dim(df)
colSums(is.na(df))

# ------------------------------------------------------------
# 3. Split columns: clinical vs. mutation vs. gene expression
# ------------------------------------------------------------
# mutation columns are suffixed "_mut"
cols_mutation <- names(df)[grepl("_mut$", names(df), ignore.case = TRUE)]

# clinical columns (schema specific to this METABRIC release)
cols_clinical <- c(
  "patient_id", "age_at_diagnosis", "type_of_breast_surgery",
  "cancer_type", "cancer_type_detailed", "cellularity", "chemotherapy",
  "pam50_._claudin.low_subtype", "cohort", "er_status_measured_by_ihc",
  "er_status", "neoplasm_histologic_grade", "her2_status_measured_by_snp6",
  "her2_status", "tumor_other_histologic_subtype", "hormone_therapy",
  "inferred_menopausal_state", "integrative_cluster",
  "primary_tumor_laterality", "lymph_nodes_examined_positive",
  "mutation_count", "nottingham_prognostic_index", "oncotree_code",
  "overall_survival_months", "overall_survival", "pr_status",
  "radio_therapy", "X3.gene_classifier_subtype", "tumor_size",
  "tumor_stage", "death_from_cancer"
)
cols_clinical <- intersect(cols_clinical, names(df))

# remaining numeric columns = gene expression
is_numeric_col <- sapply(df, is.numeric)
cols_gene <- setdiff(names(df)[is_numeric_col], c(cols_clinical, cols_mutation))

cat("Clinical columns:", length(cols_clinical), "\n")
cat("Mutation columns:", length(cols_mutation), "(unused)\n")
cat("Gene columns    :", length(cols_gene), "\n")

# ------------------------------------------------------------
# 4. Build expression (samples x genes) and trait (samples x clinical) matrices
# ------------------------------------------------------------
id_sample <- if ("patient_id" %in% names(df)) "patient_id" else names(df)[1]

df_expr <- df[, cols_gene]
rownames(df_expr) <- df[[id_sample]]

df_clinical <- df[, cols_clinical]
rownames(df_clinical) <- df[[id_sample]]
df_clinical <- df_clinical[rownames(df_expr), , drop = FALSE]  # keep sample order aligned

dim(df_expr)
dim(df_clinical)

# ------------------------------------------------------------
# 5. Data quality check (required before network construction)
# ------------------------------------------------------------
gsg <- goodSamplesGenes(df_expr, verbose = 3)
gsg$allOK

if (!gsg$allOK) {
  df_expr <- df_expr[gsg$goodSamples, gsg$goodGenes]
  df_clinical <- df_clinical[gsg$goodSamples, , drop = FALSE]
  cat("Flagged genes/samples removed. New dimensions:\n")
  print(dim(df_expr))
}

# ------------------------------------------------------------
# 6. Soft-thresholding power selection
# ------------------------------------------------------------
power_range <- 1:20
sft <- pickSoftThreshold(df_expr, powerVector = power_range, verbose = 5)

par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft threshold (power)", ylab = "Scale-free topology fit (R^2)",
     main = "Scale independence")
abline(h = 0.90, col = "red")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft threshold (power)", ylab = "Mean connectivity",
     main = "Mean connectivity")

soft_power <- sft$powerEstimate
cat("Selected soft-thresholding power:", soft_power, "\n")

# ------------------------------------------------------------
# 7. Network construction & module detection
# ------------------------------------------------------------
net <- blockwiseModules(
  df_expr,
  power = soft_power,
  TOMType = "unsigned",
  minModuleSize = 15,          # lowered from the default 30: this panel has only 489 genes
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "output/metabric_TOM",
  verbose = 3
)

table(net$colors)
moduleColors <- labels2colors(net$colors)

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

png("output/figures/dendrogram_modul.png", width = 1200, height = 800, res = 150)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                     "Module colors",
                     dendroLabels = FALSE, hang = 0.05,
                     addGuide = TRUE, guideHang = 0.08)
dev.off()

# ------------------------------------------------------------
# 8. Module membership
# ------------------------------------------------------------
gene_module_df <- data.frame(gene = colnames(df_expr), module = moduleColors)
table(gene_module_df$module)

# per-module gene lists
split(gene_module_df$gene, gene_module_df$module)

# ------------------------------------------------------------
# 9. Module-trait relationships
# ------------------------------------------------------------
# cancer_type is excluded: 1903/1904 samples fall in a single category
table(df_clinical$cancer_type)
table(df_clinical$pam50_._claudin.low_subtype)

# drop the "NC" (not classified) PAM50 category
df_clinical_filtered <- df_clinical[df_clinical$pam50_._claudin.low_subtype != "NC", ]
df_expr_filtered <- df_expr[rownames(df_clinical_filtered), ]

MEs0 <- moduleEigengenes(df_expr_filtered, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

# PAM50 subtype is nominal (unordered) -> one-hot encode, not label-encode
trait_categorical <- dummy_cols(
  df_clinical_filtered["pam50_._claudin.low_subtype"],
  select_columns = "pam50_._claudin.low_subtype",
  remove_first_dummy = FALSE,
  remove_selected_columns = TRUE
)

chosen_traits <- c("tumor_size", "tumor_stage")
trait_numeric <- df_clinical_filtered[, chosen_traits]

traits_combined <- cbind(trait_numeric, trait_categorical)

stopifnot(nrow(MEs) == nrow(traits_combined))

moduleTraitCor <- cor(MEs, traits_combined, use = "pairwise.complete.obs")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(traits_combined))

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

png("output/figures/module_trait_relationship.png", width = 1800, height = 1100, res = 150)
par(mar = c(14, 8.5, 3, 3))
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(traits_combined),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  main = "Module-trait relationships"
)
dev.off()

# ------------------------------------------------------------
# 10. Export results
# ------------------------------------------------------------
write.csv(gene_module_df, "output/tables/metabric_gene_modules.csv", row.names = FALSE)

for (mod in unique(gene_module_df$module)) {
  genes_in_module <- gene_module_df$gene[gene_module_df$module == mod]
  write.csv(
    data.frame(gene = genes_in_module),
    paste0("output/tables/module_", mod, "_genes.csv"),
    row.names = FALSE
  )
}

moduleTraitCor_df <- as.data.frame(moduleTraitCor)
moduleTraitCor_df$module <- rownames(moduleTraitCor_df)
write.csv(moduleTraitCor_df, "output/tables/module_trait_correlation.csv", row.names = FALSE)

moduleTraitPvalue_df <- as.data.frame(moduleTraitPvalue)
moduleTraitPvalue_df$module <- rownames(moduleTraitPvalue_df)
write.csv(moduleTraitPvalue_df, "output/tables/module_trait_pvalue.csv", row.names = FALSE)

save(net, moduleColors, MEs, df_expr, df_clinical,
     moduleTraitCor, moduleTraitPvalue,
     file = "output/metabric_WGCNA_result.RData")

cat("Done. Results written to output/figures/ and output/tables/\n")
