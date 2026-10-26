# Seurat v5 objects from the Pancreas sample, for build-format-fixtures.py.
#
#   Rscript build-format-fixtures.R EXPORT_DIR OUT_DIR
#
# EXPORT_DIR holds what the Python half exported (X.mtx genes x cells, the
# cell and gene names, obs, UMAP, PCA). Writes OUT_DIR/pancreas.rds (one
# Seurat object, as saveRDS() writes it) and OUT_DIR/pancreas-objects.rda
# (two Seurat objects: the whole sample and its ductal cells, as save()
# writes several objects).
args <- commandArgs(TRUE)
export <- args[[1]]
out <- args[[2]]
suppressPackageStartupMessages({
  library(Matrix)
  library(SeuratObject)
})
counts <- as(readMM(file.path(export, "X.mtx")), "CsparseMatrix")
cells <- readLines(file.path(export, "cells.txt"))
genes <- make.unique(gsub("_", "-", readLines(file.path(export, "genes.txt"))))
dimnames(counts) <- list(genes, cells)
obs <- read.csv(file.path(export, "obs.csv"), row.names = 1, check.names = FALSE)
obs$clusters <- factor(obs$clusters)
obs$clusters_coarse <- factor(obs$clusters_coarse)
umap <- as.matrix(read.csv(file.path(export, "umap.csv"), row.names = 1))
colnames(umap) <- c("umap_1", "umap_2")
pca <- as.matrix(read.csv(file.path(export, "pca.csv"), row.names = 1))
colnames(pca) <- paste0("PC_", seq_len(ncol(pca)))

options(Seurat.object.assay.version = "v5")
pancreas <- CreateSeuratObject(counts = counts, meta.data = obs, project = "pancreas")
LayerData(pancreas, layer = "data") <- log1p(counts)
pancreas[["pca"]] <- CreateDimReducObject(embeddings = pca, key = "PC_", assay = "RNA")
pancreas[["umap"]] <- CreateDimReducObject(embeddings = umap, key = "umap_", assay = "RNA")
Idents(pancreas) <- "clusters"
saveRDS(pancreas, file.path(out, "pancreas.rds"))

pancreas_ductal <- subset(pancreas, idents = "Ductal")
save(pancreas, pancreas_ductal, file = file.path(out, "pancreas-objects.rda"))
cat("R fixtures written\n")
