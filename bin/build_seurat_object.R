args <- commandArgs(trailingOnly=TRUE)
sample <- args[1]
gene_matrix_fn <- args[2]
gene_features_fn <- args[3]
gene_barcodes_fn <- args[4]
iso_matrix_fn <- args[5]
iso_features_fn <- args[6]
iso_barcodes_fn <- args[7]

library(Matrix)
library(Seurat)

read_sparse_matrix <- function(matrix_fn, features_fn, barcodes_fn) {
  mtx <- readMM(matrix_fn)
  features <- readLines(features_fn)
  barcodes <- readLines(barcodes_fn)
  # rdmap2counts.py/htseq-count-barcodes both write (cells x features); Seurat
  # wants (features x cells) -- transpose whichever way matches the features length.
  if (dim(mtx)[2] == length(features)) {
    mtx <- t(mtx)
  }
  rownames(mtx) <- make.unique(features)
  colnames(mtx) <- barcodes
  return(mtx)
}

gene_mtx <- read_sparse_matrix(gene_matrix_fn, gene_features_fn, gene_barcodes_fn)
iso_mtx <- read_sparse_matrix(iso_matrix_fn, iso_features_fn, iso_barcodes_fn)

# htseq-count-barcodes' barcodes already carry a "-1" suffix (added in
# runHtseq); rdmap2counts.py's do not (its own comment defers this to "a
# subsequent step, when needed for Seurat load" -- this is that step).
if (!any(grepl("-1$", colnames(iso_mtx)))) {
  colnames(iso_mtx) <- paste0(colnames(iso_mtx), "-1")
}

common_barcodes <- intersect(colnames(gene_mtx), colnames(iso_mtx))
cat(sprintf("Gene-level cells: %d, isoform-level cells: %d, shared: %d\n",
            ncol(gene_mtx), ncol(iso_mtx), length(common_barcodes)))

so <- CreateSeuratObject(counts = gene_mtx[, common_barcodes], project = sample, assay = "RNA")
so[["ISO"]] <- CreateAssayObject(counts = iso_mtx[, common_barcodes])

# Tag every cell with its sample of origin (orig.ident) and make barcodes
# globally unique across samples for when per-sample objects are later merged.
so$orig.ident <- sample
so <- RenameCells(so, add.cell.id = sample)

saveRDS(so, file = paste0(sample, ".seurat.rds"))
