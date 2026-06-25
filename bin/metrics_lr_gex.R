args <- commandArgs(trailingOnly=TRUE)
sample <- args[1]

library(Matrix)
library(dplyr)
library(tidyr)
library(stringr)

#lr_metrics <- read.table(paste0("../", sample, ".metrics.txt"), sep='\t', header=F)
lr_metrics <- read.table(paste0("../", sample, ".counts.txt"), sep='\t', header=F)
colnames(lr_metrics) <- c('disease', 'sample', 'label_nr_reads')

lr_metrics <- lr_metrics %>% separate(label_nr_reads, c('text_label', 'nr_reads'), sep=' : ')
lr_metrics$nr_reads <- as.integer(lr_metrics$nr_reads)
dedup_row <- lr_metrics %>% filter(str_detect(text_label, 'de-duplicated reads'))
dedup_reads <- dedup_row[1, 'nr_reads']
bc_row <- lr_metrics %>% filter(str_detect(text_label, 'barcode reads'))
bc_reads <- bc_row[1, 'nr_reads']

fn_suffix = ifelse(file.exists(paste0(sample, '.matrix.mtx.gz')), '.gz', '')
lr_gex_mtx <- readMM(paste0(sample, '.matrix.mtx', fn_suffix))
lr_barcodes <- readLines(paste0(sample, '.barcodes.tsv', fn_suffix))
features <- readLines(paste0(sample, '.features.tsv', fn_suffix))
print(paste("Input matrix dimensions:", dim(lr_gex_mtx)[1], dim(lr_gex_mtx)[2]))

if (dim(lr_gex_mtx)[2] == length(features)) {
lr_gex_mtx <- t(lr_gex_mtx)
}
rownames(lr_gex_mtx) <- features
colnames(lr_gex_mtx) <- lr_barcodes

#lr_gex_mtx[1:10, 1:10]
nr_cells <- dim(lr_gex_mtx)[2]
#reads_per_cell <- dedup_reads / nr_cells
reads_per_cell <- bc_reads / nr_cells
counts_per_cell <- Matrix::colSums(lr_gex_mtx)
genes_per_cell <- Matrix::colSums(lr_gex_mtx>0)
cells_per_gene <- Matrix::rowSums(lr_gex_mtx)
nr_genes <- length(cells_per_gene[cells_per_gene > 0])

metrics <- c(paste("Cells detected :", nr_cells),
             paste("Median genes/cell :", median(genes_per_cell)),
             paste("Mean reads/cell :", round(reads_per_cell, 1)),
             paste("Mean transcripts/cell :", round(mean(counts_per_cell),1)),
             paste("Median transcripts/cell :", median(counts_per_cell)),
             paste("Genes detected : ", nr_genes)
			 )

metrics_df <- data.frame(label_text_nr=metrics)
metrics_df$disease <- lr_metrics[1, 'disease']
metrics_df$sample <- lr_metrics[1, 'sample']

write.table(metrics_df[, c('disease', 'sample', 'label_text_nr')], 
              paste0(sample, ".cell_metrics.txt"), sep='\t', row.names=F, col.names=F, quote=F)
                     