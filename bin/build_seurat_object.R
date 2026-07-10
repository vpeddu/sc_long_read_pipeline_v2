args <- commandArgs(trailingOnly=TRUE)
sample <- args[1]
gene_matrix_fn <- args[2]
gene_features_fn <- args[3]
gene_barcodes_fn <- args[4]
iso_matrix_fn <- args[5]
iso_features_fn <- args[6]
iso_barcodes_fn <- args[7]
iso_ranges_fn <- args[8]
gene_ranges_fn <- args[9]
sqanti_classification_fn <- args[10]
sample_metrics_fn <- args[11]
transcript_xref_fn <- args[12]
# Only present when --predict_orfs is true (args[13:14] are NA otherwise --
# indexing an R character vector past its length returns NA, not an error).
orf_pep_fn <- args[13]
pfam_domtblout_fn <- args[14]

library(Matrix)
library(Seurat)

# Seurat replaces underscores with dashes in feature names on assay creation
# (eg "GIM3802_000001" -> "GIM3802-000001", "novel_intergenic_000001" ->
# "novel-intergenic-000001"); apply the same substitution to any join key
# sourced from files that still have the original underscore-based ids.
to_seurat_name <- function(x) gsub("_", "-", x)

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

# chrom/start/end/strand/id TSV (no header), as produced by the awk extraction
# in modules/A06_seurat.nf from the transcript/gene rows of a GTF. A handful of
# gene symbols in the reference GTF map to more than one Ensembl gene id (the
# same collision rename_gtf_txid.py already documents/works around for gene
# names) -- keep only the first occurrence of each id so rownames stay unique;
# genomic ranges are an enrichment, not core data, so under-covering the rare
# duplicate is preferable to crashing or guessing which entry is "right".
read_ranges <- function(ranges_fn, id_col_name) {
  df <- read.table(ranges_fn, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                    col.names = c("chrom", "start", "end", "strand", id_col_name))
  df <- df[!duplicated(df[[id_col_name]]), ]
  rownames(df) <- to_seurat_name(df[[id_col_name]])
  df[[id_col_name]] <- NULL
  return(df)
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
so[["ISO"]] <- CreateAssay5Object(counts = iso_mtx[, common_barcodes])

# Tag every cell with its sample of origin (orig.ident) and make barcodes
# globally unique across samples for when per-sample objects are later merged.
so$orig.ident <- sample
so <- RenameCells(so, add.cell.id = sample)

# --- Sample-level metrics (raw/pre-dedup/post-dedup read counts) -----------
# Broadcast as per-cell columns (standard Seurat convention, eg orig.ident)
# so they're usable directly in VlnPlot/FetchData, and also stashed as a
# compact list in @misc for programmatic access.
metrics <- read.table(sample_metrics_fn, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
metrics_list <- setNames(as.list(metrics$value), metrics$metric)
for (name in names(metrics_list)) {
  so[[name]] <- as.numeric(metrics_list[[name]])
}
so$dedup_rate <- so$postdup_reads / so$predup_reads
so@misc$sample_metrics <- c(metrics_list, list(dedup_rate = unique(so$dedup_rate)))

# --- Isoform-level metadata: genomic ranges + SQANTI3 classification -------
# Both are keyed by isoform/transcript id and joined onto the ISO assay's
# feature metadata; Seurat aligns by rowname and NA-fills anything that
# doesn't match (eg isoforms SQANTI3's rules filter dropped), so no manual
# reindexing to the assay's feature order is needed.
iso_ranges <- read_ranges(iso_ranges_fn, "transcript_id")
so[["ISO"]][[colnames(iso_ranges)]] <- iso_ranges

sqanti <- read.table(sqanti_classification_fn, sep = "\t", header = TRUE,
                      stringsAsFactors = FALSE, quote = "", comment.char = "")
sqanti_cols <- c("structural_category", "associated_gene", "associated_transcript",
                  "length", "exons", "all_canonical", "min_cov", "FL", "coding",
                  "predicted_NMD", "filter_result", "perc_A_downstream_TTS", "RTS_stage")
sqanti_meta <- sqanti[, sqanti_cols]
rownames(sqanti_meta) <- to_seurat_name(sqanti$isoform)
so[["ISO"]][[colnames(sqanti_meta)]] <- sqanti_meta

# --- Isoform-level metadata: flair's own per-transcript read support -------
# transcript_xref.tsv's "score" column carries flair's supporting-read count
# for that transcript (SQANTI3's own corrected gtf hardcodes score as "." and
# loses this, which is why it's read from the xref instead). Keyed on
# whichever id ended up as this transcript's final name -- novel_tx for novel
# transcripts, else the unchanged transcript_id -- the same fallback
# rename_gtf_txid.py itself uses when it fills tx_name.
xref <- read.table(transcript_xref_fn, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                    quote = "", comment.char = "")
final_id <- ifelse(!is.na(xref$novel_tx) & xref$novel_tx != "", xref$novel_tx, xref$transcript_id)
read_support <- data.frame(read_support = xref$score, id = to_seurat_name(final_id), stringsAsFactors = FALSE)
read_support <- read_support[!duplicated(read_support$id), ]
rownames(read_support) <- read_support$id
read_support$id <- NULL
so[["ISO"]][[colnames(read_support)]] <- read_support

# --- Isoform-level metadata: ORF prediction + Pfam domain classification ---
# TD2 (+ its internal PSAURON coding-likelihood scorer) and a Pfam hmmsearch,
# run fully independently of the isoSeQL-pinned SQANTI3 branch (modules/
# A08_orf.nf) since that fork's own ORF prediction needs the license-gated
# GeneMarkS-T. Transcript ids get a "trailing .p<N>" ORF suffix from TD2
# (eg "ISO_REAL_ORF.p2"); strip it to recover the id used elsewhere.
strip_orf_suffix <- function(x) to_seurat_name(sub("\\.p[0-9]+$", "", x))

if (!is.na(orf_pep_fn) && file.exists(orf_pep_fn)) {
  headers <- grep("^>", readLines(orf_pep_fn), value = TRUE)
  if (length(headers) > 0) {
    orf_meta <- data.frame(
      id = strip_orf_suffix(sub("^>(\\S+).*", "\\1", headers)),
      orf_type = sub(".*\\bORF type:(\\S+).*", "\\1", headers),
      orf_len = as.numeric(sub(".*\\blen:([0-9]+).*", "\\1", headers)),
      psauron_score = as.numeric(sub(".*psauron_score=([0-9.eE+-]+).*", "\\1", headers)),
      stringsAsFactors = FALSE
    )
    orf_meta <- orf_meta[!duplicated(orf_meta$id), ]
    rownames(orf_meta) <- orf_meta$id
    orf_meta$id <- NULL
    so[["ISO"]][[colnames(orf_meta)]] <- orf_meta
  }
}

if (!is.na(pfam_domtblout_fn) && file.exists(pfam_domtblout_fn)) {
  dom_lines <- readLines(pfam_domtblout_fn)
  dom_lines <- dom_lines[!grepl("^#", dom_lines)]
  if (length(dom_lines) > 0) {
    parts <- strsplit(trimws(dom_lines), "\\s+")
    dom <- data.frame(
      id = strip_orf_suffix(sapply(parts, `[`, 1)),
      domain = sapply(parts, `[`, 4),
      evalue = as.numeric(sapply(parts, `[`, 7)),
      stringsAsFactors = FALSE
    )
    dom_meta <- do.call(rbind, lapply(split(dom, dom$id), function(d) {
      best <- d[which.min(d$evalue), ]
      data.frame(pfam_domain_count = nrow(d), pfam_best_domain = best$domain,
                 pfam_best_evalue = best$evalue, stringsAsFactors = FALSE)
    }))
    so[["ISO"]][[colnames(dom_meta)]] <- dom_meta
  }
}

# --- Gene-level metadata: genomic ranges ------------------------------------
gene_ranges <- read_ranges(gene_ranges_fn, "gene_name")
so[["RNA"]][[colnames(gene_ranges)]] <- gene_ranges

saveRDS(so, file = paste0(sample, ".seurat.rds"))
