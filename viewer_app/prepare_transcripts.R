# ============================================================
# prepare_transcripts.R
# ------------------------------------------------------------
# Pre-fetches exon/transcript structure for a fixed list of genes
# from Ensembl (GRCh37/hg19) and caches it to transcripts.rds.
#
# app.R loads transcripts.rds automatically if present.
#
# Edit GENES below, then run:
#   Rscript prepare_transcripts.R
# or:
#   source("prepare_transcripts.R")
#
# Re-run this whenever you want to add more genes to the cache, or
# refresh it (e.g. if Ensembl's annotation changes).
#
# NOTE ON FAILED GENES: this script batches all genes into just two
# large queries (rather than one pair of queries per gene) to avoid
# tripping rate limits / transient outages on Ensembl's GRCh37
# archive server — that was the most common cause of "Could NOT
# fetch" for genes that are definitely valid symbols. If a gene
# still fails after this, it's more likely a real symbol mismatch:
# GRCh37's annotation was frozen around 2013 (Ensembl 75), so genes
# renamed by HGNC since then are only findable under their OLD name.
# A few known examples — try the old name if the current one fails:
#   AFG2B -> SPATA5L1      GBA1 -> GBA
#   MMUT  -> MUT           PRKN -> PARK2
# (worth double-checking against HGNC's symbol history if others
# turn up missing too)
# ============================================================

library(biomaRt)

# Candidate NON-ESTABLISHED genes from the pipeline (06_additional_filtering_save.Rmd).
# Path is relative to viewer_app/, so run this script from that folder.
INPUT_FILE <- "../output/results/Non_ESTABLISHED_genes_with_hits_SKSC_3_MAF1perrank1_orabove1_remove_varsPhenoScore0.5_no_monogen_old_gene_Exomiser14.0.0.txt"
GENES <- read.delim(INPUT_FILE)$GENE_SYMBOL   # <- put your gene list here
OUTPUT_FILE <- "transcripts.rds"

# Swap current HGNC symbols for the old GRCh37-era names Ensembl's
# frozen archive still uses
old_symbol_map <- c(
  AFG2B = "SPATA5L1",
  GBA1  = "GBA",
  MMUT  = "MUT",
  PRKN  = "PARK2"
)
GENES <- ifelse(GENES %in% names(old_symbol_map), old_symbol_map[GENES], GENES)

# Same connect-with-retry logic as app.R (duplicated here rather than
# shared, so this script stays a standalone, runnable-on-its-own file).
connect_ensembl <- function() {
  attempts <- list(
    function() useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", GRCh = 37),
    function() useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", GRCh = 37),
    function() useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                        dataset = "hsapiens_gene_ensembl",
                        host = "https://grch37.ensembl.org")
  )
  for (i in seq_along(attempts)) {
    result <- tryCatch(attempts[[i]](), error = function(e) e)
    if (!inherits(result, "error")) return(result)
    if (i < length(attempts)) Sys.sleep(2)
  }
  stop("Ensembl (GRCh37) appears to be unavailable right now. Original error: ",
       conditionMessage(result))
}

# Retries a single getBM() call a couple of times before giving up —
# large batched queries are more reliable than many small ones, but
# still worth retrying individually if one transiently fails.
getBM_retry <- function(..., attempts = 3, wait = 3) {
  for (i in seq_len(attempts)) {
    result <- tryCatch(getBM(...), error = function(e) e)
    if (!inherits(result, "error")) return(result)
    if (i < attempts) Sys.sleep(wait)
  }
  stop("getBM() failed after ", attempts, " attempts: ", conditionMessage(result))
}

mart <- connect_ensembl()
genes_upper <- toupper(GENES)

message("Fetching transcripts for ", length(genes_upper), " gene(s) in one batch ...")

# ---- Step 1: all transcripts for all genes, in a single call ----
tx_info <- getBM_retry(
  attributes = c("hgnc_symbol", "ensembl_transcript_id", "transcript_biotype"),
  filters = "hgnc_symbol",
  values = genes_upper,
  mart = mart
)

found_genes <- unique(toupper(tx_info$hgnc_symbol))
not_found_at_all <- setdiff(genes_upper, found_genes)

# ---- Step 2: pick candidate transcripts per gene, then fetch exon
#              structure for ALL of them in a single call ----
tx_by_gene <- split(tx_info, toupper(tx_info$hgnc_symbol))

candidate_ids <- unique(unlist(lapply(tx_by_gene, function(df) {
  cand <- df$ensembl_transcript_id[df$transcript_biotype == "protein_coding"]
  if (length(cand) == 0) cand <- df$ensembl_transcript_id
  cand
})))

message("Fetching exon structure for ", length(candidate_ids), " candidate transcript(s) ...")

exons_all <- if (length(candidate_ids) > 0) {
  getBM_retry(
    attributes = c("ensembl_transcript_id", "ensembl_exon_id",
                   "exon_chrom_start", "exon_chrom_end",
                   "chromosome_name", "strand", "rank"),
    filters = "ensembl_transcript_id",
    values = candidate_ids,
    mart = mart
  )
} else {
  data.frame()
}

# ---- Step 3: assemble the per-gene cache entries ----
transcripts_cache <- list()
no_exon_structure <- character(0)

for (symbol in names(tx_by_gene)) {
  df <- tx_by_gene[[symbol]]
  candidates <- df$ensembl_transcript_id[df$transcript_biotype == "protein_coding"]
  if (length(candidates) == 0) candidates <- df$ensembl_transcript_id

  gene_exons <- exons_all[exons_all$ensembl_transcript_id %in% candidates, ]
  if (nrow(gene_exons) == 0) {
    no_exon_structure <- c(no_exon_structure, symbol)
    next
  }

  # Heuristic stand-in for "canonical": the transcript with the most
  # exons — same approach as app.R's live-fetch fallback.
  exon_counts <- table(gene_exons$ensembl_transcript_id)
  tx_id <- names(exon_counts)[which.max(exon_counts)]
  exons <- gene_exons[gene_exons$ensembl_transcript_id == tx_id, ]
  exons <- exons[order(exons$exon_chrom_start), ]

  transcripts_cache[[symbol]] <- list(symbol = symbol, transcript_id = tx_id, exons = exons)
}

# Store renamed genes under their CURRENT symbol (e.g. MUT -> MMUT), so the
# app finds them in the variant file, which uses current symbols
current_symbol <- setNames(names(old_symbol_map), old_symbol_map)
for (old in intersect(names(transcripts_cache), names(current_symbol))) {
  new <- current_symbol[[old]]
  transcripts_cache[[new]] <- transcripts_cache[[old]]
  transcripts_cache[[new]]$symbol <- new
  transcripts_cache[[old]] <- NULL
}
transcripts_cache <- transcripts_cache[sort(names(transcripts_cache))]

saveRDS(transcripts_cache, OUTPUT_FILE)

cat("\nWrote", length(transcripts_cache), "gene(s) to", OUTPUT_FILE, "\n")
cat("Genes cached:", paste(names(transcripts_cache), collapse = ", "), "\n")

if (length(not_found_at_all) > 0) {
  cat("\nNo transcript found at all for:", paste(not_found_at_all, collapse = ", "), "\n")
  cat("(likely a renamed/outdated HGNC symbol on this frozen GRCh37 build — see the\n",
      " note at the top of this script for a few known renames)\n", sep = "")
}
if (length(no_exon_structure) > 0) {
  cat("\nTranscript(s) found but no exon structure returned for:",
      paste(no_exon_structure, collapse = ", "), "\n")
}
