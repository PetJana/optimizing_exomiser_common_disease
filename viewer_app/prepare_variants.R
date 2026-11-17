# ============================================================
# prepare_variants.R
# ------------------------------------------------------------
# One-time (or re-run-whenever-your-data-updates) transform:
# converts a per-carrier variant dataset (one row per variant per
# carrier) into the per-gene list of per-variant data.frames that
# app.R expects — see get_variants() in app.R for the target shape.
#
# Input assumptions — EDIT THESE to match your actual file:
#   - one row per (variant, carrier) pair
#   - columns: gene_symbol, pos, label, type, carrier_id
#   - pos is GRCh37/hg19 (must match app.R's pinned assembly —
#     mismatched positions silently disappear from the plot rather
#     than erroring, so double-check this before trusting the app)
#
# Output: written next to this script, using a filename that matches
# one of the entries in VARIANT_FILE_OPTIONS at the top of app.R:
#   variants_exomiser_14.0.0.rds  or  variants_exomiser_15.1.0.rds
# app.R loads whichever file matches the dropdown selection, falling
# back to illustrative example data if that version's file isn't found.
# Run this script once per Exomiser version (change INPUT_FILE and
# OUTPUT_FILE below each time).
#
# Run with:
#   Rscript prepare_variants.R
# or from an R session:
#   source("prepare_variants.R")
# ============================================================

Sys.setenv(LANG = "en")

library(tidyverse)

# Candidate variants from the pipeline (06_additional_filtering_save.Rmd).
# Path is relative to viewer_app/, so run this script from that folder.
INPUT_FILE  <- "../output/results/Non_ESTABLISHED_genes_with_hits_SKSC_3_MAF1perrank1_potential_monogen_new_remove_varsPhenoScore0.5_no_monogen_old_gene_Exomiser14.0.0.txt"
OUTPUT_FILE <- "variants_exomiser_14.0.0.rds"        # <- must match a name in VARIANT_FILE_OPTIONS in app.R


# ---- load ----
raw <- read_tsv(INPUT_FILE, show_col_types = FALSE)

# rename columns 
lookup <- c(gene_symbol = "GENE_SYMBOL", pos = "START", label = "HGVS", type = "FUNCTIONAL_CLASS", carrier_id = "Participant_ID")
raw <- rename(raw, all_of(lookup))

# optimize "type" column
raw$type[raw$type %in% c("stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant", "transcript_ablation", "frameshift_elongation", "frameshift_truncation", "splice_donor_variant")] <- "lof"
raw$type <- gsub("_variant", "", raw$type)
raw$type <- gsub("5", "five", raw$type)

# extracting c-position and p.position from label
raw$c.position  <- str_split_i(raw$label, ":", 3)
raw$p.position <- str_split_i(raw$label, ":", 4)

# format acmg
acmg_map <- c(
  PATHOGENIC = "P",
  LIKELY_PATHOGENIC  = "LP",
  UNCERTAIN_SIGNIFICANCE  = "US"
)
raw$acmg <- raw$EXOMISER_ACMG_CLASSIFICATION
raw$acmg <- ifelse(raw$acmg %in% names(acmg_map), acmg_map[raw$acmg], raw$acmg)

# format MAF
raw$max_freq_source <- raw$MAX_FREQ_SOURCE
raw$max_freq <- signif(raw$MAX_FREQ, 2)

# format path
raw$max_path_source <- raw$MAX_PATH_SOURCE
raw$max_path <- round(raw$MAX_PATH, 2)

# format rsid
raw$rs_id <- raw$RS_ID

raw$ref <- raw$REF
raw$alt <- raw$ALT

# ---- sanity check expected columns ----
required_cols <- c("gene_symbol", "pos", "label", "ref", "alt", "type", "carrier_id", "c.position", "p.position", "rs_id", "max_freq_source", "max_freq", "max_path_source", "max_path", "acmg")
missing_cols <- setdiff(required_cols, names(raw))
if (length(missing_cols) > 0) {
  stop("Input file is missing expected column(s): ",
       paste(missing_cols, collapse = ", "),
       "\nEither rename your columns to match, or edit the ",
       "group_by()/select() calls below to use your actual column names: ",
       paste(names(raw), collapse = ", "))
}

# ---- aggregate: one row per unique variant, across all its carriers ----
variants_agg <- raw %>%
  group_by(gene_symbol, pos, label, type, ref, alt, c.position, p.position, rs_id, max_freq_source, max_freq, max_path_source, max_path, acmg) %>%
  summarise(
    n_carriers = n_distinct(carrier_id),
    n_controls = n_distinct(grepv("^G", carrier_id)),
    n_cases = n_carriers - n_controls, 
    # `score` drives the lollipop stem height in the plot. Using
    # carrier count is a reasonable default — swap in a different
    # column/formula here if you have a better severity/frequency
    # measure to plot instead.
    score      = n_carriers,
    carriers   = paste0("Observed in ", n_cases, " case(s) and ", n_controls, " control(s)"),
    .groups = "drop"
  )

if (nrow(variants_agg) == 0) {
  stop("No rows after aggregation — check that ", INPUT_FILE, " loaded correctly.")
}

# ---- split into the per-gene list structure app.R expects ----
example_variants <- split(
  variants_agg %>%
    ungroup() %>%
    select(pos, label, type, score, carriers, ref, alt, c.position, p.position, rs_id, max_freq_source, max_freq, max_path_source, max_path, acmg) %>%
    as.data.frame(stringsAsFactors = FALSE),
  toupper(variants_agg$gene_symbol)
)

saveRDS(example_variants, OUTPUT_FILE)

cat("Wrote", length(example_variants), "gene(s) to", OUTPUT_FILE, "\n")
cat("Genes:", paste(names(example_variants), collapse = ", "), "\n")
cat("Total unique variants:", nrow(variants_agg), "\n")
