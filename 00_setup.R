# ---- Shared setup, sourced by every .Rmd --------------------------------------
# Inputs are read in place; everything the pipeline writes goes to output/
# next to these scripts.

library(tidyverse)

# ---- Run settings -------------------------------------------------------------

exomiser_version <- "14.0.0"
my_run           <- "SKSC_3_MAF1per"            # ESTABLISHED run
nonhit_run       <- paste0(my_run, "_nonHIT")   # NON-ESTABLISHED run
remove_vars      <- "_remove_vars"              # "_remove_vars" or ""
genes_to_remove  <- c("PAPSS2", "FAM20C")       # dropped when remove_vars is set

# Score cutoffs
var_score_cutoff             <- 0.8
pheno_score_cutoff_benchmark <- 0.2  # benchmarking (Figure 2A)
pheno_score_cutoff_candidate <- 0.5  # candidate gene search (step 6)

# ---- Input paths --------------------------------------------------------------

exomiser_dir      <- file.path(paste0("Exomiser_", exomiser_version))  # one folder per run
manual_vars_file  <- "Vars_from_Munch_2026_Suppl_table_green_only.txt"
study_filter_file <- "ids_filtering/IDs_for_study_filtering.txt"
ks_genes_file     <- "KS_genes.txt"

# Participants whose manually annotated variant is intronic (not scored by Exomiser);
# one SKSC ID per line, kept out of version control
intronic_carriers_file <- "ids_filtering/intronic_carriers.txt"
intronic_carriers      <- readLines(intronic_carriers_file)

# ---- Output paths -------------------------------------------------------------

output_dir  <- "output"
results_dir <- file.path(output_dir, "results")
counts_dir  <- file.path(output_dir, "counts_per_rank")
figures_dir <- file.path(output_dir, "figures")
images_dir  <- file.path(figures_dir, "images")

walk(c(results_dir, counts_dir, images_dir), dir.create, recursive = TRUE, showWarnings = FALSE)

file_suffix <- paste0(remove_vars, "_", exomiser_version)

established_file    <- file.path(results_dir, paste0("All_vars_ESTABLISHED_", my_run, file_suffix, ".tsv"))
nonestablished_file <- file.path(results_dir, paste0("All_vars_NONESTABLISHED_", my_run, file_suffix, ".tsv"))

# Step 6 outputs (inputs to step 7)
candidate_max_rank  <- 1  # only top-ranked genes
candidate_tail      <- paste0(remove_vars, "PhenoScore", pheno_score_cutoff_candidate,
                              "_no_monogen_old_gene_Exomiser", exomiser_version, ".txt")
candidate_prefix    <- paste0("Non_ESTABLISHED_genes_with_hits_", my_run, "rank", candidate_max_rank)
all_vars_file       <- file.path(results_dir, paste0("All_vars_", my_run, file_suffix, ".tsv"))
candidates_file     <- file.path(results_dir, paste0(candidate_prefix, "_orabove1", candidate_tail))
candidate_vars_file <- file.path(results_dir, paste0(candidate_prefix, "_potential_monogen_new", candidate_tail))

counts_file <- function(run) {
  file.path(counts_dir, paste0(run, "_counts_per_rank", remove_vars, ".txt"))
}

# ---- Plot colours -------------------------------------------------------------

colorBlindGrey8 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

group_colors <- c("ESTABLISHED"     = colorBlindGrey8[[7]],
                  "NON-ESTABLISHED" = colorBlindGrey8[[6]])

# ---- Helpers ------------------------------------------------------------------

# Combine the Exomiser variant tables of all participants in one run
read_exomiser_run <- function(run) {
  res_files <- list.files(file.path(exomiser_dir, run),
                          pattern = "\\.variants\\.tsv$", full.names = TRUE)
  if (length(res_files) == 0) stop("No *.variants.tsv files found for run ", run)

  res_files %>%
    set_names(~ str_remove(basename(.x), "\\.hard-filtered-exomiser\\.variants\\.tsv$")) %>%
    map(\(f) read_tsv(f, show_col_types = FALSE) %>% rename(RANK = 1)) %>%  # "#RANK"
    list_rbind(names_to = "Participant_ID")
}

# Drop all variants in `genes` and close the rank gaps they leave: each
# remaining gene moves up by the number of distinct removed ranks above it
# for the same participant.
remove_genes_and_rerank <- function(vars, genes = genes_to_remove) {
  removed <- vars %>% filter(GENE_SYMBOL %in% genes)
  message("Removing ", nrow(removed), " rows in ", paste(genes, collapse = ", "))
  print(table(removed$MOI))

  removed_ranks <- split(removed$RANK, removed$Participant_ID) %>% map(\(r) sort(unique(r)))

  vars %>%
    filter(!GENE_SYMBOL %in% genes) %>%
    group_by(Participant_ID) %>%
    mutate(RANK = RANK - findInterval(RANK - 0.5,
                                      removed_ranks[[first(Participant_ID)]] %||% numeric())) %>%
    ungroup()
}

maybe_remove_genes <- function(vars) {
  if (remove_vars == "_remove_vars") remove_genes_and_rerank(vars) else vars
}

# Participants passing study QC
read_study_filter <- function() {
  read_tsv(study_filter_file, show_col_types = FALSE) %>%
    filter(!FILTERED)
}

# Manually confirmed KS gene / participant pairs, keyed by VCF file ID
read_manual_genes <- function(data_filter = read_study_filter()) {
  vars_john <- read_tsv(manual_vars_file, show_col_types = FALSE)

  message("vars_john: ", nrow(vars_john), " rows, ",
          n_distinct(vars_john$ID), " unique IDs, ",
          n_distinct(vars_john$VAR), " unique variants")

  vars_john <- vars_john %>%
    filter(!ID %in% intronic_carriers) %>%
    inner_join(data_filter, by = c("ID" = "SKSC_ID")) %>%
    rename(GENE_SYMBOL = GENE, Participant_ID = VCF_FILE_ID)

  message("vars_john after merge: ", nrow(vars_john), " rows, ",
          n_distinct(vars_john$VAR), " unique variants")
  print(table(vars_john$INH))

  genes_john <- vars_john %>% distinct(Participant_ID, GENE_SYMBOL)

  message("Carriers: ", n_distinct(genes_john$Participant_ID),
          " | Genes: ", n_distinct(genes_john$GENE_SYMBOL),
          " | Pairs: ", nrow(genes_john))

  genes_john
}

# Keep the best-ranked row per gene / participant (collapses inheritance modes).
# Same result as slice_min(RANK, with_ties = FALSE) per group (arrange() is
# stable, so ties keep the first row), but fast on a full raw run.
best_rank_per_gene <- function(vars) {
  vars %>%
    arrange(RANK) %>%
    distinct(GENE_SYMBOL, Participant_ID, .keep_all = TRUE)
}

read_all_vars <- function() {
  bind_rows(
    read_tsv(established_file, show_col_types = FALSE),
    read_tsv(nonestablished_file, show_col_types = FALSE)
  )
}
