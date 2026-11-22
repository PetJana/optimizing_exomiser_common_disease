# ---- Render the whole pipeline in order ------------------------------------------
# Works from any working directory: Rscript run_all.R, source("run_all.R"),
# or RStudio's Source button.
# HTML reports go to output/reports/, data and figures to output/ (see 00_setup.R).

# Folder containing this script (and the .Rmd files)
script_dir <- function() {
  file_arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(file_arg) == 1) return(dirname(normalizePath(file_arg)))

  for (i in rev(seq_len(sys.nframe()))) {
    f <- sys.frame(i)$ofile
    if (!is.null(f)) return(dirname(normalizePath(f)))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }

  stop("Could not locate the script folder; setwd() to the repository folder (containing run_all.R) and try again.")
}

pipeline_dir <- script_dir()

steps <- c(
  "01_benchmarking_combine_exomiser_results_multirun.Rmd",
  "02_read_filter_save.Rmd",
  "03_benchmarking_figure_2A.Rmd",
  "04_make_figure_2B_violin.Rmd",
  "05_make_figure_2C_acmg.Rmd",
  "06_additional_filtering_save.Rmd",
  "07_go_kegg_enrichment.Rmd",
  "08_top_ranked_per_participant.Rmd",
  "09_gene_level_tests_overlap.Rmd",
  "10_combine_figures.Rmd",
  "11_magma_manhattan_plots.Rmd"
)

for (step in steps) {
  message("==== ", step)
  rmarkdown::render(file.path(pipeline_dir, step),
                    output_dir = file.path(pipeline_dir, "output", "reports"),
                    knit_root_dir = pipeline_dir,
                    envir = new.env(), quiet = TRUE)
}
