# ============================================================
# Gene Transcript Variant Viewer — R Shiny App
# ------------------------------------------------------------
# Displays genetic variants as a lollipop plot on top of a
# gene's transcript (exon/intron) structure.
#
# - Transcript structure: served ONLY from the pre-fetched cache
#   (transcripts.rds, built by prepare_transcripts.R). The gene
#   dropdown only offers genes present in that cache — there is no
#   live Ensembl lookup.
#
# Required packages (run once):
#   install.packages(c("shiny", "bslib", "DT", "ggplot2"))
#
# The lollipop plot is drawn with ggplot2 (lollipop_plot() below).
# It used to use trackViewer::lolliplot(), but trackViewer pulls in
# ~140 Bioconductor packages (~900 MB RAM), which crashed the app on
# shinyapps.io's ~1 GB free instance. The old version is kept in
# app_trackViewer.R.
#
# Run with:
#   shiny::runApp("app.R")
# ============================================================

library(shiny)
library(bslib)
library(DT)
library(ggplot2)

# ---------------------------------------------------------------
# Variant data.
# get_variants(gene_symbol, variants_list) must return a data.frame
# with columns: pos (genomic, 1-based), label, type, score, carriers.
#
# IMPORTANT: positions must be GRCh37/hg19, to match the assembly the
# transcript cache (transcripts.rds) was built from. Mixing assemblies
# silently drops variants — they get plotted outside the fetched exon
# range instead of throwing an error.
#
# Real data, per Exomiser version: one file holds ALL variants for
# that version. Generate the file by running prepare_variants.R (edit
# INPUT_FILE/OUTPUT_FILE at the top of that script). If it doesn't
# exist yet, the app falls back to the illustrative example data
# below rather than erroring.
#
# The version dropdown is commented out below (UI + here) — only one
# version is in use right now, so it's fixed via DEFAULT_VERSION_FILE.
# To bring back multiple versions later: uncomment EXOMISER_VERSION_FILES,
# uncomment the selectInput("exomiser_version", ...) in the UI, and swap
# DEFAULT_VERSION_FILE for EXOMISER_VERSION_FILES[[input$exomiser_version]]
# inside variants_source() in the server.
# EXOMISER_VERSION_FILES <- c(
#   "Exomiser 14.0.0" = "variants_exomiser_14.0.0.rds",
#   "Exomiser 15.1.0" = "variants_exomiser_15.1.0.rds"
# )
DEFAULT_VERSION_FILE <- "variants_exomiser_14.0.0.rds"

# Pre-fetched exon/transcript structure for a fixed gene list, generated
# by prepare_transcripts.R. This is the ONLY source of transcript data —
# the gene dropdown below only offers genes present here, and there is
# no live Ensembl fallback in this app.
TRANSCRIPTS_RDS <- "transcripts.rds"
transcripts_cache <- if (file.exists(TRANSCRIPTS_RDS)) readRDS(TRANSCRIPTS_RDS) else list()

# The gene dropdown only ever offers genes present in the cache — this
# is the fixed list from your prepare_transcripts.R run.
GENE_CHOICES <- sort(names(transcripts_cache))

default_example_variants <- list(
  TP53 = data.frame(
    pos   = c(7577121, 7577538, 7578190, 7578406, 7579472),
    label = c("R175H", "R248Q", "R273H", "Y220C", "R282W"),
    type  = c("missense", "missense", "missense", "missense", "missense"),
    score = c(5, 8, 9, 4, 6),
    carriers = c("Observed in 4 unrelated individuals (illustrative)",
                 "Observed in 7 unrelated individuals (illustrative)",
                 "Observed in 9 unrelated individuals (illustrative)",
                 "Observed in 3 unrelated individuals (illustrative)",
                 "Observed in 5 unrelated individuals (illustrative)"),
    stringsAsFactors = FALSE
  ),
  BRCA1 = data.frame(
    pos   = c(41245475, 41244000, 41234593, 41228577),
    label = c("5382insC", "185delAG", "Q1756fs", "R1699Q"),
    type  = c("frameshift", "frameshift", "frameshift", "missense"),
    score = c(9, 9, 7, 4),
    carriers = c("Observed in 6 unrelated individuals (illustrative)",
                 "Observed in 5 unrelated individuals (illustrative)",
                 "Observed in 2 unrelated individuals (illustrative)",
                 "Observed in 3 unrelated individuals (illustrative)"),
    stringsAsFactors = FALSE
  )
)

load_variants_file <- function(path) {
  if (!is.null(path) && file.exists(path)) readRDS(path) else default_example_variants
}

get_variants <- function(gene_symbol, variants_list) {
  v <- variants_list[[toupper(gene_symbol)]]
  if (is.null(v)) {
    return(data.frame(pos = numeric(0), label = character(0),
                       type = character(0), score = numeric(0),
                       carriers = character(0)))
  }
  v
}

# Plot colours by variant category: missense and lof (incl. start_lost
# and stop_lost) get their own colour, every other Exomiser type
# (splice_region, inframe indels, intergenic, ...) is grouped as
# "other". The table keeps the exact type. The legend lists
# categories in this order.
# Hex codes only: the legend is HTML/CSS, which doesn't know R names like "grey50".
type_colors <- c(missense = "#2E86AB", lof = "#C73E1D", other = "#7F7F7F")

LOF_TYPES <- c("lof", "start_lost", "stop_lost")

type_category <- function(type) {
  ifelse(type %in% LOF_TYPES, "lof",
         ifelse(type %in% names(type_colors), type, "other"))
}

# ---------------------------------------------------------------
# Lollipop plot (ggplot2 replacement for trackViewer::lolliplot).
# Gene model: backbone line + one box per exon. Each variant gets a
# stem from its true genomic position up to a stack of circles, one
# circle per carrier (score). When variants are too close to read,
# their circles are spread apart and the stem bends, like lolliplot.
# ---------------------------------------------------------------

# Spread sorted positions so neighbours are at least `gap` apart,
# staying inside [lo, hi] (as long as they fit at all).
spread_positions <- function(x, gap, lo, hi) {
  n <- length(x)
  if (n < 2) return(x)
  gap <- min(gap, (hi - lo) / (n - 1))
  d <- x
  for (i in 2:n) d[i] <- max(d[i], d[i - 1] + gap)
  d[n] <- min(d[n], hi)
  for (i in (n - 1):1) d[i] <- min(d[i], d[i + 1] - gap)
  d[1] <- max(d[1], lo)
  for (i in 2:n) d[i] <- max(d[i], d[i - 1] + gap)
  d
}

# Display label: the Exomiser label without the leading "GENE:", i.e.
# "ENST00000205557.7:c.3421C>T:p.(Arg1141*)". Variants without a cDNA
# change (intergenic/upstream, label "ENST...::") show the genomic change.
variant_short_label <- function(variants) {
  col <- function(name) if (is.null(variants[[name]])) rep(NA_character_, nrow(variants)) else as.character(variants[[name]])
  ok  <- function(a) !is.na(a) & nzchar(a)
  tx  <- sub("^[^:]+:", "", col("label"))
  g   <- ifelse(ok(col("ref")) & ok(col("alt")),
                paste0("g.", variants$pos, col("ref"), ">", col("alt")), col("label"))
  out <- ifelse(ok(col("c.position")) & ok(tx), tx, g)
  # Long insertions can be 100+ characters; the full notation is in the table.
  ifelse(nchar(out) > 50, paste0(substr(out, 1, 48), "…"), out)
}

lollipop_plot <- function(exons, variants, type_colors) {
  chrom <- exons$chromosome_name[1]

  # x range covers the exons and every variant (upstream/intergenic
  # variants can lie outside the transcript), plus a small margin.
  xr  <- range(c(exons$exon_chrom_start, exons$exon_chrom_end, variants$pos))
  pad <- max(diff(xr) * 0.03, 50)
  lo  <- xr[1] - pad
  hi  <- xr[2] + pad

  # y layout in (approximate) screen pixels of the 500px-high plot,
  # so circles and labels keep their size whatever the gene: exon boxes
  # around 0, stem rises to y_bend, bends over to y_top, circles stack
  # from y_circle upwards, one per carrier.
  exon_h   <- 7
  y_bend   <- 25
  y_top    <- 45
  y_circle <- 58
  step     <- 12   # circles of one variant just touch
  char_px  <- 6.5   # width of one label character at text size 4 (generous)
  plot_px  <- 430   # vertical drawing area of plotOutput(height = "500px")

  variants <- variants[order(variants$pos), , drop = FALSE]
  variants$x_disp <- spread_positions(variants$pos, gap = (hi - lo) * 0.035, lo, hi)
  variants$color  <- unname(type_colors[type_category(variants$type)])
  n_circles <- pmax(1L, as.integer(round(ifelse(is.na(variants$score), 1, variants$score))))
  variants$y_label <- y_circle + (n_circles - 1) * step + 13
  variants$short   <- variant_short_label(variants)

  circles <- data.frame(
    x     = rep(variants$x_disp, n_circles),
    y     = y_circle + unlist(lapply(n_circles, function(k) seq_len(k) - 1)) * step,
    color = rep(variants$color, n_circles)
  )
  stems <- data.frame(
    x = variants$pos, y0 = exon_h,
    xb = variants$pos, yb = y_bend,
    xt = variants$x_disp, yt = y_top,
    yc = y_circle
  )

  # Room for the rotated labels; only grows past the drawing area (and
  # so shrinks everything slightly) for very tall stacks + long labels.
  y_max <- max(variants$y_label + char_px * nchar(variants$short) + 5, plot_px)

  ggplot() +
    geom_segment(aes(x = lo, xend = hi, y = 0, yend = 0), linewidth = 0.4) +
    geom_rect(data = exons,
              aes(xmin = exon_chrom_start, xmax = exon_chrom_end,
                  ymin = -exon_h, ymax = exon_h),
              fill = "#8D99AE", colour = "black", linewidth = 0.3) +
    geom_segment(data = stems, aes(x = x, xend = xb, y = y0, yend = yb), linewidth = 0.4) +
    geom_segment(data = stems, aes(x = xb, xend = xt, y = yb, yend = yt), linewidth = 0.4) +
    geom_segment(data = stems, aes(x = xt, xend = xt, y = yt, yend = yc), linewidth = 0.4) +
    geom_point(data = circles, aes(x = x, y = y, fill = I(color)),
               shape = 21, size = 5.5, colour = "black", stroke = 0.5) +
    geom_text(data = variants, aes(x = x_disp, y = y_label, label = short),
              angle = 90, hjust = 0, size = 4) +
    scale_x_continuous(name = paste0("Chromosome ", chrom, " position (GRCh37)"),
                       labels = function(b) format(b, big.mark = ",", scientific = FALSE, trim = TRUE),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(-15, y_max), expand = c(0, 0)) +
    coord_cartesian(xlim = c(lo, hi), clip = "off") +
    theme_void(base_size = 13) +
    theme(axis.line.x  = element_line(linewidth = 0.4),
          axis.ticks.x = element_line(linewidth = 0.4),
          axis.ticks.length.x = unit(4, "pt"),
          axis.text.x  = element_text(margin = margin(t = 3)),
          axis.title.x = element_text(margin = margin(t = 8)),
          plot.margin  = margin(10, 20, 10, 20))
}

# ---------------------------------------------------------------
# Per-gene figures, browsed from a local folder.
# Expected file naming: GENE_SOMETHING.png (or .jpg/.jpeg/.gif),
# e.g. TP53_domains.png, TP53_expression.png — "SOMETHING" becomes
# the dropdown label. Create this folder next to app.R and drop
# figures into it; nothing else needs to change.
# ---------------------------------------------------------------
FIGURES_DIR <- "figures"

# Drop-down names for the biochemistry figures: the figures' y-axis
# titles without units (plasma calcium shortened: "(albumin-corrected)"
# dropped).
# Keys are the file-name labels (GENE_<label>.png). A label missing here
# is shown with underscores replaced by spaces.
FIGURE_LABELS <- c(
  "Calcium_URINE"            = "Calcium, urine",
  "Phosphate_URINE"          = "Phosphate, urine",
  "Uric_acid_URINE"          = "Uric acid, urine",
  "Oxalate_URINE"            = "Oxalate, urine",
  "Magnesium_URINE"          = "Magnesium, urine",
  "Citrate_URINE"            = "Citrate, urine",
  "pH_URINE"                 = "pH, urine",
  "Volume_URINE"             = "Volume, urine",
  "Phosphate..Pi._BLOOD"     = "Phosphate (Pi), plasma",
  "Uric.acid_BLOOD"          = "Uric acid, plasma",
  "Parathormone..PTH._BLOOD" = "Parathyroid hormone (PTH), plasma",
  "X25.OH.Vit.D_BLOOD"       = "25-OH vitamin D, plasma",
  "X1.25..OH.2.Vit.D_BLOOD"  = "1,25-(OH)2 vitamin D, plasma",
  "Hb1Ac..NGSP._BLOOD"       = "HbA1c (NGSP), plasma",
  "FGF.23_BLOOD"             = "FGF-23, plasma",
  "eGFR..CKD.EP._BLOOD"      = "eGFR (CKD-EPI)",
  "Calcium_BLOOD"            = "Calcium, plasma",
  "Urea_BLOOD"               = "Urea, plasma",
  "TmP_GFR"                  = "TmP/GFR"
)

figure_display_name <- function(label) {
  ifelse(label %in% names(FIGURE_LABELS), FIGURE_LABELS[label], gsub("_", " ", label))
}

list_gene_figures <- function(gene_symbol) {
  if (!dir.exists(FIGURES_DIR)) return(character(0))
  files <- list.files(FIGURES_DIR)
  pattern <- paste0("^", toupper(gene_symbol), "_(.+)\\.(png|jpg|jpeg|gif)$")
  matches <- files[grepl(pattern, files, ignore.case = TRUE)]
  if (length(matches) == 0) return(character(0))
  labels <- sub(pattern, "\\1", matches, ignore.case = TRUE)
  paths <- file.path(FIGURES_DIR, matches)
  stats::setNames(paths, labels)
}

# ---------------------------------------------------------------
# UI
# ---------------------------------------------------------------
# Two top-level pages via page_navbar (rather than tabs inside one
# page): "Viewer" carries the gene sidebar; "Documentation" is a
# separate page with no sidebar and no gene selector on it at all.
# Link that opens in a new tab (used in the Documentation page).
ext_link <- function(href, text) {
  tags$a(href = href, target = "_blank", rel = "noopener noreferrer", text)
}

ui <- bslib::page_navbar(
  title = tags$span(
    "Swiss Kidney Stone Cohort (SKSC):", tags$br(),
    "Gene viewer for novel nominated genes and top 1 ranking variants in Exomiser",
    style = "font-size: 1.5em; font-weight: 600; line-height: 1.3; display: inline-block;"
  ),
  theme = bslib::bs_theme(version = 5, bootswatch = "cosmo"),

  bslib::nav_panel(
    title = "Viewer",
    # Two columns: gene selection and biochemistry figure on the left
    # (1/3), lollipop plot and variant table on the right (2/3). On
    # narrow screens the columns stack, left column first.
    bslib::layout_columns(
      col_widths = c(4, 8),

      # Left column
      div(
        # 1. Gene name / selection box.
        bslib::card(
          bslib::card_body(
            # Commented out along with EXOMISER_VERSION_FILES above — only one
            # version in use right now. Uncomment both to bring back the picker.
            # selectInput("exomiser_version", "Exomiser version",
            #             choices = names(EXOMISER_VERSION_FILES)),
            selectInput("gene_symbol", "Gene",
                        choices = if (length(GENE_CHOICES) > 0) GENE_CHOICES else character(0)),
            if (length(GENE_CHOICES) == 0) {
              helpText("No genes found in transcripts.rds — run prepare_transcripts.R first.")
            },
            hr(),
            uiOutput("transcript_info")
          )
        ),

        # 2. Biochemistry figure, below the gene box. The image scales to
        # the column width (height follows the figure's aspect ratio).
        bslib::card(
          bslib::card_header("Biochemistry figure"),
          bslib::card_body(
            selectInput("figure_select", NULL, choices = NULL),
            imageOutput("gene_figure", height = "auto")
          )
        )
      ),

      # Right column: gene viewer — plot, legend, variant table.
      bslib::card(
        bslib::card_body(
          # fill = FALSE: keep the full 500px instead of letting the bslib
          # card squash the plot on short screens (lollipop_plot() lays out
          # circles and labels for 500px).
          plotOutput("lolliplot", height = "500px", fill = FALSE),
          uiOutput("color_legend"),
          h4("Variant details"),
          helpText("Hover over a row to see carrier information. ",
                    "Click a column header to sort."),
          DT::dataTableOutput("variant_table")
        )
      )
    )
  ),

  bslib::nav_panel(
    title = "Documentation",
    br(),
    h4("About this app"),
    p("This app displays relevant genetic variants from the ",
      ext_link("https://doi.org/10.1159/000529094", "Swiss Kidney Stone Cohort (SKSC)"),
      " identified by Exomiser v.14.0.0, in genes that were ", strong("not"),
      " targeted by a genetic panel in ",
      ext_link("https://doi.org/10.1172/jci196277", "Münch et al., 2026"), ". "),

    h4("How to use the viewer"),

    h5("Gene"),
    p("Please select a gene from the drop-down menu. It includes 62 newly nominated genes, to be ",
      "further investigated for involvement in (monogenic) kidney stone disease (KSD). The list was ",
      "created based on automated variant prioritization with the tool Exomiser, benchmarked and ",
      "optimized based on whole exome sequencing (WES) data from the Swiss Kidney Stone Cohort (SKSC)."),

    h5("Lollipop plot"),
    p("The lollipop plot provides localization of the variants prioritized by Exomiser (Rank 1), ",
      "based on the optimized pipeline, for the gene selected from the drop-down menu. Each circle ",
      "represents a participant."),

    h5("Variant details"),
    p("The variant details are taken directly from Exomiser’s output. The ACMG column is produced ",
      "by an automated ACMG assignment, with P = pathogenic; LP = likely pathogenic; US = uncertain ",
      "significance."),
    p("The max_freq_in column reports the population with highest minor allele frequency (MAF), and ",
      "the max_freq the respective MAF from that population. This is the MAF Exomiser considers for ",
      "scoring and ranking."),
    p("The max_path_in column reports the in-silico pathogenicity prediction algorithm that predicted ",
      "the highest pathogenicity for the variant. max_path reports the respective pathogenicity score ",
      "for that algorithm. This is the pathogenicity Exomiser considers for scoring and ranking."),
    p("Please see ",
      ext_link("https://exomiser.readthedocs.io/en/stable/", "Exomiser’s documentation"),
      " for more details."),

    h5("Biochemistry figure"),
    p("Using the drop-down menu, plasma and urine biochemical measurements which are relevant for KSD ",
      "can be selected. The respective measurement (if available) will be depicted for variant ",
      "carriers depicted in the lollipop plot on the right, in orange. To provide context, non-carriers are ",
      "depicted in gray, in two groups: individuals with KSD (left box) and controls (right box). ",
      "Only individuals without LP/P variants in KSD established genes were included in this gray ",
      "“background” group."),
    p("The middle line in each boxplot represents the median of the non-carriers in that group. The box edges are the ",
      "25th and 75th percentiles, which span the interquartile range (IQR). The whiskers extend to ",
      "the furthest point within 1.5 × IQR of the box."),

    h4("Citation"),
    p("For more information, please see < our paper>.
      We are happy if the provided information helps your work. You are free to use it, just please cite us. "),

    h4("Genome assembly"),
    p("Variant positions are given in ",
      strong("GRCh37/hg19"), " genomic coordinates. "),

    h4("Developer"),
    p("This app was developed by ",
      ext_link("https://www.linkedin.com/in/jana-petrovska", "Jana Petrovska, PhD"), ".")
  )
)

# ---------------------------------------------------------------
# Server
# ---------------------------------------------------------------
server <- function(input, output, session) {

  gene_data <- reactive({
    req(input$gene_symbol)
    symbol <- toupper(input$gene_symbol)
    cached <- transcripts_cache[[symbol]]
    validate(need(!is.null(cached),
                   paste0("No cached transcript for ", symbol,
                          " — add it via prepare_transcripts.R.")))
    cached
  })

  variants_source <- reactive({
    # path <- EXOMISER_VERSION_FILES[[input$exomiser_version]]  # when re-enabled
    path <- DEFAULT_VERSION_FILE
    load_variants_file(path)
  })

  variant_data <- reactive({
    req(gene_data())
    all_variants <- get_variants(gene_data()$symbol, variants_source())
    # Defensive schema check: fill in any missing expected column rather
    # than letting a downstream data.frame() call crash with something
    # like "arguments imply differing number of rows" (e.g. if a
    # variant file was generated by hand, or predates a column being
    # added). Recycled NA fills every row.
    expected_cols <- list(pos = NA_real_, label = NA_character_,
                           type = NA_character_, score = NA_real_,
                           carriers = NA_character_)
    for (col in names(expected_cols)) {
      if (!(col %in% names(all_variants))) {
        all_variants[[col]] <- expected_cols[[col]]
      }
    }
    all_variants
  })

  output$color_legend <- renderUI({
    variants <- variant_data()
    if (nrow(variants) == 0) return(NULL)
    # Categories present for this gene, in fixed order: missense, lof, other.
    used <- names(type_colors)[names(type_colors) %in% type_category(variants$type)]

    swatches <- lapply(used, function(t) {
      col <- unname(type_colors[t])
      tags$span(
        style = "display:inline-flex; align-items:center; margin-right:18px;",
        tags$span(style = sprintf(
          "display:inline-block; width:12px; height:12px; background-color:%s; border-radius:2px; margin-right:6px;",
          col)),
        t
      )
    })

    tags$div(style = "margin: 4px 0 14px 0;", swatches)
  })

  figure_choices <- reactive({
    req(gene_data())
    list_gene_figures(gene_data()$symbol)
  })

  observeEvent(figure_choices(), {
    choices <- figure_choices()
    # Shown with the figure's y-axis title ("Calcium, urine"), sorted
    # alphabetically; the selected value stays the original file label
    # ("Calcium_URINE") used to look up the file.
    if (length(choices) == 0) {
      updateSelectInput(session, "figure_select", choices = character(0))
    } else {
      labels  <- names(choices)
      display <- unname(figure_display_name(labels))
      ord     <- order(tolower(display))
      updateSelectInput(session, "figure_select",
                         choices = stats::setNames(labels[ord], display[ord]))
    }
  })

  output$gene_figure <- renderImage({
    choices <- figure_choices()
    validate(need(length(choices) > 0,
                   paste0("No figures found in '", FIGURES_DIR, "/' for ",
                          gene_data()$symbol, ".")))
    req(input$figure_select)
    path <- choices[[input$figure_select]]
    validate(need(!is.null(path) && file.exists(path), "Selected figure file not found."))

    list(src = path,
         contentType = paste0("image/", tolower(tools::file_ext(path))),
         width = "100%",
         style = "width: 100%; height: auto;",
         alt = paste(gene_data()$symbol, figure_display_name(input$figure_select)))
  }, deleteFile = FALSE)

  output$transcript_info <- renderUI({
    gd <- tryCatch(gene_data(), error = function(e) NULL)
    if (is.null(gd)) return(NULL)
    genecards_url <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", gd$symbol)
    omim_url <- paste0("https://www.omim.org/search/?search=", gd$symbol)
    tagList(
      p(strong("Gene: "), gd$symbol),
      p(strong("Transcript: "), gd$transcript_id),
      p(strong("Exons: "), nrow(gd$exons)),
      p(tags$a(href = genecards_url, target = "_blank", rel = "noopener noreferrer",
               "View on GeneCards ", HTML("&#8599;"))),
      p(tags$a(href = omim_url, target = "_blank", rel = "noopener noreferrer",
               "View on OMIM ", HTML("&#8599;")))
    )
  })

  output$lolliplot <- renderPlot({
    gd <- gene_data()
    variants <- variant_data()
    validate(need(nrow(variants) > 0,
                   paste0("No example variants bundled for ", gd$symbol,
                          ". Try TP53 or BRCA1, or add data in get_variants().")))
    lollipop_plot(gd$exons, variants, type_colors)
  })

  output$variant_table <- DT::renderDataTable({
    variants <- variant_data()
    validate(need(nrow(variants) > 0, "No variants to list for this gene."))
    variants <- variants[order(variants$pos), ]

    df <- data.frame(
      Position = variants$pos,
      #Label    = variants$label,
      #Score    = variants$score,
      Ref = variants$ref, 
      Alt = variants$alt, 
      Type     = variants$type,
      c.pos = variants$c.position, 
      p.pos = variants$p.position, 
      rsid = variants$rs_id, 
      max_freq_in = variants$max_freq_source, 
      max_freq = variants$max_freq, 
      max_path_in = variants$max_path_source, 
      max_path = variants$max_path, 
      ACMG = variants$acmg,
      Carriers = variants$carriers,   # kept but hidden — source for tooltips
      stringsAsFactors = FALSE
    )

    tooltip_render <- DT::JS(
      "function(data, type, row) {",
      "  if (type !== 'display') return data;",
      "  var tip = String(row[12]).replace(/\"/g, '&quot;');",
      "  return '<span data-bs-toggle=\"tooltip\" data-bs-placement=\"top\" title=\"' +",
      "         tip + '\">' + data + '</span>';",
      "}"
    )

    DT::datatable(
      df,
      rownames = FALSE,
      selection = "none",
      options = list(
        pageLength = 10,
        columnDefs = list(
          list(targets = 12, visible = FALSE),      # hide raw Carriers column
          list(targets = 0:11, render = tooltip_render)
        ),
        # Re-attach Bootstrap tooltips every time DT (re)draws the table
        # body — this fires on initial load, sorting, and pagination.
        drawCallback = DT::JS(
          "function(settings) {",
          "  var api = this.api();",
          "  $(api.table().node()).find('[data-bs-toggle=\"tooltip\"]').each(function() {",
          "    if (!this._tooltipInstance) {",
          "      this._tooltipInstance = new bootstrap.Tooltip(this);",
          "    }",
          "  });",
          "}"
        )
      )
    )
  })
}

shinyApp(ui, server)
