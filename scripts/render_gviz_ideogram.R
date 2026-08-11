#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Gviz))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("usage: render_gviz_ideogram.R IDEOGRAM_PNG")
}

ideogram_png <- args[[1]]
dir.create(dirname(ideogram_png), recursive = TRUE, showWarnings = FALSE)

chromosome <- "chr3"
genome <- "hg38"
locus_start <- 4803502
locus_end <- 5144387

# Gviz's IdeogramTrack marks the `from`/`to` interval on the chromosome
# ideogram. Supplying the exact paper coordinates therefore produces the
# requested cytoband ideogram with the locus boxed in red.
ideogram <- IdeogramTrack(genome = genome, chromosome = chromosome, name = "chr3")

png(ideogram_png, width = 2200, height = 430, res = 220, bg = "white")
plotTracks(
  ideogram,
  chromosome = chromosome,
  from = locus_start,
  to = locus_end,
  showTitle = FALSE,
  main = "hg38 chromosome 3 - Figure 3D locus boxed in red"
)
grid::grid.text(
  "chr3:4,803,502-5,144,387 (340,886 bp)",
  x = grid::unit(0.5, "npc"), y = grid::unit(0.055, "npc"),
  gp = grid::gpar(col = "#b00000", fontsize = 11, fontface = "bold")
)
dev.off()
