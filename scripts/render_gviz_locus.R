#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Gviz))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("usage: render_gviz_locus.R LOCUS_PNG")
}

locus_png <- args[[1]]
dir.create(dirname(locus_png), recursive = TRUE, showWarnings = FALSE)

chromosome <- "chr3"
genome <- "hg38"
locus_start <- 4803502
locus_end <- 5144387
deletion_start <- 4976068
deletion_end <- 4976790
subtad_start <- 4977047
subtad_end <- 5056047

axis_track <- GenomeAxisTrack(name = "hg38", add53 = TRUE)

gene_track <- AnnotationTrack(
  start = c(4493345, 4896806, 4978715, 5122231),
  end = c(4847518, 4980414, 4985323, 5180916),
  chromosome = chromosome,
  genome = genome,
  strand = c("+", "-", "+", "+"),
  id = c("ITPR1", "BHLHE40-AS1", "BHLHE40", "ARL8B"),
  group = c("ITPR1", "BHLHE40-AS1", "BHLHE40", "ARL8B"),
  name = "Genes",
  fill = c("#3182bd", "#7a3db8", "#3182bd", "#3182bd"),
  col = "#333333",
  fontcolor = "#202020",
  stacking = "squish",
  showFeatureId = TRUE,
  featureAnnotation = "id"
)

subtad_track <- AnnotationTrack(
  start = subtad_start,
  end = subtad_end,
  chromosome = chromosome,
  genome = genome,
  id = "K562 active sub-TAD",
  group = "K562 active sub-TAD",
  name = "Published sub-TAD",
  fill = "#3182bd55",
  col = "#1565c0",
  lwd = 2
)

edit_track <- AnnotationTrack(
  start = c(4976065, deletion_start, 4976788),
  end = c(4976084, deletion_end, 4976807),
  chromosome = chromosome,
  genome = genome,
  id = c("sgRNA 1", "723-bp cut-to-cut deletion", "sgRNA 2"),
  group = c("sgRNA 1", "deletion", "sgRNA 2"),
  name = "CRISPR audit",
  fill = c("#fdae6b", "#d00000", "#fdae6b"),
  col = c("#e6550d", "#8b0000", "#e6550d"),
  lwd = 2,
  stacking = "dense"
)

png(locus_png, width = 2200, height = 900, res = 220, bg = "white")
plotTracks(
  list(axis_track, gene_track, subtad_track, edit_track),
  chromosome = chromosome,
  from = locus_start,
  to = locus_end,
  sizes = c(0.55, 1.25, 0.7, 0.85),
  rotation.title = 0,
  cex.title = 0.7,
  cex.axis = 0.8,
  background.title = "#f2f2f2",
  col.title = "#202020",
  main = "Gviz annotation of the exact NAR Figure 3D interval"
)
dev.off()
