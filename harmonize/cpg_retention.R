#!/usr/bin/env Rscript
# Share of the CpG sites covered in the BAM files that the methylation caller kept.
#
# retention = CpG sites in the methylation-call union / CpG sites in the BAM union
#
# both at the same minimum depth, and the BAM union counted at MAPQ 0, because that is
# the BAM the methylation caller reads. Retention therefore measures the filtering by
# the methylation caller, not read loss during alignment.
#
# Usage:
#   Rscript cpg_retention.R --calls CALL_SHEET.tsv --bam-summary bam_cpg_coverage_summary.tsv \
#       --out retention.tsv [--min-depth 1]
#
# CALL_SHEET.tsv: aligner, sample, file. Each file holds one workflow's CpG calls for one
# sample, tab-separated, with a header and at least the columns chr, pos, and coverage,
# after the two strands of each CpG have been combined on the forward-strand C. Files may
# be gzipped. bam_cpg_coverage_summary.tsv comes from bam_cpg_coverage.R run with
# --mapq 0, and its aligner names must match the call sheet.
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
opt <- list(calls = NULL, bam_summary = NULL, out = NULL, min_depth = "1")
for (i in seq(1, length(args), by = 2)) {
  key <- gsub("-", "_", sub("^--", "", args[[i]]))
  if (!key %in% names(opt) || i == length(args)) stop("bad argument: ", args[[i]])
  opt[[key]] <- args[[i + 1]]
}
if (is.null(opt$calls) || is.null(opt$bam_summary) || is.null(opt$out)) {
  stop("usage: cpg_retention.R --calls F --bam-summary F --out F [--min-depth N]")
}
depth_cut <- as.integer(opt$min_depth)

sheet <- fread(opt$calls)
bam <- fread(opt$bam_summary)
bam <- bam[mapq == 0 & min_depth == depth_cut]
if (!nrow(bam)) stop("the BAM summary has no MAPQ 0 rows at minimum depth ", depth_cut)

out <- rbindlist(lapply(unique(sheet$aligner), function(a) {
  files <- sheet[aligner == a, file]
  sites <- lapply(files, function(f) {
    x <- fread(f, select = c("chr", "pos", "coverage"))
    unique(x[coverage >= depth_cut, .(chr, pos)])
  })
  union <- unique(rbindlist(sites))
  inter <- Reduce(function(x, y) fintersect(x, y), sites)
  b <- bam[aligner == a]
  if (nrow(b) != 1) stop("no BAM summary row for aligner ", a)
  data.table(aligner = a, min_depth = depth_cut, samples = length(files),
             call_union_sites = nrow(union), call_intersection_sites = nrow(inter),
             bam_union_sites_mapq0 = b$union_sites,
             retention_pct = 100 * nrow(union) / b$union_sites)
}))
fwrite(out, opt$out, sep = "\t")
message("wrote ", opt$out)
