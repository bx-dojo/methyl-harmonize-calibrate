#!/usr/bin/env Rscript
# CpG coverage counted directly from deduplicated BAM files, with the same samtools depth
# procedure for every aligner.
#
# Methylation callers apply their own read filters, so the CpG sites they report mix the
# coverage of the alignments with the filtering of the caller. Counting read depth at
# every reference CpG site from the BAM files separates the two. For each aligner the
# script reports, at each minimum depth:
#   union         CpG sites reaching the depth in at least one sample
#   intersection  CpG sites reaching it in every sample
# with their share of all reference CpG sites and the mean depth at those sites.
#
# samtools depth excludes unmapped, secondary, QC-failed, and duplicate records by
# default. Supplementary records are also excluded (-G 0x800), the overlapping part of
# the two reads of a pair is counted once (-s), and zero-depth sites are kept (-aa).
# Run at --mapq 0 and --mapq 1 to see how much coverage comes from multi-mapping reads,
# which BWA-MEM based aligners write with MAPQ 0.
#
# Usage:
#   Rscript bam_cpg_coverage.R --bams BAM_SHEET.tsv --cpg-bed reference_cpg.bed \
#       --out-dir DIR [--mapq 1] [--baseq 0] [--thresholds 1,5,10] [--threads 4]
#
# BAM_SHEET.tsv: aligner, sample, bam   (deduplicated BAMs; the same samples for every
# aligner). reference_cpg.bed comes from reference_cpg_bed.sh. BAMs must be coordinate
# sorted and indexed; unsorted ones are sorted into DIR/sorted_bams/ first.
#
# Writes bam_cpg_coverage_per_sample.tsv and bam_cpg_coverage_summary.tsv.

options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
opt <- list(bams = NULL, cpg_bed = NULL, out_dir = NULL, mapq = "1", baseq = "0",
            thresholds = "1,5,10", threads = "4")
for (i in seq(1, length(args), by = 2)) {
  key <- gsub("-", "_", sub("^--", "", args[[i]]))
  if (!key %in% names(opt) || i == length(args)) stop("bad argument: ", args[[i]])
  opt[[key]] <- args[[i + 1]]
}
if (is.null(opt$bams) || is.null(opt$cpg_bed) || is.null(opt$out_dir)) {
  stop("usage: bam_cpg_coverage.R --bams F --cpg-bed F --out-dir D [options]")
}
mapq <- as.integer(opt$mapq); baseq <- as.integer(opt$baseq)
threads <- as.integer(opt$threads)
thresholds <- sort(unique(as.integer(strsplit(opt$thresholds, ",")[[1]])))
if (anyNA(c(mapq, baseq, threads, thresholds)) || any(thresholds < 1)) stop("invalid number")
if (!nzchar(Sys.which("samtools"))) stop("samtools not found in PATH")
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

sheet <- read.delim(opt$bams)
samples <- unique(sheet$sample)
for (a in unique(sheet$aligner)) {
  if (!setequal(sheet$sample[sheet$aligner == a], samples)) {
    stop("aligner ", a, " does not list the same samples as the others")
  }
}
missing <- sheet$bam[!file.exists(sheet$bam)]
if (length(missing)) stop("missing BAM(s):\n", paste(missing, collapse = "\n"))

sq <- function(x) shQuote(x, type = "sh")
run <- function(cmd) {
  status <- system2("bash", c("-o", "pipefail", "-c", sq(cmd)))
  if (!identical(status, 0L)) stop("command failed:\n", cmd)
}

cpg_bed <- normalizePath(opt$cpg_bed, mustWork = TRUE)
n_ref <- as.numeric(system2("awk", c(sq("END{print NR}"), sq(cpg_bed)), stdout = TRUE))
if (!is.finite(n_ref) || n_ref < 1) stop("empty CpG BED: ", cpg_bed)

# samtools depth with -b needs coordinate-sorted, indexed BAMs.
is_sorted <- function(bam) {
  hd <- system2("samtools", c("view", "-H", sq(bam)), stdout = TRUE)
  any(grepl("^@HD.*\tSO:coordinate", hd))
}
has_index <- function(bam) {
  any(file.exists(paste0(bam, c(".bai", ".csi")), sub("\\.bam$", ".bai", bam)))
}
for (i in seq_len(nrow(sheet))) {
  bam <- sheet$bam[i]
  if (!is_sorted(bam)) {
    sorted <- file.path(opt$out_dir, "sorted_bams",
                        paste0(sheet$aligner[i], "_", sheet$sample[i], ".bam"))
    dir.create(dirname(sorted), recursive = TRUE, showWarnings = FALSE)
    if (!file.exists(sorted)) {
      message("sorting ", bam)
      run(paste("samtools sort -@", threads, "-o", sq(sorted), sq(bam)))
    }
    bam <- sorted
  }
  if (!has_index(bam)) run(paste("samtools index -@", threads, sq(bam)))
  sheet$bam[i] <- bam
}

# One pass per aligner. awk tallies, per threshold, each sample's covered sites and the
# union and intersection across samples, without keeping per-site data.
awk_prog <- '
BEGIN { OFS = "\t"; nt = split(tlist, t, ",") }
{
  sites++
  for (j = 1; j <= ns; j++) { d[j] = $(j + 2); raw[j] += d[j] }
  for (i = 1; i <= nt; i++) {
    n = 0; s = 0
    for (j = 1; j <= ns; j++) if (d[j] >= t[i]) { n++; s += d[j]; cnt[i, j]++; dep[i, j] += d[j] }
    if (n > 0)  { un[i]++; uns[i] += s }
    if (n == ns) { inter[i]++; ints[i] += s }
  }
}
END {
  print "N", sites
  for (i = 1; i <= nt; i++) {
    for (j = 1; j <= ns; j++) print "S", t[i], j, cnt[i, j] + 0, dep[i, j] + 0
    print "A", t[i], un[i] + 0, inter[i] + 0, uns[i] + 0, ints[i] + 0
  }
}'
awk_file <- file.path(opt$out_dir, ".depth.awk")
writeLines(awk_prog, awk_file)
on.exit(unlink(awk_file), add = TRUE)

per_sample <- list(); summary <- list()
for (a in unique(sheet$aligner)) {
  bams <- sheet$bam[sheet$aligner == a][match(samples, sheet$sample[sheet$aligner == a])]
  raw <- file.path(opt$out_dir, ".depth.tsv")
  message("counting depth: ", a)
  run(paste("samtools depth -aa -s -G 0x800 -q", baseq, "-Q", mapq, "-@", threads,
            "-b", sq(cpg_bed), paste(sq(bams), collapse = " "),
            "| awk", paste0("-v tlist=", paste(thresholds, collapse = ",")),
            paste0("-v ns=", length(bams)),
            "-f", sq(awk_file), ">", sq(raw)))
  x <- read.delim(raw, header = FALSE, fill = TRUE, col.names = paste0("V", 1:6),
                  colClasses = c("character", rep("numeric", 5)))
  unlink(raw)
  if (x$V2[x$V1 == "N"] != n_ref) {
    stop("samtools depth returned ", x$V2[x$V1 == "N"], " sites, expected ", n_ref,
         ". Check that the BED and BAM sequence names match.")
  }
  for (th in thresholds) {
    for (j in seq_along(samples)) {
      z <- x[x$V1 == "S" & x$V2 == th & x$V3 == j, ]
      per_sample[[length(per_sample) + 1]] <- data.frame(
        aligner = a, sample = samples[j], min_depth = th, mapq = mapq,
        covered_sites = z$V4, covered_pct = 100 * z$V4 / n_ref,
        mean_depth = if (z$V4 > 0) z$V5 / z$V4 else NA)
    }
    z <- x[x$V1 == "A" & x$V2 == th, ]
    summary[[length(summary) + 1]] <- data.frame(
      aligner = a, min_depth = th, mapq = mapq, reference_cpg_sites = n_ref,
      union_sites = z$V3, union_pct = 100 * z$V3 / n_ref,
      intersection_sites = z$V4, intersection_pct = 100 * z$V4 / n_ref,
      union_mean_depth = z$V5 / (z$V3 * length(samples)),
      intersection_mean_depth = z$V6 / (z$V4 * length(samples)))
  }
}
w <- function(d, f) write.table(d, file.path(opt$out_dir, f), sep = "\t", quote = FALSE,
                                row.names = FALSE, na = "")
w(do.call(rbind, per_sample), "bam_cpg_coverage_per_sample.tsv")
w(do.call(rbind, summary), "bam_cpg_coverage_summary.tsv")
message("wrote ", opt$out_dir, "/bam_cpg_coverage_{per_sample,summary}.tsv")
