#!/usr/bin/env Rscript
# Count methylKit calls for every label assignment: the observed labels and each
# permuted set from make_assignments.R. Everything except the labels stays fixed, so the
# permuted counts show how many calls the same pipeline returns by chance.
#
# Steps, as in a standard methylKit analysis:
#   1. per sample, drop CpG sites with fewer than --min-coverage reads or above the
#      --max-percentile of the sample's coverage (filterByCoverage)
#   2. --level dml: keep the sites that pass in every sample (unite)
#      --level dmr: count reads in adjacent, non-overlapping windows of --window bp
#      (tileMethylCounts), then unite
#   3. test with calculateDiffMeth at its defaults (logistic regression), or with
#      --overdispersion MN, its correction for variation between samples, which scales
#      each site's test by an estimated overdispersion factor and uses the F test
#   4. count the sites or windows with q-value < each --q and absolute methylation
#      difference >= each --diff (percentage points). For --level dmr, runs of adjacent
#      significant windows that change in the same direction are also merged into
#      regions and counted.
#
# Usage:
#   Rscript run_methylkit.R --samples samples.tsv --assignments assignments.tsv \
#       --out-dir DIR [--level dml] [--min-coverage 10] [--max-percentile 99.9] \
#       [--q 0.05,0.01,0.001] [--diff 0,5,10,15,20,25,30,40,50] [--window 100] \
#       [--overdispersion none] [--ids observed,perm_001,...] [--cores 1]
#
# samples.tsv: sample, group, file (see read_calls.R for the call file format). One
# checkpoint per assignment is written to DIR/counts/, finished ones are skipped, and
# all are combined into DIR/counts.tsv for summarize_calibration.R. --ids runs a subset,
# so a long calibration can be split across jobs into the same DIR.
suppressPackageStartupMessages({
  library(methylKit)
  library(data.table)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))))
source(file.path(script_dir, "read_calls.R"))

opt <- parse_opts(commandArgs(trailingOnly = TRUE), list(
  samples = NULL, assignments = NULL, out_dir = NULL, level = "dml",
  min_coverage = "10", max_percentile = "99.9", q = "0.05,0.01,0.001",
  diff = "0,5,10,15,20,25,30,40,50", window = "100", overdispersion = "none",
  ids = NULL, cores = "1"))
if (is.null(opt$samples) || is.null(opt$assignments) || is.null(opt$out_dir)) {
  stop("usage: run_methylkit.R --samples F --assignments F --out-dir D [options]")
}
if (!opt$level %in% c("dml", "dmr")) stop("--level must be dml or dmr")
if (!opt$overdispersion %in% c("none", "MN")) stop("--overdispersion must be none or MN")
q_cut <- num_list(opt$q); d_cut <- num_list(opt$diff)
design <- read_design(opt$samples, opt$assignments, opt$ids)
samples <- design$samples; assignments <- design$assignments
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

todo <- assignments[!file.exists(checkpoint_path(opt$out_dir, assignments$assignment_id)), ]
if (!nrow(todo)) {
  message("all assignments are done")
  combine_checkpoints(opt$out_dir)
  quit(save = "no")
}

# ---- build the united object once; only the treatment vector changes afterwards
observed <- as.integer(samples$group == unique(samples$group)[2])
raw <- lapply(seq_len(nrow(samples)), function(i) {
  x <- read_calls(samples$file[i])
  new("methylRaw",
      data.frame(chr = x$chr, start = x$pos, end = x$pos, strand = "+",
                 coverage = x$coverage, numCs = x$methylated,
                 numTs = x$coverage - x$methylated),
      sample.id = samples$sample[i], assembly = "reference", context = "CpG",
      resolution = "base")
})
obj <- methylRawList(raw, treatment = observed)
obj <- filterByCoverage(obj, lo.count = as.numeric(opt$min_coverage), lo.perc = NULL,
                        hi.count = NULL, hi.perc = as.numeric(opt$max_percentile))
if (opt$level == "dmr") {
  w <- as.integer(opt$window)
  obj <- tileMethylCounts(obj, win.size = w, step.size = w)
}
united <- unite(obj, destrand = FALSE)
unit <- if (opt$level == "dml") "sites" else "windows"
caller <- if (opt$overdispersion == "MN") "methylKit (MN)" else "methylKit"
message("tested ", unit, ": ", nrow(united))

run_one <- function(i) {
  a <- todo[i, ]
  group_b <- strsplit(a$group_b, ",")[[1]]
  m <- united
  m@treatment <- as.integer(m@sample.ids %in% group_b)
  res <- getData(if (opt$overdispersion == "MN") {
    calculateDiffMeth(m, overdispersion = "MN", test = "F")
  } else {
    calculateDiffMeth(m)
  })
  valid <- is.finite(res$qvalue) & is.finite(res$meth.diff)
  rows <- list()
  for (qc in q_cut) for (dc in d_cut) {
    pass <- valid & res$qvalue < qc & abs(res$meth.diff) >= dc
    rows[[length(rows) + 1]] <- data.table(unit = unit, calls = sum(pass),
                                           sig_cutoff = qc, diff_cutoff = dc)
    if (opt$level == "dmr") {
      rows[[length(rows) + 1]] <- data.table(
        unit = "regions", sig_cutoff = qc, diff_cutoff = dc,
        calls = merge_windows(res$chr[pass], res$start[pass], res$end[pass], res$meth.diff[pass]))
    }
  }
  d <- rbindlist(rows, use.names = TRUE)
  d <- cbind(data.table(assignment_id = a$assignment_id, observed = a$observed,
                        caller = caller, level = opt$level, tested = nrow(res)), d)
  write_checkpoint(d, checkpoint_path(opt$out_dir, a$assignment_id))
  message("done ", a$assignment_id)
  invisible(TRUE)
}
cores <- as.integer(opt$cores)
if (cores > 1) {
  res <- parallel::mclapply(seq_len(nrow(todo)), run_one, mc.cores = cores,
                            mc.preschedule = FALSE)
  if (any(vapply(res, inherits, logical(1), "try-error"))) stop("an assignment failed")
} else {
  invisible(lapply(seq_len(nrow(todo)), run_one))
}
combine_checkpoints(opt$out_dir)
