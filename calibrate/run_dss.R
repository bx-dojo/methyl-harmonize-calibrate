#!/usr/bin/env Rscript
# Count DSS calls for every label assignment: the observed labels and each permuted set
# from make_assignments.R. Everything except the labels stays fixed, so the permuted
# counts show how many calls the same pipeline returns by chance.
#
# Steps:
#   1. per sample, drop CpG sites with fewer than --min-coverage reads or above the
#      --max-percentile of the sample's coverage (the same filter as run_methylkit.R)
#   2. combine the samples. A site is kept if it passed in at least one sample, and a
#      sample that failed there enters with no reads, which DSS treats as missing data.
#      --all-samples keeps only the sites that passed in every sample, as methylKit does.
#   3. test with DMLtest (two groups, smoothing off unless --smoothing)
#   4. count the CpG sites with FDR < each --q and absolute methylation difference >=
#      each --diff (percentage points)
#   5. with --level dmr, also call regions with callDMR at each --dmr-p, using the
#      region rules given below, and count regions whose absolute mean difference is
#      >= each --diff. callDMR applies no region-level FDR, which is why a permutation
#      null is the only way to judge its region count.
#
# Usage:
#   Rscript run_dss.R --samples samples.tsv --assignments assignments.tsv --out-dir DIR \
#       [--level dml] [--min-coverage 10] [--max-percentile 99.9] [--all-samples] \
#       [--q 0.05,0.01,0.001] [--diff 0,5,10,15,20,25,30,40,50] [--smoothing] \
#       [--dmr-p 0.01] [--dmr-minlen 50] [--dmr-mincg 3] [--dmr-dis-merge 100] \
#       [--dmr-pct-sig 0.5] [--ids observed,perm_001,...] [--ncores 1]
#
# samples.tsv: sample, group, file. Outputs as for run_methylkit.R. DMLtest is the slow,
# memory-hungry step; --ncores sets its worker processes.
suppressPackageStartupMessages({
  library(parallel)    # some DSS releases call detectCores() unqualified
  library(DSS)
  library(data.table)
})
script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))))
source(file.path(script_dir, "read_calls.R"))

opt <- parse_opts(commandArgs(trailingOnly = TRUE), list(
  samples = NULL, assignments = NULL, out_dir = NULL, level = "dml",
  min_coverage = "10", max_percentile = "99.9", all_samples = FALSE,
  q = "0.05,0.01,0.001", diff = "0,5,10,15,20,25,30,40,50", smoothing = FALSE,
  dmr_p = "0.01", dmr_minlen = "50", dmr_mincg = "3", dmr_dis_merge = "100",
  dmr_pct_sig = "0.5", ids = NULL, ncores = "1"))
if (is.null(opt$samples) || is.null(opt$assignments) || is.null(opt$out_dir)) {
  stop("usage: run_dss.R --samples F --assignments F --out-dir D [options]")
}
if (!opt$level %in% c("dml", "dmr")) stop("--level must be dml or dmr")
q_cut <- num_list(opt$q); d_cut <- num_list(opt$diff); p_dmr <- num_list(opt$dmr_p)
design <- read_design(opt$samples, opt$assignments, opt$ids)
samples <- design$samples; assignments <- design$assignments
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

todo <- assignments[!file.exists(checkpoint_path(opt$out_dir, assignments$assignment_id)), ]
if (!nrow(todo)) {
  message("all assignments are done")
  combine_checkpoints(opt$out_dir)
  quit(save = "no")
}

# ---- build the BSseq object once; only the group membership changes afterwards
calls <- lapply(samples$file, function(f) {
  x <- coverage_filter(read_calls(f), as.numeric(opt$min_coverage),
                       as.numeric(opt$max_percentile))
  data.frame(chr = x$chr, pos = x$pos, N = x$coverage, X = x$methylated)
})
bs <- makeBSseqData(calls, samples$sample)
if (opt$all_samples) {
  keep <- rowSums(as.matrix(bsseq::getCoverage(bs, type = "Cov")) > 0) == ncol(bs)
  bs <- bs[keep, ]
}
message("CpG sites: ", nrow(bs))

ncores <- as.integer(opt$ncores)
for (i in seq_len(nrow(todo))) {
  a <- todo[i, ]
  g1 <- strsplit(a$group_a, ",")[[1]]
  g2 <- strsplit(a$group_b, ",")[[1]]
  # DMLtest draws a progress bar on standard output; keep the log readable
  invisible(utils::capture.output(dml <- suppressMessages(
    DMLtest(bs, group1 = g1, group2 = g2, smoothing = opt$smoothing, ncores = ncores))))
  rows <- list()
  for (qc in q_cut) for (dc in d_cut) {
    rows[[length(rows) + 1]] <- data.table(
      unit = "sites", sig_cutoff = qc, diff_cutoff = dc,
      calls = sum(dml$fdr < qc & abs(dml$diff) >= dc / 100, na.rm = TRUE))
  }
  if (opt$level == "dmr") {
    for (p in p_dmr) {
      regions <- callDMR(dml, delta = 0, p.threshold = p,
                         minlen = as.numeric(opt$dmr_minlen), minCG = as.numeric(opt$dmr_mincg),
                         dis.merge = as.numeric(opt$dmr_dis_merge),
                         pct.sig = as.numeric(opt$dmr_pct_sig))
      d <- if (is.null(regions)) numeric(0) else regions$diff.Methy
      for (dc in d_cut) {
        rows[[length(rows) + 1]] <- data.table(
          unit = "regions", sig_cutoff = p, diff_cutoff = dc,
          calls = sum(abs(d) >= dc / 100, na.rm = TRUE))
      }
    }
  }
  out <- cbind(data.table(assignment_id = a$assignment_id, observed = a$observed,
                          caller = if (opt$smoothing) "DSS (smoothed)" else "DSS",
                          level = opt$level, tested = nrow(dml)),
               rbindlist(rows, use.names = TRUE))
  write_checkpoint(out, checkpoint_path(opt$out_dir, a$assignment_id))
  message("done ", a$assignment_id, " (", i, "/", nrow(todo), ")")
  rm(dml); invisible(gc())
}
combine_checkpoints(opt$out_dir)
