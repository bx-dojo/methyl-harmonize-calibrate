#!/usr/bin/env Rscript
# Compare each observed call count with its permutation null.
#
# For every caller, level, unit, and threshold pair in the counts, the script reports
#   observed           calls with the observed labels
#   permutations       number of permuted assignments (n)
#   null_mean, null_sd, null_median, null_min, null_max
#   ratio_to_null_mean observed / null mean. Near or below 1, the observed labels returned
#                      no more calls than random ones.
#   ratio_to_null_max  observed / the largest permuted count. Above 1, no permutation
#                      returned as many calls.
#   exceeding          permutations returning at least as many calls as observed (k)
#   exceedance_p       (k + 1) / (n + 1), the permutation p-value. Its smallest possible
#                      value is 1 / (n + 1), so a few permutations cannot give a small one.
#
# A count at chance level does not mean that no true difference exists, and an observed
# count far above the null does not say which calls are true. The null describes the set
# of calls, not individual calls.
#
# Usage:
#   Rscript summarize_calibration.R --out summary.tsv counts.tsv [more_counts.tsv ...]
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3 || args[[1]] != "--out") {
  stop("usage: summarize_calibration.R --out summary.tsv counts.tsv [...]")
}
out <- args[[2]]
counts <- rbindlist(lapply(args[-(1:2)], fread), use.names = TRUE, fill = TRUE)
by <- c("caller", "level", "unit", "sig_cutoff", "diff_cutoff")
need <- c("assignment_id", "observed", "calls", by)
miss <- setdiff(need, names(counts))
if (length(miss)) stop("counts lack columns: ", paste(miss, collapse = ", "))
counts[, observed := as.logical(observed)]

summ <- counts[, {
  obs <- as.numeric(calls[observed])
  perm <- as.numeric(calls[!observed])
  if (length(obs) != 1) stop("expected one observed count for each threshold, found ", length(obs))
  k <- sum(perm >= obs)
  .(observed = obs, permutations = length(perm),
    null_mean = mean(perm), null_sd = sd(perm), null_median = median(perm),
    null_min = min(perm), null_max = max(perm),
    ratio_to_null_mean = obs / mean(perm), ratio_to_null_max = obs / max(perm),
    exceeding = k, exceedance_p = (k + 1) / (length(perm) + 1))
}, by = by]
setorderv(summ, by, order = c(1, 1, 1, -1, 1))
fwrite(summ, out, sep = "\t")
message("wrote ", out)
