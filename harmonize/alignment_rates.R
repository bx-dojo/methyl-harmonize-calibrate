#!/usr/bin/env Rscript
# Harmonized alignment rates from the counts written by count_alignments.sh.
#
# All rates are percentages of the clean read pairs given to the aligners, which should
# be the same for every aligner (for example the read pairs left after trimming).
#
#   mapped_pct        mapped reads / (2 x clean pairs). Aligners that write only uniquely
#                     aligned reads (unique_mapq 0) cannot give an overall mapping rate
#                     from their BAM, so supply their report's count of mapped pairs with
#                     --reported-mapped, or the column is left empty for them.
#   unique_pct        uniquely mapped fragments / clean pairs, before deduplication
#   dedup_unique_pct  uniquely mapped fragments / clean pairs, after deduplication. These
#                     are the reads passed on to methylation calling.
#   duplicate_pct     fragments removed by deduplication / fragments before it
#
# A fragment is a read pair, or a single read whose mate was not uniquely placed:
#   fragments = unique - proper_pairs / 2
# This also holds for a deduplicated BAM that stores one record per fragment (for
# example after dnmtools uniq), where no record is flagged as properly paired.
#
# Usage:
#   Rscript alignment_rates.R --counts counts.tsv --clean-pairs clean_pairs.tsv \
#       --out rates.tsv [--reported-mapped reported_mapped.tsv]
#
# clean_pairs.tsv:     sample, clean_pairs
# reported_mapped.tsv: aligner, sample, mapped_pairs   (optional)
# Writes rates.tsv (one row per aligner and sample) and rates_mean.tsv (mean and SD per
# aligner).

args <- commandArgs(trailingOnly = TRUE)
opt <- list(counts = NULL, clean_pairs = NULL, out = NULL, reported_mapped = NULL)
for (i in seq(1, length(args), by = 2)) {
  key <- sub("^--", "", args[[i]])
  key <- gsub("-", "_", key)
  if (!key %in% names(opt) || i == length(args)) stop("bad argument: ", args[[i]])
  opt[[key]] <- args[[i + 1]]
}
if (is.null(opt$counts) || is.null(opt$clean_pairs) || is.null(opt$out)) {
  stop("usage: alignment_rates.R --counts F --clean-pairs F --out F [--reported-mapped F]")
}

read_tsv <- function(f) read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
counts <- read_tsv(opt$counts)
pairs <- read_tsv(opt$clean_pairs)
reported <- if (is.null(opt$reported_mapped)) NULL else read_tsv(opt$reported_mapped)

counts$fragments <- counts$unique - counts$proper_pairs / 2
key <- function(a, s) paste(a, s, sep = "\r")
get <- function(stage, col) {
  x <- counts[counts$stage == stage, ]
  setNames(x[[col]], key(x$aligner, x$sample))
}

ids <- unique(counts[, c("aligner", "sample")])
k <- key(ids$aligner, ids$sample)
clean <- setNames(pairs$clean_pairs, pairs$sample)[ids$sample]
if (anyNA(clean)) stop("no clean_pairs for: ", paste(unique(ids$sample[is.na(clean)]), collapse = ", "))

unique_only <- get("align", "unique_mapq")[k] == 0
mapped_pct <- 100 * get("align", "mapped")[k] / (2 * clean)
mapped_pct[unique_only] <- NA
if (!is.null(reported)) {
  rep <- setNames(reported$mapped_pairs, key(reported$aligner, reported$sample))[k]
  mapped_pct[!is.na(rep)] <- 100 * rep[!is.na(rep)] / clean[!is.na(rep)]
}

frag_align <- get("align", "fragments")[k]
frag_dedup <- get("dedup", "fragments")[k]
rates <- data.frame(
  aligner = ids$aligner, sample = ids$sample, clean_pairs = unname(clean),
  mapped_pct = unname(mapped_pct),
  unique_pct = unname(100 * frag_align / clean),
  dedup_unique_pct = unname(100 * frag_dedup / clean),
  duplicate_pct = unname(100 * (frag_align - frag_dedup) / frag_align)
)
rates <- rates[order(rates$aligner, rates$sample), ]
write.table(rates, opt$out, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

metrics <- c("mapped_pct", "unique_pct", "dedup_unique_pct", "duplicate_pct")
summ <- do.call(rbind, lapply(split(rates, rates$aligner), function(d) {
  data.frame(aligner = d$aligner[1], samples = nrow(d),
             t(unlist(lapply(metrics, function(m) {
               setNames(c(mean(d[[m]]), sd(d[[m]])), paste0(m, c("_mean", "_sd")))
             }))))
}))
mean_out <- sub("(\\.tsv)?$", "_mean.tsv", opt$out)
write.table(summ, mean_out, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
message("wrote ", opt$out, " and ", mean_out)
