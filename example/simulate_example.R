#!/usr/bin/env Rscript
# Simulate a small data set that exercises every script in this repository, and record
# the values the harmonize scripts should return, so run_example.sh can check them.
#
#   Rscript simulate_example.R OUT_DIR
#
# Alignment part: a 40-kb reference and two samples aligned by two imaginary aligners.
#   unique_only  writes only uniquely aligned read pairs (like Bismark or ABISMAL), and its
#                deduplication tool removes duplicates. Its report counts the
#                multi-mapping pairs that its BAM leaves out.
#   bwa_like     also writes multi-mapping pairs with MAPQ 0, unmapped pairs, secondary
#                and supplementary records, and pairs with one mate not uniquely placed
#                (like BWA-meth or BISCUIT). Duplicates are only marked, so one BAM serves
#                both stages.
# Methylation part: 12 samples, six pairs (individuals) with one sample in each of
# groups A and B, with true A-B differences in 30 regions, and a second set of 12 samples
# in groups C and D with no true difference, in which each individual deviates strongly
# from the others at about 1% of sites, so a test that pools reads can still make calls.
set.seed(20261005)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: simulate_example.R OUT_DIR")
out <- args[[1]]
dir.create(file.path(out, "bam"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "calls"), recursive = TRUE, showWarnings = FALSE)
tsv <- function(d, f) write.table(d, file.path(out, f), sep = "\t", quote = FALSE,
                                  row.names = FALSE)

# ------------------------------------------------------------------ reference
chroms <- c(chr1 = 20000, chr2 = 20000)
ref <- lapply(chroms, function(n) sample(c("A", "C", "G", "T"), n, replace = TRUE))
fa <- unlist(lapply(names(ref), function(ch) {
  s <- paste(ref[[ch]], collapse = "")
  c(paste0(">", ch, " simulated"), substring(s, seq(1, nchar(s), 60), pmin(seq(60, nchar(s) + 59, 60), nchar(s))))
}))
writeLines(fa, file.path(out, "reference.fa"))
cpg <- do.call(rbind, lapply(names(ref), function(ch) {
  b <- ref[[ch]]
  data.frame(chr = ch, pos = which(b[-length(b)] == "C" & b[-1] == "G"))   # 1-based C
}))

# ------------------------------------------------------------------ alignments
L <- 50; insert <- 150
read_seq <- function(ch, p) paste(ref[[ch]][p:(p + L - 1)], collapse = "")
qual <- strrep("I", L)
sam_pair <- function(name, ch, p, f1, f2, q1 = 60, q2 = 60) {
  p2 <- p + insert - L
  c(paste(name, f1, ch, p, q1, paste0(L, "M"), "=", p2, insert, read_seq(ch, p), qual, sep = "\t"),
    paste(name, f2, ch, p2, q2, paste0(L, "M"), "=", p, -insert, read_seq(ch, p2), qual, sep = "\t"))
}
unmapped_pair <- function(name) {
  s <- strrep("N", L)
  c(paste(name, 77, "*", 0, 0, "*", "*", 0, 0, s, qual, sep = "\t"),
    paste(name, 141, "*", 0, 0, "*", "*", 0, 0, s, qual, sep = "\t"))
}
header <- c("@HD\tVN:1.6\tSO:unsorted",
            paste0("@SQ\tSN:", names(chroms), "\tLN:", chroms))
rand_pos <- function(n) {
  ch <- sample(names(chroms), n, replace = TRUE)
  data.frame(ch = ch, p = sample.int(min(chroms) - insert - 1, n, replace = TRUE))
}

bams <- list(); counts <- list(); clean <- list(); reported <- list(); reads <- list()
for (s in c("S1", "S2")) {
  n_u <- 3000; n_dup <- 300; n_multi <- 450; n_unm <- 150; n_half <- 90
  u <- rand_pos(n_u); m <- rand_pos(n_multi); h <- rand_pos(n_half)
  d <- u[sample.int(n_u, n_dup), ]
  clean[[s]] <- data.frame(sample = s, clean_pairs = n_u + n_dup + n_multi + n_unm + n_half)

  # unique_only: unique pairs and their duplicates. Multi-mapping pairs and pairs with an
  # ambiguous mate are left out of the BAM, and unmapped pairs are not written.
  uo_align <- c(unlist(Map(function(i) sam_pair(paste0(s, "_u", i), u$ch[i], u$p[i], 99, 147), seq_len(n_u))),
                unlist(Map(function(i) sam_pair(paste0(s, "_d", i), d$ch[i], d$p[i], 99, 147), seq_len(n_dup))))
  uo_dedup <- uo_align[seq_len(2 * n_u)]
  reported[[s]] <- data.frame(aligner = "unique_only", sample = s, mapped_pairs = n_u + n_dup + n_multi)

  # bwa_like: everything, duplicates marked (0x400), plus MAPQ 0 pairs, secondary and
  # supplementary records, unmapped pairs, and pairs with one mate at MAPQ 0 that are not
  # flagged as proper pairs.
  bl <- c(uo_align[seq_len(2 * n_u)],
          unlist(Map(function(i) sam_pair(paste0(s, "_d", i), d$ch[i], d$p[i], 99 + 1024, 147 + 1024), seq_len(n_dup))),
          unlist(Map(function(i) sam_pair(paste0(s, "_m", i), m$ch[i], m$p[i], 99, 147, 0, 0), seq_len(n_multi))),
          unlist(Map(function(i) sam_pair(paste0(s, "_h", i), h$ch[i], h$p[i], 97, 145, 60, 0), seq_len(n_half))),
          unlist(lapply(seq_len(n_unm), function(i) unmapped_pair(paste0(s, "_n", i)))),
          unlist(Map(function(i) sam_pair(paste0(s, "_u", i), u$ch[i], u$p[i] + 7, 99 + 256, 147 + 256, 0, 0), 1:60))[c(TRUE, FALSE)],
          unlist(Map(function(i) sam_pair(paste0(s, "_u", i), u$ch[i], u$p[i] + 3, 99 + 2048, 147 + 2048), 61:90))[c(TRUE, FALSE)])
  for (nm in c("unique_only_align", "unique_only_dedup", "bwa_like")) {
    lines <- switch(nm, unique_only_align = uo_align, unique_only_dedup = uo_dedup, bwa_like = bl)
    writeLines(c(header, lines), file.path(out, "bam", paste0(s, "_", nm, ".sam")))
  }
  # expected harmonized counts (reads, after the flag filters of count_alignments.sh)
  counts[[s]] <- rbind(
    data.frame(aligner = "unique_only", sample = s, stage = "align", unique_mapq = 0,
               total = 2 * (n_u + n_dup), mapped = 2 * (n_u + n_dup), unique = 2 * (n_u + n_dup),
               proper_pairs = 2 * (n_u + n_dup)),
    data.frame(aligner = "unique_only", sample = s, stage = "dedup", unique_mapq = 0,
               total = 2 * n_u, mapped = 2 * n_u, unique = 2 * n_u, proper_pairs = 2 * n_u),
    data.frame(aligner = "bwa_like", sample = s, stage = "align", unique_mapq = 1,
               total = 2 * (n_u + n_dup + n_multi + n_half + n_unm),
               mapped = 2 * (n_u + n_dup + n_multi + n_half),
               unique = 2 * (n_u + n_dup) + n_half, proper_pairs = 2 * (n_u + n_dup)),
    data.frame(aligner = "bwa_like", sample = s, stage = "dedup", unique_mapq = 1,
               total = 2 * (n_u + n_multi + n_half + n_unm), mapped = 2 * (n_u + n_multi + n_half),
               unique = 2 * n_u + n_half, proper_pairs = 2 * n_u))
  # reads that count toward depth after deduplication: (chr, start, mapq)
  rd <- function(df, q1, q2) rbind(data.frame(ch = df$ch, p = df$p, q = q1),
                                   data.frame(ch = df$ch, p = df$p + insert - L, q = q2))
  reads[[s]] <- list(unique_only = rd(u, 60, 60),
                     bwa_like = rbind(rd(u, 60, 60), rd(m, 0, 0), rd(h, 60, 0)))
}
tsv(do.call(rbind, clean), "clean_pairs.tsv")
tsv(do.call(rbind, reported), "reported_mapped.tsv")
tsv(do.call(rbind, counts), "expected_alignment_counts.tsv")

# expected CpG coverage, union across the two samples, at MAPQ 0 and 1
depth_at <- function(r, minq) {
  r <- r[r$q >= minq, ]
  key <- paste(cpg$chr, cpg$pos)
  dep <- integer(nrow(cpg))
  for (ch in names(chroms)) {
    cov <- tabulate(unlist(lapply(r$p[r$ch == ch], function(p) p:(p + L - 1))), chroms[[ch]])
    i <- cpg$chr == ch
    dep[i] <- cov[cpg$pos[i]]
  }
  dep
}
expected_cov <- list(); calls_sheet <- list()
for (a in c("unique_only", "bwa_like")) for (minq in c(0, 1)) {
  d <- sapply(c("S1", "S2"), function(s) depth_at(reads[[s]][[a]], minq))
  for (th in c(1, 5, 10)) expected_cov[[length(expected_cov) + 1]] <- data.frame(
    aligner = a, mapq = minq, min_depth = th, union_sites = sum(rowSums(d >= th) > 0),
    intersection_sites = sum(rowSums(d >= th) == 2))
}
tsv(do.call(rbind, expected_cov), "expected_cpg_coverage.tsv")

# methylation calls from each workflow: the unique_only caller keeps every covered CpG,
# the bwa_like caller only reads with MAPQ >= 40, as BISCUIT does by default
expected_ret <- list()
for (a in c("unique_only", "bwa_like")) {
  d <- sapply(c("S1", "S2"), function(s) depth_at(reads[[s]][[a]], if (a == "bwa_like") 40 else 0))
  for (j in 1:2) {
    f <- file.path(out, "calls", paste0(a, "_S", j, ".tsv.gz"))
    k <- d[, j] > 0
    data.table::fwrite(data.frame(chr = cpg$chr[k], pos = cpg$pos[k], coverage = d[k, j],
                                  methylated = rbinom(sum(k), d[k, j], 0.8)), f, sep = "\t")
    calls_sheet[[length(calls_sheet) + 1]] <- data.frame(aligner = a, sample = paste0("S", j), file = f)
  }
  d0 <- sapply(c("S1", "S2"), function(s) depth_at(reads[[s]][[a]], 0))
  expected_ret[[a]] <- data.frame(aligner = a, call_union_sites = sum(rowSums(d > 0) > 0),
                                  bam_union_sites_mapq0 = sum(rowSums(d0 > 0) > 0))
}
tsv(do.call(rbind, calls_sheet), "call_sheet.tsv")
tsv(do.call(rbind, expected_ret), "expected_retention.tsv")

# ------------------------------------------------------------------ methylation calls
# 20,000 CpG sites in clusters on two sequences, 12 samples per data set
n_sites <- 20000
sites <- data.frame(chr = rep(c("chr1", "chr2"), each = n_sites / 2),
                    pos = c(sort(sample.int(400000, n_sites / 2)), sort(sample.int(400000, n_sites / 2))))
sites <- sites[!duplicated(sites), ]
base <- ifelse(runif(nrow(sites)) < 0.75, rbeta(nrow(sites), 8, 1.5), rbeta(nrow(sites), 1.5, 8))
dmr_start <- sample(which(diff(sites$pos) < 60 & diff(sites$pos) > 0), 30)
effect <- numeric(nrow(sites))
for (i in dmr_start) effect[i:min(i + 9, nrow(sites))] <- sample(c(-0.4, 0.4), 1)
logit <- function(p) log(p / (1 - p)); expit <- function(x) 1 / (1 + exp(-x))
p_base <- pmin(pmax(base, 0.02), 0.98)
write_sample <- function(name, p) {
  cov <- rpois(nrow(sites), 15)
  sparse <- runif(nrow(sites)) < 0.03                      # a few sparse sites
  cov[sparse] <- rpois(sum(sparse), 2)
  meth <- rbinom(nrow(sites), cov, p)
  f <- file.path(out, "calls", paste0(name, ".tsv.gz"))
  data.table::fwrite(data.frame(chr = sites$chr, pos = sites$pos, coverage = cov,
                                methylated = meth)[cov > 0, ], f, sep = "\t")
  f
}
paired <- list(); unpaired <- list()
for (i in 1:6) {
  ind <- rnorm(nrow(sites), 0, 0.3)                       # individual differences
  pa <- expit(logit(p_base) + ind + rnorm(nrow(sites), 0, 0.2))
  pb <- pmin(pmax(pa + effect, 0.01), 0.99)
  paired[[length(paired) + 1]] <- data.frame(sample = paste0("A", i), group = "A", pair = i,
                                             file = write_sample(paste0("A", i), pa))
  paired[[length(paired) + 1]] <- data.frame(sample = paste0("B", i), group = "B", pair = i,
                                             file = write_sample(paste0("B", i), pb))
}
for (g in c("C", "D")) for (i in 1:6) {
  p <- expit(logit(p_base) + rnorm(nrow(sites), 0, 0.35))
  odd <- runif(nrow(sites)) < 0.01                       # sites where this individual
  p[odd] <- 1 - p[odd]                                   # deviates strongly from the rest
  unpaired[[length(unpaired) + 1]] <- data.frame(sample = paste0(g, i), group = g,
                                                 file = write_sample(paste0(g, i), p))
}
tsv(do.call(rbind, paired), "samples_paired.tsv")
tsv(do.call(rbind, unpaired), "samples_free.tsv")
message("simulated example written to ", out)
