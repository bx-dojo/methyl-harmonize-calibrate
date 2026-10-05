#!/usr/bin/env Rscript
# Label assignments for permutation calibration: the observed labels plus permuted ones.
#
# Rerunning a differential methylation analysis on permuted sample labels gives the
# number of calls the same pipeline returns when the groups differ only by chance (the
# null distribution). Two designs are supported.
#
#   free    Labels are shuffled freely among all samples, keeping the group sizes. When
#           the groups are the same size, an assignment and its mirror image (the two
#           groups swapped) give the same two-sided test, so only one of each pair is
#           kept. The observed assignment is excluded, and --n assignments are drawn at
#           random (with --seed) from the rest, or all of them if there are fewer.
#   paired  For paired samples (for example two tissues from each individual), labels
#           are exchanged only within pairs. The first pair is never exchanged, which
#           removes the mirror images, so P pairs give 2^(P-1) - 1 assignments besides
#           the observed one. With --n, a random subset is drawn.
#
# Usage:
#   Rscript make_assignments.R --samples samples.tsv --design free|paired \
#       --out assignments.tsv [--n 100] [--seed 1]
#
# samples.tsv: tab-separated with a header and the columns sample and group (exactly two
# group labels), plus pair for the paired design. Other columns are ignored. The order of
# the rows matters for the random draw, so keep it fixed to reproduce an analysis.
#
# assignments.tsv has one row per assignment: assignment_id ("observed" or perm_001,
# perm_002, ...), observed (TRUE/FALSE), and group_a and group_b, the comma-separated
# samples given the first and the second group label.

args <- commandArgs(trailingOnly = TRUE)
opt <- list(samples = NULL, design = NULL, out = NULL, n = NA, seed = "1")
for (i in seq(1, length(args), by = 2)) {
  key <- gsub("-", "_", sub("^--", "", args[[i]]))
  if (!key %in% names(opt) || i == length(args)) stop("bad argument: ", args[[i]])
  opt[[key]] <- args[[i + 1]]
}
if (is.null(opt$samples) || is.null(opt$design) || is.null(opt$out)) {
  stop("usage: make_assignments.R --samples F --design free|paired --out F [--n N] [--seed S]")
}
seed <- as.integer(opt$seed)
n_wanted <- as.integer(opt$n)

sheet <- read.delim(opt$samples, colClasses = "character")
ids <- sheet$sample
labels <- unique(sheet$group)
if (length(labels) != 2) stop("the group column must hold exactly two labels")
if (anyDuplicated(ids)) stop("duplicated sample names")
observed <- as.integer(sheet$group == labels[2])     # 0 = group_a, 1 = group_b
n <- length(ids)

key <- function(tr) paste(tr, collapse = "")
equal_size <- sum(observed) * 2 == n
canonical <- function(tr) if (equal_size) min(key(tr), key(1L - tr)) else key(tr)

draw <- function(treatments, n_wanted) {
  total <- length(treatments)
  if (is.na(n_wanted) || n_wanted >= total) {
    if (!is.na(n_wanted) && n_wanted > total) {
      message("only ", total, " distinct assignments exist; using all of them")
    }
    return(treatments)
  }
  set.seed(seed)
  treatments[sample.int(total, n_wanted)]
}

if (opt$design == "free") {
  seen <- new.env(parent = emptyenv())
  treatments <- list()
  for (idx in utils::combn(n, sum(observed), simplify = FALSE)) {
    tr <- integer(n)
    tr[idx] <- 1L
    k <- canonical(tr)
    if (k == canonical(observed) || exists(k, envir = seen, inherits = FALSE)) next
    assign(k, TRUE, envir = seen)
    if (equal_size && tr[1] == 1L) tr <- 1L - tr   # one fixed orientation per mirror pair
    treatments[[length(treatments) + 1]] <- tr
  }
  treatments <- draw(treatments, n_wanted)
} else if (opt$design == "paired") {
  if (is.null(sheet$pair)) stop("the paired design needs a pair column")
  pairs <- sort(unique(sheet$pair))
  members <- lapply(pairs, function(p) which(sheet$pair == p))
  ok <- vapply(members, function(m) length(m) == 2 && sum(observed[m]) == 1, logical(1))
  if (!all(ok)) stop("every pair must hold one sample of each group: ",
                     paste(pairs[!ok], collapse = ", "))
  swaps <- cbind(0L, as.matrix(expand.grid(rep(list(c(0L, 1L)), length(pairs) - 1))))
  swaps <- swaps[rowSums(swaps) > 0, , drop = FALSE]
  treatments <- lapply(seq_len(nrow(swaps)), function(i) {
    tr <- observed
    for (m in members[swaps[i, ] == 1L]) tr[m] <- 1L - tr[m]
    tr
  })
  treatments <- draw(treatments, n_wanted)
} else {
  stop("--design must be free or paired")
}

row <- function(id, obs, tr) {
  data.frame(assignment_id = id, observed = obs,
             group_a = paste(ids[tr == 0L], collapse = ","),
             group_b = paste(ids[tr == 1L], collapse = ","))
}
out <- rbind(row("observed", TRUE, observed),
             do.call(rbind, lapply(seq_along(treatments), function(i) {
               row(sprintf("perm_%03d", i), FALSE, treatments[[i]])
             })))
write.table(out, opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote ", opt$out, ": the observed assignment and ", length(treatments),
        " permuted ones (", opt$design, " design, groups ", labels[1], " and ", labels[2], ")")
