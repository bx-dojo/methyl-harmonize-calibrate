# Shared helpers for run_methylkit.R and run_dss.R. Sourced, not run.

parse_opts <- function(args, defaults) {
  opt <- defaults
  i <- 1
  while (i <= length(args)) {
    key <- gsub("-", "_", sub("^--", "", args[[i]]))
    if (!key %in% names(opt)) stop("unknown argument: ", args[[i]])
    if (is.logical(defaults[[key]])) {          # flags take no value
      opt[[key]] <- TRUE
      i <- i + 1
      next
    }
    if (i == length(args)) stop("missing value after ", args[[i]])
    opt[[key]] <- args[[i + 1]]
    i <- i + 2
  }
  opt
}

num_list <- function(x) as.numeric(strsplit(as.character(x), ",")[[1]])

# Per-sample CpG calls: tab-separated with a header and the columns chr, pos (1-based
# position of the forward-strand C, both strands combined), coverage, and methylated.
read_calls <- function(file) {
  x <- data.table::fread(file, select = c("chr", "pos", "coverage", "methylated"))
  x[coverage > 0]
}

# The same per-sample coverage filter as methylKit::filterByCoverage: drop sites below
# min_coverage reads and above the max_percentile-th percentile of the sample's coverage.
coverage_filter <- function(x, min_coverage, max_percentile) {
  ceiling <- stats::quantile(x$coverage, max_percentile / 100, names = FALSE)
  x[coverage >= min_coverage & coverage <= ceiling]
}

read_design <- function(samples_file, assignments_file, ids = NULL) {
  samples <- read.delim(samples_file, colClasses = "character")
  assignments <- read.delim(assignments_file, colClasses = "character")
  assignments$observed <- as.logical(assignments$observed)
  if (!is.null(ids)) {
    keep <- strsplit(ids, ",")[[1]]
    bad <- setdiff(keep, assignments$assignment_id)
    if (length(bad)) stop("unknown assignment IDs: ", paste(bad, collapse = ", "))
    assignments <- assignments[assignments$assignment_id %in% keep, ]
  }
  for (i in seq_len(nrow(assignments))) {
    s <- c(strsplit(assignments$group_a[i], ",")[[1]], strsplit(assignments$group_b[i], ",")[[1]])
    if (!setequal(s, samples$sample) || anyDuplicated(s)) {
      stop("assignment ", assignments$assignment_id[i], " does not use every sample once")
    }
  }
  list(samples = samples, assignments = assignments)
}

# One checkpoint file per assignment, so an interrupted run can be resumed and a long
# run can be split across jobs with --ids.
checkpoint_path <- function(out_dir, id) file.path(out_dir, "counts", paste0(id, ".tsv"))

write_checkpoint <- function(d, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp")
  data.table::fwrite(d, tmp, sep = "\t")
  file.rename(tmp, path)
}

combine_checkpoints <- function(out_dir) {
  files <- list.files(file.path(out_dir, "counts"), pattern = "\\.tsv$", full.names = TRUE)
  all <- data.table::rbindlist(lapply(files, data.table::fread), use.names = TRUE, fill = TRUE)
  data.table::fwrite(all, file.path(out_dir, "counts.tsv"), sep = "\t")
  message("wrote ", file.path(out_dir, "counts.tsv"), " (", length(files), " assignments)")
}

# Merge runs of adjacent significant windows that change in the same direction into
# regions. Windows are adjacent when one starts right after the previous one ends.
merge_windows <- function(chr, start, end, diff) {
  if (!length(chr)) return(0L)
  o <- order(chr, start)
  chr <- chr[o]; start <- start[o]; end <- end[o]; up <- diff[o] > 0
  n <- length(chr)
  new_region <- c(TRUE, chr[-1] != chr[-n] | start[-1] != end[-n] + 1 | up[-1] != up[-n])
  sum(new_region)
}
