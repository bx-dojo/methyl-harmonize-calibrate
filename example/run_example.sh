#!/usr/bin/env bash
# Run every script on the simulated example and check the harmonized metrics against
# the values the simulation knows to be true. Takes about 12 minutes.
#
#   bash example/run_example.sh [OUT_DIR]      (default: example_output)
set -euo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
out=${1:-example_output}
mkdir -p "$out"
out=$(cd "$out" && pwd)

echo "== simulating data"
Rscript "$repo/example/simulate_example.R" "$out"

# ---------------------------------------------------------------- harmonized metrics
echo "== alignment counts"
printf 'aligner\tsample\tstage\tbam\tunique_mapq\n' > "$out/bam_sheet.tsv"
printf 'aligner\tsample\tbam\n' > "$out/dedup_bam_sheet.tsv"
for s in S1 S2; do
    for f in unique_only_align unique_only_dedup bwa_like; do
        samtools sort -o "$out/bam/${s}_$f.bam" "$out/bam/${s}_$f.sam" 2>/dev/null
        samtools index "$out/bam/${s}_$f.bam"
    done
    printf 'unique_only\t%s\talign\t%s\t0\n' $s "$out/bam/${s}_unique_only_align.bam" >> "$out/bam_sheet.tsv"
    printf 'unique_only\t%s\tdedup\t%s\t0\n' $s "$out/bam/${s}_unique_only_dedup.bam" >> "$out/bam_sheet.tsv"
    printf 'bwa_like\t%s\talign\t%s\t1\n'    $s "$out/bam/${s}_bwa_like.bam" >> "$out/bam_sheet.tsv"
    printf 'bwa_like\t%s\tdedup\t%s\t1\n'    $s "$out/bam/${s}_bwa_like.bam" >> "$out/bam_sheet.tsv"
    printf 'unique_only\t%s\t%s\n' $s "$out/bam/${s}_unique_only_dedup.bam" >> "$out/dedup_bam_sheet.tsv"
    printf 'bwa_like\t%s\t%s\n'    $s "$out/bam/${s}_bwa_like.bam" >> "$out/dedup_bam_sheet.tsv"
done
bash "$repo/harmonize/count_alignments.sh" "$out/bam_sheet.tsv" "$out/alignment_counts.tsv"
Rscript "$repo/harmonize/alignment_rates.R" --counts "$out/alignment_counts.tsv" \
    --clean-pairs "$out/clean_pairs.tsv" --reported-mapped "$out/reported_mapped.tsv" \
    --out "$out/alignment_rates.tsv"

echo "== CpG coverage and retention"
bash "$repo/harmonize/reference_cpg_bed.sh" "$out/reference.fa" "$out/reference_cpg.bed"
for q in 0 1; do
    Rscript "$repo/harmonize/bam_cpg_coverage.R" --bams "$out/dedup_bam_sheet.tsv" \
        --cpg-bed "$out/reference_cpg.bed" --mapq $q --threads 1 --out-dir "$out/coverage_mapq$q"
done
Rscript "$repo/harmonize/cpg_retention.R" --calls "$out/call_sheet.tsv" \
    --bam-summary "$out/coverage_mapq0/bam_cpg_coverage_summary.tsv" --out "$out/retention.tsv"

echo "== checking harmonized metrics against the simulation"
Rscript - "$out" <<'EOF'
out <- commandArgs(trailingOnly = TRUE)[1]
r <- function(f) read.delim(file.path(out, f), stringsAsFactors = FALSE)
k <- c("aligner", "sample", "stage")
a <- merge(r("expected_alignment_counts.tsv"), r("alignment_counts.tsv"), by = k)
cols <- c("total", "mapped", "unique", "proper_pairs")
ok1 <- nrow(a) == 8 && all(a[paste0(cols, ".x")] == a[paste0(cols, ".y")])
cov <- rbind(r("coverage_mapq0/bam_cpg_coverage_summary.tsv"), r("coverage_mapq1/bam_cpg_coverage_summary.tsv"))
c2 <- merge(r("expected_cpg_coverage.tsv"), cov, by = c("aligner", "mapq", "min_depth"))
ok2 <- nrow(c2) == 12 && all(c2$union_sites.x == c2$union_sites.y) &&
       all(c2$intersection_sites.x == c2$intersection_sites.y)
re <- merge(r("expected_retention.tsv"), r("retention.tsv"), by = "aligner")
ok3 <- nrow(re) == 2 && all(re$call_union_sites.x == re$call_union_sites.y) &&
       all(re$bam_union_sites_mapq0.x == re$bam_union_sites_mapq0.y)
cat(sprintf("  alignment counts: %s\n  CpG coverage:     %s\n  retention:        %s\n",
            if (ok1) "match" else "MISMATCH", if (ok2) "match" else "MISMATCH",
            if (ok3) "match" else "MISMATCH"))
if (!(ok1 && ok2 && ok3)) quit(status = 1)
EOF

# ---------------------------------------------------------------- permutation calibration
echo "== permutation calibration (about 12 minutes; progress in $out/calibration.log)"
cal="$repo/calibrate"
log="$out/calibration.log"
: > "$log"
run() { "$@" >> "$log" 2>&1 || { echo "failed: $*" >&2; tail -n 20 "$log" >&2; exit 1; }; }
# A vs B: two samples from each of six individuals, with true differences (paired design)
run Rscript "$cal/make_assignments.R" --samples "$out/samples_paired.tsv" --design paired \
    --out "$out/assignments_paired.tsv"
# C vs D: twelve unrelated samples with no true difference (free design, 20 permutations)
run Rscript "$cal/make_assignments.R" --samples "$out/samples_free.tsv" --design free \
    --n 20 --seed 1 --out "$out/assignments_free.tsv"

for design in paired free; do
    run Rscript "$cal/run_methylkit.R" --samples "$out/samples_$design.tsv" \
        --assignments "$out/assignments_$design.tsv" --out-dir "$out/methylkit_dml_$design" --cores 2
    run Rscript "$cal/run_methylkit.R" --samples "$out/samples_$design.tsv" --overdispersion MN \
        --assignments "$out/assignments_$design.tsv" --out-dir "$out/methylkit_mn_dml_$design" --cores 2
    run Rscript "$cal/run_dss.R" --samples "$out/samples_$design.tsv" --level dmr \
        --assignments "$out/assignments_$design.tsv" --out-dir "$out/dss_$design" --ncores 2
done
run Rscript "$cal/run_methylkit.R" --samples "$out/samples_paired.tsv" --level dmr \
    --assignments "$out/assignments_paired.tsv" --out-dir "$out/methylkit_dmr_paired" --cores 2

for design in paired free; do
    files=("$out"/methylkit_dml_$design/counts.tsv "$out"/methylkit_mn_dml_$design/counts.tsv
           "$out"/dss_$design/counts.tsv)
    [ "$design" = paired ] && files+=("$out/methylkit_dmr_paired/counts.tsv")
    run Rscript "$cal/summarize_calibration.R" --out "$out/calibration_$design.tsv" "${files[@]}"
done

Rscript - "$out" <<'EOF'
out <- commandArgs(trailingOnly = TRUE)[1]
for (d in c("paired", "free")) {
  x <- read.delim(file.path(out, paste0("calibration_", d, ".tsv")))
  x <- x[x$sig_cutoff == 0.01 & x$diff_cutoff %in% c(0, 25), ]
  cat("\n", if (d == "paired") "A vs B (true differences, paired)" else "C vs D (no true difference, free)",
      " at FDR < 0.01\n", sep = "")
  print(x[, c("caller", "level", "unit", "diff_cutoff", "observed", "permutations", "null_mean",
              "null_max", "ratio_to_null_mean", "exceedance_p")], row.names = FALSE, digits = 3)
}
EOF
echo
echo "All outputs are in $out"
