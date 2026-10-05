#!/usr/bin/env bash
# Count reads in aligner BAM files with the same SAMtools rules for every aligner.
#
# Aligners report mapping and duplication rates on different definitions, so their own
# reports cannot be compared directly. This script recounts every BAM with one set of
# FLAG and MAPQ filters:
#
#   total        primary, QC-passed records        -F 0xB00 (align) / -F 0xF00 (dedup)
#   mapped       ... that are mapped               -F 0xB04 (align) / -F 0xF04 (dedup)
#   unique       ... with MAPQ >= unique_mapq      -q unique_mapq
#   proper_pairs ... unique records in proper pairs -f 0x2
#
# 0xB00 = secondary (0x100) + QC-failed (0x200) + supplementary (0x800). 0x004 is
# unmapped, and 0x400 (duplicate) is added after deduplication, so marked duplicates are
# excluded whether the deduplication tool marks or removes them.
#
# unique_mapq depends on how the aligner writes multi-mapping reads:
#   0  the BAM holds only uniquely aligned reads (for example Bismark, ABISMAL)
#   1  multi-mapping reads are written with MAPQ 0 (BWA-MEM based: BWA-meth, BISCUIT)
#
# Usage:
#   count_alignments.sh BAM_SHEET.tsv OUT.tsv
#
# BAM_SHEET.tsv is tab-separated with a header line and the columns
#   aligner  sample  stage  bam  unique_mapq
# where stage is "align" (before deduplication) or "dedup" (after). A BAM in which
# duplicates are only marked can be listed for both stages.
#
# Environment: THREADS (samtools threads per count, default 2) and JOBS (BAMs counted
# at once, default 4). Needs samtools and bash >= 4.3.
set -euo pipefail

if [ "$#" -ne 2 ]; then
    sed -n '2,30p' "$0" | sed 's/^# \{0,1\}//' >&2
    exit 2
fi
sheet=$1
out=$2
THREADS=${THREADS:-2}
JOBS=${JOBS:-4}
command -v samtools >/dev/null || { echo "samtools not found in PATH" >&2; exit 1; }

tmp=$(mktemp -d)
trap 'rm -r "$tmp"' EXIT

count_one() {  # index aligner sample stage bam unique_mapq
    local i=$1 aligner=$2 sample=$3 stage=$4 bam=$5 q=$6 base map
    case $stage in
        align) base=0xB00; map=0xB04 ;;
        dedup) base=0xF00; map=0xF04 ;;
        *) echo "stage must be align or dedup: $stage" >&2; return 1 ;;
    esac
    local total mapped unique pp
    total=$( samtools view -@ "$THREADS" -c -F "$base"               "$bam")
    mapped=$(samtools view -@ "$THREADS" -c -F "$map"                "$bam")
    unique=$(samtools view -@ "$THREADS" -c -F "$map" -q "$q"        "$bam")
    pp=$(    samtools view -@ "$THREADS" -c -F "$map" -q "$q" -f 0x2 "$bam")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$aligner" "$sample" "$stage" "$q" "$total" "$mapped" "$unique" "$pp" > "$tmp/$i"
}

throttle() { while (( $(jobs -rp | wc -l) >= JOBS )); do wait -n; done; }

i=0
while IFS=$'\t' read -r aligner sample stage bam q; do
    [ "$aligner" = aligner ] && continue          # header
    [ -z "$aligner" ] && continue
    [ -r "$bam" ] || { echo "cannot read $bam" >&2; exit 1; }
    i=$((i + 1))
    throttle
    count_one "$i" "$aligner" "$sample" "$stage" "$bam" "$q" &
done < "$sheet"
wait

{
    printf 'aligner\tsample\tstage\tunique_mapq\ttotal\tmapped\tunique\tproper_pairs\n'
    cat "$tmp"/* | sort
} > "$out"
echo "wrote $out ($i BAM counts)" >&2
