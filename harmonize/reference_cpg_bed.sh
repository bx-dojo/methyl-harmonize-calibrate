#!/usr/bin/env bash
# List every CpG site of a reference genome as BED.
#
# One line per CpG site: sequence name, zero-based start of the forward-strand C, and
# end (start + 1). Calls from both strands of a CpG are combined on this position, and
# the number of lines is the denominator for "percentage of reference CpG sites".
#
# Usage:
#   reference_cpg_bed.sh REFERENCE.fa[.gz] OUT.bed
#
# The FASTA is streamed once. A CpG split across two sequence lines is found by
# carrying the last base of the previous line.
set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "usage: $0 REFERENCE.fa[.gz] OUT.bed" >&2
    exit 2
fi
reference=$1
out=$2
[ -r "$reference" ] || { echo "$0: cannot read $reference" >&2; exit 1; }

reader=cat
case "$reference" in
    *.gz) reader="gzip -dc" ;;
esac

$reader "$reference" | awk '
    BEGIN { OFS = "\t"; offset = 0; previous = ""; chrom = "" }
    /^>/ {
        chrom = substr($1, 2)      # the identifier is the first word of the header
        offset = 0
        previous = ""
        next
    }
    {
        line = toupper($0)
        n = length(line)
        for (i = 1; i <= n; i++) {
            base = substr(line, i, 1)
            if (previous == "C" && base == "G") {
                start = offset + i - 2     # the C is one position back
                print chrom, start, start + 1
            }
            previous = base
        }
        offset += n
    }
' > "$out"

echo "$0: wrote $(wc -l < "$out") CpG sites to $out" >&2
