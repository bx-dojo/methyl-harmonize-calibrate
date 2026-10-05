# methyl-harmonize-calibrate

Scripts for two checks that make DNA methylation sequencing analyses easier to compare and
to trust. Neither depends on the conversion chemistry (bisulfite or enzymatic) or on the
genome.

- **Harmonized metrics.** Aligners define a "mapped" or "uniquely mapped" read
  differently, so the rates they report cannot be compared directly. These scripts recount
  every aligner's BAM files with the same SAMtools rules, and count CpG coverage directly
  from the BAM files, which separates what the aligner covered from what the methylation
  caller kept.
- **Permutation calibration.** A differential methylation caller can return thousands of
  significant sites even when the groups do not differ. Rerunning the whole analysis with
  shuffled sample labels shows how many calls the same pipeline returns by chance, and
  the observed count can then be compared with that null distribution.

The scripts are generalized from the analysis in Saito et al. (in preparation),
*Benchmarking DNA methylation sequencing analysis with harmonized metrics and permutation
calibration in a non-model genome*. They implement the same procedures, but they are not
the pipeline used for the paper. Given the same counts, `alignment_rates.R` and
`summarize_calibration.R` reproduce the values reported there, and `make_assignments.R`
reproduces its permuted labels. Paths, computing resources, and the steps before these
scripts must be adapted to your own data.

## Contents

| Script | What it does |
|---|---|
| `harmonize/count_alignments.sh` | Counts total, mapped, uniquely mapped, and properly paired reads in each BAM with the same FLAG and MAPQ filters, before and after deduplication |
| `harmonize/alignment_rates.R` | Turns those counts into mapping, unique, deduplicated-unique, and duplicate rates over the same denominator |
| `harmonize/reference_cpg_bed.sh` | Lists every CpG site of a reference genome as BED |
| `harmonize/bam_cpg_coverage.R` | Counts read depth at every reference CpG site from the BAM files, with the same `samtools depth` options for every aligner |
| `harmonize/cpg_retention.R` | Share of the CpG sites covered in the BAM files that the methylation caller kept |
| `calibrate/make_assignments.R` | Observed and permuted label assignments, for unrelated (free) or paired samples |
| `calibrate/run_methylkit.R` | methylKit call counts for every assignment, for single CpG sites or windows and regions |
| `calibrate/run_dss.R` | DSS call counts for every assignment, for single CpG sites and, optionally, regions from `callDMR` |
| `calibrate/summarize_calibration.R` | Observed count against the null: ratio to the null mean and maximum, and the permutation p-value |
| `example/` | A small simulated data set that runs every script and checks the harmonized metrics against known values |

## Requirements

- bash 4.3 or later, and SAMtools 1.10 or later (tested with 1.19)
- R 4.3 or later with data.table
- For calibration: methylKit and DSS (with bsseq) from Bioconductor. Tested with
  methylKit 1.30 and DSS 2.52.

methylKit 1.30 fails to load with data.table 1.17 or later (`object 'key<-' is not
exported`). If you meet this error, use a newer methylKit or install data.table 1.16.4.

## Quick start

```bash
bash example/run_example.sh example_output
```

This simulates small sets of BAM files and methylation calls, runs every script, and
checks the alignment counts, CpG coverage, and retention against the values the
simulation knows to be true. It then calibrates methylKit and DSS on two comparisons.
Groups A and B are paired samples that differ in 30 regions, and every caller returns
far more calls than with any permuted labels. Groups C and D do not differ, but each
individual deviates strongly from the others at a few sites. There methylKit, whose
default test pools the reads of a group, still returns dozens of calls at FDR < 0.01,
yet no more than with shuffled labels (a ratio near 1). The example takes about 10
minutes, most of it for the permutations.

## Harmonized metrics

### 1. Alignment counts and rates

List every BAM in a tab-separated sheet with a header:

```
aligner   sample  stage  bam                         unique_mapq
bismark   S1      align  /path/S1_bismark.bam        0
bismark   S1      dedup  /path/S1_bismark.dedup.bam  0
bwameth   S1      align  /path/S1_bwameth.bam        1
bwameth   S1      dedup  /path/S1_bwameth.dedup.bam  1
```

`stage` is `align` before deduplication and `dedup` after it. A BAM in which duplicates
are only marked (as by BISCUIT's dupsifter) can be listed for both stages.
`unique_mapq` says how a uniquely mapped read is recognized. Use `0` when the BAM holds
only uniquely aligned reads (Bismark, ABISMAL), and `1` when multi-mapping reads are
written with MAPQ 0 (BWA-MEM based aligners such as BWA-meth and BISCUIT).

```bash
THREADS=2 JOBS=8 bash harmonize/count_alignments.sh bam_sheet.tsv counts.tsv
Rscript harmonize/alignment_rates.R --counts counts.tsv --clean-pairs clean_pairs.tsv \
    --reported-mapped reported_mapped.tsv --out rates.tsv
```

`clean_pairs.tsv` (columns `sample`, `clean_pairs`) gives the read pairs given to the
aligners, the same for every aligner. An aligner that writes only unique reads cannot
give its overall mapping rate from the BAM, so `reported_mapped.tsv` (columns `aligner`,
`sample`, `mapped_pairs`) can supply the count from its own report.

The counting rules:

| Quantity | SAMtools filter |
|---|---|
| Primary, QC-passed records | `-F 0xB00` (secondary, QC-failed, and supplementary excluded) |
| Mapped | `-F 0xB04` |
| After deduplication | `-F 0xF04` (duplicates also excluded) |
| Uniquely mapped | as mapped, with `-q unique_mapq` |
| Fragments | uniquely mapped reads minus half the uniquely mapped reads in proper pairs (`-f 0x2`) |
| Duplicate fraction | fragments removed by deduplication, as a share of fragments before it |

### 2. CpG coverage from the BAM files

```bash
bash harmonize/reference_cpg_bed.sh reference.fa reference_cpg.bed
for q in 0 1; do
    Rscript harmonize/bam_cpg_coverage.R --bams dedup_bam_sheet.tsv \
        --cpg-bed reference_cpg.bed --mapq $q --out-dir coverage_mapq$q
done
```

`dedup_bam_sheet.tsv` lists the deduplicated BAMs (columns `aligner`, `sample`, `bam`).
For each aligner and minimum depth, the summary gives the CpG sites covered in at least
one sample (union) and in every sample (intersection), as counts and as percentages of
all reference CpG sites, with their mean depth. Comparing MAPQ 0 with MAPQ 1 shows how
much coverage comes from multi-mapping reads.

### 3. Retention by the methylation caller

```bash
Rscript harmonize/cpg_retention.R --calls call_sheet.tsv \
    --bam-summary coverage_mapq0/bam_cpg_coverage_summary.tsv --out retention.tsv
```

`call_sheet.tsv` lists each workflow's methylation call files (columns `aligner`,
`sample`, `file`) in the format described below. Retention is the methylation-call union
divided by the BAM union at MAPQ 0, the BAM the caller reads, so it measures the caller's
filtering rather than read loss during alignment.

## Permutation calibration

### Input

A sample sheet (tab-separated, with a header) and one methylation call file per sample:

```
sample  group  pair  file
A1      A      1     calls/A1.tsv.gz
B1      B      1     calls/B1.tsv.gz
...
```

`pair` is needed only for paired samples, such as two tissues from each individual.
Each call file holds the CpG calls of one sample, with the two strands of each CpG
combined on the forward-strand C, as tab-separated columns with a header:

```
chr   pos     coverage  methylated
chr1  10468   14        12
```

`pos` is the 1-based position of the C, `coverage` the reads covering the site, and
`methylated` the reads showing it methylated. Most methylation callers can write this
format directly or after a small conversion.

### Steps

```bash
# 1. label assignments: the observed labels and the permuted ones
Rscript calibrate/make_assignments.R --samples samples.tsv --design free \
    --n 100 --seed 1 --out assignments.tsv          # or --design paired

# 2. rerun the analysis for every assignment
Rscript calibrate/run_methylkit.R --samples samples.tsv --assignments assignments.tsv \
    --out-dir methylkit_dml
Rscript calibrate/run_dss.R --samples samples.tsv --assignments assignments.tsv \
    --level dmr --out-dir dss --ncores 8

# 3. compare the observed counts with the null
Rscript calibrate/summarize_calibration.R --out calibration.tsv \
    methylkit_dml/counts.tsv dss/counts.tsv
```

**Designs.** With `free`, labels are shuffled among all samples, keeping the group sizes.
When the groups are the same size, an assignment and its mirror image (the two groups
swapped) give the same two-sided test, so only one of each is kept, and `--n` of the
remaining assignments are drawn with `--seed`. With `paired`, labels are exchanged only
within pairs, and the first pair is never exchanged, so P pairs give 2^(P-1) - 1
assignments. Some of these exchange only one pair and still differ mostly by the true
contrast, so they can return many calls. For a paired contrast, compare the observed
count with the largest permuted count as well as with the null mean.

**Analysis settings.** Both runners apply the same per-sample filter, at least
`--min-coverage` reads (10) and no more than the `--max-percentile` of the sample's
coverage (99.9). methylKit then tests only the sites that pass in every sample, at its
defaults. DSS tests every site that passes in at least one sample, and a sample that
failed there enters with no reads, which DSS treats as missing. `--all-samples` makes DSS
test the same sites as methylKit. Both count calls over a grid of significance cutoffs
(`--q`) and minimum absolute differences in percentage points (`--diff`). Settings other
than the labels are fixed across assignments, so differences in the counts come from
the labels alone.

**Run time.** Each assignment reruns the full test, so a calibration costs as much as the
original analysis times the number of assignments. Every assignment is written to its
own checkpoint in `OUT_DIR/counts/`, finished ones are skipped on a rerun, and `--ids`
runs a subset, so the work can be split across jobs that write to the same directory.

### Reading the output

`summarize_calibration.R` reports, for each caller, level, and threshold pair:

| Column | Meaning |
|---|---|
| `observed` | calls with the observed labels |
| `null_mean`, `null_max` | mean and largest count over the permuted labels |
| `ratio_to_null_mean` | observed / null mean. At or below 1, the observed labels returned no more calls than random ones |
| `ratio_to_null_max` | observed / null maximum. Above 1, no permutation returned as many calls |
| `exceedance_p` | (k + 1) / (n + 1), where k of n permutations returned at least as many calls as observed. The smallest possible value is 1 / (n + 1) |

The null describes the set of calls, not individual calls. An observed count at chance
level means that the count did not exceed chance, not that no true difference exists,
and a count far above the null does not say which calls are true. Inspecting the
methylation of each sample at called sites, and confirming selected sites with an
independent method, remain useful alongside calibration.

## Citation

If you use these scripts, please cite: Saito T, et al. Benchmarking DNA methylation
sequencing analysis with harmonized metrics and permutation calibration in a non-model
genome. In preparation.

## License

MIT. See [LICENSE](LICENSE).
