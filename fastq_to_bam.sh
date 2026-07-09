#!/usr/bin/env bash
#
# fastq_to_bam.sh
# Paired-end FASTQ -> aligned/sorted BAM -> single-contig extraction ->
# coverage bedGraph.
#
# Pipeline:
#   Trimmomatic (adapter/quality trim)
#   BWA-MEM     (align paired reads, and surviving unpaired reads single-end)
#   samtools    (sort, merge, index, extract one contig, rename it)
#   bedtools    (genomecov -> bedGraph)
#
#
# Requires on PATH: java, bwa, samtools, bedtools  (+ a Trimmomatic jar)
# The reference FASTA must already be BWA-indexed:  bwa index reference.fna
# ---------------------------------------------------------------------------

set -euo pipefail

# --------------------------- defaults --------------------------------------
threads=8
adapters="TruSeq3-PE.fa"
contig="NC_004354.4"                              # dm6 X chromosome
rename_to="X"
trimmomatic_jar="${HOME}/genomics/trimmomatic-0.39.jar"
prefix=""
keep_intermediates=0

# --------------------------- helpers ---------------------------------------
log() { printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*" >&2; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

usage() {
    cat >&2 <<EOF
Usage: $(basename "$0") -1 R1.fastq.gz -2 R2.fastq.gz -r reference.fna [options]

Required:
  -1 FILE   Forward (R1) reads, fastq(.gz)
  -2 FILE   Reverse (R2) reads, fastq(.gz)
  -r FILE   Reference FASTA (must be BWA-indexed)

Options:
  -o STR    Output prefix        (default: derived from the R1 filename)
  -t INT    Threads              (default: ${threads})
  -a FILE   Adapter FASTA        (default: ${adapters})
  -c STR    Contig to extract    (default: ${contig})
  -n STR    Rename contig to     (default: ${rename_to})
  -j FILE   Trimmomatic jar      (default: ${trimmomatic_jar})
  -k        Keep intermediate files
  -h        Show this help

Outputs: <prefix>.bam (+ .bai), <prefix><renamed>.bam, <prefix>.bedgraph
EOF
    exit "${1:-1}"
}

# --------------------------- parse args ------------------------------------
r1="" r2="" ref=""
while getopts "1:2:r:o:t:a:c:n:j:kh" opt; do
    case "$opt" in
        1) r1="$OPTARG" ;;
        2) r2="$OPTARG" ;;
        r) ref="$OPTARG" ;;
        o) prefix="$OPTARG" ;;
        t) threads="$OPTARG" ;;
        a) adapters="$OPTARG" ;;
        c) contig="$OPTARG" ;;
        n) rename_to="$OPTARG" ;;
        j) trimmomatic_jar="$OPTARG" ;;
        k) keep_intermediates=1 ;;
        h) usage 0 ;;
        *) usage 1 ;;
    esac
done

[[ -n "$r1" && -n "$r2" && -n "$ref" ]] || { echo "Missing required argument." >&2; usage 1; }

# Derive an output prefix from R1 if none was given
# (strips .gz / .fastq / .fq and a trailing _f, _R1 or _1)
if [[ -z "$prefix" ]]; then
    prefix="${r1%.gz}"; prefix="${prefix%.fastq}"; prefix="${prefix%.fq}"
    prefix="${prefix%_f}"; prefix="${prefix%_R1}"; prefix="${prefix%_1}"
fi

# --------------------------- sanity checks ---------------------------------
for tool in java bwa samtools bedtools; do
    command -v "$tool" >/dev/null 2>&1 || die "'$tool' not found on PATH."
done
[[ -f "$trimmomatic_jar" ]] || die "Trimmomatic jar not found: $trimmomatic_jar"
[[ -f "$r1"        ]] || die "R1 not found: $r1"
[[ -f "$r2"        ]] || die "R2 not found: $r2"
[[ -f "$ref"       ]] || die "Reference not found: $ref"
[[ -f "${ref}.bwt" ]] || die "Reference not BWA-indexed (missing ${ref}.bwt). Run: bwa index '$ref'"
[[ -f "$adapters"  ]] || die "Adapter FASTA not found: $adapters"

mkdir -p "$(dirname "$prefix")"

# --------------------------- names -----------------------------------------
r1p="${prefix}_R1_paired.fastq.gz"
r1u="${prefix}_R1_unpaired.fastq.gz"
r2p="${prefix}_R2_paired.fastq.gz"
r2u="${prefix}_R2_unpaired.fastq.gz"
paired_sorted="${prefix}_paired_sorted.bam"
singles_sorted="${prefix}_singles_sorted.bam"
full_bam="${prefix}.bam"
region_bam="${prefix}${rename_to}.bam"
bedgraph="${prefix}.bedgraph"

# escape dots in the contig name so sed matches it literally
esc_contig="$(printf '%s' "$contig" | sed 's/\./\\./g')"

# --------------------------- 1. trim ---------------------------------------
log "Trimmomatic (${threads} threads)"
java -jar "$trimmomatic_jar" PE -threads "$threads" -phred33 \
    "$r1" "$r2" \
    "$r1p" "$r1u" "$r2p" "$r2u" \
    ILLUMINACLIP:"${adapters}":2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36

# --------------------------- 2. align (piped straight to sort) -------------
# Piping bwa -> samtools sort avoids writing multi-GB intermediate SAM/BAM files.
log "BWA-MEM: paired reads -> sort"
bwa mem -t "$threads" "$ref" "$r1p" "$r2p" \
    | samtools sort -@ "$threads" -o "$paired_sorted" -

# Surviving unpaired reads are NOT mates of each other, so they are aligned
# single-end (concatenated), not as a fake pair.
log "BWA-MEM: unpaired reads (single-end) -> sort"
cat "$r1u" "$r2u" \
    | bwa mem -t "$threads" "$ref" - \
    | samtools sort -@ "$threads" -o "$singles_sorted" -

# --------------------------- 3. merge + index ------------------------------
log "Merge + index"
samtools merge -f -@ "$threads" "$full_bam" "$singles_sorted" "$paired_sorted"
samtools index "$full_bam"

# --------- 4. extract + rename X chromosome from NCBI to flybase -----------
log "Extract ${contig} and rename to ${rename_to}"
samtools view -@ "$threads" -b "$full_bam" "$contig" \
    | samtools reheader <(samtools view -H "$full_bam" | sed "s/${esc_contig}/${rename_to}/") - \
    > "$region_bam"

# --------------------------- 5. coverage bedGraph --------------------------
log "bedGraph coverage"
bedtools genomecov -ibam "$region_bam" -bg > "$bedgraph"

# --------------------------- 6. cleanup ------------------------------------
if [[ "$keep_intermediates" -eq 0 ]]; then
    log "Removing intermediates"
    rm -f "$r1p" "$r1u" "$r2p" "$r2u" "$paired_sorted" "$singles_sorted"
fi

log "Done."
log "  Whole-genome BAM : $full_bam (+ .bai)"
log "  ${rename_to} BAM  : $region_bam"
log "  bedGraph         : $bedgraph"
