#!/bin/zsh

#takes .bam files from fastq_to_bam script and generates normalized .bed coverage files with coverage normalized to unbroken ring
#requires covdiv and winavg (https://github.com/jgbaldwinbrown/fastats)
#unbroken ring bedgraph files were generated using the bedtools genomecov -ibam command without follow-up steps

# reference (input unbroken ring .bed coverage file)
ref="unbroken_ring.bedgraph"

# 1. Check if a filename prefix was provided
if [ -z "$1" ]; then
    echo "Usage: process_genomics <filename_prefix>"
    echo "Example: process_genomics /Volumes/JTRv2/ring_breakage/fastqs/19787R/25F/25F"
    exit 1
fi

# 2. Assign variables
filename="$1"

# 3. Verify the BAM file exists
if [[ ! -f "${filename}.bam" ]]; then
    echo "Error: ${filename}.bam not found."
    exit 1
fi

echo "--- Starting Pipeline for ${filename} ---"

# 4. Index BAM file
echo "[1/5] Indexing BAM..."
samtools index "${filename}.bam"

# 5. Graph coverage
echo "[2/5] Calculating genome coverage (bedtools)..."
bedtools genomecov -ibam "${filename}.bam" -bg > "${filename}_autosomes.bedgraph"

# 6. Covdiv
echo "[3/5] Running covdiv against reference..."
covdiv \
    "${filename}_autosomes.bedgraph" \
    "${ref}" \
> "${filename}_autosomes_div.bedgraph"

# 7. Sort bedgraph
echo "[4/5] Sorting bedgraph..."
sort -t $'\t' -k1,1 -k2,2n -k3,3n "${filename}_autosomes_div.bedgraph" > "${filename}_autosomes_divsorted.bedgraph"

# 8. Windowed averages
echo "[5/5] Generating winavg (1k and 10k)..."
winavg -w 1000 -s 100 "${filename}_autosomes_divsorted.bedgraph" > "${filename}_autosomes_divsorted_win1k.bedgraph" 2>/dev/null
winavg -w 10000 -s 1000 "${filename}_autosomes_divsorted.bedgraph" > "${filename}_autosomes_divsorted_win10k.bedgraph" 2>/dev/null

echo "--- Pipeline Complete ---"
