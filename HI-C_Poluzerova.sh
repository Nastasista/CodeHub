#!/bin/bash

echo
echo ------------------------------ START -----------------------------------
echo " -> Number of CPU cores to use: "
read CPU

# Check if the entered value is a number
if ! [[ $CPU =~ ^[0-9]+$ ]]; then
    echo "Error: Please enter an integer."
    exit 1
fi

echo
echo ------------------------------ STEP 0 ----------------------------------
echo " -> Creating REF_GENOM_INDEX directory"
mkdir REF_GENOM_INDEX

echo " -> Add REFERENCE GENOME and INDEX (if any) here: REF_GENOM_INDEX"
read -rp " -> If added, enter 'y': " continue_step1

if [[ $continue_step1 != "y" ]]; then
    echo " -> Script execution terminated."
    exit 0
fi

echo
echo ------------------------------ STEP 1 -----------------------------------
echo "      >>>>> PREPARING REFERENCE GENOME AND INDEX <<<<<"
echo

# Search for .fa files in the REF_GENOM_INDEX folder
ref_genome=$(find REF_GENOM_INDEX -name "*.fa" -type f)
echo " -> REF_GENOM_INDEX:"
ls -l REF_GENOM_INDEX

echo " -> Do you need to create a reference genome index?"
read -rn2 -p " -> If yes, enter 'y'; if not, enter any other character: " continue_step2

if [[ $continue_step2 == "y" ]]; then
    bwa index "$ref_genome"
else
    echo " -> Continuing script execution without creating an index..."
fi

echo
echo ----------------------------- STEP 2 -----------------------------------
echo "                  >>>>> QUALITY CONTROL <<<<<"
echo
echo " -> Creating HTML QC reports..."

mkdir QC_html
fast_html="QC_html"

files=(*.fastq)

for file in "${files[@]}"; do
    fastqc "$file" -o "$fast_html"
done

# Open HTML files
outdir="QC_html"
for file in "$outdir"/*.html; do
    xdg-open "$file"
done
echo
echo "                 >>>>> QUALITY CHECK <<<<<"
echo
echo " -> Check the HTML reports..."
read -rn2 -p " -> Quality check passed? If yes, enter 'y': " continue_step3

if [[ $continue_step3 != "y" ]]; then
    echo " -> Script execution terminated."
    exit 0
fi

echo
echo ----------------------------- STEP 3 ------------------------------------
echo "            >>>>> ALIGNMENT TO INDEX -> .BAM <<<<<"
echo
echo " -> Creating RAW_ALN directory..."
mkdir -p RAW_ALN

find . -name "*.fastq" -type f | while read -r fastq_file; do
    base_name=$(basename "${fastq_file}" .fastq)
    bwa mem -t "$CPU" "$ref_genome" "$fastq_file" | samtools view -Sb > "RAW_ALN/${base_name}_aln.bam"
done

echo "RAW_ALN :"
ls -l RAW_ALN

echo
echo ----------------------------- STEP 4 ------------------------------------
echo "                   >>>>> '5' FILTERING <<<<<"
echo
echo " -> Creating FILT_ALN directory..."
mkdir FILT_ALN

find RAW_ALN -name "*.bam" -type f | while read -r bam_file; do
    base_name=$(basename "${bam_file}" .bam)
    samtools view -@ "$CPU" -h "$bam_file" | perl filter_5end.pl | samtools view -@ "$CPU" -Sb > "FILT_ALN/${base_name}_filt.bam"
done

echo "FILT_ALN :"
ls -l FILT_ALN

echo
echo ------------------------------ STEP 5 ------------------------------------
echo "                 >>>>> MERGING ALIGNMENTS <<<<<"
echo
echo " -> Creating COMB_ALN directory..."
mkdir COMB_ALN
perl two_read_bam_combiner.pl \
    FILT_ALN/*.bam |\
    samtools view -@ "$CPU" -Sb > COMB_ALN/aln.bam

echo "COMB_ALN :"
ls -l COMB_ALN

echo
echo ------------------------------ STEP 6 -----------------------------------
echo "                     >>>>> DUPLICATION REMOVAL <<<<<"
echo
echo " -> Creating DEDUP_ALN directory..."
mkdir DEDUP_ALN
mkdir DEDUP_ALN/tmp   # What is this for?

echo " -> Sorting alignments by coordinates..."
samtools sort \
    -@ "$CPU" \
    -T DEDUP_ALN/sort_aln.bam.tmp \
    -m2G \
    -O bam \
    -o DEDUP_ALN/sort_aln.bam \
    COMB_ALN/aln.bam

echo " -> Deduplicating using Picard..."
java -jar -Xmx6g -Djava.io.tmpdir=DEDUP_ALN/tmp picard.jar MarkDuplicates \
    -REMOVE_DUPLICATES true \
    -I DEDUP_ALN/sort_aln.bam \
    -O DEDUP_ALN/dedup_aln.bam \
    -M DEDUP_ALN/dedup_aln.metrics.txt \
    -ASSUME_SORT_ORDER coordinate \
    -MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1024

echo "DEDUP_ALN :"
ls -l DEDUP_ALN

echo
echo --------------------------------- STEP 7 -------------------------------------
echo "                    >>>>> STATISTICS EXTRACTION <<<<<"
echo
perl get_stats.pl DEDUP_ALN/dedup_aln.bam > DEDUP_ALN/dedups_bam.stats

echo " -> dedups_bam.stats"
cat DEDUP_ALN/dedups_bam.stats
echo
echo --------------------------------- STEP 8 --------------------------------------
echo "             >>>>> CHROMOSOME SIZE FILE PREPARATION <<<<<"
echo
echo " -> Creating COOLER directory..."

mkdir COOLER
python create_chrom_sizes.py "$ref_genome" COOLER
echo " -> Chromosome sizes text file generated"

echo "COOLER :"
ls -l COOLER

echo
echo --------------------------------- STEP 9 ------------------------------------------
echo  "   >>>>> CONVERT ALIGNMENT FILE TO CONTACT MATRIX TABULAR FORMAT <<<<<"
echo
label='N'
bam2pairs -c COOLER/chrom.sizes DEDUP_ALN/dedup_aln.bam COOLER/${label}

echo "COOLER :"
ls -l COOLER
echo
echo ---------------------------------- STEP 10 ---------------------------------------
echo "                     >>>>> BINNING USING COOLER <<<<<"
echo
echo " -> Creating contact map with a single resolution..."
cooler cload pairix \
    -p "$CPU" \
    COOLER/chrom.sizes:100000 \
    COOLER/${label}.bsorted.pairs.gz \
    COOLER/${label}_100k.cool
 echo " -> .cool file generated "
 echo "COOLER :"
ls -l COOLER

echo
echo ---------------------------------- STEP 11 --------------------------------------
echo "       >>>>> BALANCING AND CREATING ADDITIONAL MATRIX LAYERS <<<<<"
echo
cooler zoomify \
    -n "$CPU"  \
    -r 100000N, \
```bash
    -r 100000N \
    --balance \
    --balance-args "--nproc ${CPU}" \
    -o COOLER/${label}_multires.mcool \
    COOLER/N_100k.cool

echo " -> .mcool file generated"
echo " -> Contents of COOLER directory:"

ls -l COOLER
echo
echo "                             >>>>> DONE <<<<<"
