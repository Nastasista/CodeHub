#!/bin/bash

echo
echo ------------------------------ START -----------------------------------
echo " -> Number of CPU cores to use: "
read CPU

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

echo " -> Create a reference genome index?"
read -rn2 -p " -> If yes, enter 'y'; if not, enter any other character: " continue_step2

if [[ $continue_step2 == "y" ]]; then
    hisat2-build "$ref_genome" ref_index
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
echo " -> Check the HTML reports, improve if necessary..."
read -rn2 -p " -> Quality check passed? If yes, enter 'y': " continue_step3

if [[ $continue_step3 != "y" ]]; then
    echo " -> Script execution terminated."
    exit 0
fi

echo ----------------------------- STEP 3 ------------------------------------
echo "            >>>>> ALIGNMENT TO INDEX -> .BAM <<<<<"
echo
echo " -> Creating RAW_ALN directory..."
mkdir -p RAW_ALN

find . -name "*.fastq" -type f | sort | while read -r fastq_file; do
    base_name=$(basename "${fastq_file}" .fastq)
    file_number=$(echo "${base_name}" | cut -d "_" -f 1)
    file_extension=$(echo "${base_name}" | cut -d "_" -f 2)

    echo "fastq_file: ${fastq_file}"
    echo "base_name: ${base_name}"

    if [ "${file_extension}" == "1" ]; then
        paired_file="${file_number}_2.fastq"

        echo "paired_file: ${paired_file}"

        # Check if the pair file exists
        if [ -f "${paired_file}" ]; then
            hisat2 -p "$CPU" -x REF_GENOM_INDEX/ref_index --dta -1 "${fastq_file}" -2 "${paired_file}" -S SRR.sam | samtools view -Sb > "RAW_ALN/$(basename ${base_name} _1)_aln.bam"
        fi
    fi
done

echo "RAW_ALN :"
ls -l RAW_ALN

echo
echo ----------------------------- STEP 4 ------------------------------------
echo "                   >>>>> SORTING .BAM <<<<<"
echo

# Create SORTED directory
mkdir -p SORTED

# Loop through all .bam files in RAW_ALN folder
for bam_file in RAW_ALN/*.bam; do
    # Get the base name of the file without the extension
    base_name=$(basename "$bam_file" .bam)

    # Perform samtools sort command
    samtools sort "$bam_file" -o "SORTED/${base_name}_sorted.bam"
done

echo "SORTED :"
ls -l SORTED

echo ----------------------------- STEP 5 ----------------------------------
echo "                >>>>> Choose Further Action <<<<<"
echo
select step in "HTSeq count" "featureCounts"; do
    case $step in
        "HTSeq count")
            echo
            echo "            >>>>> HTSeq count <<<<<"
            echo

            # Perform htseq-count on all .bam files in SORTED folder
            htseq-count -f bam SORTED/*.bam "$ref_genome" > matrix_htseq

            echo "HTSeq count results saved in matrix_htseq file"
            break
            ;;
        "featureCounts")
            echo
            echo "          >>>>> featureCounts <<<<<"
            echo

            # Perform featureCounts on all .bam files in SORTED folder
            featureCounts \
                -T 4 \
                -p --countReadPairs \
                -a *gtf -s 2 -g gene_id \
                -t exon \
                -o featCounts.txt \
                SORTED/*.bam

            echo "featureCounts results saved in featCounts.txt file"
            break
            ;;
        *)
            echo "Invalid choice. Enter 1 or 2."
            ;;
    esac
done

ls -l

echo
echo "                             >>>>> DONE <<<<<"

