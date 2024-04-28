params.CPU = 20

 STEP 0 Creating REF_GENOM_INDEX directory
process createRefGenomeIndex {
    output
    directory REF_GENOM_INDEX

    script
    
    mkdir -p REF_GENOM_INDEX
    
}

 STEP 1 PREPARING REFERENCE GENOME AND INDEX
process prepareRefGenomeAndIndex {
    input
    directory refGenomeIndex

    output
    file(refGenome) into refGenomeChannel

    script
    
    refGenome=$(find REF_GENOM_INDEX -name .fa -type f)
    echo $refGenome  refGenome.txt
    

    errorStrategy ignore
}

 STEP 2 QUALITY CONTROL
process qualityControl {
    input
    file refGenome from refGenomeChannel
    path '.fastq' from rawFastq

    output
    directory 'QC_html' into qcHtml
    file 'QC_html.html' into htmlReports

    script
    
    mkdir -p QC_html
    fast_html=QC_html

    for file in .fastq; do
        fastqc $file -o $fast_html
    done
    

    errorStrategy ignore
}

 STEP 3 ALIGNMENT TO INDEX - .BAM
process alignmentToBam {
    input
    file refGenome from refGenomeChannel
    path '.fastq' from rawFastq

    output
    directory 'RAW_ALN' into rawAlignment

    script
    
    mkdir -p RAW_ALN

    for file in .fastq; do
        base_name=$(basename $file .fastq)
        file_number=$(echo $base_name  cut -d _ -f 1)
        file_extension=$(echo $base_name  cut -d _ -f 2)

        if [ $file_extension == 1 ]; then
            paired_file=${file_number}_2.fastq

            if [ -f $paired_file ]; then
                hisat2 -p ${params.CPU} -x REF_GENOM_INDEXref_index --dta -1 $file -2 $paired_file -S SRR.sam  samtools view -Sb  RAW_ALN$(basename $base_name _1)_aln.bam
            fi
        fi
    done
    

    errorStrategy ignore
}

 STEP 4 SORTING .BAM
process sortBam {
    input
    path '.bam' from rawAlignment

    output
    directory 'SORTED' into sortedBam

    script
    
    mkdir -p SORTED

    for bam_file in RAW_ALN.bam; do
        base_name=$(basename $bam_file .bam)
        samtools sort $bam_file -o SORTED${base_name}_sorted.bam
    done
    

    errorStrategy ignore
}

 Choose Further Action
process furtherAction {
    input
    file refGenome from refGenomeChannel
    path '.bam' from sortedBam

    script
    
    select step in HTSeq count featureCounts; do
        case $step in
            HTSeq count)
                htseq-count -f bam SORTED.bam $refGenome  matrix_htseq
                break
                ;;
            featureCounts)
                featureCounts 
                    -T ${params.CPU} 
                    -p --countReadPairs 
                    -a gtf -s 2 -g gene_id 
                    -t exon 
                    -o featCounts.txt 
                    SORTED.bam
                break
                ;;
            )
                echo Invalid choice.
                ;;
        esac
    done
    

    errorStrategy ignore
}

workflow {
     Execute processes in order
    createRefGenomeIndex
    prepareRefGenomeAndIndex(refGenomeIndex)
    qualityControl(refGenome)
    alignmentToBam(refGenome, rawFastq)
    sortBam(rawAlignment)
    furtherAction(refGenome, sortedBam)
}
