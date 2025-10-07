#!/bin/bash
# STAR alignment 
#Define directories

TRIM_DIR="plants/data/trimmed"
STAR_INDEX="plants/references/star_index"
ALIGNMENTS_DIR="plants/alignments"
IGV_DIR="$ALIGNMENTS_DIR/IGV"

mkdir -p $ALIGNMENTS_DIR $IGV_DIR

echo "Running STAR alignment for single-end reads "
for f in ${TRIM_DIR}/*_trim.fastq.gz; do
    base=$(basename $f _trim.fastq.gz)
    outPrefix=${ALIGNMENTS_DIR}/${base}_

    STAR --runThreadN 8 \
         --genomeDir $STAR_INDEX \
         --readFilesIn $f \
         --readFilesCommand zcat \
         --outFileNamePrefix $outPrefix \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes All

    cp ${outPrefix}Aligned.sortedByCoord.out.bam $IGV_DIR/
done

# Index BAMs
for bam in $IGV_DIR/*.bam; do
    samtools index $bam
done

echo "STAR alignment complete "
