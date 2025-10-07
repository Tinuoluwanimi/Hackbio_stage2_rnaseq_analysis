#!/bin/bash
# Run featureCounts on STAR-aligned BAM files 

MAPPED_DIR="plants/alignments/IGV"
OUT_DIR="plants/counts"
GTF="plants/references/Arabidopsis_thaliana.TAIR10.59.gtf"

mkdir -p $OUT_DIR
OUT_FILE="$OUT_DIR/counts.txt"

echo "Running featureCounts "
featureCounts -T 8 \
  -a $GTF \
  -o $OUT_FILE \
  $MAPPED_DIR/*Aligned.sortedByCoord.out.bam

cd $OUT_DIR
echo "FeatureCounts complete "


